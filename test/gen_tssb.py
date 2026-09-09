#!/usr/bin/env python3
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - gen_tssb.py - - - - - - - - - - - - - - - - -
#
# writes an InvestmentBlock over a TwoStageStochasticBlock, i.e., the instance
# the tssb test reads, with as many scenarios and time steps as asked: one
# asset, whose design sizes an intermittent unit, and a slack unit that makes
# the subproblem always feasible, so that what is being compared is how the
# two Benders forms handle the optimality cuts alone
#
# usage: gen_tssb.py <file> <scenarios> <time steps> [ seed ]
#
#   Donato Meoli
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

import sys

import numpy as np
from netCDF4 import Dataset

fn = sys.argv[ 1 ]
S = int( sys.argv[ 2 ] )
T = int( sys.argv[ 3 ] )
seed = int( sys.argv[ 4 ] ) if len( sys.argv ) > 4 else 1

rng = np.random.default_rng( seed )

# a daily profile the scenarios perturb, so that they are neither equal nor
# unrelated: an instance whose scenarios say the same thing is no test of a
# decomposition by scenario
base = 6000 + 3000 * np.sin( np.linspace( 0 , 2 * np.pi , T , endpoint=False ) )
demand = np.maximum( 500 , base + rng.normal( 0 , 1200 , ( S , T ) ) )

peak = float( demand.max() )

f = Dataset( fn , "w", format="NETCDF4" )
f.setncattr( "SMS++_file_type" , np.int64( 1 ) )

b = f.createGroup( "Block_0" )
b.type = "InvestmentBlock"
b.id = "0"
b.createDimension( "NumAssets" , 1 )
b.createVariable( "Cost" , "f8" , ( "NumAssets" , ) )[ : ] = [ 55.0 ]
b.createVariable( "LowerBound" , "f8" , ( "NumAssets" , ) )[ : ] = [ 1200.0 ]
b.createVariable( "UpperBound" , "f8" , ( "NumAssets" , ) )[ : ] = [ 100 * peak ]
b.createVariable( "InstalledQuantity" , "f8" , ( "NumAssets" , ) )[ : ] = [ 0.0 ]
b.createVariable( "Assets" , "u8" , ( "NumAssets" , ) )[ : ] = [ 0 ]
b.createVariable( "AssetType" , "i8" , ( "NumAssets" , ) )[ : ] = [ 0 ]

i = b.createGroup( "InnerBlock" )
i.type = "TwoStageStochasticBlock"
i.id = "0"
i.createDimension( "NumberScenarios" , S )

d = i.createGroup( "DiscreteScenarioSet" )
d.type = "DiscreteScenarioSet"
d.createDimension( "NumberScenarios" , S )
d.createDimension( "ScenarioSize" , T )
d.createVariable( "Scenarios" , "f8" ,
                  ( "NumberScenarios" , "ScenarioSize" ) )[ : ] = demand
d.createVariable( "PoolWeights" , "f8" ,
                  ( "NumberScenarios" , ) )[ : ] = np.full( S , 1.0 / S )

s = i.createGroup( "StochasticBlock" )
s.type = "StochasticBlock"
s.createDimension( "NumberDataMappings" , 1 )
s.createDimension( "SetSize_dim" , 2 )
s.createDimension( "SetElements_dim" , 4 )
s.createVariable( "FunctionName" , str , ( "NumberDataMappings" , )
                  )[ 0 ] = "UCBlock::set_active_power_demand"
s.createVariable( "Caller" , "S1" , ( "NumberDataMappings" , ) )[ : ] = b"B"
s.createVariable( "DataType" , "S1" , ( "NumberDataMappings" , ) )[ : ] = b"D"
s.createVariable( "SetSize" , "u4" , ( "SetSize_dim" , ) )[ : ] = [ 0 , 0 ]
s.createVariable( "SetElements" , "u4" ,
                  ( "SetElements_dim" , ) )[ : ] = [ 0 , T , 0 , T ]

p = s.createGroup( "AbstractPath" )
p.type = "AbstractPath"
p.createDimension( "PathDim" , 1 )
p.createDimension( "TotalLength" , None )
p.createVariable( "PathStart" , "u4" , ( "PathDim" , ) )[ : ] = [ 0 ]
p.createVariable( "PathNodeTypes" , "S1" , ( "TotalLength" , ) )
p.createVariable( "PathGroupIndices" , "u4" , ( "TotalLength" , ) )
p.createVariable( "PathElementIndices" , "u4" , ( "TotalLength" , ) )
p.createVariable( "PathRangeIndices" , "u4" , ( "TotalLength" , ) )

u = s.createGroup( "Block" )
u.type = "UCBlock"
u.id = "0"
u.createDimension( "TimeHorizon" , T )
u.createDimension( "NumberUnits" , 2 )
u.createDimension( "NumberElectricalGenerators" , 2 )
u.createDimension( "NumberNodes" , 1 )
u.createDimension( "NumberLines" , None )
u.createVariable( "ActivePowerDemand" , "f8" ,
                  ( "NumberNodes" , "TimeHorizon" ) )[ : ] = demand[ 0 ]
u.createVariable( "NodeName" , str , ( "NumberNodes" , ) )[ 0 ] = "b"
u.createVariable( "LineName" , str , ( "NumberLines" , ) )
u.createVariable( "GeneratorNode" , "i8" ,
                  ( "NumberElectricalGenerators" , ) )[ : ] = [ 0 , 0 ]

# the unit the design sizes: one unit of power per unit of design
g = u.createGroup( "UnitBlock_0" )
g.type = "IntermittentUnitBlock"
g.setncattr( "name" , "diesel" )
g.createVariable( "MaxPower" , "f8" , ( "TimeHorizon" , ) )[ : ] = np.ones( T )
g.createVariable( "MinPower" , "f8" , ( "TimeHorizon" , ) )[ : ] = np.zeros( T )
g.createVariable( "ActivePowerCost" , "f8" ,
                  ( "TimeHorizon" , ) )[ : ] = np.full( T , 3.0 )
g.createVariable( "MinGeneration" , "f8" )[ : ] = -np.inf
g.createVariable( "MaxGeneration" , "f8" )[ : ] = np.inf

# the slack, which serves whatever the design cannot: it is what makes the
# subproblem feasible for every design
g = u.createGroup( "UnitBlock_1" )
g.type = "SlackUnitBlock"
g.setncattr( "name" , "slack_unit b" )
g.createVariable( "ActivePowerCost" , "f8" ,
                  ( "TimeHorizon" , ) )[ : ] = np.full( T , 10000.0 )
g.createVariable( "MaxPower" , "f8" ,
                  ( "TimeHorizon" , ) )[ : ] = np.full( T , 50 * peak )

f.close()
print( "%s: %d scenarios, %d time steps, peak %.1f" % ( fn , S , T , peak ) )
