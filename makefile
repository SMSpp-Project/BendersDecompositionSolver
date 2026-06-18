##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of BendersDecompositionSolver                                   #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by BendersDecompositionSolver, i.e., core SMS++ (which   #
#   provides BendersBFunction). The master and the subproblem Solver are      #
#   instantiated through the Solver factory at run time, hence they are not   #
#   compile-time dependencies of this module.                                #
#                                                                            #
#   Note that $(SMS++INC) is assumed to include any -I directive             #
#   corresponding to external libraries needed by SMS++, at least to the      #
#   extent in which they are needed by the parts of SMS++ used by             #
#   BendersDecompositionSolver.                                              #
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(SMS++INC)    = the -I$( core SMS++ directory )                 #
#           $(SMS++OBJ)    = the core SMS++ library                          #
#           $(BnDSLVSDR)   = the directory where the source is               #
#                                                                            #
#   Output: $(BnDSLVOBJ)   = the final object(s) / library                   #
#           $(BnDSLVH)     = the .h files to include                         #
#           $(BnDSLVINC)   = the -I$( source directory )                     #
#                                                                            #
#                              Antonio Frangioni                             #
#                                Donato Meoli                                #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

BnDSLVOBJ = $(BnDSLVSDR)/obj/BendersDecompositionSolver.o

BnDSLVINC = -I$(BnDSLVSDR)/include

BnDSLVH   = $(BnDSLVSDR)/include/BendersDecompositionSolver.h

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(BnDSLVOBJ) $(BnDSLVSDR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(BnDSLVSDR)/obj/BendersDecompositionSolver.o: \
	$(BnDSLVSDR)/src/BendersDecompositionSolver.cpp \
	$(BnDSLVSDR)/include/BendersDecompositionSolver.h $(SMS++OBJ)
	$(CC) -c $(BnDSLVSDR)/src/BendersDecompositionSolver.cpp -o $@ \
	$(BnDSLVINC) $(SMS++INC) $(SW)

########################## End of makefile ###################################
