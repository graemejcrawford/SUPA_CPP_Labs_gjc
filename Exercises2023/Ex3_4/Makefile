CC=g++ #Name of compiler
FLAGS=-std=c++20 -w #Compiler flags (the s makes it silent)
TARGET=FunctionScriptFinal #Executable name
OBJECTS=FunctionScriptFinal.o FiniteFunctionsFinal.o 
LIBS=-I ../../GNUplot/ -lboost_iostreams

#First target in Makefile is default
${TARGET}:${OBJECTS} #Make target from objects
	@echo "Linking..."
	${CC} ${FLAGS} ${OBJECTS} ${LIBS} -o ${TARGET}
	@make clean

FunctionScriptFinal.o : FunctionScriptFinal.cxx FiniteFunctionsFinal.h 
	${CC} ${FLAGS} ${LIBS} -c FunctionScriptFinal.cxx

FiniteFunctionsFinal.o : FiniteFunctionsFinal.cxx FiniteFunctionsFinal.h
	${CC} ${FLAGS} ${LIBS} -c FiniteFunctionsFinal.cxx

clean: #No targets just run shell command to rm object files and emacs backups
	@rm -f *.o *~

cleantarget: #Delete the exectuable
	@rm -f ${TARGET}
