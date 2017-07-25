# location of the Python header files
PYTHON_VERSION = 2.7
PYTHON_INCLUDE = /p/home/sblair/anaconda2/include/python$(PYTHON_VERSION)
PYTHON_LIB=/p/home/sblair/anaconda2/lib
BOOST_PYLIB = py27
 
# location of the Boost Python include files and library
 
BOOST_INC = /app/COST/boost/1.58.0/gnu/include
BOOST_LIB = /app/COST/boost/1.58.0/gnu/lib

# compile mesh classes
TARGET = LBM_Interface
FILE=PyLBM_Interface
EXT=cpp

MPI_CC=CC
ifeq ($(PE_ENV),PGI)
	MPI_FLAGS= -std=c++11 -fast -acc -ta=tesla:cc35 -fPIC
else
	MPI_FLAGS=-O3 -hacc -hlist=m
endif
MY_LIBS= -lboost_python-$(BOOST_PYLIB) -lpython$(PYTHON_VERSION)

SOURCES= Lattice.cpp D3Q15Lattice.cpp D3Q19Lattice.cpp D3Q27Lattice.cpp LBM_DataHandler.cpp \
	LBM_HaloData.cpp LBM_HaloDataOrganizer.cpp  
OBJECTS= Lattice.o  D3Q15Lattice.o D3Q19Lattice.o D3Q27Lattice.o LBM_DataHandler.o  \
	LBM_HaloData.o LBM_HaloDataOrganizer.o
	 

$(FILE).so: $(FILE).o $(OBJECTS)
	$(MPI_CC) $(MPI_FLAGS) -shared -Wl,--export-dynamic $(FILE).o -L$(BOOST_LIB) -lboost_python -L$(PYTHON_LIB) -lpython$(PYTHON_VERSION) -o $(TARGET).so $(OBJECTS)
 
$(FILE).o: $(FILE).cpp
	$(MPI_CC) $(MPI_FLAGS) -I$(PYTHON_INCLUDE) -I$(BOOST_INC)  -c $(FILE).$(EXT)

%.o:%.cpp
	$(MPI_CC) $(MPI_FLAGS)  -c $^

clean:
	rm -f *.o *.so $(TARGET) *~



