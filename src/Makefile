SWIG=swig
CC=g++
libraries = -lpython2.7
hostname = $(shell hostname)

ifeq ($(hostname),rneher-iMac)
	#RICHARD_linux
	pythonpath := /usr/local/EPD/include/python2.7/
	numpypath := /usr/local/EPD/lib/python2.7/site-packages/numpy/core/include/
	LFLAGS = -L/usr/local/EPD/lib/python2.7/config   -lpython2.7 -O2 -shared
	SWIG = /usr/bin/swig
else
	#CLUSTER
	pythonpath := /ebio/ag-neher/share/epd_free-7.1-2-rh5-x86_64/include/python2.7/
	numpypath := /ebio/ag-neher/share/epd_free-7.1-2-rh5-x86_64/lib/python2.7/site-packages/numpy/core/include/
	LFLAGS = -L/ebio/ag-neher/share/epd_free-7.1-2-rh5-x86_64/lib/python2.7/config -L/ebio/ag-neher/share/programs/lib  -lpython2.7 -O2 -shared
endif

ifeq ($(hostname),ag-neher-xps13)
	#RICHARD_laptop
	pythonpath := /usr/local/EPD/include/python2.7/
	numpypath := /usr/local/EPD/lib/python2.7/site-packages/numpy/core/include/
	LFLAGS = -L/usr/local/EPD/lib/python2.7/config   -lpython2.7 -O2 -shared
endif

#RICHARD
#pythonpath := /Library/Frameworks/EPD64.framework/Versions/7.0/include/python2.7/
#numpypath := /Library/Frameworks/EPD64.framework/Versions/7.0/lib/python2.7/site-packages/numpy/core/include/
#LFLAGS = -L/Library/Frameworks/EPD64.framework/Versions/7.0/lib/python2.7/config -ldl -framework CoreFoundation -lpython2.7 -O2 -shared

#LFLAGS = -L/Library/Frameworks/Python.framework/Versions/Current/lib/python2.7/config -ldl -framework CoreFoundation -lpython2.7 -O2 -shared -arch i386
#LFLAGS = -L/Library/Frameworks/EPD64.framework/Versions/7.0/lib/python2.7/config -ldl -framework CoreFoundation -lpython2.7 -O2 -shared
#LFLAGS = -L/Library/Frameworks/EPD64.framework/Versions/7.0/lib/python2.7/config -ldl
#CFLAGS = -I$(pythonpath) -I$(numpypath) -O2 -fPIC -c -arch i386


CFLAGS = -I$(pythonpath) -I$(numpypath) -O2 -fPIC -c


SWIG_INTERFACE := cfit.i
SWIG_SOURCE = $(SWIG_INTERFACE:.i=_wrap.cpp)
SWIG_OBJECT = $(SWIG_SOURCE:.cpp=.o)

HEADERS = cfit.h
SOURCES := cfit.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = _$(SOURCES:.cpp=.so)

all: $(SOURCES) $(EXECUTABLE)

#$(EXECUTABLE): $(SOURCES) $(HEADERS)
#	./setup.py build_ext --inplace
#

$(EXECUTABLE): $(OBJECTS) $(SWIG_OBJECT)
	$(CC) $(LFLAGS) $^ $(libraries) -o $@

$(OBJECTS) $(SWIG_OBJECT) : $(HEADERS) $(SOURCES) $(SWIG_SOURCE)
	$(CC) $(CFLAGS) $(HEADERS) $(SOURCES) $(SWIG_SOURCE)

$(SWIG_SOURCE): $(SWIG_INTERFACE) $(SOURCES) $(HEADERS)
	echo $(hostname)
	echo $(pythonpath)
	$(SWIG) -c++ -module cfit -o $(SWIG_SOURCE)  -python $(SWIG_INTERFACE)


clean:
	rm -rf *.o $(EXECUTABLE) $(SWIG_SOURCE) cfit*.py{,c,o} *.gch
	echo $(hostname)
