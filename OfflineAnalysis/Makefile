CC=g++
CFLAGS=-c -Wall `root-config --cflags`
LDFLAGS=`root-config --libs` `root-config --glibs` -lgmp -lm
SOURCES=noWaveformSensitivities.cc
OBJECTS=$(SOURCES:.cc=.o)
EXECUTABLE=noWaveformSensitivities

all: $(SOURCES) $(EXECUTABLE)
    
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cc.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm ./*.o ./MyMain

