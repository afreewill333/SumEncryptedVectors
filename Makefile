COMPILER = $(CXX)
CPPFLAGS=-Wall -DNDEBUG -std=c++11 -g2 -I$(HELIB) -pthread -DFHE_THREADS -DFHE_DCRT_THREADS -DFHE_BOOT_THREADS # -fopenmp
LIBS=$(HELIB)/fhe.a -lntl -lgmp -lm -lrt
OBJS = paillier.o
TARGET=vectorSum

.PHONY: all clean run

all: $(TARGET)

$(TARGET): main.cpp $(OBJS)
	$(COMPILER) $(CPPFLAGS) $(OBJS) $< -o $@ $(LIBS)

paillier.o: paillier.c paillier.h 
	$(COMPILER) $(CPPFLAGS) $< -c

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)
