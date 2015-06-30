


all: zeroMix.exe

zeroMix.exe: zeroMix.o segyIO_class.o
	gcc -O2 -o  zeroMix.exe zeroMix.o segyIO_class.o -lm -fopenmp

zeroMix.o: zeroMix.c
	gcc -O2 -c zeroMix.c -lm  -fopenmp

segyIO_class.o: segyIO_class.c
	gcc -O2 -c segyIO_class.c 


