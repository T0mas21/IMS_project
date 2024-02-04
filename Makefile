## IMS 
# T9: Spojitý model z oblasti fyziky a biologie
# Téma: Hemolýza erytrocytů v důsledku ozvučení ultrazvukem
# Brief: Makefile
# Author: Tomáš Janečka(xjanec35)
# Author: Jan Novák (xnovak3i)
# Date: 10.12.2022


CC = g++
CFLAGS = -Wall -std=c++11

all: main

main: main.o
	$(CC) $(CFLAGS) -o main main.o -lsimlib -lm

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp -lsimlib -lm

run: main
	./main 76 285 1 0.1 > pokus_1.txt 
	./main 125 0 1 0.1 > pokus_2.txt
	
	./main 30 100 1 0.1 > pokus_3.txt
	
	./main 80 100 1 0.1 > pokus_4.txt

	./main 40 50 5 2 > pokus_5.txt
	./main 40 50 10 10 > pokus_6.txt


clean:
	rm -f main main.o pokus_*.txt