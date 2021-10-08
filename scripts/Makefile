objects = main.o 
photonTreeProducer.exe : $(objects)
	g++ -o photonTreeProducer.exe $(objects) `root-config --libs`
main.o : main.C photonTreeProducer.h photonTreeProducer.C
	g++ -Wall -Wno-deprecated `root-config --cflags` main.C -c
clean :
	rm photonTreeProducer.exe $(objects)

