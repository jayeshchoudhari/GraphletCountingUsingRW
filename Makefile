main-organized-make: mainCount.o count3Graphlets.o count4Graphlets.o graphIO.o utilities.o
	g++ -O3 mainCount.o count3Graphlets.o count4Graphlets.o  graphIO.o utilities.o -o main-organized-make

mainCount.o: mainCount.cpp  ./include/namespace.h
	g++ -O3 -c mainCount.cpp

utilities.o: ./include/utilities.cpp  ./include/utilities.h  ./include/namespace.h
	g++ -O3 -c ./include/utilities.cpp

count3Graphlets.o: ./include/count3Graphlets.cpp ./include/count3Graphlets.h ./include/graphIO.h ./include/namespace.h 
	g++ -O3 -c ./include/count3Graphlets.cpp

count4Graphlets.o: ./include/count4Graphlets.cpp  ./include/count4Graphlets.h  ./include/graphIO.h ./include/namespace.h
	g++ -O3 -c ./include/count4Graphlets.cpp

graphIO.o: ./include/graphIO.cpp  ./include/graphIO.h  ./include/namespace.h
	g++ -O3 -c ./include/graphIO.cpp

clean:
	rm *.o main-organized-make