all: main-DiffStartPoint

main-DiffStartPoint: mainCount_DiffStartPoint.o countCliques.o count3Graphlets.o count4Graphlets.o count5Graphlets.o SRWcount5cliques.o graphIO.o utilities.o
	g++ -O3 mainCount_DiffStartPoint.o countCliques.o count3Graphlets.o count4Graphlets.o count5Graphlets.o SRWcount5cliques.o graphIO.o utilities.o -o main-DiffStartPoint

main-SingleRW: mainCount_SingleRW.o countCliques.o count3Graphlets.o count4Graphlets.o count5Graphlets.o graphIO.o utilities.o
	g++ -O3 mainCount_SingleRW.o countCliques.o count3Graphlets.o count4Graphlets.o count5Graphlets.o graphIO.o utilities.o -o main-SingleRW

main-organized-make: mainCount.o countCliques.o count3Graphlets.o count4Graphlets.o count5Graphlets.o graphIO.o utilities.o
	g++ -O3 mainCount.o countCliques.o count3Graphlets.o count4Graphlets.o count5Graphlets.o graphIO.o utilities.o -o main-organized-make

g2Degree: computeSRWG2Degree.o graphIO.o utilities.o
	g++ -O3 computeSRWG2Degree.o graphIO.o utilities.o -o g2Degree

mainCount.o: mainCount.cpp  ./include/namespace.h
	g++ -O3 -c mainCount.cpp

mainCount_DiffStartPoint.o: mainCount_DiffStartPoint.cpp  ./include/namespace.h
	g++ -O3 -c mainCount_DiffStartPoint.cpp

mainCount_SingleRW.o: mainCount_SingleRW.cpp  ./include/namespace.h
	g++ -O3 -c mainCount_SingleRW.cpp



SRWcount5cliques.o: ./include/SRWcount5cliques.cpp  ./include/count5Graphlets.h  ./include/utilities.h ./include/utilities.cpp ./include/graphIO.h ./include/namespace.h
	g++ -O3 -c ./include/SRWcount5cliques.cpp	

computeSRWG2Degree.o: ./computeSRWG2Degree.cpp  ./include/namespace.h ./include/graphIO.h  ./include/utilities.h
	g++ -O3 -c computeSRWG2Degree.cpp

count5Graphlets.o: ./include/count5Graphlets.cpp  ./include/count5Graphlets.h  ./include/graphIO.h ./include/namespace.h
	g++ -O3 -c ./include/count5Graphlets.cpp

count4Graphlets.o: ./include/count4Graphlets.cpp  ./include/count4Graphlets.h  ./include/graphIO.h ./include/namespace.h
	g++ -O3 -c ./include/count4Graphlets.cpp

count3Graphlets.o: ./include/count3Graphlets.cpp ./include/count3Graphlets.h ./include/graphIO.h ./include/namespace.h 
	g++ -O3 -c ./include/count3Graphlets.cpp

countCliques.o: ./include/countCliques.cpp ./include/countCliques.h ./include/graphIO.h ./include/namespace.h 
	g++ -O3 -c ./include/countCliques.cpp

graphIO.o: ./include/graphIO.cpp  ./include/graphIO.h  ./include/namespace.h
	g++ -O3 -c ./include/graphIO.cpp

utilities.o: ./include/utilities.cpp  ./include/utilities.h  ./include/namespace.h
	g++ -O3 -c ./include/utilities.cpp

clean:
	rm *.o main-organized-make main-SingleRW main-DiffStartPoint
