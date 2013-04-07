objects = main.o

atom : $(objects)
	g++ $(objects) -o atom

main.o : main.cpp
	g++ -c main.cpp

clean :
	rm $(objects) atom