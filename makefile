all:	propagacao.o	main.o
	g++ -o exe propagacao.o main.o
propagacao.o: propagacao.cpp propagacao.h
	g++ -c propagacao.cpp
main.o: propagacao.h main.cpp
	g++ -c main.cpp
go:
	./exe
clean:
	rm *.o
cleanall:
	rm exe && rm *.o