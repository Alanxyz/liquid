liquid: main.c
	gcc -std=c99 -o main main.c

test:
	gcc -std=c99 -o main main.c
	./main
