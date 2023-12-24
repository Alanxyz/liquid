liquid: main.c
	gcc -std=c99 -o liquid main.c

test:
	gcc -std=c99 -g -o liquid main.c
	./liquid
