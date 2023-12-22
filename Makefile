liquid: main.c
	gcc -std=c99 -o liquid main.c

test:
	gcc -std=c99 -o liquid main.c
	./main
