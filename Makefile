all: cylmain

cylmain: cylmain.c
	gcc -I. cylmain.c -lOpenCL -o cylmain

add: main.c
	gcc -I. main.c -lOpenCL -o add
