all : bsprmtxex.exe
	./bsprmtxex.exe


bsprmtxex.exe : ./bsprmtxex.c
	gcc -o ./bsprmtxex.exe ./bsprmtxex.c -lm -O
