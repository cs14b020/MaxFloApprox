main: kernels.o main.o
	nvcc -Wno-deprecated-gpu-targets kernels.o main.o -o main

kernels.o: kernels.cu kernels.h
	nvcc -Wno-deprecated-gpu-targets -std=c++11 -c kernels.cu

main.o: main.cu
	nvcc -Wno-deprecated-gpu-targets -std=c++11 -c main.cu

clean: 
	rm -rf *.o main