CUDA_HOME ?= /usr/local/cuda
EX_HOME ?= /home/deepl/sooho/test/K-com/test

all: revise4t

revise4t: main.cu include/revise4t.cu
	$(CUDA_HOME)/bin/nvcc -g -lineinfo -gencode arch=compute_70,code=sm_70 include/revise4t.cu main.cu -o /home/deepl/sooho/test/K-com/test/exec/revise4t -I$(CUDA_HOME)/include,$(EX_HOME)/include,$(EX_HOME)/.. -L$(CUDA_HOME)/lib64,pthread -std=c++11