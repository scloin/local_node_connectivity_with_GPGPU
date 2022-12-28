CUDA_HOME ?= /usr/local/cuda
EX_HOME ?= /home/deepl/sooho/test/K-com/test
CODE_HOME ?= /home/deepl/sooho/test/K-com/test/code

all: revise4t revise4

revise4t: $(CODE_HOME)/CUDA/main.cu $(CODE_HOME)/CUDA/include/revise4t.cu
	$(CUDA_HOME)/bin/nvcc -g -lineinfo -gencode arch=compute_70,code=sm_70 $(CODE_HOME)/CUDA/include/revise4t.cu $(CODE_HOME)/CUDA/main.cu -o $(CODE_HOME)/exec/revise4t -I$(CUDA_HOME)/include,$(CODE_HOME)/CUDA/include,$(EX_HOME)/.. -L$(CUDA_HOME)/lib64,pthread -std=c++11

revise4: $(CODE_HOME)/CUDA/revise4.cu
	$(CUDA_HOME)/bin/nvcc -g -lineinfo -gencode arch=compute_70,code=sm_70 $(CODE_HOME)/CUDA/revise4.cu -o $(CODE_HOME)/exec/revise4 -I$(CUDA_HOME)/include,$(EX_HOME)/.. -L$(CUDA_HOME)/lib64 -std=c++11