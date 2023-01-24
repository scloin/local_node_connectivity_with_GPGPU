CUDA_HOME ?= /usr/local/cuda
EX_HOME ?= /home/deepl/sooho/test/K-com/test
CODE_HOME ?= /home/deepl/sooho/test/K-com/test/code

all: frontier_thd frontier no_frontier

no_frontier: $(CODE_HOME)/CUDA/no_frontier/main.cu $(CODE_HOME)/CUDA/no_frontier/include/revise5t.cu
	$(CUDA_HOME)/bin/nvcc -g -lineinfo $(CODE_HOME)/CUDA/no_frontier/include/revise5t.cu $(CODE_HOME)/CUDA/no_frontier/main.cu -o $(CODE_HOME)/exec/no_frontier -I$(CUDA_HOME)/include,$(CODE_HOME)/CUDA/no_frontier/include,$(EX_HOME)/.. -L$(CUDA_HOME)/lib64,pthread -std=c++11

frontier_thd: $(CODE_HOME)/CUDA/frontier/main.cu $(CODE_HOME)/CUDA/frontier/include/revise4t.cu
	$(CUDA_HOME)/bin/nvcc -g -lineinfo $(CODE_HOME)/CUDA/frontier/include/revise4t.cu $(CODE_HOME)/CUDA/frontier/main.cu -o $(CODE_HOME)/exec/frontier_thd -I$(CUDA_HOME)/include,$(CODE_HOME)/CUDA/frontier/include,$(EX_HOME)/.. -L$(CUDA_HOME)/lib64,pthread -std=c++11

frontier: $(CODE_HOME)/CUDA/frontier/revise4.cu
	$(CUDA_HOME)/bin/nvcc -g -lineinfo $(CODE_HOME)/CUDA/revise4.cu -o $(CODE_HOME)/exec/frontier -I$(CUDA_HOME)/include,$(EX_HOME)/.. -L$(CUDA_HOME)/lib64 -std=c++11