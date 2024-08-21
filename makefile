CUDA_HOME ?= /usr/local/cuda
EX_HOME ?= /home/deepl/sooho/local_node_connectivity_with_GPGPU
CODE_HOME ?= /home/deepl/sooho/local_node_connectivity_with_GPGPU/code

all: frontier_thd no_frontier cpptest frontier_thd_test no_frontier_test frontier justgpu justgpu_test

test: frontier_thd_test no_frontier_test

cpp: cpptest

just: justgpu justgpu_test

no_frontier: $(CODE_HOME)/CUDA/no_frontier/main.cu $(CODE_HOME)/CUDA/no_frontier/include/revise5t.cu
	$(CUDA_HOME)/bin/nvcc -g -lineinfo $(CODE_HOME)/CUDA/no_frontier/include/revise5t.cu $(CODE_HOME)/CUDA/no_frontier/main.cu -o $(CODE_HOME)/exec/no_frontier -I$(CUDA_HOME)/include,$(CODE_HOME)/CUDA/no_frontier/include,$(EX_HOME)/.. -L$(CUDA_HOME)/lib64,pthread -std=c++11

frontier_thd: $(CODE_HOME)/CUDA/frontier/main.cu $(CODE_HOME)/CUDA/frontier/include/revise4t.cu
	$(CUDA_HOME)/bin/nvcc -g -lineinfo $(CODE_HOME)/CUDA/frontier/include/revise4t.cu $(CODE_HOME)/CUDA/frontier/main.cu -o $(CODE_HOME)/exec/frontier_thd -I$(CUDA_HOME)/include,$(CODE_HOME)/CUDA/frontier/include,$(EX_HOME)/.. -L$(CUDA_HOME)/lib64,pthread -std=c++11

frontier: $(CODE_HOME)/CUDA/frontier/revise4.cu
	$(CUDA_HOME)/bin/nvcc -g -lineinfo $(CODE_HOME)/CUDA/frontier/revise4.cu -o $(CODE_HOME)/exec/frontier -I$(CUDA_HOME)/include,$(EX_HOME)/.. -L$(CUDA_HOME)/lib64 -std=c++11

frontier_thd_test: $(CODE_HOME)/CUDA/frontier/test.cu $(CODE_HOME)/CUDA/frontier/include/revise4t.cu
	$(CUDA_HOME)/bin/nvcc -g -lineinfo $(CODE_HOME)/CUDA/frontier/include/revise4t.cu $(CODE_HOME)/CUDA/frontier/test.cu -o $(CODE_HOME)/exec/frontier_thd_test -I$(CUDA_HOME)/include,$(CODE_HOME)/CUDA/frontier/include,$(EX_HOME)/.. -L$(CUDA_HOME)/lib64,pthread -std=c++11

no_frontier_test: $(CODE_HOME)/CUDA/no_frontier/test.cu $(CODE_HOME)/CUDA/no_frontier/include/revise5t.cu
	$(CUDA_HOME)/bin/nvcc -g -lineinfo $(CODE_HOME)/CUDA/no_frontier/include/revise5t.cu $(CODE_HOME)/CUDA/no_frontier/test.cu -o $(CODE_HOME)/exec/no_frontier_test -I$(CUDA_HOME)/include,$(CODE_HOME)/CUDA/no_frontier/include,$(EX_HOME)/.. -L$(CUDA_HOME)/lib64,pthread -std=c++11

cpptest: $(CODE_HOME)/cpp/main.cpp $(CODE_HOME)/cpp/include/cyver_t.cpp
	g++ -g -Wall $(CODE_HOME)/cpp/include/cyver_t.cpp $(CODE_HOME)/cpp/main.cpp -o $(CODE_HOME)/exec/cyver_t -I$(CODE_HOME)/cpp/include -std=c++11


justgpu: $(CODE_HOME)/CUDA/justgpu/main.cu $(CODE_HOME)/CUDA/justgpu/include/justgpu.cu
	$(CUDA_HOME)/bin/nvcc -g -lineinfo $(CODE_HOME)/CUDA/justgpu/include/justgpu.cu $(CODE_HOME)/CUDA/justgpu/main.cu -o $(CODE_HOME)/exec/justgpu -I$(CUDA_HOME)/include,$(CODE_HOME)/CUDA/justgpu/include,$(EX_HOME)/.. -L$(CUDA_HOME)/lib64,pthread -std=c++11

justgpu_test: $(CODE_HOME)/CUDA/justgpu/test.cu $(CODE_HOME)/CUDA/justgpu/include/justgpu.cu
	$(CUDA_HOME)/bin/nvcc -g -lineinfo $(CODE_HOME)/CUDA/justgpu/include/justgpu.cu $(CODE_HOME)/CUDA/justgpu/test.cu -o $(CODE_HOME)/exec/justgpu_test -I$(CUDA_HOME)/include,$(CODE_HOME)/CUDA/justgpu/include,$(EX_HOME)/.. -L$(CUDA_HOME)/lib64,pthread -std=c++11
