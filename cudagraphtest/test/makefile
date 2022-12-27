# Copyright (c) 2022, NVIDIA Corporation
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

CUDA_HOME ?= /usr/local/cuda

all: run.x

run.x: revise5a.cu gpu_graph.hpp gpu_graph.cpp cuda_helper.hpp common_string.h common.h
	$(CUDA_HOME)/bin/nvcc -g -lineinfo -gencode arch=compute_70,code=sm_70 revise5a.cu gpu_graph.cpp -o /home/deepl/sooho/test/K-com/test/exec/run.x -I$(CUDA_HOME)/include -L$(CUDA_HOME)/lib64 -lcudart -lpthread -std=c++11

# all: graph1.x

# graph1.x: revise6.cu gpu_graph.hpp gpu_graph.cpp cuda_helper.hpp
# 	$(CUDA_HOME)/bin/nvcc -g -lineinfo -gencode arch=compute_70,code=sm_70 revise6.cu gpu_graph.cpp -o /home/deepl/sooho/test/K-com/test/exec/graph1.x -I$(CUDA_HOME)/include -I/home/deepl/sooho/test/K-com/test/.. -L$(CUDA_HOME)/lib64 -lcudart -std=c++11
