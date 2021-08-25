INCLUDE=-I/home/vicini/NVIDIA_GPU_Computing_SDK/CUDALibraries/common/inc/ -I/home/vicini/NVIDIA_GPU_Computing_SDK/shared/inc/ -I/home/vicini/mcf1819/common/ -I/usr/local/cuda/include/
OBJS=
NVCC=nvcc  -Wno-deprecated-gpu-targets

ECHO=/bin/echo



default: GaussJordan_GPU_manyBlocks.x
#default: dot2gpu.x
#default: test_time.x
 
%.x: %.cu
	$(NVCC) $(INCLUDE) $(OBJS) $(LIBS) $< -o $@
	./GaussJordan_GPU_manyBlocks.x
	#./test_time.x

checks: GaussJordan_checks.C
	g++ GaussJordan_checks.C -o gaussJordan_checks.x
	./gaussJordan_checks.x

clean:
	rm -f *.x *.o