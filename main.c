// Evolved from https://github.com/arkanis/minimal-opencl-on-windows.git
// 2017-06-01 Paul Kienzle modified for testing double precision kernels

#include <stdio.h>
#include <string.h>
#include <CL/cl.h>


int main() {
	// Find the first GPU device
	cl_device_id device = 0;
	
	cl_uint platform_count = 0;
	clGetPlatformIDs(0, NULL, &platform_count);
	cl_platform_id platform_ids[platform_count];
	clGetPlatformIDs(platform_count, platform_ids, &platform_count);
	
	size_t i;
	for(i = 0; i < platform_count; i++) {
		cl_platform_id platform_id = platform_ids[i];
		
		cl_uint device_count = 0;
		clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, 0, NULL, &device_count);
		cl_device_id device_ids[device_count];
		clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, device_count, device_ids, &device_count);
		
		if (device_count > 0) {
			device = device_ids[0];
			break;
		}
	}
	
	if (!device) {
		fprintf(stderr, "Failed to find any OpenCL GPU device. Sorry.\n");
		return 1;
	}
	
	
	// Setup OpenCL
	cl_context context = clCreateContext(NULL, 1, &device, NULL, NULL, NULL);
	cl_command_queue command_queue = clCreateCommandQueue(context, device, 0, NULL);
	
	// Compile the kernel stored in hello.c
	const char* program_code = "#include <kernel.c>\n";
	cl_program program = clCreateProgramWithSource(context, 1, (const char*[]){program_code}, NULL, NULL);
	cl_int error = clBuildProgram(program, 0, NULL, "-I.", NULL, NULL);
	if (error) {
  	    char compiler_log[4096];
	    clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(compiler_log), compiler_log, NULL);
	    printf("Program build info:\n%s", compiler_log);
		printf("OpenCL compiler failed: %d\n", error);
		return 2;
	}
	cl_kernel kernel = clCreateKernel(program, "hello", NULL);
	
	// Setup GPU buffers
	double transformed[8] = {0.,1.,2.,3.,4.,5.,6.,7.};
	size_t transformed_length = 8;
	cl_mem buffer_in = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, transformed_length*sizeof(*transformed), transformed, NULL);
	cl_mem buffer_out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, transformed_length*sizeof(*transformed), NULL, NULL);
	
	// Execute kernel
	clSetKernelArg(kernel, 0, sizeof(buffer_in), &buffer_in);
	clSetKernelArg(kernel, 1, sizeof(buffer_out), &buffer_out);
	size_t global_work_size[1];
	global_work_size[0] = transformed_length;
	error = clEnqueueNDRangeKernel(command_queue, kernel, (cl_uint)1, NULL, global_work_size, NULL, (cl_uint)0, NULL, NULL);
	if (error) {
		printf("Enqueue kernel failed: %d\n", error);
		return 2;
	}
	
	// Output result
	double* result = clEnqueueMapBuffer(command_queue, buffer_out, CL_TRUE, CL_MAP_READ, 0, transformed_length*sizeof(*transformed), 0, NULL, NULL, NULL);
	printf("in: "); for (int k=0; k<transformed_length; k++) printf(" %8.5f", transformed[k]); printf("\n");
	printf("out:"); for (int k=0; k<transformed_length; k++) printf(" %8.5f", result[k]); printf("\n");
	clEnqueueUnmapMemObject(command_queue, buffer_out, result, 0, NULL, NULL);
	
	// Clean up
	clReleaseMemObject(buffer_in);
	clReleaseMemObject(buffer_out);
	clReleaseProgram(program);
	clReleaseCommandQueue(command_queue);
	clReleaseContext(context);
	
	return 0;
}