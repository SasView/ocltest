// Evolved from https://github.com/arkanis/minimal-opencl-on-windows.git
// 2017-06-01 Paul Kienzle modified for testing double precision kernels

//#define SOURCE "#include <add.c>\n"
#define SOURCE "#include <cylinder.c>\n"

#include <stdio.h>
#include <string.h>
#include <CL/cl.h>

#define MAX_PD 4
typedef struct {
    int32_t pd_par[MAX_PD];
    int32_t pd_length[MAX_PD];
    int32_t pd_offset[MAX_PD];
    int32_t pd_stride[MAX_PD];
    int32_t num_eval;
    int32_t num_weights;
    int32_t num_active;
    int32_t theta_par;
} ProblemDetails;

// Parameters: 
//     scale, background, sld, sld_solvent, radius, length, theta, phi, 
//     up:frac_i, up:frac_f, up:angle, M0:sld, mtheta:sld, mphi:sld,
//     M0:sld_solvent, mtheta:sld_solvent, mphi:sld_solvent
const double call_values[] = {
	1., 0.,                               //  2 (scale, background)
        3., 0., 100., 300., 30., 60.,         //  8 (parameter slots)
        0., 0., 0., 0., 0., 0., 0., 0., 0.,   // 17 (magnetism)
	3., 0., 100., 300., 30., 60.,         // 23 (parameter value vectors)
        0., 0., 0., 0., 0., 0., 0., 0., 0.,   // 32 (magnetism value vectors)
	1., 1., 1., 1., 1., 1.,               // 38 (parameter weight vectors)
        1., 1., 1., 1., 1., 1., 1., 1., 1.,   // 47 (magnetism weight vectors)
        0.,                             // 48 padding
        0., 0., 0., 0., 0., 0., 0., 0., // 56 padding
        0., 0., 0., 0., 0., 0., 0., 0.  // 64 padding
};
const ProblemDetails call_details = {
    {2, 3, 4, 5},
    {1, 1, 1, 1},
    {2, 3, 4, 5},
    {1, 1, 1, 1},
    1, 15, 0, 4
};
const double call_q[] = { // qx, qy pairs
    0.001, 0.001,
    0.002, 0.002, 
    0.005, 0.003, 
    0.01, 0.004,
    0.02, -0.01, 
    0.05, 0.0, 
    0.1, 0.0,
    0.2, 0.0, 
    0.5, 0.0, 
    1., 0.0,
    2., 0.0, 
    5., 0.0,
    0., 0.,  // padding
    0., 0.,  // padding
    0., 0.,  // padding
    0., 0.,  // padding
};


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
	const char* program_code = SOURCE;
	cl_program program = clCreateProgramWithSource(context, 1, (const char*[]){program_code}, NULL, NULL);
        printf("building...\n%s", program_code);
	cl_int error = clBuildProgram(program, 0, NULL, "-I.", NULL, NULL);
	if (error) {
  	    char compiler_log[4096];
	    clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(compiler_log), compiler_log, NULL);
	    printf("Program build info:\n%s", compiler_log);
		printf("OpenCL compiler failed: %d\n", error);
		return 2;
	}
        printf("loading...\n");
	cl_kernel kernel = clCreateKernel(program, "cylinder_Iqxy", NULL);
	
        // place holder for the results
        double call_result[16];

	// Setup GPU buffers
        int32_t nq = 12;
        int32_t pd_start = 0;
        int32_t pd_stop = 1;
	cl_mem details = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(call_details), (void *)&call_details, NULL);
        cl_mem values = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(call_values), (void *)call_values, NULL);
        cl_mem q = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(call_q), (void *)call_q, NULL);
	cl_mem buffer_out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(call_result), NULL, NULL);
        double cutoff = 0.;
	
	// Execute kernel
        printf("running...\n");
	clSetKernelArg(kernel, 0, sizeof(nq), &nq);
	clSetKernelArg(kernel, 1, sizeof(pd_start), &pd_start);
	clSetKernelArg(kernel, 2, sizeof(pd_stop), &pd_stop);
	clSetKernelArg(kernel, 3, sizeof(details), &details);
	clSetKernelArg(kernel, 4, sizeof(values), &values);
	clSetKernelArg(kernel, 5, sizeof(q), &q);
	clSetKernelArg(kernel, 6, sizeof(buffer_out), &buffer_out);
	clSetKernelArg(kernel, 7, sizeof(cutoff), &cutoff);
	size_t global_work_size[1];
	global_work_size[0] = 16;
	error = clEnqueueNDRangeKernel(command_queue, kernel, (cl_uint)1, NULL, global_work_size, NULL, (cl_uint)0, NULL, NULL);
	if (error) {
		printf("Enqueue kernel failed: %d\n", error);
		return 2;
	}
	
	// Output result
	double* result = clEnqueueMapBuffer(command_queue, buffer_out, CL_TRUE, CL_MAP_READ, 0, sizeof(call_result), 0, NULL, NULL, NULL);
	printf("in: "); for (int k=0; k<nq; k++) printf(" %8.5f", call_q[2*k]); printf("\n");
	printf("out:"); for (int k=0; k<nq; k++) printf(" %8.5f", result[k]); printf("\n");
	clEnqueueUnmapMemObject(command_queue, buffer_out, result, 0, NULL, NULL);
	
	// Clean up
	clReleaseMemObject(details);
	clReleaseMemObject(values);
	clReleaseMemObject(q);
	clReleaseMemObject(buffer_out);
	clReleaseProgram(program);
	clReleaseCommandQueue(command_queue);
	clReleaseContext(context);
	
	return 0;
}
