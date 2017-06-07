kernel void hello(global double* in, global double* out)
{
	size_t i = get_global_id(0);
	out[i] = in[i] + 1.235;
}