__kernel void local_max(__global float* in_mem, __global float* out_mem, float mean_sigma_std)
{

    int id = get_global_id(0);
    
    if(in_mem[id] > mean_sigma_std)
    {
        out_mem[id] = 1.0f;
    }
    else
    {
        out_mem[id] = 0.0f;
    }



}

