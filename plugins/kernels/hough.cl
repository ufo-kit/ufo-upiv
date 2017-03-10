const sampler_t smp = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;
#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable

kernel void likelihood (read_only image3d_t input, 
                        global float *output,
                        constant int *mask,
                        private int maskSizeH, 
                        private int maskNumOnes,
                        local float *local_mem,
                        global int *count,
                        private float threshold)
{
    int shift = 6;

    unsigned glb_size_x = get_global_size(0);
    unsigned glb_size_y = get_global_size(1);

    unsigned loc_size_x = get_local_size(0);
    unsigned loc_size_y = get_local_size(1);
    unsigned loc_mem_size_x = loc_size_x + 2 * shift;
    
    unsigned glb_x = get_global_id(0);
    unsigned glb_y = get_global_id(1);
    unsigned glb_z = get_global_id(2);
    const int4 glb_pos = (int4) (glb_x, glb_y, glb_z, 0);

    unsigned loc_x = get_local_id(0);
    unsigned loc_y = get_local_id(1);
    unsigned local_tmp_id;

    local_tmp_id = (shift + loc_y) * loc_mem_size_x + (shift + loc_x);
    local_mem[local_tmp_id] = read_imagef(input, smp, glb_pos).x; 

    if (loc_x < shift) {
        local_tmp_id = (shift + loc_y) * loc_mem_size_x + loc_x;
        local_mem[local_tmp_id] = read_imagef(input, smp, (int4)(glb_x - shift, glb_y, glb_z, 0)).x;
    }

    if (loc_y < shift) {
        local_tmp_id = loc_y * loc_mem_size_x + (shift + loc_x);
        local_mem[local_tmp_id] = read_imagef(input, smp, (int4)(glb_x, glb_y - shift, glb_z, 0)).x;
    }

    if (loc_x + shift >= loc_size_x) {
        local_tmp_id = (loc_y + shift) * loc_mem_size_x + (loc_x + 2*shift);
        local_mem[local_tmp_id] = read_imagef(input, smp, (int4)(glb_x + shift, glb_y, glb_z, 0)).x;
    }

    if (loc_y + shift >= loc_size_y) {
        local_tmp_id = (loc_y + 2*shift) * loc_mem_size_x + (loc_x + shift);
        local_mem[local_tmp_id] = read_imagef(input, smp, (int4)(glb_x, glb_y + shift, glb_z, 0)).x;
    }

    if (loc_x < shift && loc_y < shift) {
        local_tmp_id = loc_y * loc_mem_size_x + loc_x;
        local_mem[local_tmp_id] = read_imagef(input, smp, (int4)(glb_x - shift, glb_y-shift, glb_z, 0)).x;
    }

    if (loc_x + shift >= loc_size_x && loc_y < shift) {
        local_tmp_id = loc_y * loc_mem_size_x + (loc_x + 2*shift);
        local_mem[local_tmp_id] = read_imagef(input, smp, (int4)(glb_x + shift, glb_y - shift, glb_z, 0)).x;
    }

    if (loc_x < shift && loc_y + shift >= loc_size_y) {
        local_tmp_id = (loc_y + 2*shift) * loc_mem_size_x + loc_x;
        local_mem[local_tmp_id] = read_imagef(input, smp, (int4)(glb_x - shift, glb_y + shift, glb_z, 0)).x;
    }

    if (loc_x + shift >= loc_size_x && loc_y + shift >= loc_size_y) {
        local_tmp_id = (loc_y + 2*shift) * loc_mem_size_x + loc_x + 2*shift;
        local_mem[local_tmp_id] = read_imagef(input, smp, (int4)(glb_x + shift, glb_y + shift, glb_z, 0)).x;
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    float res, mean = 0.0f, std = 0.0f, peak = 0.0f, f;

    for (int a = -maskSizeH; a < maskSizeH + 1; a++) {
        for (int b = -maskSizeH; b < maskSizeH + 1; b++) {
            local_tmp_id = (loc_y + shift + b) * loc_mem_size_x + (loc_x + shift +a);
            mean += local_mem[local_tmp_id] * mask[a+maskSizeH+(b+maskSizeH)*(maskSizeH*2+1)];
        }
    }
    mean = mean / maskNumOnes;

    for (int a = -maskSizeH; a < maskSizeH + 1; a++) {
        for (int b = -maskSizeH; b < maskSizeH + 1; b++) {
            local_tmp_id = (loc_y + shift + b) * loc_mem_size_x + (loc_x + shift +a);
            f = local_mem[local_tmp_id];
            std += (f - mean) * (f -mean) * mask[a+maskSizeH+(b+maskSizeH)*(maskSizeH*2+1)];
        }
    }
    std = sqrt(std / maskNumOnes);

#if 1
    for (int a = -1; a <= 1; a++) for (int b = -1; b <= 1; b++) {
        local_tmp_id = (loc_y + shift + b) * loc_mem_size_x + (loc_x + shift +a);
        peak += local_mem[local_tmp_id];
    }
    peak /= 9.0f;
#else
    local_tmp_id = (loc_y + shift) * loc_mem_size_x + (loc_x + shift);
    peak = local_mem[local_tmp_id];
#endif

    res = exp((peak - mean)/std);
/*
 *    unsigned idx = glb_x 
 *                   + glb_y * get_global_size(0) 
 *                   + glb_z * get_global_size(0) * get_global_size(1);
 *
 *    output[idx] = res;
 */
    if (res > threshold) {
        int old = atomic_inc(count);
        output[5*old + 2] = (float)glb_x; //save x coordinate
        output[5*old + 3] = (float)glb_y; //save y coordinate
        output[5*old + 4] = (float)glb_z;
        output[5*old + 5] = res;
        output[5*old + 6] = (float)1;
    }
}

kernel void likelihood_image (read_only image3d_t input, 
                        global float *output,
                        constant int *mask,
                        private int maskSizeH, 
                        private int maskNumOnes,
                        local float *local_mem)
{
    int shift = 6;

    unsigned glb_size_x = get_global_size(0);
    unsigned glb_size_y = get_global_size(1);

    unsigned loc_size_x = get_local_size(0);
    unsigned loc_size_y = get_local_size(1);
    unsigned loc_mem_size_x = loc_size_x + 2 * shift;
    
    unsigned glb_x = get_global_id(0);
    unsigned glb_y = get_global_id(1);
    unsigned glb_z = get_global_id(2);
    const int4 glb_pos = (int4) (glb_x, glb_y, glb_z, 0);

    unsigned loc_x = get_local_id(0);
    unsigned loc_y = get_local_id(1);
    unsigned local_tmp_id;

    local_tmp_id = (shift + loc_y) * loc_mem_size_x + (shift + loc_x);
    local_mem[local_tmp_id] = read_imagef(input, smp, glb_pos).x; 

    if (loc_x < shift) {
        local_tmp_id = (shift + loc_y) * loc_mem_size_x + loc_x;
        local_mem[local_tmp_id] = read_imagef(input, smp, (int4)(glb_x - shift, glb_y, glb_z, 0)).x;
    }

    if (loc_y < shift) {
        local_tmp_id = loc_y * loc_mem_size_x + (shift + loc_x);
        local_mem[local_tmp_id] = read_imagef(input, smp, (int4)(glb_x, glb_y - shift, glb_z, 0)).x;
    }

    if (loc_x + shift >= loc_size_x) {
        local_tmp_id = (loc_y + shift) * loc_mem_size_x + (loc_x + 2*shift);
        local_mem[local_tmp_id] = read_imagef(input, smp, (int4)(glb_x + shift, glb_y, glb_z, 0)).x;
    }

    if (loc_y + shift >= loc_size_y) {
        local_tmp_id = (loc_y + 2*shift) * loc_mem_size_x + (loc_x + shift);
        local_mem[local_tmp_id] = read_imagef(input, smp, (int4)(glb_x, glb_y + shift, glb_z, 0)).x;
    }

    if (loc_x < shift && loc_y < shift) {
        local_tmp_id = loc_y * loc_mem_size_x + loc_x;
        local_mem[local_tmp_id] = read_imagef(input, smp, (int4)(glb_x - shift, glb_y-shift, glb_z, 0)).x;
    }

    if (loc_x + shift >= loc_size_x && loc_y < shift) {
        local_tmp_id = loc_y * loc_mem_size_x + (loc_x + 2*shift);
        local_mem[local_tmp_id] = read_imagef(input, smp, (int4)(glb_x + shift, glb_y - shift, glb_z, 0)).x;
    }

    if (loc_x < shift && loc_y + shift >= loc_size_y) {
        local_tmp_id = (loc_y + 2*shift) * loc_mem_size_x + loc_x;
        local_mem[local_tmp_id] = read_imagef(input, smp, (int4)(glb_x - shift, glb_y + shift, glb_z, 0)).x;
    }

    if (loc_x + shift >= loc_size_x && loc_y + shift >= loc_size_y) {
        local_tmp_id = (loc_y + 2*shift) * loc_mem_size_x + loc_x + 2*shift;
        local_mem[local_tmp_id] = read_imagef(input, smp, (int4)(glb_x + shift, glb_y + shift, glb_z, 0)).x;
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    float res, mean = 0.0f, std = 0.0f, peak = 0.0f, f;

    for (int a = -maskSizeH; a < maskSizeH + 1; a++) {
        for (int b = -maskSizeH; b < maskSizeH + 1; b++) {
            local_tmp_id = (loc_y + shift + b) * loc_mem_size_x + (loc_x + shift +a);
            mean += local_mem[local_tmp_id] * mask[a+maskSizeH+(b+maskSizeH)*(maskSizeH*2+1)];
        }
    }
    mean = mean / maskNumOnes;

    for (int a = -maskSizeH; a < maskSizeH + 1; a++) {
        for (int b = -maskSizeH; b < maskSizeH + 1; b++) {
            local_tmp_id = (loc_y + shift + b) * loc_mem_size_x + (loc_x + shift +a);
            f = local_mem[local_tmp_id];
            std += (f - mean) * (f -mean) * mask[a+maskSizeH+(b+maskSizeH)*(maskSizeH*2+1)];
        }
    }
    std = sqrt(std / maskNumOnes);

#if 1
    for (int a = -1; a <= 1; a++) for (int b = -1; b <= 1; b++) {
        local_tmp_id = (loc_y + shift + b) * loc_mem_size_x + (loc_x + shift +a);
        peak += local_mem[local_tmp_id];
    }
    peak /= 9.0f;
#else
    local_tmp_id = (loc_y + shift) * loc_mem_size_x + (loc_x + shift);
    peak = local_mem[local_tmp_id];
#endif

    res = exp((peak - mean)/std);
    unsigned idx = glb_x 
                   + glb_y * get_global_size(0) 
                   + glb_z * get_global_size(0) * get_global_size(1);

    output[idx] = res;
}

kernel void zero(global float *mem) {
    unsigned glb_x = get_global_id(0);
    mem[glb_x] = -100.0f;
}

kernel void likelihood_old (read_only image3d_t input, 
                            global float *output, 
                            constant int *mask,
                            private int maskSizeH,
                            private int maskNumOnes)
{
   const int4 pos = (int4) (get_global_id (0), get_global_id (1), get_global_id (2), 0);

   float f;
   float mean = 0.0;
   float std = 0.0;

   for (int a = -maskSizeH; a < maskSizeH + 1; a++) {
       for (int b = -maskSizeH; b < maskSizeH + 1; b++) {
           if (mask[a+maskSizeH+(b+maskSizeH)*(maskSizeH*2+1)] == 1) {
               mean += read_imagef(input, smp, pos + (int4)(a,b,0,0)).x;
           }
       }
   }

   int2 here = (int2) (get_global_id (0), get_global_id (1));
   mean /= maskNumOnes;

   for (int a = -maskSizeH; a < maskSizeH + 1; a++) {
       for (int b = -maskSizeH; b < maskSizeH + 1; b++) {
           if (mask[a+maskSizeH+(b+maskSizeH)*(maskSizeH*2+1)] == 1) {
               f = read_imagef(input, smp, pos + (int4)(a,b,0,0)).x;
               std += (f - mean) * (f - mean);
           }
       }
   }

   std = sqrt(std/maskNumOnes);

   f = read_imagef(input, smp, pos).x;
   output[pos.x + pos.y*get_global_size(0)] = exp((f - mean)/std);
}
