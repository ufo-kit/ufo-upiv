 __kernel void
hessian_det (__global float *out,
             __global float *in)
{
    const int idx = get_global_id(0);
    const int idy = get_global_id(1);
    
    const int width = get_global_size(0);
    const int slice = get_global_size(0) * get_global_size(1);
    const int pitch = width * idy + idx;

    float hxx = 1000 * in[pitch];
    float hyy = 1000 * in[slice + pitch];
    float hxy = 1000 * in[slice*2 +  pitch];

    out[idy * width + idx] = hxx*hyy - hxy*hxy;
}


__kernel void
hessian_eigval (__global float *out, __global float *in)
{
    float res;

    const int idx = get_global_id(0);
    const int idy = get_global_id(1);
    
    const int width = get_global_size(0);
    const int slice = get_global_size(0) * get_global_size(1);
    const int pitch = width * idy + idx;

    float hxx = in[pitch];
    float hyy = in[slice + pitch];
    float hxy = in[slice*2 +  pitch];

    float tmp1 = 0.5 * (hxx  + hyy);
    float tmp2 = 0.5 * sqrt(4 * hxy * hxy + (hxx - hyy)*(hxx - hyy));

    float lmd1 = tmp1 + tmp2;
    float lmd2 = tmp1 - tmp2;

    if (lmd1 < 0 && lmd2 < 0) {
        res = (lmd2 > lmd1) ? lmd2 / lmd1: lmd1 / lmd2;
    } else {
        res = 0.0f;
    }

    out[idy * width + idx] = res;
}
