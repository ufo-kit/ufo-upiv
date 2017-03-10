__kernel
void mult(global float* a, global float* b, global float* res)
{
    int width = get_global_size(0) * 2;
    int x = get_global_id(0) * 2;
    int y = get_global_id(1);
    int idx_r = x + y * width;
    int idx_i = 1 + x + y * width;

    float ra = a[idx_r];
    float rb = b[idx_r];
    float ia = a[idx_i];
    float ib = b[idx_i];

    res[idx_r] = ra * rb - ia * ib;
    res[idx_i] = ra * ib + rb * ia;
}

__kernel
void mult3d(global float* a, global float* b, global float* res)
{
    const int width = get_global_size(0);
    const int height = get_global_size(1);

    const int stride_x = width * 2;
    const int stride_y = stride_x * height;

    const int x = get_global_id(0);
    const int y = get_global_id(1);
    const int z = get_global_id(2);

    int idx_r = 2*x + y*stride_x;
    int idx_i = 1 + idx_r;

    float ra = a[idx_r];
    float ia = a[idx_i];
    float rb = b[idx_r + z*stride_y];
    float ib = b[idx_i + z*stride_y];

    res[idx_r + z*stride_y] = ra * rb - ia * ib;
    res[idx_i + z*stride_y] = ra * ib + rb * ia;
}
