__kernel
void fillzero(global float* a)
{
    int width = get_global_size(0);
    int x = get_global_id(0);
    int y = get_global_id(1);

    int idx = x + y * width;
    
    a[idx] = 0.0f;
}
