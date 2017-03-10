__kernel void 
sigmoid (__global float* res, __global float* data, float mean, float sigma, 
                  float max, float gamma, float c1, float c2)
{
    float f0, f1, f2; 

    int id = get_global_id(0);
    float low = mean + c1*sigma;
    float high = mean + c2*sigma;
    high = (high > max) ? max : high;

    // data < low, ==> 0
    // data > high, ==> 1
    // otherwise, pow(d, gamma), where d in (0,1)

    if (data[id] > high) {
        f0 = 1.0;
        f1 = 1.0;
    } else if (data[id] < low) {
        f0 = 1.0;
        f1 = 0.0f;
    } else {
        f0 = (data[id] - low) / (high - low);
        f1 = 1.0;
    }

    res[id] = f1 * pow(f0, gamma);
}
//    ;float c1 = 0.25; // lower threshold
//    ;float c2 = 10;  // higher threshold
//    ;float c3 = 4.0; // shift of sigmoid (in unit of sigma)
//    ;float c4 = 2.0; // width of sigmoid (in unit of sigma)
//        f0 = data[id];
//        f1 = exp( (mean + c3 * sigma - data[id]) / (c4 *sigma) );
//        f1 = 1.0f / (1.0f + f1);
