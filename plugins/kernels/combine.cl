__kernel void combine(__global float* in_mem, __global float* out_mem,unsigned  int counter, unsigned int mem_size_p)
{
	int id = get_global_id(0);
	float tmp1 = 0.0f; //Due to this initialization the else part can be vanished.
	float tmp2 = 0.0f;
	int i = 1;

	//store variables private, thus the global memoy won't be used 
	unsigned int loc_mem_p = (unsigned int)mem_size_p;
	unsigned int loc_counter = (unsigned int) counter;

	//first and operation
	if(in_mem[id] > 0 && in_mem[id] > 0)
	{
		tmp1 = 1.0f;
	}
	
	out_mem[id] = tmp1; //first picture does not have a left partner
	for(;i<loc_counter;i++)
	{
		if(in_mem[id + i*loc_mem_p] > 0 && in_mem[id + (1+i)*loc_mem_p])
		{
			tmp2 = 1.0f;
		}

		if(tmp1 > 0 || tmp2 > 0)
		{
			out_mem[id + i*loc_mem_p] = 1.0f;
		}
		else
		{
			out_mem[id + i*loc_mem_p] = 0.0f;
		}
		
		tmp1 = tmp2;
		tmp2 = 0.0f;
				



	}



}
