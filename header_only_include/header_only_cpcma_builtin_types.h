#ifndef __cplusplus
#ifndef Included_header_only_cpcma_builtin_types_h
#define Included_header_only_cpcma_builtin_types_h
#include<math.h>
#include<cpcma_builtin_types.h>
/**
 * Uses trial division with 6n+-1 optimization
 * Returns zero if prime and non-zero if not prime
 */
int cpcma_check_prime(cpcma____uint64 x)
{
	// get rid of obvious cases
	// cases that break the 3x speed algorithm
	if(x == 2 || x == 3)
		return 0;
	else if(x == 1 || x % 2 == 0 || x % 3 == 0)
		return 1;
	else
	{
		// the factor
		int fac = 0;
		cpcma____uint64 y = sqrt(x);
		// loop up to the square root
		// i+=2 to keep divisor odd
		for(cpcma____uint64 i = 5; i <= y ; i+=2)
		{
			// trial division
			if(x % i == 0)
			{
				fac = i;
				i = y;
			}
			// actually, we do not want odd numbers that are multiples of three
			// to skip 3 mod 6 we add 2 to 1 mod 6 twice
			if(i % 6 == 1)
				i+=2;
		}
		return fac;
	}
}
#endif
#endif
