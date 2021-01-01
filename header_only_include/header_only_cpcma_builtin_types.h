#ifndef __cplusplus
#ifndef Included_header_only_cpcma_builtin_types_h
#define Included_header_only_cpcma_builtin_types_h
#include<math.h>
#include<cpcma_builtin_types.h>
// number of primes under one thousand
#define CPCMA____NPUOT 168

// prime checking base
#define CPCMA____PCB 3

// maximum exponent of prime number
#define CPCMA____MEP 70

// primes up to one thousand
static int cpcma____putot[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
		41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
		109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
		181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251,
		257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331,
		337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
		419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487,
		491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577,
		587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653,
		659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743,
		751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829,
		839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929,
		937, 941, 947, 953, 967, 971, 977, 983, 991, 997};

/**
 * Faster algorithm for checking if a number is prime
 * Works for numbers less than 10 to the 19
 */
int cpcma_probably_prime(cpcma____uint64 x)
{
	// get rid of the obvious cases
	if(x == 2 || x == 3)
		return 0;
	else if(x == 1 || x % 2 == 0 || x % 3 == 0)
		return 1;
	else
	{
		// the factor
		int fact = 0;
		// loop through all the primes less than one thousand
		for(int i = 2; i < CPCMA____NPUOT; ++i)
		{
			// trial division
			if(x != cpcma____putot[i] && x % cpcma____putot[i] == 0)
			{
				fact = i;
				i = CPCMA____NPUOT;
			}
		}
		// if still prime and there are more factors to check
		// switch algorithms
		if(fact == 0 && 997 * 997 < x)
		{
			// stores the most significant 19 digits
			// and the least significant 19 digits
			cpcma____uint64 msig19d = 0, lsig19d = CPCMA____PCB;
			// cache for modular exponentiation
			cpcma____uint64 cache[CPCMA____MEP];
			cache[0] = lsig19d;
			// the most significant digit
			cpcma____uint64 msd;
			// compute the base raised to powers of two
			// and cache them
			cpcma____uint64 pow = 1, ind = 0;
			cpcma____uint64 tmp;
			while(pow && pow < x)
			{
				lsig19d = cache[ind] % 1000000000;
				msig19d = cache[ind] / 1000000000 % 1000000000;
				msd = cache[ind] / 1000000000000000000ull % 10;
				tmp = 2 * lsig19d * msig19d;
				lsig19d *= lsig19d;
				msig19d *= msig19d;
				lsig19d += tmp % 1000000000 * 1000000000;
				msig19d += tmp / 1000000000;
				tmp = (cache[ind] % 1000000000) * msd * 2;
				msig19d += tmp;
				tmp = (cache[ind] / 1000000000 % 1000000000) * msd * 2;
				msig19d += tmp % 1000000000 * 1000000000;
				msd *= msd;
				msd += tmp / 1000000000 % 100;
				// double the power, exponent is increased by one
				// using exponent laws
				pow *= 2;
			}
		}
		return fact;
	}
}
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
