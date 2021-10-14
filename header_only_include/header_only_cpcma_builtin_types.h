#ifndef __cplusplus
#ifndef Included_header_only_cpcma_builtin_types_h
#define Included_header_only_cpcma_builtin_types_h
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
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
 * Finds all primes up to a number, note size is one less than the actual size of buf
 * Specifically, sieves the array buf, so that only prime indices are non-zero
 */
void cpcma_sieve_eratosthenes(size_t size,char buf[])
{
	// set everything to be prime at first
	memset(buf, 1, size);
	for(size_t i = 2; i <= size; ++i)
	{
		if(buf[i])
		{
			// set multiples to be composite
			for(size_t j = i * 2; j <= size; j += i)
				buf[j] = 0;
		}
	}
	// explicitly set 0 and 1 to be not prime
	buf[0] = 0;
	buf[1] = 0;
}

/**
 * Correct modular arithmetic for ints
 */
int cpcma_correct_mod(int x, int y)
{
	int n = x % y + y;
	return n % y;
}

/**
 * Correct modular arithmetic for longs
 */
long cpcma_correct_modl(long x, long y)
{
	long n = x % y + y;
	return n % y;
}

/**
 * Correct modular arithmetic for ints
 */
long long cpcma_correct_modll(long long x, long long y)
{
	long long n = x % y + y;
	return n % y;
}

/**
 * Modular exponentiation
 */
int cpcma_mod_pow(int base, int exp, int mod)
{
	long long x = 1;
	long long cache[35];
	cache[0] = base;
	for(int i = 1; i < 35; ++i)
		cache[i] = cpcma_correct_modll(cache[i - 1] * cache[i - 1], mod);
	for(int i = 0; i < 31; ++i)
	{
		if((exp >> i) % 2)
			x = cpcma_correct_modll(x * cache[i], mod);
	}
	return x;
}

/**
 * Modular exponentiation with long exponent
 */
int cpcma_mod_pow_long_exp(int base, long exp, int mod)
{
	return cpcma_mod_pow_llong_exp(base, exp, mod);
}

/**
 * Greatest common divisor for ints
 */
int cpcma_gcd_int32(int x, int y)
{
	int r = x;
	while(y > 0)
	{
		r = x % y;
		x = y;
		y = r;
	}
	return x;
}

/**
 * Greatest common divisor for 64-bit ints
 */
cpcma____int64 cpcma_gcd_int64(cpcma____int64 x, cpcma____int64 y)
{
	cpcma____int64 r = x;
	while(y > 0)
	{
		r = x % y;
		x = y;
		y = r;
	}
	return x;
}

/**
 * Least common multiple for ints
 */
int cpcma_lcm_int32(int x, int y)
{
	return x * y / cpcma_gcd_int32(x, y);
}

/**
 * Least common multiple for 64-bit ints
 */
cpcma____int64 cpcma_lcm_int64(cpcma____int64 x, cpcma____int64 y)
{
	return x * y / cpcma_gcd_int64(x, y);
}

/**
 * Factors a number and stores them in an array pointed to by factorp, numfac points to the size
 */
void cpcma_factor_uint64(cpcma____uint64 x, cpcma____uint64 *factorp[], size_t *numfac)
{
	// figure out number of factors
	cpcma____uint64 y = sqrt(x);
	*numfac = 0;
	for(cpcma____uint64 i = 1; i <= y; ++i)
	{
		if(x % i == 0)
			*numfac += i * i == x ? 1 : 2;
	}

	// write all the factors
	*factorp = malloc(*numfac * sizeof(cpcma____uint64));
	size_t size = 0;
	for(cpcma____uint64 i = 1; i <= y; ++i)
	{
		if(x % i == 0)
		{
			if(i * i == x)
			{
				(*factorp)[size] = i;
				size++;
			}
			else
			{
				(*factorp)[size] = i;
				(*factorp)[size + 1] = x / i;
				size += 2;
			}
		}
	}
}

/**
 * Modular exponentiation with long long exponent
 */
int cpcma_mod_pow_llong_exp(int base, long long exp, int mod)
{
	long long x = 1;
	long long cache[70];
	cache[0] = base;
	for(int i = 1; i < 70; ++i)
		cache[i] = cpcma_correct_modll(cache[i - 1] * cache[i - 1], mod);
	for(int i = 0; i < 63; ++i)
	{
		if((exp >> i) % 2)
			x = cpcma_correct_modll(x * cache[i], mod);
	}
	return x;
}

/**
 * 64-bit modular exponentiation
 */
cpcma____int64 cpcma_mod_pow64(cpcma____int64 base, cpcma____int64 exp, cpcma____int64 mod)
{
	base = (base % mod + mod) % mod;
	return cpcma_mod_pow64u(base, exp, mod);
}

/**
 * 64-bit unsigned modular exponentiation
 * For numbers up to 18 digits
 */
cpcma____uint64 cpcma_mod_pow64u(cpcma____uint64 base, cpcma____uint64 exp, cpcma____uint64 mod)
{
	// caches for mod pows
	cpcma____uint64 cache[65];
	base %= mod;
	cache[0] = base;
	// most significant and least significant eighteen digits, and the 2xy in the binomial expansion
	cpcma____uint64 msds, lsds, mid;
	// the amount of times the dividend can be divided by 10 during long division
	int digs;
	for(int i = 1; i < sizeof(cpcma____uint64) * 8; ++i)
	{
		msds = cache[i - 1] / 1000000000;
		lsds = cache[i - 1] % 1000000000;
		mid = msds * lsds << 1;
		lsds *= lsds;
		msds *= msds;
		lsds += mid % 1000000000 * 1000000000;
		msds += mid / 1000000000;
		if(lsds > 999999999999999999)
		{
			lsds -= 1000000000000000000;
			++msds;
		}
		//printf("%llu %llu\n", msds, lsds);
		// perform long division
		digs = 18;
		while(digs > 0)
		{
			while(digs > 0 && msds < mod)
			{
				msds *= 10;
				msds += lsds / 100000000000000000;
				lsds %= 100000000000000000;
				lsds *= 10;
				--digs;
			}
			msds %= mod;
		}
		cache[i] = msds;
		//printf("%llu\n", msds);
	}
	cpcma____uint64 res = 1;
	for(int i = 0; i < sizeof(cpcma____uint64) * 8; ++i)
	{
		if((exp >> i & 1) == 1)
		{
			lsds = res % 1000000000 * (cache[i] % 1000000000);
			msds = res / 1000000000 * (cache[i] / 1000000000);
			mid = res % 1000000000 * (cache[i] / 1000000000) + res / 1000000000 * (cache[i] % 1000000000);
			lsds += mid % 1000000000 * 1000000000;
			msds += mid / 1000000000;
			if(lsds > 999999999999999999)
			{
				lsds -= 1000000000000000000;
				++msds;
			}
			// perform long division
			digs = 18;
			while(digs > 0)
			{
				while(digs > 0 && msds < mod)
				{
					msds *= 10;
					msds += lsds / 100000000000000000;
					lsds %= 100000000000000000;
					lsds *= 10;
					--digs;
				}
				msds %= mod;
			}
			res = msds;
		}
	}
	return res;
}

/**
 * Faster algorithm for checking if a number is prime
 * Works for numbers less than 10 to the 18
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
			if(cpcma_mod_pow64u(CPCMA____PCB, x - 1, x) != 1)
				fact = 1;
		}
		return fact;
	}
}

cpcma____uint64 cpcma_get_fib(int x)
{
	if(x == 1 || x == 2)
		return 1;
	else
	{
		cpcma____uint64 m = 1, n = 1, o = 2;
		while(x > 3)
		{
			m = n;
			n = o;
			o = m + n;
			--x;
		}
		return o;
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
