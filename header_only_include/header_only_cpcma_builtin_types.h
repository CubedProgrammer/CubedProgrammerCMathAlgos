#ifndef __cplusplus
#ifndef Included_header_only_cpcma_builtin_types_h
#define Included_header_only_cpcma_builtin_types_h
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<cpcma_builtin_types.h>
// number of primes under one thousand
#define CPCMA____NPUOT 168

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
 * Finds a number x such that pow(x, x) == num
 * For values between 0 and 1, returns the greater of the two
 */
double cpcma_inv_power_tower(double num)
{
	if(num < exp(-1))
		return NAN;
	else
	{
		double guess = num * log(num) + 2;
		double check = pow(guess, guess) - num;
		while(check >= 0.000000000001)
		{
			guess -= (1 - num / (check + num)) / (1 + log(guess));
			check = pow(guess, guess) - num;
		}
		return guess;
	}
}

/**
 * Calculates the factorial of a number
 * Numbers higher than 20 causes overflow
 */
cpcma____uint64 cpcma_fact(unsigned u)
{
	cpcma____uint64 prod = 1;
	for(unsigned i = 1; i <= u; ++i)
		prod *= i;
	return prod;
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
		// split into nine-digit parts
		msds = cache[i - 1] / 1000000000;
		lsds = cache[i - 1] % 1000000000;
		// 2xy of the binomial expansion
		mid = msds * lsds << 1;
		// square the number
		lsds *= lsds;
		msds *= msds;
		lsds += mid % 1000000000 * 1000000000;
		msds += mid / 1000000000;
		// check overflow
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
			// multiply by 10 until at least mod
			while(digs > 0 && msds < mod)
			{
				msds *= 10;
				msds += lsds / 100000000000000000;
				lsds %= 100000000000000000;
				lsds *= 10;
				--digs;
			}
			// perform modular arithmetic
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
			// break res and cache into binomials and multiply
			// last
			lsds = res % 1000000000 * (cache[i] % 1000000000);
			// first
			msds = res / 1000000000 * (cache[i] / 1000000000);
			// inner + outer
			mid = res % 1000000000 * (cache[i] / 1000000000) + res / 1000000000 * (cache[i] % 1000000000);
			lsds += mid % 1000000000 * 1000000000;
			msds += mid / 1000000000;
			// check for overflow
			if(lsds > 999999999999999999)
			{
				lsds -= 1000000000000000000;
				++msds;
			}
			// perform long division
			digs = 18;
			while(digs > 0)
			{
				// multiply by 10 until at least mod
				while(digs > 0 && msds < mod)
				{
					msds *= 10;
					msds += lsds / 100000000000000000;
					lsds %= 100000000000000000;
					lsds *= 10;
					--digs;
				}
				// perform modular arithmetic
				msds %= mod;
			}
			res = msds;
		}
	}
	return res;
}

/**
 * Cross product
 */
void cpcma_crossf(float x1, float y1, float z1, float x2, float y2, float z2, float *restrict x, float *restrict y, float *restrict z)
{
	*x = y1 * z2 - z1 * y2;
	*y = z1 * x2 - x1 * z2;
	*z = x1 * y2 - y1 * x2;
}

/**
 * Cross product
 */
void cpcma_cross(double x1, double y1, double z1, double x2, double y2, double z2, double *restrict x, double *restrict y, double *restrict z)
{
	*x = y1 * z2 - z1 * y2;
	*y = z1 * x2 - x1 * z2;
	*z = x1 * y2 - y1 * x2;
}

/**
 * Normalizes a vector, returns zero on success
 */
int cpcma_normalizef(float *restrict x, float *restrict y, float *restrict z)
{
	float len = sqrtf(*x * *x + *y * *y + *z * *z);
	if(len == 0)
		return-1;
	else
	{
		// do one division instead of three
		float div = 1.0 / len;
		*x *= div;
		*y *= div;
		*z *= div;
		return 0;
	}
}

/**
 * Normalizes a vector, returns zero on success
 */
int cpcma_normalize(double *restrict x, double *restrict y, double *restrict z)
{
	double len = sqrt(*x * *x + *y * *y + *z * *z);
	if(len == 0)
		return-1;
	else
	{
		// do one division instead of three
		double div = 1.0 / len;
		*x *= div;
		*y *= div;
		*z *= div;
		return 0;
	}
}

/**
 * Rotate the vector (x, y) by an angle counterclockwise, in radians
 */
void cpcma_rotate_vectorf(float *restrict x, float *restrict y, float angle)
{
	float re = cosf(angle), im = sinf(angle);
	float nx = re * *x - im * *y, ny = re * *y + im * *x;
	*x = nx;
	*y = ny;
}

/**
 * Rotate the vector (x, y) by an angle counterclockwise, in radians
 */
void cpcma_rotate_vector(double *restrict x, double *restrict y, double angle)
{
	double re = cos(angle), im = sin(angle);
	double nx = re * *x - im * *y, ny = re * *y + im * *x;
	*x = nx;
	*y = ny;
}

/**
 * Multiplies two matrices
 * Returns zero is successful, nonzero if matrices cannot be multiplied
 */
int cpcma_matmulf(float *restrict dest, const float *restrict x, size_t w1, size_t h1, const float *restrict y, size_t w2, size_t h2)
{
	if(w1 == h2)
	{
		size_t rs = h1, cs = w2;
		for(size_t i = 0; i < rs; ++i)
		{
			for(size_t j = 0; j < cs; ++j)
			{
				dest[i * cs + j] = 0;
				for(size_t k = 0; k < h2; ++k)
					dest[i * cs + j] += x[i * w1 + k] * y[k * w2 + j];
			}
		}
		return 0;
	}
	else
		return 13;
}

int cpcma_matmul(double *restrict dest, const double *restrict x, size_t w1, size_t h1, const double *restrict y, size_t w2, size_t h2)
{
	if(w1 == h2)
	{
		size_t rs = h1, cs = w2;
		for(size_t i = 0; i < rs; ++i)
		{
			for(size_t j = 0; j < cs; ++j)
			{
				dest[i * cs + j] = 0;
				for(size_t k = 0; k < h2; ++k)
					dest[i * cs + j] += x[i * w1 + k] * y[k * w2 + j];
			}
		}
		return 0;
	}
	else
		return 13;
}

/**
 * Quadratic formula, number of real roots is returned
 */
int cpcma_quadratic_formulaf(float a, float b, float c, float *restrict r1, float *restrict r2)
{
	float crit = b * b - 4 * a * c;
	int roots;
	if(crit < 0)
		roots = 0;
	else if(crit == 0)
	{
		roots = 1;
		*r1 = *r2 = -b / (2 * a);
	}
	else
	{
		roots = 2;
		crit = sqrt(crit);
		float div = 0.5 / a;
		*r1 = (-b + crit) * div;
		*r2 = (-b - crit) * div;
	}
	return roots;
}

int cpcma_quadratic_formula(double a, double b, double c, double *restrict r1, double *restrict r2)
{
	double crit = b * b - 4 * a * c;
	int roots;
	if(crit < 0)
		roots = 0;
	else if(crit == 0)
	{
		roots = 1;
		*r1 = *r2 = -b / (2 * a);
	}
	else
	{
		roots = 2;
		crit = sqrt(crit);
		double div = 0.5 / a;
		*r1 = (-b + crit) * div;
		*r2 = (-b - crit) * div;
	}
	return roots;
}

/**
 * Evaluates a polynomial, coefficients assumed to be low-order to high-order
 */
float cpcma_eval_polyf(const float *coefs, int sz, float arg)
{
	float res = 0;
	for(int i = sz - 1; i >= 0; --i)
	{
		res *= arg;
		res += coefs[i];
	}
	return res;
}

double cpcma_eval_poly(const double *coefs, int sz, double arg)
{
	double res = 0;
	for(int i = sz - 1; i >= 0; --i)
	{
		res *= arg;
		res += coefs[i];
	}
	return res;
}

/**
 * Probabilistic prime checking algorithm
 * Returns zero if x is strong-pseudoprime or prime
 * Works for numbers less than 10 to the 18
 */
int cpcma_strong_pseudoprime_base(cpcma____uint64 x, cpcma____uint64 base)
{
	cpcma____uint64 exp = x - 1;
	while(exp % 2 == 0)
		exp >>= 1;
	cpcma____uint64 m = cpcma_mod_pow64u(base, exp, x);
	while(m != 1 && exp < x - 1)
	{
		m = cpcma_mod_pow64u(m, 2, x);
		exp <<= 1;
	}
	if(m == 1)
		return 0;
	else
		return 1;
}

/**
 * Faster algorithm for checking if a number is prime
 * Custom base for Fermat's little theorem
 * Works for numbers less than 10 to the 18
 */
int cpcma_probably_prime_base(cpcma____uint64 x, cpcma____uint64 base)
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
			if(cpcma_mod_pow64u(base, x - 1, x) != 1)
				fact = 1;
		}
		return fact;
	}
}

/**
 * Gets the xth fibonacci number
 */
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
