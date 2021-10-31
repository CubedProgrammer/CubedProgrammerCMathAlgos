#ifndef __cplusplus
#ifndef Included_cpcma_builtin_types_h
#define Included_cpcma_builtin_types_h
#include<math.h>

// prime checking base
#define CPCMA____PCB 3

typedef long long unsigned cpcma____uint64;
typedef long long cpcma____int64;

/**
 * Finds all primes up to a number
 * Specifically, sieves the array buf, so that only prime indices are non-zerp
 */
void cpcma_sieve_eratosthenes(size_t size,char buf[]);

/**
 * Finds a number x such that pow(x, x) == num
 * For values between 0 and 1, returns the greater of the two
 */
double cpcma_inv_power_tower(double num);

/**
 * num raised to the power of num
 * Zero to the zero is one in this case
 */
static inline double cpcma_power_tower(double num)
{
	return pow(num, num);
}

/**
 * Uses trial division with 6n+-1 optimization
 * Returns zero if prime and non-zero if not prime
 */
int cpcma_check_prime(cpcma____uint64 x);

/**
 * Correct modular arithmetic for ints
 */
int cpcma_correct_mod(int x, int y);

/**
 * Correct modular arithmetic for longs
 */
long cpcma_correct_modl(long x, long y);

/**
 * Correct modular arithmetic for ints
 */
cpcma____int64 cpcma_correct_modll(cpcma____int64 x, cpcma____int64 y);

/**
 * Modular exponentiation
 */
int cpcma_mod_pow(int base, int exp, int mod);

/**
 * Modular exponentiation with long exponent
 */
int cpcma_mod_pow_long_exp(int base, long exp, int mod);

/**
 * Modular exponentiation with cpcma____int64 exponent
 */
int cpcma_mod_pow_llong_exp(int base, cpcma____int64 exp, int mod);

/**
 * 64-bit modular exponentiation
 */
cpcma____int64 cpcma_mod_pow64(cpcma____int64 base, cpcma____int64 exp, cpcma____int64 mod);

/**
 * 64-bit unsigned modular exponentiation
 */
cpcma____uint64 cpcma_mod_pow64u(cpcma____uint64 base, cpcma____uint64 exp, cpcma____uint64 mod);

/**
 * Least common multiple for ints
 */
int cpcma_lcm_int32(int x, int y);

/**
 * Least common multiple for 64-bit ints
 */
cpcma____int64 cpcma_lcm_int64(cpcma____int64 x, cpcma____int64 y);

/**
 * Greatest common divisor for ints
 */
int cpcma_gcd_int32(int x, int y);

/**
 * Greatest common divisor for 64-bit ints
 */
cpcma____int64 cpcma_gcd_int64(cpcma____int64 x, cpcma____int64 y);

/**
 * Factors a number and stores them in an array pointed to by factorp, numfac points to the size
 */
void cpcma_factor_uint64(cpcma____uint64 x, cpcma____uint64 *factorp[], size_t *numfac);

/**
 * Gets the xth fibonacci number
 */
cpcma____uint64 cpcma_get_fib(int x);

/**
 * Faster algorithm for checking if a number is prime
 * Custom base for Fermat's little theorem
 * Works for numbers less than 10 to the 18
 */
int cpcma_probably_prime_base(cpcma____uint64 x, cpcma____uint64 base);

/**
 * Faster algorithm for checking if a number is prime
 * Works for numbers less than 10 to the 18
 */
static inline int cpcma_probably_prime(cpcma____uint64 x)
{
	return cpcma_probably_prime_base(x, CPCMA____PCB);
}

#endif
#endif
