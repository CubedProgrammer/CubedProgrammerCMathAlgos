#ifndef __cplusplus
#ifndef Included_cpcma_builtin_types_h
#define Included_cpcma_builtin_types_h
typedef long long unsigned cpcma____uint64;
typedef long long cpcma____int64;

/**
 * Finds all primes up to a number
 * Specifically, sieves the array buf, so that only prime indices are non-zerp
 */
void cpcma_sieve_eratosthenes(size_t size,char buf[]);

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
cpcma____uint64 cpcma_get_fib(int x);
int cpcma_probably_prime(cpcma____uint64 x);
#endif
#endif
