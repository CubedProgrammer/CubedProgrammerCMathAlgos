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
long long cpcma_correct_modll(long long x, long long y);

/**
 * Modular exponentiation
 */
int cpcma_mod_pow(int base, int exp, int mod);
int cpcma_probably_prime(cpcma____uint64 x);
#endif
#endif
