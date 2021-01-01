#ifndef __cplusplus
#ifndef Included_cpcma_builtin_types_h
#define Included_cpcma_builtin_types_h
typedef long long unsigned cpcma____uint64;
typedef long long cpcma____int64;
/**
 * Uses trial division with 6n+-1 optimization
 * Returns zero if prime and non-zero if not prime
 */
int cpcma_check_prime(cpcma____uint64 x);
int cpcma_probably_prime(cpcma____uint64 x);
#endif
#endif
