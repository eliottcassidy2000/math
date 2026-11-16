# Exact extensions of 1024 = 1000 + 24

* [Exact radix classifications, cubic norm orbits and coordinate returns, elliptic and Pell deformations](EXACT.md).
* [Fermat-divisor mechanism, decimal and ternary clocks, MRB interpolation, and the totient fiber](FERMAT_AND_PRECISION.md).
* [Previous exponential classification and the tournament-controller context](../chanpak_collatz_mrb_fibers_20261116/README.md).

The point is to distinguish what is being preserved. The literal affine arrow is rigid; the cubic norm can be preserved indefinitely; the vanishing coordinate that makes it a two-term equation returns only twice; and a single clock can fit arbitrarily many decimal digits without ever returning to its original absolute value.

## An exhaustive coupled sum-product extension

The small pair {3,5} has sum 8 and product plus one 16. For b+c=2^m and bc+1=2^ℓ, directly expanding gives
`(b·2^m+1)(c·2^m+1)=2^(2m+ℓ)+1`.
In particular 1025=(3·8+1)(5·8+1)=25·41. The entire nontrivial family of such b,c is two away from one another. To see this, order 1<b≤c. Now bc+1>b+c, so ℓ>m and b²≡1 mod 2^m. For m≥3 the four roots are ±1 and 2^(m−1)±1: exactly one of b−1 and b+1 is 2 mod 4 and the other must be 0 mod 2^(m−1). In the eligible range 1<b≤2^(m−1), only 2^(m−1)−1 remains. Smaller m has no nontrivial solution. Thus b,c=2^j−1,2^j+1 for j≥2 is exhaustive apart from b=1; it yields the exact factorization
`2^(4j+2)+1=(2^(2j+1)−2^(j+1)+1)(2^(2j+1)+2^(j+1)+1)`.
This is an instance of the Sophie Germain factorization, with 25·41 its first nontrivial occurrence in the coupled sum-product family.

## The positive transfer exists in every dimension

Write α_d=2^(1/d). The inverse to α_d−1 is U_d=1+α_d+⋯+α_d^(d−1). On power-basis coefficients multiplication by U_d has entry 1 on and below the diagonal, 2 above it: terms whose exponent wraps around d acquire the factor 2. For odd d this preserves the norm of 1+α_d, namely 3, since N(α_d−1)=1. Hence the cubic positive recurrence belongs to an all-odd-dimension family of explicit integer transfers, although higher-dimensional classification of norm fibers is not asserted. For odd prime d=p, Frobenius in the ring modulo p gives (α_p−1)^p≡1, and binomial lifting gives U_p^(p^s)≡1 mod p^s. This supplies exact simultaneous coordinate resets *modulo* p^s at unbounded precision, without any claim of an actual zero-coordinate return.

For the cubic d=3, the complex embeddings of U have modulus U^(-1/2). If t_n is the real embedding of (1+α)U^n, the Fourier bound in the exact note gives x_n,y_n,z_n by rounding t_n/3, t_n/(3α), t_n/(3α²) to the nearest integer: for n≥1 the unscaled error is at most (2/3)√(3/t_n)<1/2; n=0 is direct. Here exponentially small error actually implies exact recovery because an explicit threshold is crossed, unlike the decimal relative-error records.
