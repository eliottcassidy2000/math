# Joint binary controller: all moduli, not just powers of two and three

This strengthens the [short two-prime construction](README.md):
for **any** positive modulus D and odd integer target N, a literal
strong tournament from the same clone core has H≡N mod D. The
proof constructs its binary tag. It is not an absolute selector.

## Control (u,F(u)) at an arbitrary odd modulus

Keep u odd after each macro word, and write U=u+1. For appending
the words A=`1`, B=`01`, C=`001`, the binary recurrence gives
the following affine updates on (U,F):

    A : (2U,    F)
    B : (4U−2,  F+U)
    C : (8U−6,  F+3U−2).

Modulo any odd M each is a bijection of the finite state space.
Consequently inverses are realized by sufficiently many positive
repetitions of the corresponding word. Compose right to left.
The normalized word B'=A^(-2)B adds −1/2 to U and U to F.
Conjugating it by A gives D'=A^(-1)B'A, with
U→U−1/4,F→F+2U. Compare three D' moves to
E=A^(-3)C: they both shift U by −3/4, and their F increments
are 6U−3/2 and 3U−2 respectively. Therefore the positive-word
realization of (D')^(-3)E fixes U and adds −3U−1/2 to F.
Translations can take U to 0 and to any other residue, since
1/4 is a unit; at U=0 this F-only shift is a unit mod M.
Choose U=0, set the needed F after accounting for the final
U-translation, then move U to its target. This controls every
pair of residues modulo M, from any odd seed.

## Set the dyadic pair without disturbing the odd pair

Given L and odd residues u₀,F₀ modulo 2^L, append a fixed
guard `1^g` with g≥L, then the prescribed L-bit suffix for u₀.
The guard sets U≡0 mod 2^L and leaves F unchanged. The suffix
sets u≡u₀ and adds a fixed known amount to F. Indeed the
starting dyadic U-class is fixed; a zero appends by
u→2u,F→F+u+[u odd], and a one doubles and adds one to u
without changing F. That amount is even. Modulo odd M this
guard/suffix is bijective on the two coordinates: parity is
fixed by the preceding bit at each step and the affine
triangular maps have unit determinant up to a power of 2.
So first use the odd controller to reach the inverse image
of the requested odd state under the guard/suffix.

To adjust F before the guard, while fixing the odd state,
use W=`1^r001` with r≥L and r+3 divisible by ord_M(2).
(For M=1 use r+3 arbitrary.) Modulo M this word gives

    U→U−6,  F→F+3·2^r U−2.

W^M acts as identity: the summed drift is a multiple of M,
using that M is odd. Modulo 2^L, W always changes F by −2,
since 2^r≡0. Thus W^M changes F by −2M; repeating it chooses
every odd F residue. Finish with the guard/suffix. The final
positive, odd address achieves the independently chosen states
mod 2^L and M. This is a literal concatenation, not a formal
combination of counts.

## Solve the Hamilton count congruence

Choose L≥5 at least the target's dyadic exponent, and apply
the [16-residue table and power-of-three lift](README.md)
to get k and a base dyadic pair (u₀,F₀). For the odd part
M of the target modulus, the same exact formula for odd u is

    H(k,u)=R_k u+h_k F(u)+2h_k−5.

The coefficients have gcd 1:
R_k+h_k=36·3^k while h_k=12·2^k−1 is coprime to 6.
An integer Bézout pair therefore chooses u and F residues
satisfying the odd-side congruence, including when M contains
primes other than 3. In tandem with the previously selected
dyadic pair, the controller gives the address of a strong
host with the desired combined residue. For an odd modulus
there is no restriction on its target residue, since an odd
representative can always be chosen.

Thus the Hamilton spectrum has no congruence obstruction at
any finite modulus on the odd layer, even in this specific
clone family. This does not ensure an answer at a fixed
absolute integer: words realizing inverses, the odd register,
and the length needed for the dyadic register can all grow
with the modulus. Nor does a local congruence certificate
distinguish 7 or 21 from realizable odd integers. The remaining
absolute-selection task needs a construction with numerical
height control, not merely more residue classes.
