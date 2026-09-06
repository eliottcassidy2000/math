---
id: THM-4417
title: "Width-two Laurent first-return parabolic critical bound"
status: >
  PROVED relative to CITED Leau-Fatou and petal-critical-point theorems;
  algebra independently audited and FINITE-EXACT controls reproduced.
  For f=z^(-2)R, R(0)!=0 and deg R=d>=3, the first nonzero constant-term
  moment is bounded by the number of distinct nonpole critical points of
  -R(0)z/R(z), hence by the number of distinct roots of R and by d.
  The sharp M+N bound for min(M,N)>=3 remains OPEN.
source: synthesis-sep05 / two-small-root involution and parabolic critical budget
depends_on:
  - THM-2111-effective-compound-root-bound-for-one-variable-constant-terms
proof: 05-knowledge/results/nc2-width-two-parabolic-seed-synthesis-sep05.md
script: 04-computation/nc2_width_two_parabolic_seed_synthesis_sep05.py
output: 05-knowledge/results/nc2_width_two_parabolic_seed_synthesis_sep05.out
script_sha256: 7eaa6ae69467252271e224ff6e5081d81e8ef97b8ba81860252f32aa1284a783
output_sha256: 89b7efc2c582aa7420653e7f14680e7d12bbe6a999f119b9d1ed11836d4b4085
hash_basis: raw LF bytes
---

# THM-4417 -- a linear width-two first-return bound

**PROVED relative to the cited dynamical input.** See the
[complete proof, primary-source pins, exact controls, and concurrent-work comparison](../../05-knowledge/results/nc2-width-two-parabolic-seed-synthesis-sep05.md).

Let R in C[z], r=R(0)!=0, and d=deg R>=3. Put f=z^-2 R(z), let
m_* be the first m>=1 with CT(f^m)!=0, and set

```text
T(z)=-r z/R(z),  D_R=(R-zR')/gcd(R,R').
```

Then

```text
m_* <= number of distinct zeros of D_R
    <= number of distinct roots of R <= d.
```

No sign, reality, squarefreeness, or genericity hypothesis is imposed.
Reflection z->z^-1 gives m_*<=M+N whenever min(M,N)=2; THM-2111
already supplies width one. Odd-degree binomials attain m_*=d, and
z^-2+z^-1+z-z^2 attains d=4.

Here is the transfer responsible for the sharpening. The two small roots
of X^2-tR(X) define a local involution iota at zero. THM-2111's logarithmic
root-product identity gives

```text
ord_0(T^2-id)=2m_*+1,
[z^(2m_*+1)](T^2-z)=-2 CT(f^m_*)/(m_* r^m_*).
```

Since T'(0)=-1, this is exactly m_* cycles of parabolic petals.
Milnor's [primary lecture notes](https://abel.math.harvard.edu/archive/118r_spring_05/docs/milnor.pdf),
Sections 7.2, 7.5, and 10.3 with its following paragraph, require at least
m_* distinct critical points in their basin. Poles and infinity land
exactly on zero and lie in its Julia backward orbit, so cannot supply those
critical points. The eligible critical points are precisely the distinct
zeros of D_R. Its degree equals the number of distinct roots of R.

The full proof has two independent algebra reviews and a checked primary
citation. Exact Laurent convolution versus rational composition passes
480 coefficient-box cases, sharp binomials, and repeated-root hostiles.
No numerical dynamics, external priority, or new Lean build is claimed.

For three or more small roots, the root-product quotient is not a local
root-swap map; the argument stops there. This improves only the named
lowest-face seed in NC2/GMC(2), not a uniform coefficient-independent final
Gaussian-moment bound.
