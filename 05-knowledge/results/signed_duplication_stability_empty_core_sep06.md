# A sharp quartic duplication margin and stability of its root pairs

**Status: PROVED ANALYTICALLY + INDEPENDENTLY AUDITED + FINITE-EXACT.**
Bounded wildcard continuation, 2026-09-06. No theorem ID or external-priority
claim. The theorem is about the actual ordinary carrier polynomial; no
positivity transfer to LRC is asserted.

## Inheritance and the changed question

The closest proved mechanism is
[THM-4440 — signed duplication SOS and real-rooted Laurent return](../../01-canon/theorems/THM-4440-signed-duplication-sos-and-real-rooted-laurent-return.md),
with its [full proof and Hermite comparison](nc2_signed_duplication_overnight_hexagon_sep05.md).
For `H(s)=H(0) prod(1+r_i s)`, the terminal subset square at `e_k(r)=0`
gives

```text
-[s^(2k)]H^2 >= H(0)^2 e_k(r_1^2,...,r_n^2)/(2k-1).
```

That theorem already has an exact SOS, its equality classification, and a
strict negative Hermite limit. Merely retaining its final square would not
be a new result. The question here is whether the multiplicative structure
of the subset-product vector strengthens its normalized margin and controls
the original roots. For `n=4,k=2`, the answer is a sharp `7/9` in place of
`1/3`, with an exact stability remainder. The dimension restriction matters.

The inherited hostile with positive individual opposite-coefficient product
is `H=(s+1)(s+2)(s-10)(s-4)`. A new degree-five hostile below blocks a
dimension-free extension of `7/9`. The least-used retained sidecar is the
sign partition and the two within-sign pair discriminants.

The live board is: full subset-product energy; the actual square coefficient;
same-sign root pairs; gauge normalization; the quartic boundary; the existing
Hermite limit. The map sends a real-rooted quartic with vanishing middle
coefficient to two root sums and products, retaining their pair variances.
It preserves the target coefficient and all root-product weights. Keeping
only the final subset norm loses the pair variances; the identity below
restores them. Arbitrary subset vectors and arbitrary complex cores are not
in the target domain.

## 1. Sharp normalized margin

Let `G(s)=sum_(j=0)^4 a_j s^j` be real-rooted, of exact degree four, with
`a_0!=0` and `a_2=0`. Write

```text
G(s)=a_0 prod_(i=1)^4(1+r_i s),       r_i real and nonzero,
E=a_0^2 sum_(i<j) r_i^2 r_j^2 >0,
D=-[s^4]G(s)^2.
```

Then

```text
D >= (7/9)E >0.                                           (1)
```

Equality holds exactly when, up to permutation, the reciprocal roots are

```text
(p,p,-q,-q),   p,q>0,   p/q=2+sqrt(3) or 2-sqrt(3).        (2)
```

Equivalently, up to nonzero real scalar and variable scaling,
`G(s)=(1+sqrt(2)s-s^2)^2`. Its coefficients are
`(1,2sqrt(2),0,-2sqrt(2),1)`, so `D=14`, `E=18`. The quadratic factor has
discriminant six and two real roots of opposite signs. This proves sharpness
inside the stated domain, including repeated real roots.

Here is an elementary proof independent of any spectral estimate. Temporarily
divide out `a_0` and write `e_j=e_j(r)`, `T=-e_1e_3`, `V=e_4`. Expanding
`e_2^2`, and using `e_2=0`, gives

```text
E/a_0^2=2(T+V),              D/a_0^2=2(T-V).               (3)
```

If `V<0`, positivity of `E` immediately gives `D/E>1`. If `V>0`, the
four reciprocal roots cannot have the same sign, since that would make
`e_2>0`. Thus there are two of each sign. Label them
`p_1,p_2,-q_1,-q_2`, with all `p_i,q_i>0`, and put

```text
x=p_1+p_2, y=q_1+q_2, P=p_1p_2, U=q_1q_2,
X=x/sqrt(P)>=2, Y=y/sqrt(U)>=2.
```

The cancellation equation is `xy=P+U`. Also
`e_1=x-y`, `e_3=xU-yP`, and `V=PU`, whence

```text
T/V=X^2Y^2-X^2-Y^2 >= 8.                                (4)
```

The last inequality follows by setting `alpha=X^2-4>=0`,
`beta=Y^2-4>=0`: its left side is
`8+3(alpha+beta)+alpha*beta`. Combining (3)-(4) gives (1).
Equality requires `alpha=beta=0`, so both same-sign pairs coincide.
Then `e_2=0` is `p^2+q^2=4pq`, yielding exactly (2).

## 2. Exact stability and root location

In the two-positive/two-negative case define the dimensionless pair spreads

```text
alpha=(p_1-p_2)^2/(p_1p_2),
beta =(q_1-q_2)^2/(q_1q_2),
z=3(alpha+beta)+alpha*beta.
```

The proof gives the full identities

```text
D-(7/9)E
 =(4a_0^2/9)[3U(p_1-p_2)^2+3P(q_1-q_2)^2
                         +(p_1-p_2)^2(q_1-q_2)^2],        (5)
D/E=(7+z)/(9+z).                                         (6)
```

Consequently, if `0<=epsilon<2/9` and `D/E<=7/9+epsilon`, then the root
list necessarily has two roots of each sign and

```text
z <=81epsilon/(2-9epsilon),
alpha+beta <=27epsilon/(2-9epsilon).                      (7)
```

Thus proximity to the optimal normalized margin quantitatively forces
each same-sign root pair to approach a repeated pair. These spreads are
unchanged when a pair is replaced by its reciprocal magnitudes, so (7)
also describes the actual roots of `G`, not only the `r_i`.

There is also a location constraint between the two pairs. Put
`t=max(sqrt(P/U),sqrt(U/P))>=1`. The cancellation equation gives

```text
t+1/t=XY=sqrt((4+alpha)(4+beta)) >=4,
and therefore t>=2+sqrt(3).                              (8)
```

When `z<=Z`, one additionally has `XY<=sqrt(16+4Z/3)`: indeed
`X^2Y^2=16+z+alpha+beta` and `alpha+beta<=z/3`. Setting
`L=sqrt(16+4Z/3)` yields the explicit interval

```text
2+sqrt(3) <=t<= (L+sqrt(L^2-4))/2.                       (9)
```

With `Z=81epsilon/(2-9epsilon)`, (7) and (9) characterize quantitative
closeness to the extremal root configuration modulo scaling, permutation,
and simultaneous sign reversal. No absolute root scale can be controlled
by a gauge-invariant quantity.

## 3. Correct monomial and scalar normalization

For `H(s)=s^ell G(s)`, with `G` as above, extract at `k=ell+2`.
Then `[s^k]H=0` and `[s^(2k)]H^2=[s^4]G^2`. The index must be shifted
together with the monomial factor; otherwise this is a different problem.

For a core admitting a real gauge define its normalized magnitude by

```text
J(H,k)= |[s^(2k)]H^2| /
       (|G(0)|^2 sum_(i<j)|r_i r_j|^2),   k=ord_0(H)+2.  (10)
```

In a real gauge at the zero coefficient, `J=D/E`. Under

```text
H(s) -> lambda*s^b*H(sigma*s),       k -> k+b,
lambda,sigma nonzero complex,
```

take integer `b` with `ell+b>=0` to keep `H` polynomial. Both numerator and
denominator in (10) acquire the factor
`|lambda|^2 |sigma|^(2k)`. In the denominator, this follows from
`G(0)->lambda*sigma^ell G(0)` and `r_i->sigma*r_i`.
Thus the magnitude margin and its equality/stability statements are
invariant. A real negative sign is meaningful only in a declared real
gauge. The denominator is an explicitly defined root-product energy; it
is not being identified with a different Hermitian Laurent coefficient.

For the Laurent consumer `f(u)=u^(-a)R(u^d)`, at a congruent mass `m`,
THM-4440 retains the exact map

```text
H=R^m,    k=am/d,
CT(f^m)=[s^k]H,    CT(f^(2m))=[s^(2k)]H^2.
```

Whenever the effective ordinary degree is four and the effective index is
two, (1), (5), and (10) therefore give the sharper return margin. The simplest
actual first-return sharp example is

```text
f(u)=u^(-2)(1+sqrt(2)u-u^2)^2,
CT(f)=0, CT(f^2)=-14, E=18.
```

Scalar or variable gauges of this Laurent example preserve `J=7/9`.
This result does not enlarge the set of ordinary cores known to be
real-rooted; in particular it does not repair the compressed-degree-four
trinomial obstruction in THM-4440.

## 4. Hostile boundaries and the inherited Hermite comparison

The new constant is already false in degree five. Take reciprocal roots

```text
(1,1,-1/6,-1/6,-13/60).
```

Their polynomial coefficients are

```text
(1,29/20,0,-769/2160,19/216,-13/2160).
```

It has the required vanishing coefficient and only real roots, but
`D/E=18501/26101<7/9`. The coefficient identities (3) still hold; the first
failed step is identifying `V=e_4` with the product of the entire root list
and then with `PU` for two exhaustive same-sign pairs. In degree five,
`e_4` is a sum of five products, so the pair argument giving `T/V>=8` is
unavailable. This preserves THM-4440's general terminal bound and sharply
limits this improvement.

The existing Hermite limit is also not a proof of the quartic constant.
At a zero `x=+/-1` of `He_2`, the limiting carrier is
`exp(x s-s^2/2)`. Its doubled fourth coefficient is `-5/6`, while the
limiting energy `e_2(r_i^2)` is `1/2`, so its normalized margin is `5/3`.
The finite degree goes to infinity in that limit; the root-pair equality
configuration above is a different regime. No new Hermite identity is
claimed here.

The inherited positive-pair-product hostile has normalized margin
`233/273>7/9`. The sign-count hostile with reciprocal roots `(1,1,1,-1)`
has margin `5/3>1`. Dropping real-rootedness still fails at `1+s^2`, and
an unshifted index outside the effective interior still fails at `s^2`.

## 5. Exact controls and audit manifest

[Source](../../04-computation/signed_duplication_stability_empty_core_sep06.py)
and [output](signed_duplication_stability_empty_core_sep06.out).

```bash
python3 -B 04-computation/signed_duplication_stability_empty_core_sep06.py
python3 -B -O 04-computation/signed_duplication_stability_empty_core_sep06.py
```

The bounded universe consists of the 56 multisets of three numbers from
`{-3,-2,-1,1,2,3}`; solve the cancellation equation for the fourth root,
reject six undefined/zero-root cases, and retain all 44 distinct sorted
quartic lists. The checker independently multiplies the literal factors
and their squares, compares root-product energy, checks the exact remainder
and inverse stability, tests three scalar/variable gauges and two monomial
shifts per row, and retains the named hostile controls. It verifies the
sharp carrier in exact `Q(sqrt(2))` arithmetic. All **543 explicit gates**
pass; normal and optimized output are identical. These controls challenge
the analytic argument, not justify its universal quantifiers.

```text
source SHA256 54b5b6ceffc26697f966d56ffeb57dbaa00ebd121416b0f9c020af6de35ed489
output SHA256 088ba7367328bdfbf5a99a0c8292036ab1b06fe8aed694c0b68ae45781e5cf7d
trace  SHA256 9837c6d36b6b69b8073cbea3a8c9a47a2ca2f949d9ee1cad37fe6791ee8c913a
```

Independent final audits by the root and orthogonal-return agents: **PASS**. Both checked the full analytic inequality, equality and stability statements, gauge factors, and degree-five failure boundary; the orthogonal-return audit also replayed all 543 explicit gates.
