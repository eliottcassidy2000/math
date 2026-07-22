---
id: THM-2097
title: "Cubic-fiber degree-thirteen Faber-tail gate"
status: >
  PROVED. For every n prime to three, the two adjacent omitted coefficients
  Phi_n=3[z^-1](z^3+pz+q)^(n/3) and
  R_n=3[z^-2](z^3+pz+q)^(n/3) vanish simultaneously only at (p,q)=(0,0).
  The proof is all-degree: their simultaneous vanishing makes the binomial
  series (1+pw^2+qw^3)^(n/3) terminate, and unique factorization then forces
  1+pw^2+qw^3 to be a cube. Applying this lemma together with an exact
  centering-pole resultant excludes reduced complementary fiber degree
  thirteen in the cubic source-fiber stratum of THM-2084. Thus a non-tame pair
  in that stratum has reduced complementary fiber degree at least fourteen.
  This is not a generic-cover, Jelonek-curve, or full JC(2) theorem.
source: codex-2026-07-22-JC2-degree13
depends_on:
  - THM-2084
related:
  - THM-2063
  - THM-2071
  - MISTAKE-229
script: 04-computation/jc2_cubic_fiber_degree13_tail_gate_codex_20260722.py
output: 05-knowledge/results/jc2_cubic_fiber_degree13_tail_gate_codex_20260722.out
script_sha256: d1221581ec75cab4ad18e8606bff8a1c8da785d2ed8b6ab1b5ad5bed81293c45
output_sha256: aafc72942a01526bedabbe68c9aed886dc9bc8bb54693e31db8527c7d90188b1
hash_basis: repository blobs with LF line endings
---

# THM-2097 -- cubic-fiber degree-thirteen Faber-tail gate

## 1. Statement

Retain the notation and hypotheses of
`THM-2084-cubic-fiber-low-complement-gauss-manin-gate.md`.  Thus one member
`P` of the output pencil of a planar Keller pair has degree three along a
fixed linear source fiber, and

```text
mu=min_H deg_y(Q-H(P)),                    H in C[T].       (1)
```

Then

```text
mu != 13.                                                   (2)
```

Consequently, if the pair is non-tame, then

```text
mu >= 14.                                                   (3)
```

The new structural input is the following all-degree fact.  For `3` not
dividing `n`, put

```text
P_0=z^3+pz+q,
Phi_n=3[z^-1]P_0^(n/3),
R_n  =3[z^-2]P_0^(n/3).                                    (4)
```

Then, as equations in `(p,q) in C^2`,

```text
Phi_n=R_n=0                 iff                 p=q=0.      (5)
```

This adjacent-tail noncollision is genuinely all-degree.  It is stronger
than the degree-thirteen calculation used below and isolates the primitive
sidecar that is lost if one keeps only the first integral `Phi_n`.

## 2. The all-degree adjacent-tail lemma

Set `w=z^-1` and write

```text
S(w)=1+pw^2+qw^3,
S(w)^(n/3)=sum_(k>=0) c_k w^k.                              (6)
```

Thus `Phi_n=3c_(n+1)` and `R_n=3c_(n+2)`.  Logarithmic
differentiation gives

```text
S F'=(n/3)S'F,

k c_k +(k-2-2n/3)p c_(k-2)+(k-3-n)q c_(k-3)=0,             (7)
```

with `c_j=0` for `j<0`.  Suppose that the two coefficients in (4) vanish.
At `k=n+3`, the coefficient of `q c_n` in (7) is zero, so `c_(n+3)=0`.
At `k=n+4`, both preceding terms already vanish, so `c_(n+4)=0`.
Equation (7) now inducts forward and gives

```text
c_k=0                         for every k>=n+1.             (8)
```

Hence `F=S^(n/3)` is a polynomial `G` and

```text
G^3=S^n.                                                     (9)
```

Because `gcd(n,3)=1`, unique factorization in `C[w]` forces every irreducible
factor of `S` to occur with multiplicity divisible by three.  Therefore
`S=T^3`.  Normalize `T(0)=1`.  Since `deg S<=3`, write `T=1+lambda w`.
The missing linear coefficient of `S` gives `3lambda=0`, whence `lambda=0`
and `S=1`.  Thus `p=q=0`.  The converse in (5) is immediate.  This proves
the lemma.

The mechanism is a manufactured valuation at infinity in its most elementary
form: two missing consecutive tail orders force *every* later order to be
missing, after which the exponent denominator is detected by factor
multiplicities in the UFD.  No monodromy, convergence, or genericity is used.

## 3. Degree-thirteen normal form

After the cube gate and depressed-cubic change of THM-2084, work over
`K=C(x)` with

```text
P=z^3+pz+q,
L=(p'z+q') partial_z-(3z^2+p) partial_x,
L(Q)=kappa/U.                                               (10)
```

Let

```text
E_m=Pol_z P^(m/3).                                          (11)
```

The recursive normal form from THM-2084 says that a reduced mate of degree
thirteen is, modulo `C[P]`,

```text
Q=cE_13+c_11E_11+c_10E_10+c_8E_8+c_7E_7
    +c_5E_5+c_4E_4+c_2E_2+c_1E_1,          c != 0,         (12)
```

with constant coefficients.  If `Phi,R` are the same linear combinations
of the columns in (4), then

```text
Phi'=0,                         R'=kappa/U.                 (13)
```

The top pair is

```text
Phi_13=(65/6561)p(p^6-63p^3q^2+189q^4),
R_13  =(91/6561)q(5p^6-60p^3q^2+27q^4).                  (14)
```

## 4. The centering coordinate has no finite pole

Let `v` be a finite-place valuation and suppose the centering coordinate
`h` has `v(h)=-H<0`.  Put `rho=-v(p)`.  Polynomiality of the original cubic
coefficient

```text
D=h^3+ph+q                                                   (15)
```

gives three regimes.

If `rho>2H`, then `q~-ph`.  The pole orders of the three monomials in
`Phi_13` are

```text
7rho,                    6rho+2H,                    5rho+4H.
```

They are strictly decreasing, and every lower `Phi_m` in (12) has smaller
order.  The `p^7` term therefore contradicts `Phi'=0`.

If `rho<2H`, then `q~-h^3`.  Directly from (11),

```text
E_13(1;0,-1)=35/243 != 0.                                  (16)
```

Thus `cE_13(h)` has the unique pole order `13H`; every lower representative
in (12) has order at most `mH<13H`.  This contradicts polynomiality of
`Q(x,0)`.

It remains to treat `rho=2H`.  Write

```text
p~a h^2,                         q~-(1+a)h^3.               (17)
```

Constancy of `Phi` and polynomiality of `Q(x,0)` respectively force the two
primitive polynomials

```text
f(a)=a(a^6-63a^5+63a^4+693a^3+1134a^2+756a+189),

e(a)=13a^6-78a^5-351a^4-468a^3-234a^2+27                  (18)
```

to vanish.  Indeed the actual leading coefficients are `(65/6561)f(a)` and
`(35/6561)e(a)`.  Exact Euclidean elimination gives

```text
Res_a(e,f)=-7856816824690792659
          =-3^21 17^5 23^2 != 0.                           (19)
```

So the balanced regime is empty as well, and

```text
h in C[x].                                                  (20)
```

## 5. The depressed coefficients have no finite poles

As in THM-2084, (20) and the original polynomial coefficients give

```text
Up in C[x],                         q+hp in C[x].           (21)
```

If `p` had a finite pole of order `rho>0`, then `q=O(p)` there.  In
`Phi_13` the `p^7` term has order `7rho`, whereas its other two terms have
orders at most `6rho` and `5rho`; every lower `Phi_m` has ordinary
`(p,q)`-degree at most six.  Again `Phi` could not be constant.  Hence `p`
has no finite pole, and (21) makes `q` polynomial.

Now `R` is polynomial, so (13) says that `kappa/U=R'` is polynomial.  Since
`U` is a polynomial and `kappa` is a nonzero constant, `U` is constant.
The depressed change is therefore an honest triangular polynomial coordinate,
and (13) reduces to

```text
p,q in C[x],                  Phi in C,       R affine nonconstant. (22)
```

## 6. The degree-thirteen branch at infinity is empty

Suppose first that both `p` and `q` are nonconstant, of degrees `a,b>0`.
The support of `Phi_13` is

```text
(7,0), (4,2), (1,4).                                      (23)
```

Every support monomial of every lower `Phi_m` in (12) is strictly
componentwise dominated by one of these three.  Lower representatives
therefore cannot join the leading cancellation.  The three degrees in (23)
are

```text
7a,                       4a+2b,                       a+4b.
```

Their maximum can be attained more than once only when `3a=2b`.  Write
`a=2s,b=3s`, and let `p_0,q_0` be the nonzero leading coefficients.
Constancy of `Phi` forces

```text
Phi_13(p_0,q_0)=0.                                          (24)
```

Under the same `2:3` weight, every lower `R_m` has degree at most `(m+2)s`,
strictly less than the degree `(13+2)s` of `R_13`.  The all-degree lemma (5), applied
to the nonzero pair `(p_0,q_0)`, now gives

```text
R_13(p_0,q_0) != 0.                                        (25)
```

Hence `R` has degree `15s`, contradicting (22).

For a completely explicit degree-thirteen certificate, put
`Y=q_0^2/p_0^3`.  The two factors in (24)--(25) are

```text
F=1-63Y+189Y^2,                  G=5-60Y+27Y^2,

(63Y-134)F+(105-441Y)G=391.                                (26)
```

Thus they cannot vanish together.

If `q` is constant and `p` is nonconstant, the `p^7` term of `Phi_13` is
uniquely highest.  If `p` is constant and `q` is nonconstant, the `q^5` term
of `R_13` is uniquely highest among all terms in (12).  If both are constant,
then `R'=0`.  Each case contradicts (22).  This proves (2).

Finally THM-2084 already leaves only the tame degrees one and two below
thirteen.  Combining that theorem with (2) proves (3).

## 7. Scope and surviving frontier

This theorem closes one more **source-fiber degree**, not a generic mapping
degree.  The cubic source-fiber stratum remains open, beginning at reduced
degree fourteen.  The adjacent-tail lemma settles the top primitive
noncollision for every degree prime to three, but two other effects remain:

1. at a centering pole one needs the different pair
   `E_n(1;a,-1-a), Phi_n(a,-1-a)`;
2. for even top degree, a lower Faber representative can join a different
   upper Newton edge, so the top pair alone does not control every branch.

Those losses explain exactly why (5) is structural progress without being an
all-degree cubic closure.  No claim about generic cover degree, Jelonek curves,
DC(2), or full `JC(2)` is made.

## 8. Exact referee

The companion script checks over exact rational polynomial rings:

- `L(E_13)=z Phi_13'+R_13'` and the formulas (14);
- the control (16) and exact resultant (19);
- strict upper-Newton separation from every lower reduced `Phi_m`;
- the Bezout identity (26);
- a projective-gcd sweep for odd `n<=61` and a discriminant-zero control in
  both residue classes; and
- the necessary excluded boundary `(p,q)=(0,0)`.

Both normal and `python -O` runs end in `RESULT: PASS`.  The recurrence/UFD
argument in Section 2, rather than the finite sweep, supplies the all-degree
quantifier.  QED.
