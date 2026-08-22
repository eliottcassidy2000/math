---
id: THM-3404
title: "Factorized Danielewski principal parts and finite-cover obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY PROOF-AUDITED.  For every field k
  and nonzero F=c product_i f_i^(e_i) in k[v], with distinct irreducible
  factors, the normal surface R=k[v,X,Y]/(XY-F) has an exact graded
  principal-part module R_X/R=direct_sum_(q>=1) X^(-q)k[v]/(F^q), split by
  CRT into one jet tower for each factor.  Hence X^(-q)h(v) is regular iff
  F^q divides h.  More generally, every single homogeneous pole
  w=b(v)X^(-a) has a degree-one-generated exact clearing filtration
  I_n=J(w)^n.  For nonconstant F, the height-one orders and Nagata
  presentation give Cl(R)=Z^r/<(e_i)>=Z^(r-1) direct_sum Z/gcd(e_i); the
  unit case is Laurent factorial.  Finite dominant covers
  preserve the principal-part obstruction, and Weil conorm plus norm shows
  that only torsion divisor classes can die.  Consequently R has a finite
  dominant normal factorial cover iff r<=1.  The canonical common-root
  cover kills exactly its cyclic torsion subgroup but preserves every old
  principal-part packet in the corresponding multiplied degree.  Over C,
  for nonconstant F a finite polynomial-plane cover therefore forces the
  original one-factor pure-power model; no nonmonomial JC branch or JC(2)
  follows without a physical finite embedding.
source: root-2608-q8-multiblocker-2026-08-15
audit: self-contained Serre/Nagata/Laurent/CRT/conorm proof; independent arbitrary-field, grading, factor-duplication, characteristic, nonnormal, nonfinite, and cancellation audit; exact rational-polynomial companion
depends_on:
  - THM-2700-danielewski-s3-resolvent-standard-plane-exclusion
  - THM-3397-torsor-killing-versus-effective-boundary-valuations
related:
  - THM-3383-terminal-monomial-cone-polynomiality-fork
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
script: 04-computation/jc_factorized_danielewski_principal_parts_thm3404.py
output: 05-knowledge/results/jc_factorized_danielewski_principal_parts_thm3404.out
script_sha256: 0a7a0c0372b1195bb61c4a9018169ef38977f97a4a78e045bee4fbae84b3cf10
output_sha256: d87195850315223f794ad3776301d3449f731fbdd11cbb0f19d64a3e48ce67ec
semantic_sha256: 6459529a71a434aa0856d38d3d21bd11314a05dec449773331aac907543690c9
hash_basis: LF-normalized bytes
---

# THM-3404 -- factorized Danielewski principal parts and finite-cover obstruction

**PROVED + VERIFIED-EXACT + INDEPENDENTLY PROOF-AUDITED.**

## 1. Inheritance and exact connection

[THM-3383](THM-3383-terminal-monomial-cone-polynomiality-fork.md)
computes one Laurent boundary cone when the invariant equation is
`XY=(1+cv)^e`.
[THM-3397](THM-3397-torsor-killing-versus-effective-boundary-valuations.md)
turns its two boundary orders into an exact denominator filtration and proves
that the same rational function cannot lose a divisorial pole on a finite
dominant normal pullback.  The canonical hostile is still a killed torsor
with a surviving pole; the corrected near miss is localization, which simply
deletes the boundary.

[THM-2700](THM-2700-danielewski-s3-resolvent-standard-plane-exclusion.md)
already proves, over `C`, normality and the class-group formula for a general
Danielewski polynomial.  Those facts are inherited, not recanonized as new.
The new payload is the complete factorwise principal-part tower, its exact
finite-cover behavior, and the resulting factorial-cover obstruction.

| field | exact connection |
|---|---|
| source | a factorized one-variable suspension `XY=F(v)` and Laurent target views |
| target | polynomial regularity in `R`, or after a finite dominant cover |
| map | put every view in unique `X`-Laurent normal form and pull its boundary orders/principal parts |
| preserved | factor arms, multiplicities, Laurent degree, coefficient jets, and poles under finite pullback |
| destroyed by `Cl(R)` | the sign, magnitude, location, and leading coefficient of a principal divisor |
| required sidecar | the oriented graded module `R_X/R`, with its CRT arm decomposition and any symmetry action |
| cheapest decisive tests | `F=v(v-1)` principal-part cancellation and `F=v^2(v-1)^2` torsion-killing double cover |

## 2. Setup, normality, and every height-one valuation

Let `k` be any field, put `A=k[v]`, and write the unique aggregated
factorization

```text
0 != F=c product_(i=1)^r f_i^(e_i),
c in k^*,  f_i monic pairwise distinct irreducibles,  e_i>=1.   (1)
```

Set

```text
R=A[X,Y]/(XY-F).                                             (2)
```

If `r=0`, then `F` is a unit and `R=A[X,X^(-1)]`; all boundary
statements below are empty.  Assume `r>=1` unless stated otherwise.

The polynomial `XY-F` is prime, so `R` is a two-dimensional hypersurface
domain and hence satisfies Serre's `S_2`.  We prove `R_1` without derivatives,
so the argument works over imperfect fields as well.

If a height-one prime avoids `X`, localizing places it in the regular Laurent
UFD

```text
R_X=A[X,X^(-1)],             Y=F/X.                       (3)
```

If a height-one prime contains `X`, it contains `F`, hence one unique `f_i`,
and dimension forces it to be

```text
D_i^X=(X,f_i).                                              (4)
```

At its generic point, `Y` and `F/f_i^(e_i)` are units.  Thus `f_i` is a
uniformizer and

```text
ord_(D_i^X)(f_i)=1,   ord_(D_i^X)(X)=e_i,
ord_(D_i^X)(Y)=0.                                           (5)
```

Every other height-one prime is uniquely

```text
P_pi=(pi)A[X,X^(-1)] intersect R                           (6)
```

for one Laurent-associate class of irreducibles
`pi in A[X,X^(-1)]`, and its valuation is the `pi`-adic order in `(3)`.
In particular,

```text
P_(f_i)=D_i^Y=(Y,f_i),
ord_(D_i^Y)(f_i)=1,   ord_(D_i^Y)(Y)=e_i,
ord_(D_i^Y)(X)=0.                                           (7)
```

All height-one local rings are DVRs.  Serre's criterion proves that `R` is
normal.  Refactoring after an arbitrary field extension gives the same
proof, so no squarefreeness, separability, or characteristic-zero hypothesis
is needed for normality.  The divisor identities are

```text
div(f_i)=D_i^X+D_i^Y,
div(X)=sum_i e_i D_i^X,
div(Y)=sum_i e_i D_i^Y.                                   (8)
```

## 3. Exact Laurent normal form and principal-part tower

Inside `(3)`, relation `(2)` gives

```text
R=A[X,F X^(-1)]
 = direct_sum_(m>=0) X^m A
   direct_sum_(q>=1) X^(-q)F^q A.                         (9)
```

Indeed a monomial `X^aY^b` of Laurent degree `a-b=-q<0` has
coefficient `F^b`, with `b>=q`; conversely `X^(-q)F^q h=Y^q h`.
Distinct Laurent degrees cannot cancel.  Therefore

```text
R_X/R = direct_sum_(q>=1) X^(-q) M_q,
M_q=A/(F^q),                                               (10)
```

as a graded `A`-module.  It is a graded `R`-module with cross-degree action

```text
X: M_q -> M_(q-1) by reduction,       M_1 -> 0,
Y: M_q -> M_(q+1),                    [h] -> [Fh].         (11)
```

The individual fixed-`q` summands are not separate `R`-submodules.  The
pairwise coprime factorization `(1)` gives the canonical CRT refinement

```text
M_q = direct_sum_i A/(f_i^(q e_i)),
R_X/R = direct_sum_i direct_sum_(q>=1)
        X^(-q) A/(f_i^(q e_i)).                           (12)
```

Each full factor arm across all degrees is stable under `(11)`.  Thus the
sidecar retains not only the pole orders `q e_i`, but every truncated local
coefficient jet.

For `0!=h in A`, let `a_i=ord_(f_i)(h)`.  Equations `(5)--(8)` give

```text
ord_(D_i^X)(X^(-q)h)=a_i-qe_i,
ord_(D_i^Y)(X^(-q)h)=a_i.                                (13)
```

Every factor of `h` not dividing `F` contributes only a nonnegative vertical
order, and every other Laurent prime has order zero.  Hence, for `q>=1`,

```text
X^(-q)h(v) in R
iff F^q divides h,
{h in A:X^(-q)h in R}=(F^q).                             (14)
```

The obstruction is exactly the class `[h] in M_q`.  This is stronger than a
list of valuations when several terms will be added.  For example, with
`F=v(v-1)`, both

```text
X^(-1),             X^(-1)(F-1)                         (15)
```

have order `-1` on both `X`-boundary arms, but their leading principal parts
cancel and

```text
X^(-1)+X^(-1)(F-1)=X^(-1)F=Y.                           (16)
```

In `A/(F)`, this is the exact cancellation `[1]+[-1]=0`.  Signed orders of
the summands alone would miss it.

## 4. Every single homogeneous Laurent pole

The THM-3397 filtration also generalizes exactly.  Let

```text
w=b(v)X^(-a),       a>=0,       0!=b in k(v),             (17)
```

and, for `n>=0`, define

```text
I_n(w)={h in A:h w^n in R}.                               (18)
```

For every monic irreducible `p in A`, write `nu_p` for its order and put

```text
J(w)=product_p p^max(0,a nu_p(F)-nu_p(b)).                (19)
```

Only finitely many exponents in `(19)` are nonzero.  Then

```text
I_0(w)=A,
I_n(w)=(F^(an):b^n)=J(w)^n                 for n>=1.       (20)
```

Here the colon is fractional notation for
`{h in A:h b^n in F^(an)A}`.  By `(9)`, membership is equivalent to

```text
nu_p(h)+n nu_p(b) >= an nu_p(F)             for every p.  (21)
```

The least nonnegative exponent of `p` in `h` is exactly `n` times `(19)`,
which proves `(20)`.  If `b=0`, separately `I_n=A`.

For THM-3397 take `F=(1+cv)^e`, `a=1`, and `b=v^(-m)`.  Formula `(19)` gives

```text
J=v^m(1+cv)^e,
I_n=(v^(mn)(1+cv)^(en)),                                  (22)
```

exactly its two-boundary filtration.  The power law `(20)` is special to
one homogeneous Laurent pole.  It is not asserted for a sum of different
Laurent degrees, whose powers mix coefficient convolutions and principal
parts.

## 5. Class group: the coarser quotient

Nagata localization at `X` applies because `(3)` is a UFD.  The only
height-one primes removed are the `D_i^X`, and

```text
(R_X)^* = k^* X^Z.                                        (23)
```

Thus `(8)` supplies the sole relation among their classes:

```text
Cl(R)=Z^r / Z(e_1,...,e_r)
     =Z^(r-1) direct_sum Z/d,
d=gcd(e_1,...,e_r).                                      (24)
```

Also `[D_i^Y]=-[D_i^X]`.  Formula `(24)` is the arbitrary-field version of
the class computation already proved over `C` in THM-2700.

The class group remembers the boundary lattice modulo the one principal
relation.  It cannot decide `(14)`: the divisor of every rational function is
already principal, while regularity asks whether its actual divisor is
effective.  Module `(10)` retains the lost signs and principal coefficients.

The factors in `(1)` must be irreducible, pairwise nonassociate, and
aggregated.  Writing `v^2=v*v` as two arms would falsely replace

```text
Cl(k[v,X,Y]/(XY-v^2))=Z/2                                (25)
```

by `Z^2/<(1,1)>=Z`; it would also invoke a false CRT split between two copies
of the non-comaximal ideal `(v)`.

## 6. Finite covers preserve principal parts

### 6.1 The same rational function cannot improve

Let `R subset B` be a finite injective inclusion of domains and put
`K=Frac(R) subset L=Frac(B)`.  Normality downstairs gives

```text
B intersect K=R.                                         (26)
```

Indeed an element of the intersection is integral over `R` and lies in `K`,
so it belongs to the integrally closed domain `R`.  Consequently pullback
induces injections

```text
R_X/R -> B_X/B -> L/B.                                   (27)
```

Thus the whole principal-part module, not merely nonregularity of one chosen
function, survives every finite dominant affine cover.  Equivalently, for
every such `B`,

```text
{h in A:X^(-q)h is regular in B}=(F^q).                  (28)
```

If `B` is normal and `Q` lies over a height-one prime `P`, the normalized
valuations satisfy

```text
ord_Q(z)=e(Q/P) ord_P(z)                    for z in K^*, (29)
```

so every negative order remains negative.

### 6.2 Weil conorm and the factorial-cover wall

Assume now that `B` is normal and let `n=[L:K]`.  The valuation-theoretic
Weil conorm is

```text
pi^*P=sum_(Q over P) e(Q/P) Q.                            (30)
```

No flatness is needed.  With

```text
pi_*Q=[kappa(Q):kappa(P)]P,                               (31)
```

the DVR degree identity and field norm give

```text
pi_* pi^* = n                 on Cl(R).                   (32)
```

Indeed pushforward takes a principal divisor upstairs to the divisor of its
field norm.  Hence

```text
ker(Cl(R)->Cl(B)) subset Cl(R)[n],
Cl(R)/torsion -> Cl(B)/torsion is injective.              (33)
```

Combining `(24)` and `(33)` yields the sharp classification

```text
R has a finite dominant normal factorial cover
iff r<=1.                                                  (34)
```

Necessity follows because a factorial target has zero class group, whereas
the free rank `r-1` cannot die.  For `r=0`, the identity cover works.  For
`r=1`, the explicit cover in the next subsection with `N=e_1` has class
group zero and is therefore factorial.

### 6.3 The canonical common-root cover

Suppose `r>=1`, let `N` divide `d` in `(24)`, and put

```text
d_i=e_i/N,               G=product_i f_i^(d_i),
S=k[v,x,y]/(xy-G).                                         (35)
```

There is a finite degree-`N` map

```text
R -> S,        X -> x^N,        Y -> c y^N.               (36)
```

Both `x` and `y` are integral over `R`, and
`Frac(S)=k(v,x)` has degree `N` over `k(v,x^N)=Frac(R)`.
The normality proof of Section 2 applies to `S`.

Pullback sends `D_i^X` to `(x,f_i)` with ramification index one.  Therefore
the class map is induced by the identity on `Z^r`:

```text
Z^r/Z(e_i) -> Z^r/Z(d_i),
ker = Z(d_i)/Z(Nd_i) = Z/N.                               (37)
```

For `N=d`, it kills exactly the torsion subgroup of `(24)` and leaves the
free part.  At the principal-part level, the old degree `q` packet maps
isomorphically to the cover's degree `Nq` packet:

```text
X^(-q) A/(F^q) -> x^(-Nq) A/(G^(Nq)),
G^(Nq)=c^(-q)F^q.                                         (38)
```

When `N>1`, the full cover module is larger: it also has negative degrees
not divisible by `N`.  For example, `F=v^2,N=2` sends the old degree-one
packet to cover degree two, while `x^(-1)A/(v)` is a new degree-one packet.

If `char(k)` does not divide `N` and `k` contains `mu_N`, `(36)` is the
ordinary cyclic quotient for
`(x,y)->(zeta x,zeta^(-1)y)` and is a torsor over the regular locus.
Without these hypotheses it remains the exact finite degree-`N` algebraic
cover above, but need not be an etale or cyclic Galois cover.

The mixed control

```text
F=v^2(v-1)^2,        N=2,        G=v(v-1)                (39)
```

has

```text
Cl(R)=Z direct_sum Z/2  ->  Cl(S)=Z.                     (40)
```

The torsion dies, the free arm difference survives, and every old
coefficient modulus `(F^q)` is unchanged in degree `2q`.

## 7. The nonmonomial JC sidecar and stopping boundary

The exact sidecar for a factorized terminal model is

```text
factorization and X/Y orientation
+ graded principal-part tower direct_sum_q A/(F^q)
+ its CRT arm jets and symmetry action
+ any additional base or infinity divisors.              (41)
```

It is strictly stronger than depth, residue degree, torsor class, or
`Cl(R)`.  For a single homogeneous decoded coordinate, `(19)--(20)` reduce
the tower to one exact degree-one debt polynomial.  For sums or orbit
averages, the residue classes in `(10)--(12)` decide tied cancellation as in
`(15)--(16)`.

There is also a sharp structural stopping obstruction.  Over `C`, every
irreducible `f_i` is linear.  If `F` is nonconstant and a factorized
Danielewski initial ring admits a finite dominant polynomial-plane cover, its
pullback class group is zero; `(33)--(34)` force `r=1`.  After an affine
change of `v`,

```text
F=c(v-alpha)^e,                                           (42)
```

which is precisely the one-factor pure-power geometry of THM-3383/3397.
Thus a genuinely multi-factor `F` cannot extend that finite polynomial-cover
mechanism.  It can enter a JC terminal analysis only if the proposed map is
nonfinite, merely local/punctured, nonnormal downstairs, or a differently
typed embedding.

To use `(41)` on an actual JC branch one must still prove that the physical
normalized initial ring is this `R`, that the decoded targets have the stated
Laurent form, that the map to the polynomial source is finite and dominant on
the whole affine chart, and that the factor arms match the actual divisors at
infinity.  None of those identifications follows from an abstract
Danielewski isomorphism, a killed torsor, or rational field decoding.
Accordingly this theorem closes no nonmonomial terminal family, `A4/S4`
branch, `JC(2)`, or `DC(2)` by itself.

## 8. Sharp hypothesis hostiles

- **Zero polynomial.**  If `F=0`, then `XY=0` is reducible and nonnormal;
  localization `R->R_X` kills `Y`, so `(9)--(10)` are not an injection.
- **Unit polynomial.**  If `F in k^*`, then `r=0`, `R=R_X` is a Laurent UFD,
  and the principal-part module is zero.  The gcd and cover-kernel formulas
  are void.
- **Duplicate factor.**  Equation `(25)` is the minimal witness that factors
  must be irreducible and aggregated.
- **Base field.**  Over `R`, `v^2+1` is one squarefree irreducible and gives
  `Cl=0`; after base change to `C` it splits into two factors and gives
  `Cl=Z`.  The factorization field must be named.
- **Downstairs normality.**  The finite normalization
  `k[t^2,t^3] subset k[t]` makes `t` regular although it was absent
  downstairs, so `(26)` fails without normality.
- **Finiteness.**  The localization `R->R_X` makes `X^(-1)` regular only by
  deleting every `D_i^X`; it is not finite.
- **Principal parts.**  Equations `(15)--(16)` show why valuation vectors of
  separate summands do not control cancellation.
- **Cover grading.**  The `F=v^2,N=2` packet after `(38)` forbids calling the
  entire cover principal-part module unchanged.

## 9. Exact companion

The standard-library companion uses only integer and rational polynomial
arithmetic.  It checks `372` factor/exponent/degree packets, `900` CRT
idempotents, `372` exact coefficient clearings, `900` one-factor-short
hostiles, `124` cancellation controls, `164` canonical covers with `225`
kernel representatives, and principal-part lengths through `36`.  A separate
bank checks `400` homogeneous-pole filtration cells, including `320` positive
powers and `636` one-step local hostiles.  It freezes `(25)`, `(39)--(40)`,
the nonclosed-base, zero-`F`, nonnormal, nonfinite, and cover-degree controls.

Reproduce with

```text
python3 04-computation/jc_factorized_danielewski_principal_parts_thm3404.py
python3 -O 04-computation/jc_factorized_danielewski_principal_parts_thm3404.py
```

The artifact hashes are pinned in the frontmatter.

**QED.**
