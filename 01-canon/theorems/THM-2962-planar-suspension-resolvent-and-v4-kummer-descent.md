---
id: THM-2962
title: "Planar suspension resolvent and V4 Kummer descent"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For a planar
  generically finite map with A4/S4 Galois closure,
  direct affine suspension preserves the Galois group and turns the full
  V4-resolvent normalization and its V4-cover into literal A1-products.
  On the regular locus, A1-invariance identifies H1_et(-,mu2)
  Q-equivariantly. Whenever the V4-cover is finite etale there—in
  particular in the Keller branch supplied by THM-2655—the canonical
  character plane, its unit/Cl[2] alternative, and its complete
  divisor-parity code are pulled back unchanged. THM-2465 therefore reduces
  the nonautomorphic affine
  point-cap degree-four Kummer gate exactly to the planar normal resolvent
  surface.  In the S4 branch W=Hom(V4,C2)=F2^2 has automorphism group S3;
  PSL2(Z)=C2*C3 reaches it only after the extra reflection relation, so a
  binary/ternary tree grammar does not retain free-word holonomy.
  Exact THM-2769 and THM-2871 hostiles show respectively that a genuine S4
  family can retain divisor word 110 and that the leading resolvent ranges
  over arbitrary depressed cubics.  Thus neither suspension nor shared
  discriminant supplies the missing Kummer carrier, Keller cubic owner,
  Jelonek square law, or source reconstruction.  No degree-four map is
  constructed or excluded; JC(2), G1, and the order-{1,3} conjecture remain
  OPEN.
source: codex-jc2-planar-suspension-resolvent-2026-07-29
depends_on:
  - THM-2465-g1-exclusion-package-for-degree-four-twojet-keller
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
  - THM-2685-equivariant-kummer-boundary-parity-completion-and-divisor-residue-gate
related:
  - THM-1310-conic-pair-fibers-and-design-equations
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2681-thm1310-s3-normalization-and-quartic-v4-torsor-exclusion
  - THM-2696-reflection-completed-s4-relative-different-and-coordinate-invariant-jacobian-gate
  - THM-2769-full-s4-pair-sum-affine-divisor-parity-hostile
  - THM-2871-quartic-leading-face-survivor-and-integral-depression-flag-dichotomy
script: 04-computation/jc2_planar_suspension_resolvent_thm2962.py
output: 05-knowledge/results/jc2_planar_suspension_resolvent_thm2962.out
script_sha256: ebceb860221eb60217b57fad5ab868a1cb00eb48dcfbc51a322291ac094ee097
output_sha256: 3053a6d107b47a61fe321422f09c5b3763e49d8d4fc1bc0779aa984707236113
hash_basis: working-tree bytes (LF)
---

# THM-2962 -- planar suspension resolvent and `V4` Kummer descent

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2465 proves that the nonautomorphic affine point-cap degree-four lane is
exactly a suspension of the planar degree-four problem.  THM-2655 and
THM-2685 identify the actual quartic obstruction: not a shared discriminant,
but a `Q`-standard plane in the Kummer cohomology of the **full Galois
resolvent normalization**, with every divisor residue zero.  This theorem
shows that suspension contributes no new coordinate to that obstruction.

## 1. The planar-suspension theorem

Let `k` be algebraically closed of characteristic zero and let

```text
g:A^2_k -> A^2_k
```

be a dominant generically finite polynomial map.  Write

```text
K=k(g_1,g_2) subset L=k(x,y) subset M                         (1)
```

for its target field, source field, and Galois closure.  Assume

```text
G=Gal(M/K)=A4 or S4,       V=V4,       E=M^V,
Q=G/V=C3 or S3.                                                (2)
```

Let `R` be the normalization of the planar target in `E`, and let `Z` be
its normalization in `M`.  Form the direct suspension

```text
g~=(g,t):A^3_k -> A^3_k.                                      (3)
```

Then:

```text
Galois closure(g~)=M(t),       M(t)^V=E(t);                    (4)

R~=R x A^1,                    Z~=Z x A^1;                     (5)

R~_reg=R_reg x A^1.                                               (6)
```

If `U=R_reg`, pullback is a natural `Q`-equivariant isomorphism

```text
H^1_et(U,mu_2)  ~=  H^1_et(U x A^1,mu_2).                     (7)
```

For an arbitrary `g`, `Z_U -> U` is a connected finite `V`-Galois cover.
If it is finite etale—in particular, when `g` is Keller and THM-2655
applies—then it is a connected `V`-torsor. Under (7), that torsor and all
three of its nonzero characters pull back to the corresponding torsor and
characters for the suspension. In the Keller case its canonical standard
plane

```text
W=Hom(V,C2) -> H^1_et(U,mu_2)                                 (8)
```

is present for `g~` if and only if it is present for `g`.

### Proof of (4)--(6)

The variable `t` is transcendental over `M`.  Every `K`-embedding of `M`
extends uniquely after fixing `t`; hence the normal closure of
`L(t)/K(t)` is `M(t)` with unchanged Galois group.  Taking `V`-fixed
fields commutes with the purely transcendental extension:

```text
M(t)^V=M^V(t)=E(t).                                           (9)
```

Let `A=k[g_1,g_2]` and let `B` be its integral closure in `E`.  The normal
ring `B[t]` is integral over `A[t]` and has fraction field `E(t)`.  If an
element of `E(t)` is integral over `A[t]`, the same monic equation makes it
integral over `B[t]`; normality then puts it in `B[t]`.  Thus `B[t]` is
the integral closure of `A[t]` in `E(t)`.  Applying the same argument to
the integral closure in `M` proves (5).

Over the perfect field `k`, a finite-type local ring is regular exactly
when its polynomial extension is regular at the corresponding primes:
regularity ascends along the smooth map and descends along its faithfully
flat localizations.  Hence the singular locus takes product with `A^1`,
which proves (6).

### Proof of (7) and the Kummer statement

The regular finite-type `k`-scheme `U` is smooth.  Étale homotopy
invariance for the prime-to-characteristic finite sheaf `mu_2` gives (7).
Equivalently, the Kummer sequence combines with the natural isomorphisms

```text
Gamma(U x A^1,O)^* = Gamma(U,O)^*,
Pic(U x A^1)       = Pic(U).                                (10)
```

All maps are natural under `Q`, so (7) is `Q`-equivariant. Equation (5)
makes `Z~ -> R~` the literal base change of `Z -> R`. Consequently,
whenever the restriction is etale, its character classes are precisely the
pullbacks of (8).

The etale hypothesis cannot be omitted from this character-plane sentence.
There are rational degree-four `P^1` covers with branch cycles

```text
(12)(34), (12), (34), (23), (23)
```

whose product is one and whose monodromy is `S4`. Writing such a rational
function as `P/Q`, the polynomial map

```text
(x,y) -> (yQ(x),yP(x))
```

is dominant and generically degree four, but the double-transposition branch
divisor has nontrivial inertia inside `V4`. Thus the induced `V4` cover
ramifies at a codimension-one regular point of the resolvent normalization.
Suspension still preserves this finite Galois cover and equations (4)--(7);
it does not turn the ramified cover into a torsor.

Since `R` is normal,

```text
Pic(U)=Cl(R),          Pic(U x A^1)=Cl(R[t]),                  (11)
```

so the unit-squareclass versus `Cl[2]` alternative of THM-2655 is also
unchanged, including its `Q`-module structure.

## 2. The divisor-parity code is unchanged

Let `u_1,u_2,u_3 in E^*` be the three squared pair sums of THM-2685.  They
represent the three nonzero characters of `W` and satisfy

```text
[u_1]+[u_2]+[u_3]=0.                                         (12)
```

For each height-one prime `D` of `R`, its pullback `D[t]` has

```text
v_(D[t])(u_i)=v_D(u_i).                                      (13)
```

Every other height-one prime of `R[t]` contracts to the generic point of
`R`; it is generated there by a nonconstant polynomial in `t`.  Because
`u_i` lies in the coefficient field `E`, its valuation at such a prime is
zero.  Therefore the parity matrix for the suspension consists of the
planar rows

```text
(v_D(u_1),v_D(u_2),v_D(u_3)) mod 2                           (14)
```

plus zero rows only.  In particular, its rank and its maximal unramified
character subspace are unchanged.

This proves more than discriminant invariance: the entire THM-2685
codimension-one Kummer gate survives suspension without loss or gain.

## 3. Exact point-cap consequence

THM-2465 proves, in its stated affine point-cap degree bounds, that every
nonautomorphic Keller map is polynomially source/target equivalent to a
direct suspension of a planar Keller pair.  Such automorphisms transport
the function-field tower and the normalizations isomorphically.  Sections
1--2 therefore give:

> **Point-cap descent.** The degree-four affine point-cap Kummer obstruction
> is already a normal-surface obstruction on the planar full Galois
> resolvent normalization.  The third affine coordinate supplies no unit
> squareclass, `Cl[2]` class, divisor residue, or Kummer character capable
> of repairing a missing planar standard plane.

This is a dimension reduction, not a degree-four exclusion.  Planar
degree-four Keller existence, hence this point-cap lane, remains open.

## 4. The exact `2`/`3` co-occurrence and its loss

The recurring binary/ternary grammar has a precise representation-theoretic
core:

```text
W=Hom(V4,C2)=F_2^2,
W\{0} has three elements,
Aut(W)=GL_2(F_2)=S3.                                        (15)
```

An order-three element cycles the three nonzero Kummer characters; an
order-two element reflects two of them.  In the `S4` branch this is exactly
the quotient action `S4/V4=S3`.  In the `A4` branch only the cyclic `C3`
subaction occurs.

After choosing an order-two reflection `s` and an order-three cycle `r`,
the standard presentation gives an epimorphism

```text
PSL_2(Z)=C2*C3 -> S3                                      (16)
```

obtained by imposing the additional relation

```text
(sr)^2=1.                                                  (17)
```

Thus a binary/ternary or Farey-tree model retains the orders of the two
generators but loses free-word holonomy when it passes to the actual
quartic quotient.  More importantly, (15) is only an abstract
representation.  The Keller problem still requires its embedding into the
geometric group (8) with zero residue rows.  No tournament orientation is
canonical here: the order-two quotient element reverses the three-character
cycle.

## 5. Why the grade-three anatomy does not transfer

The exact source/target contract is:

| datum | quartic resolvent operation |
|---|---|
| source | quartic Galois closure `M/K` with `V normal in G` |
| target | full Galois resolvent field `E=M^V` and normalization `R` |
| map | fixed-field quotient followed by normalization |
| preserved | branch-discriminant support up to index/different squares; the standard `V4` character representation |
| destroyed | chosen quartic sheet, affine source coordinates, Keller Jacobian, and a distinguished cubic root |
| needed sidecar | geometric Kummer embedding with zero residues, an actual Keller cubic owner, Jelonek square law, and unit source reconstruction |

THM-2598 and THM-2696 make the loss load-bearing.  In the `S4` case the
quartic root field and one non-Galois cubic-resolvent root field are
incomparable, while the reflection quotient has an intrinsic different.
Discriminant equality therefore does not identify either field with a
degree-three Keller source.

## 6. Two exact stopping hostiles

The companion verifies both boundaries with exact arithmetic.

### 6.1 A genuine `S4` packet can retain parity word `110`

For THM-2769's family

```text
f_t(Y)=Y^4-2Y^2-8tY+1-4t,
S_t(U)=U^3-4U^2+16tU-64t^2,                               (18)
```

the two discriminants are identically

```text
-4096 t^2(27t^2-14t+3).                                   (19)
```

At `t=0`, the coefficient Newton points of `S_t` are

```text
(0,2),(1,1),(2,0),(3,0),                                  (20)
```

so its root valuations are `(1,1,0)` and its divisor-parity word is
`110`.  This is a nonzero THM-2685 row.

The specialization `t=2` is irreducible mod `5`, has factorization type
`1+3` mod `3`, and has negative discriminant.  Its transitive Galois group
contains a four-cycle and a three-cycle and is not contained in `A4`;
hence it is `S4`.  The hostile is therefore not a degenerate abstract
quartic.  By Section 2, suspension preserves its `110` row.

### 6.2 The leading cubic is arbitrary

THM-2871's leading face gives, for

```text
h(X)=X^3+pX+q,
S_0(U)=U^3+16pU-64q,                                      (21)
```

the identity

```text
disc(S_0)=4096 disc(h).                                    (22)
```

As `(p,q)` vary, (21) ranges over every depressed cubic after invertible
coefficient rescaling.  Therefore (22) alone cannot imply THM-1310's
special square/Jelonek law.  The missing grade-three structure is an actual
Keller owner and reconstruction map, not another discriminant identity.

## 7. Residue-first next test

The cheapest correctly typed test for a bounded planar degree-four Keller
ansatz is:

1. compute the **full** `A4`/`S4` resolvent normalization, of degree `3` or
   `6`; do not replace the `S4` object by one cubic root field;
2. run THM-2685 first: at every prime divisor compute (14) and reject any
   nonzero word among `110,101,011`;
3. only if all rows vanish, compute the `Q`-module structures of
   `O(R)^*/O(R)^{*2}` and `Cl(R)[2]`, rejecting the ansatz if neither
   contains the standard plane;
4. only for survivors that enter THM-2871's simple `A`-leading /
   `A=B=0` flag chart, test that theorem's integral depression, square law,
   Jelonek ownership, and source-coordinate unit reconstruction.

The test is planar by Sections 1--3.  It neither constructs nor excludes a
degree-four Keller map, and it does not reopen the already closed inherited
degree-eighteen or degree-twenty-two branches.

## 8. Exact evidence

Reproduce with

```text
python 04-computation/jc2_planar_suspension_resolvent_thm2962.py
python -O 04-computation/jc2_planar_suspension_resolvent_thm2962.py
```

Both executions are byte-identical to
`05-knowledge/results/jc2_planar_suspension_resolvent_thm2962.out`.
The companion checks (18)--(22), the `S4` specialization, the nonzero
parity word, `|GL_2(F_2)|=6`, and the exact order-`2`, order-`3`,
order-`2` quotient relations in (16)--(17).  It uses explicit runtime
requirements rather than Python `assert`, so optimized execution retains
every truth-bearing check.

An independent hostile audit caught and repaired the distinction between a
finite `V4`-Galois cover and an etale `V4`-torsor, verified the inserted
ramified `S4` boundary by branch cycles and Riemann--Hurwitz, and restricted
the THM-2871 follow-up to its actual flag chart. It then rederived the
remaining base-change, Kummer, divisor-parity, and quotient claims and
independently reproduced the normal, optimized, and stored transcripts and
declared hashes.
