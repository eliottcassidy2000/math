---
id: THM-2374
title: "Binary-allocation complete-subcube Dirichlet spectrum"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. For an
  exactly two-block prime-token allocation, averaging squared energy
  differences over every complete Boolean coordinate subcube gives an
  exact 0/2 Walsh multiplier. The resulting Dirichlet bank recovers
  every squared ANOVA component. For two prime types P,Q, the defect
  D_P+D_Q-D_(P union Q) is nonnegative and equals twice the total ANOVA
  energy meeting both types, or one half of the mean squared
  conditional-rectangle difference. Hence all pair defects vanish
  exactly when THM-2348's robust prime-type rectangularity holds. The fixed-optimum
  hostile [[0,M],[M,M]] has exact defect M^2/8. The theorem requires the
  full two-block cube and labelled-uniform measure, loses ANOVA signs
  and owners, and proves no new Gordian distance or knot realization.
source: codex-2026-07-25-knot-complete-subcube-dirichlet
depends_on:
  - THM-2346-global-allocation-anova-normal-form
  - THM-2348-prime-type-rectangularity-and-target-token-conditioning
related:
  - THM-2248-higher-interaction-defect-complex-and-tropical-trace-spectrum
  - THM-2369-complete-line-target-dirichlet-and-balanced-observable-no-go
script: 04-computation/binary_allocation_subcube_dirichlet_thm2374.py
output: 05-knowledge/results/binary_allocation_subcube_dirichlet_thm2374.out
script_sha256: 225b8d386fe3536365c3dcefd8e1da70435f98a921f9bc3305dfd129ec55646f
output_sha256: b6dcc575992b6a64adfaa4a8a0f8d261d9c20329d6978637857ca75b29e995d5
hash_basis: working-tree bytes (LF)
---

# THM-2374 -- complete Boolean subcubes measure every mixed allocation mode

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2346 gives the unique labelled-token ANOVA expansion, and THM-2348
characterizes robust factorization by the vanishing of every component
meeting two prime types. When the packet has exactly two source blocks,
that qualitative condition has an exact nonnegative Dirichlet
certificate.

## 1. The genuine two-block cube

Let the labelled target-token set be partitioned by oriented prime type:

```text
S=disjoint_union_(P in Pcal) S_P.                  (1)
```

Assume the current packet has exactly two source blocks

```text
C={a,b}.
```

Choose either identification `C=F_2`. The allocation energy is then

```text
E:F_2^S->R.
```

All expectations below use uniform measure on the **labelled** cube.
THM-2346's ANOVA expansion becomes the normalized Walsh expansion

```text
E=sum_(T subset S)H_T,

H_T=Ehat(T)chi_T.                                  (2)
```

The summands are orthogonal and

```text
||H_T||_2^2=|Ehat(T)|^2.                           (3)
```

Changing the identification by swapping `a` and `b` translates the
whole cube. Every squared difference below is unchanged.

## 2. Complete-subcube Dirichlet energy

For `U subset S`, let `F_2^U` denote the subgroup of shifts supported on
`U`, including the zero shift. Define

```text
D_U(E)
 =E_(x in F_2^S, h in F_2^U)
   |E(x+h)-E(x)|^2.                                (4)
```

Let `A_U` be averaging over that subgroup:

```text
(A_U E)(x)=E_(h in F_2^U)E(x+h).
```

The exact complete-subcube identity is

```text
D_U(E)
 =2 sum_(T:T intersect U!=empty)||H_T||_2^2
 =2||E-A_U E||_2^2.                               (5)
```

Indeed, the average squared multiplier of `chi_T` is zero if
`T intersect U` is empty and two otherwise. Equivalently `A_U` is the
orthogonal projection retaining precisely the characters disjoint from
`U`.

The zero shift in (4) is part of the normalized subgroup average. Removing
it changes the constants.

## 3. Pairwise prime-type defects are squared rectangles

For disjoint coordinate sets `U,V`, put

```text
J_(U,V)
 =D_U+D_V-D_(U union V).                           (6)
```

Applying (5) frequency by frequency gives

```text
J_(U,V)
 =2 sum_(
    T:T intersect U!=empty,
      T intersect V!=empty
   ) ||H_T||_2^2
 >=0.                                              (7)
```

If

```text
Delta_h E(x)=E(x+h)-E(x),
```

then the same character calculation gives

```text
J_(U,V)
 =1/2 E_(x,h in F_2^U,k in F_2^V)
       |Delta_h Delta_k E(x)|^2.                   (8)
```

Thus `J_(U,V)=0` exactly when every conditional rectangle using one
direction from `U` and one from `V` vanishes.

Taking `U=S_P` and `V=S_Q`,

```text
J_(P,Q):=J_(S_P,S_Q)
 =2 sum_(
    T meets S_P and S_Q
   )||H_T||_2^2.                                   (9)
```

Consequently

```text
J_(P,Q)=0 for every distinct P,Q

iff

H_T=0 for every mixed-prime-type T.                (10)
```

Under THM-2348's repeated-type symmetry, (10) is also equivalent to
additive type separation, every conditional rectangle on the equal-prime
quotient, and perturbation-robust Cartesian tied optima.

## 4. Exact support-energy inversion

The complete family `(D_U)_(U subset S)` determines the squared norm of
every nonconstant Walsh/ANOVA component. For nonempty `T`,

```text
||H_T||_2^2
 =-1/2 sum_(J subset T)
   (-1)^(|T|-|J|)D_(S minus J).                    (11)
```

This is Boolean Mobius inversion applied to (5). The constant component
and every sign or phase are lost.

There are two useful aggregate formulas. Put

```text
r(T)=number of prime types P met by T.
```

Summing pair defects gives the multiplicity-weighted spectrum

```text
sum_(P<Q)J_(P,Q)
 =2 sum_T binom(r(T),2)||H_T||_2^2.                (12)
```

For three or more types this is not twice the unweighted mixed energy.
The latter is instead

```text
sum_(T:r(T)>=2)||H_T||_2^2
 =1/2[
    sum_P D_(S minus S_P)
    -( |Pcal|-1 )D_S
   ].                                              (13)
```

Equations (11)--(13) make the qualitative THM-2348 obstruction a full
squared-support atlas.

## 5. The sharp two-by-two hostile

THM-2348's accidental-product-optimum hostile is

```text
E_M=
 [0 M
  M M].                                           (14)
```

Its normalized Walsh coefficients are

```text
3M/4, -M/4, -M/4, -M/4.
```

For the two token coordinates `p,q`,

```text
D_p=D_q=M^2/4,

D_(p,q)=3M^2/8,

J_(p,q)=M^2/8.                                     (15)
```

Thus a Cartesian singleton optimum does not imply robust
rectangularity. More generally, a real `2 x 2` table with rectangle
contrast

```text
Delta=E_00-E_01-E_10+E_11
```

has

```text
J=Delta^2/8.                                       (16)
```

This is not a universal lower bound in terms of table range unless that
range is the displayed contrast.

## 6. Scope boundaries

The theorem is global for a packet with exactly two source blocks. A
chosen binary face inside a packet with at least three blocks is not
decisive. For example,

```text
E(x_p,x_q)=1_(x_p=c,x_q=c)
```

vanishes on the `{a,b}^2` face but has a nonzero `a/c` rectangle.
Selecting two colours from a larger packet is not the same as globally
relabelling the two actual blocks.

The measure is uniform on labelled assignments. On the equal-prime
quotient this becomes the induced binomial measure, not uniform measure
on quotient count cells.

The Dirichlet bank preserves:

- the squared norm and type support of every nonconstant ANOVA component;
- every mixed-type vanishing decision; and
- the exact quantitative rectangle defect.

It loses:

- the sign or phase of every component;
- the constant term;
- the minimizing owner set; and
- the knot-theoretic realization of an abstract energy table.

Accordingly (10) is a new exact diagnostic and repackaging of
THM-2346/2348. It proves no new Gordian distance, owner identity,
connected-sum law, or knot realization.

## 7. Exact companion

The dependency-free `Fraction` companion:

- checks all `16,864` complete-subcube spectral identities on `1,054`
  distinct exact rational four-token difference profiles;
- checks `3,162` pairwise prime-type defects;
- checks `15,810` squared-support inversion cells;
- verifies (13) on every table;
- verifies the exact hostile constants at `M=7/3`;
- checks `373` distinct prescribed rectangle contrasts; and
- checks the global-swap and larger-packet-face boundaries.

Run

```bash
python3 04-computation/binary_allocation_subcube_dirichlet_thm2374.py
python3 -O 04-computation/binary_allocation_subcube_dirichlet_thm2374.py
```

Both transcripts must match

```text
05-knowledge/results/binary_allocation_subcube_dirichlet_thm2374.out
```

byte-for-byte after LF normalization. Every executable check raises
explicitly under optimized Python.

Independent hostile audit rederived (5)--(13), all hostile constants,
the labelled-measure boundary, and the exact coverage claims; normal,
optimized, and stored transcripts agree. QED.
