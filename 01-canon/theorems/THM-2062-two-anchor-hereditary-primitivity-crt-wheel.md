---
id: THM-2062
title: "Two-anchor hereditary-primitivity determinantal sieve and CRT wheel"
status: >
  PROVED. For nonzero coefficient rows generating Z^2, every rational-rank-two
  deletion has determinantal index I_l, and its specialized gcd divides I_l.
  Modulo each prime the bad primitive parameters form one projective
  annihilator per bad deletion, with at most two distinct directions across
  all deletions. Hence every fixed-N THM-2058 coprime interval carries an
  exact squarefree CRT wheel, including zero local factors. A rank-one
  deletion instead forces one primitive linear form to equal plus or minus
  one and leaves an explicit one-dimensional coprimality wheel.
source: codex-2026-07-21-LRC-two-anchor-CRT-wheel
script: 04-computation/lrc_two_anchor_hereditary_crt_wheel_codex_20260721.py
result: 05-knowledge/results/lrc_two_anchor_hereditary_crt_wheel_codex_20260721.out
script_sha256: e46757c4acc6508ff199240d07f18e185435549967cbc909318d621dbdd702b5
result_sha256: 16ae242af27a18ccead4865edeb56b8e67113e5ad239dd24e612714a0f1d7b4f
hash_basis: normalized repository blobs (LF)
depends_on:
  - THM-2053
  - THM-2058
  - THM-765
related:
  - THM-2052
  - THM-2060
  - THM-2061
  - HYP-8846
  - HYP-2764
---

# THM-2062 -- Hereditary primitivity is a two-anchor CRT wheel

Let `r_1,...,r_n` (`n>=3`) be nonzero rows in `Z^2` which generate `Z^2`.
Let

```text
d=(N,M) in Hom(Z^2,Z),       gcd(N,M)=1,
v_i=d(r_i)=r_i.d.                                        (1)
```

For a deletion `ell`, put

```text
L_ell=sum_(i!=ell) Z r_i,
g_ell(d)=gcd_(i!=ell)|v_i|.                              (2)
```

The specialized row is **hereditarily primitive** exactly when
`g_ell(d)=1` for every `ell`.

This is the coefficient-column language of a saturated THM-2052/2058
two-anchor star. Positivity, distinctness, owner sectors, and phase height are
separate filters; the present theorem computes only hereditary primitivity.

## A. Determinantal divisibility and the local bad direction

A deletion has rational rank zero only if all its rows vanish, impossible
because the rows are nonzero and `n>=3`. If `rank_Q(L_ell)=2`, define

```text
I_ell=[Z^2:L_ell]
     =gcd_(i<j; i,j!=ell)|det(r_i,r_j)| >0.               (3)
```

Then

```text
g_ell(d) divides I_ell.                                  (4)
```

For a prime `p`, the following are equivalent:

```text
p|g_ell(d),
d mod p annihilates L_ell mod p,
[d mod p] is the unique projective annihilator K_(ell,p). (5)
```

Such a direction exists exactly when `p|I_ell`.

### Proof

Reduction of the primitive covector `d` modulo `g_ell` is a surjection
`Z^2 -> Z/g_ell Z` and vanishes on `L_ell`. Hence `Z^2/L_ell` surjects onto
`Z/g_ell Z`, proving (4). Equivalently, the elementary identities

```text
m_j v_i-m_i v_j=det(r_i,r_j)N,
a_i v_j-a_j v_i=det(r_i,r_j)M                            (6)
```

for `r_i=(a_i,m_i)`, together with `gcd(N,M)=1`, give the same result.

The condition `p|I_ell` says that the deletion rows have rank at most one
modulo `p`. They cannot all vanish: then the one omitted row would leave the
full row family with rank at most one modulo `p`, contradicting generation of
`Z^2`. Thus their rank is exactly one and their annihilator in the dual
projective line is unique. This proves (5). QED.

## B. At most two forbidden projective directions per prime

For a fixed prime `p`, let

```text
B_p={K_(ell,p):p|I_ell and rank_Q(L_ell)=2}.              (7)
```

Then

```text
|B_p|<=2.                                                 (8)
```

If equality holds, there are unique indices `ell,j` such that every common
row `r_i`, `i!=ell,j`, is zero modulo `p`, while `r_ell,r_j` are nonzero and
independent modulo `p`. The two forbidden directions are their two
annihilators.

### Proof

Let `W_ell` be the one-dimensional row span of a bad deletion. Two different
bad deletions cannot have the same span: `W_ell=W_j` would contain `r_j`,
`r_ell`, and all common rows, making the whole family rank one. If their spans
are different, every common row lies in their zero intersection. The full
family can span `F_p^2` only when the two omitted rows are nonzero and
independent. Any third deletion retains both and therefore has rank two.
This proves (8) and the equality description. QED.

There is a useful global guardrail. Put `I=product I_ell`, assume all
deletions have rational rank two, and let `Q=rad(I)`. The number of unimodular
vector residues modulo `Q` (nonzero modulo every `p|Q`) avoiding every bad
direction is

```text
product_(p|Q)(p-1)(p+1-|B_p|),                            (8a)
```

and the number of allowed projective classes is
`product_(p|Q)(p+1-|B_p|)`. Relative to primitive directions the local density
is

```text
product_(p|Q)(1-|B_p|/(p+1))>0.                           (8b)
```

Indeed, each projective class has `p-1` nonzero representatives. Ordinary
CRT gives (8a). More precisely, if `Omega` is a bounded Jordan-measurable set
inside a parameter cone, Mobius inversion gives the number of good primitive
points in `T Omega` as

```text
area(Omega)T^2/zeta(2)
  * product_(p|Q)(p+1-|B_p|)/(p+1) + o(T^2).              (8c)
```

Thus hereditary primitivity has positive density globally and cannot
terminate a whole rank-two cone by itself.

## C. The exact fixed-`N` CRT wheel

Assume first that every deletion has rational rank two, and fix `N>0`. Put

```text
I=product_ell I_ell,
Q_0=rad(NI).                                               (9)
```

For each prime `p|Q_0`, define the local number of admissible longitudinal
residues

```text
kappa_p(N)=#{m mod p:
  (p does not divide N or m!=0)
  and [N:m] not in B_p}.                                  (10)
```

Here `B_p` is empty when `p` does not divide `I`. More explicitly:

- if `p` does not divide `N`, let `nu_p(N)` be the number of directions in
  `B_p` whose first coordinate is nonzero; then

  ```text
  kappa_p(N)=p-nu_p(N),       0<=nu_p(N)<=2;               (11)
  ```

- if `p|N`, then every primitive residue has `m!=0` and projective direction
  `[0:1]`, so

  ```text
  kappa_p(N)=0       if [0:1] in B_p,
             p-1    otherwise.                            (12)
  ```

Then the residue classes `M mod Q_0` for which `(N,M)` is primitive and the
specialized row is hereditarily primitive form the exact CRT product of the
local allowed sets. Their number is

```text
product_(p|Q_0) kappa_p(N).                               (13)
```

Zero factors are possible in either branch. In particular, when `p=2` does
not divide `N`, two forbidden affine directions can consume both residues;
one must not describe fiber death as a `p|N` phenomenon only.

More exactly, the fixed-`N` wheel is empty if and only if either

1. some `p|N` has `[0:1] in B_p`; or
2. `N` is odd and `B_2` contains both affine directions `[1:0],[1:1]`.

For every odd `p` not dividing `N`, (11) is at least `p-2>=1`. By Part B, the
second case has a rigid form modulo `2`: two omitted rows are independent and
all other rows vanish. Thus local fiber death is exactly the parameter-space
shadow of a twelve-of-thirteen common-divisibility obstruction.

### Proof

The condition `gcd(N,M)=1` excludes exactly `M=0 mod p` for primes `p|N`.
By Parts A--B, a prime divides one deletion gcd exactly when `[N:M]` belongs
to `B_p`; every possible prime lies in `I`. Thus the local conditions (10)
are necessary and sufficient, and they involve only `M mod p`. The ordinary
Chinese remainder theorem proves the product description and (13). QED.

For an integer interval `[A,B]`, if `R_N subset Z/Q_0Z` is the allowed CRT
set, its exact population is the elementary floor sum

```text
sum_(r in R_N)(
  floor((B-r)/Q_0)-floor((A-1-r)/Q_0)
).                                                        (14)
```

Thus THM-2058's coprime longitudinal interval can be sieved symbolically
before its explicit collision walls are removed. The local factors also give
an exact density and a finite Jacobsthal-gap problem, rather than a raw scan
of every longitudinal integer.

## D. Rank-one deletion is an affine terminal plus a 1D wheel

Suppose `rank_Q(L_ell)=1`. Write

```text
r_i=lambda_i c       (i!=ell),                            (15)
```

where `c` is a primitive row. Saturation of the full row family forces

```text
gcd_(i!=ell)|lambda_i|=1,
|det(c,r_ell)|=1.                                        (16)
```

Therefore

```text
g_ell(d)=|d(c)|,                                         (17)
```

and this deletion is primitive exactly on the two affine parameter lines

```text
d(c)=+1       or       d(c)=-1.                           (18)
```

On either line, use the unimodular coordinate

```text
t=d(r_ell),
h_j=gcd_(i!=ell,j)|lambda_i|       (j!=ell).               (19)
```

The remaining deletion gcds are exactly

```text
g_j(d)=gcd(|t|,h_j),                                     (20)
```

so full hereditary primitivity is the one-dimensional wheel

```text
gcd(t,product_(j!=ell)h_j)=1.                             (21)
```

There is at most one rational-rank-one deletion when all rows are nonzero.

### Proof

Every full-rank minor uses `r_ell` and equals
`lambda_i det(c,r_ell)`. Their gcd is one because the rows generate `Z^2`,
which proves (16). Equations (17) and (20) follow by taking the displayed
gcds after `d(c)=+-1`. The equivalence (21) is prime by prime. If two
deletions had rational rank one, their `n-2` common nonzero rows would lie in
the intersection of two distinct rational lines; equal lines would make the
whole family rank one. Both alternatives are impossible. QED.

The condition (18) alone is not full heredity. For example, the rows

```text
(2,0),(3,0),(4,0),(0,1)
```

generate `Z^2`. Deleting `(0,1)` and taking `d=(1,2)` gives `|d(1,0)|=1`,
but deleting `(3,0)` leaves specialized values `2,4,2` of gcd `2`. This is
exactly the residual wheel (21).

## E. Atlas effect and scope

THM-2053/2058 reduce a surviving two-anchor plane to finitely many bad
transverse clocks and explicit coprime intervals in `M`. Parts C--D insert a
lossless hereditary-primitivity layer:

```text
rank-two deletions -> squarefree CRT wheel, <=2 bad directions/prime;
rank-one deletion  -> two affine lines + one-dimensional coprime wheel.
```

This is genuine atlas compression: the wheel depends only on the coefficient
template and fixed `N`, not on enumerating specialized speed gcds one by one.
It does not itself prove phase safety. The surviving residue classes must
still meet THM-2058's phase packets, THM-2055 owner sector, pair-sum exits,
and collision walls.

The ingredients existed separately—THM-765's hereditary-primitivity
necessity, THM-2053's saturated plane coordinates, and THM-2058's coprime
intervals—but no earlier canonical result states the deletion indices,
projective annihilators, two-direction bound, or exact CRT product.

## F. Arithmetic-matroid and higher-anchor extension

The mechanism is not a rank-two coordinate accident. For nonzero rows
spanning `Z^r`, a rational-rank-`r` deletion has index equal to the gcd of its
maximal minors, and the specialized deletion gcd divides that index. If a
prime divides the index, the deletion has rank exactly `r-1` modulo `p`:
adding its one omitted row must recover rank `r`. Its bad primitive covectors
therefore form one projective annihilator point.

These bad deletions are precisely the coloops of the rank-`r` row matroid
after reduction modulo `p`, so there are at most `r` of them. If their number
is `r_p`, the local good unimodular-vector and projective counts are

```text
p^r-1-r_p(p-1),
(p^r-1)/(p-1)-r_p,                                       (22)
```

respectively. The relative primitive-density factor is

```text
1-r_p/(1+p+...+p^(r-1)).                                 (23)
```

Thus the two-anchor wheel is the `r=2` member of a general deletion-index
wheel. HYP-2764's coloopless/cosimple zonotope lane anticipated the correct
matroid vocabulary; the new content here is the arithmetic multiplicity,
specialized speed gcd, and fixed-coordinate CRT interval.

## Assumption challenge and tournament analysis

The vertices are deletion obligations, not runners or parameter points. At a
prime, orient two obligations by their omitted-row index, using projective
annihilator as the gauge and index order as the tie Hamiltonian path. The
incidence structure has at most two distinct vertices after quotient; if two
survive, their omitted rows are the only nonzero independent rows modulo the
prime. Hence every such tournament is transitive, with score histogram `(0)`
or `(0,1)`, no directed cycle, singleton SCCs, and one Hamiltonian path.

That tiny tournament is only a fingerprint. It destroys the actual forbidden
projective directions and the distinction between accessible affine classes
and `[0:1]` on a fixed-`N` fiber. The faithful carrier is the prime-labelled
CRT wheel with those directions retained.

## Computational audit

The companion uses integer arithmetic. It exhausts the first `5,000`
saturated four-row templates from `[-2,2]^2`, checks determinantal divisibility,
local projective kernels, the two-direction equality pattern, and direct
fixed-`N` residue counts against (13). It separately checks every rank-one
template and the guardrail above. Runtime checks survive `python -O`; the
frozen output ends in `PASS`.
