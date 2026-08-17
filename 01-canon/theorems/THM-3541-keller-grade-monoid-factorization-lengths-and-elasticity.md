---
id: THM-3541
title: "Keller grade-monoid factorization lengths and elasticity"
status: >
  PROVED + VERIFIED-EXACT.  The multiplicative Keller degree monoid
  M={1} union {n>=3} has atoms exactly the odd primes, 4, 8, and twice an
  odd prime.  If N=2^a m with m odd and b=Omega(m), its numerical-atom
  factorization lengths form the complete integer interval from
  b+ceil(max(a-b,0)/3) to b+floor(a/2); the monoid elasticity is 3/2.
  Grade 24 is the first with nonunique numerical atom factorizations and 36
  the first with different lengths.  In contrast, THM-3438 supplies a
  composition-atomic Keller map in every grade N>=3.  The theorem classifies
  numerical grades and realizable expression lengths across maps, not
  factorization uniqueness of any fixed map.
source: codex/turnpike-atlas/2026-08-16
depends_on:
  - THM-1330-keller-monoid-exact-picture-inverse-jelonek-cusp-rule
  - THM-3438-weighted-lift-keller-degree-spectrum
related:
  - HYP-9030-keller-degree-semigroup
script: 04-computation/keller_grade_monoid_factorization_lengths_thm3541_20260816.py
output: 05-knowledge/results/keller_grade_monoid_factorization_lengths_thm3541_20260816.out
script_sha256: 11541108facd297a55287b1657e58eac777ef97eae36a6d0cfb7c84e5024c414
output_sha256: 46ec02241ded0a8cba53801cc1b521f9dd9802db997356e5d802114965fc6394
semantic_sha256: c31c4c64207547cf2b239fd59be8fc53925f75e84f4b6123772ef098a49432a3
hash_basis: LF-normalized bytes
---

# THM-3541 -- exact arithmetic of the Keller grade monoid

**PROVED + VERIFIED-EXACT.**

Fix any ambient dimension `r>=3`.  THM-3438 proves that the generic-degree
values of Keller endomorphisms are

```text
M=KDeg(r)={1} union {n:n>=3}.                           (1)
```

Composition makes `(1)` a multiplicative monoid.  This theorem determines
its irreducibles and complete sets of factorization lengths, then states the
exact consequence for composition expressions of Keller maps.

The central distinction is:

```text
atom of the numerical grade monoid M
    !=
composition-atomic Keller map.                         (2)
```

The two notions agree only in the forced direction.

## 1. The numerical atoms

Call `N in M\{1}` a **grade atom** if it cannot be written `N=uv` with
`u,v in M\{1}`.  Then

```text
Atom(M)
 ={odd primes} union {4,8} union {2p:p an odd prime}.  (3)
```

Proof.  An odd `N` is reducible in `M` exactly when it is composite.  If
`N=2m` with `m` odd, a factorization in `M` exists exactly when `m` is
composite: write `m=ab` and `N=a(2b)`; conversely a factorization of `2m`
must allocate the sole factor two to one side and factor `m`.  Thus the atoms
in this class are `2p`.

Finally, `4` and `8` have no factorization with both factors at least three.
Every larger multiple of four is

```text
N=4(N/4),                 N>=12,                       (4)
```

with both factors in `M\{1}`.  This proves `(3)`.

THM-3438's forced-atomic grade list is therefore exactly `Atom(M)`: if a
Keller map in one of these grades decomposed, degree multiplicativity would
factor its grade in `M`.  The converse fails maximally.  THM-3438 constructs
an `S_N` composition atom in **every** grade `N>=3`, including every reducible
grade of `M`.

## 2. Parameterization of every numerical atom factorization

Write

```text
N=2^a m,             m odd,          b=Omega(m),       (5)
```

where `Omega` counts odd prime factors with multiplicity.  In a factorization
of `N` into the atoms `(3)`, let

```text
x = number of atoms 2p,
y = number of atoms 4,
z = number of atoms 8.                                 (6)
```

The `x` atoms of type `2p` consume `x` powers of two and `x` of the `b` odd
prime factors.  Every remaining odd prime is itself an atom.  Hence the
complete constraints and length are

```text
0<=x<=min(a,b),
2y+3z=a-x,
ell=x+y+z+(b-x)=b+y+z.                                (7)
```

Conversely, every nonnegative solution of `(7)`, together with any allocation
of the chosen odd primes to the `2p` atoms, gives a factorization.  Thus `(7)`
is an exact parameterization, not only a bound.

## 3. Every set of lengths is an interval

Put `k=ell-b`.  Given a proposed `k`, equations `(7)` can be rewritten

```text
y=k-z,               x=a-2k-z.                        (8)
```

There is an integer `z` making `x,y,z` lawful exactly when

```text
max(0,a-2k-b) <= z <= min(k,a-2k).                     (9)
```

The interval `(9)` is nonempty exactly for

```text
ceil(max(a-b,0)/3) <= k <= floor(a/2).                 (10)
```

Indeed, the upper endpoint is nonnegative precisely when `2k<=a`, while the
only additional lower-versus-upper condition is `a-b<=3k`.  The exceptional
formal input `(a,b)=(1,0)` would be `N=2`, which is not in `M`; every other
endpoint construction is lawful.

Therefore the complete factorization-length set is

```text
L_M(N)
 ={b+ceil(max(a-b,0)/3), ... , b+floor(a/2)}.          (11)
```

In particular, there are no gaps.

## 4. Elasticity and first nonunique grades

Let `ell_min,ell_max` be the endpoints of `(11)`.  If `a<=b`, then

```text
ell_max <= b+a/2 <= 3b/2 =3ell_min/2.                  (12)
```

If `a>b`, put `d=a-b`.  Then

```text
ell_max <= 3b/2+d/2,
3ell_min/2 >= 3b/2+d/2.                               (13)
```

Thus

```text
ell_max/ell_min <= 3/2.                                (14)
```

At `N=36=2^2*3^2`, the two factorizations

```text
36=6*6=4*3*3                                           (15)
```

have lengths two and three, so the bound is attained.  The elasticity of
`M` is exactly `3/2`.

The first three boundaries are

```text
first reducible grade:                    9=3*3;
first nonunique atom factorization:      24=3*8=4*6;
first different factorization lengths:  36=6*6=4*3*3. (16)
```

The first statement is immediate from `(3)`.  Below `24`, the reducible
grades are `9,12,15,16,18,20,21`, each with its unique unordered atom
factorization.  Formula `(11)` gives a singleton length set below `36`, and
`(15)` realizes its first nonsingleton.  The exact companion independently
enumerates all unordered factorizations for these boundary checks.

## 5. Exact Keller-map expression-length spectrum across a grade

For a grade `N`, define `ExprLen_r(N)` to be the set of integers `k` for
which **some** Keller endomorphism of `A^r` of degree `N` is displayed as a
composition of `k` composition atoms.  This varies the map inside the grade;
it is not the set of factorizations of one fixed map.

Every nonunit factor has degree at least three, so any such expression gives
a numerical factorization of `N` and

```text
k <= b+floor(a/2).                                     (17)
```

Conversely, take a numerical atom factorization of maximum length in `(11)`.
Grouping its factors produces a factorization of `N` into exactly `k`
integers at least three for every

```text
1<=k<=b+floor(a/2).                                    (18)
```

For each grouped degree `d`, THM-3438 supplies an `S_d` composition atom
`F_d`.  Their composition realizes `(18)`.  Therefore

```text
ExprLen_r(N)={1,2,...,b+floor(a/2)}.                    (19)
```

This is an existence spectrum across maps.  Equations such as `(15)` do not
prove that one Keller map has two inequivalent decompositions, and `(19)`
does not classify all maps in a grade.

## 6. The precise failure of the strong-atom-floor analogy

The old proposal treated degree four as a Keller analogue of an exceptional
forbidden value.  Equations `(3)` and THM-3438 show the opposite anatomy:

- grades `3,4,5,6,7,8` are all forced atomic;
- grade `9` is the first **mixed** grade, containing both the `S_9` atom and
  a composite `F_3 o F_3`;
- degree four is only the first ordinary composite integer made irreducible
  in `M` by the missing Keller grade two; and
- at every reducible grade, degree alone cannot distinguish the primitive
  `S_N` cover from an imprimitive composition.

The lawful survivor of the H-spectrum grammar is therefore a necessary
factorization floor: composition forces a nontrivial monodromy block system
whose block and quotient sizes lie in `M\{1}`.  The block lattice, not the
integer grade, is the restoring sidecar.  Tournament strong components and
ternary branch trees may visualize a supplied decomposition, but neither is
recovered from the scalar degree.

## 7. Harmonic subset sidecar

Equation `(3)` is also an exact subset identity:

```text
Atom(M)=P_odd disjoint-union 2P_odd disjoint-union {4,8}. (20)
```

Hence, for every cutoff `X>=8`,

```text
sum_(n<=X,n in Atom(M)) 1/n
 =sum_(p<=X,p odd prime)1/p
  +(1/2)sum_(p<=X/2,p odd prime)1/p+3/8.               (21)
```

Euler's divergence of the prime reciprocal series makes the grade-atom
support harmonic-divergent.  The reducible grades are harmonic-divergent too,
because they contain every multiple of nine.  Thus reciprocal mass cannot
separate the two factorization phases: both complementary arithmetic supports
are harmonic-thick.

This is the same controlled-forgetting lesson as the Boolean `V4` atlas.
The support is primary; one scalar reciprocal mass loses the factorization
word and the monodromy block sidecar.

## 8. Exact boundary

THM-3541 does not prove:

1. unique factorization, half-factoriality, or a length theorem for one fixed
   Keller map;
2. that the two numerical factorizations of `24` or `36` give equivalent
   polynomial maps;
3. classification up to tame conjugacy, W1/W2 moves, or stable equivalence;
4. any new `JC(2)` or `DC(2)` statement; or
5. a tournament identification, physical current, or LRC consequence.

It closes the arithmetic of the **degree grading** and the realizable
composition-atom expression lengths across each grade.  Monodromy and
intermediate fields remain indispensable for classifying an individual map.

## 9. Exact companion

Reproduce with

```text
python -B 04-computation/keller_grade_monoid_factorization_lengths_thm3541_20260816.py
python -B -O 04-computation/keller_grade_monoid_factorization_lengths_thm3541_20260816.py
```

The standard-library companion compares the closed atom classification with
an independent divisor search through `100000`, compares `(11)` with a
recursive factorization engine through `5000`, checks every interval and the
elasticity bound, enumerates the first nonunique factorizations, and verifies
the exact finite-prefix form of `(21)`.  Ordinary and optimized transcripts
match the stored output exactly, and the script has no executable `assert`.

**QED.**
