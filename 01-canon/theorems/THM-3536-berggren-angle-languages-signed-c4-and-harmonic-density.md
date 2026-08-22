---
id: THM-3536
title: "Berggren angle languages, signed C4 holonomy, and harmonic density"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED, after the
  Fibonacci sidecar repair in MISTAKE-418.  The two irrational descendant-
  angle walls have inverse itineraries (BABB)^omega and (BCBB)^omega.  The
  exact U-obtuse/acute/D-obtuse counts satisfy four-step anti-recurrences,
  have densities 16/41,9/41,16/41, and have bounded period-eight residuals.
  Their ancestry languages are regular; shortlex counting has discrepancy
  O(log N), giving the same coefficients in the three harmonic subseries.
  The phase carrier is a directed C4 with two missing antipodal pairs and the
  nonzero class in H^1(|C4|;F2), not a canonical tournament.  The normalized
  Fibonacci three-ray sidecar has chamber word D,D,U, not all D.
source: codex/berggren-angle-language/2026-08-16
audit: >
  A read-only independent audit rederived both walls, every inverse-orbit
  point, CDF composition order, all count recurrences, densities, residuals,
  generating functions, the finite automata, the uniform-cylinder and Abel-
  summation arguments, and the signed-C4/tournament boundary.  It found and
  repaired only the omitted odd/odd Fibonacci normalization (MISTAKE-418).
depends_on:
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-3334-berggren-parabolic-spine-gaussian-collision-torsor
  - THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction
related:
  - MISTAKE-418
script: 04-computation/berggren_descendant_angle_language_density_20260816.py
output: 05-knowledge/results/berggren_descendant_angle_language_density_20260816.out
script_sha256: f7d3eb39285751f7be0db083ceff050f8347dddfefe43c0efe8b3ab6d6108f45
output_sha256: 662e7ccd78722b27fc27653239151b25dac28f35af2c5b8acba60b1bef6b5a54
semantic_sha256: acc0170636f7ce47b074214c116930ea6728fa3ea96314cc622e41f83966b2ec
hash_basis: LF-normalized bytes
---

# THM-3536 -- exact angle languages on the Berggren tree

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Use THM-2596's reduced Euclid slope `x=m/n in (0,1)` and branch maps

```text
A(x)=1/(2-x),       B(x)=1/(2+x),       C(x)=x/(2x+1).   (1)
```

Their images are the disjoint intervals `(1/2,1)`, `(1/3,1/2)`, and
`(0,1/3)`.  Starting at `x_0=1/2`, words in `{A,B,C}` enumerate the primitive
Pythagorean/Berggren tree.

## 1. The two angle walls

For the parent triple

```text
(a,b,c)=(n^2-m^2,2mn,n^2+m^2),
```

THM-3334 classifies its descendant triangle by the ratio `b/a`.  Since

```text
b/a=2x/(1-x^2)                                          (2)
```

is strictly increasing on `(0,1)`, the U-obtuse/acute/D-obtuse walls are

```text
alpha=(sqrt(145)-9)/8,       beta=(sqrt(145)-8)/9.      (3)
```

They satisfy

```text
4alpha^2+9alpha-4=0,       9beta^2+16beta-9=0,
0<alpha<beta<1/2.                                           (4)
```

Thus a rational tree slope is U-obtuse below `alpha`, acute between the
walls, and D-obtuse above `beta`.  No equality can occur because both walls
are irrational.

The inverse reduction map selected by the three branch intervals is

```text
T(t)=t/(1-2t) on (0,1/3),
T(t)=1/t-2     on (1/3,1/2),
T(t)=2-1/t     on (1/2,1).                              (5)
```

Exact arithmetic in `Q(sqrt(145))` gives

```text
alpha --B--> (sqrt(145)-7)/8
      --A--> (17-sqrt(145))/12
      --B--> (sqrt(145)-7)/12
      --B--> alpha,

beta  --B--> (sqrt(145)-10)/9
      --C--> sqrt(145)/29
      --B--> -2+sqrt(145)/5
      --B--> beta.                                      (6)
```

All displayed points are irrational and avoid both seams.  The inverse wall
words are therefore exactly `(BABB)^omega` and `(BCBB)^omega`.

## 2. Four-step CDF returns

For any rational starting slope `x`, define

```text
F_n(t;x)=3^(-n) #{w in {A,B,C}^n:w(x)<t}.               (7)
```

The interval partition and the fact that `B` reverses orientation give

```text
t in (0,1/3): F_(n+1)(t)=      F_n(Tt)/3,
t in (1/3,1/2): F_(n+1)(t)=2/3-F_n(Tt)/3,
t in (1/2,1): F_(n+1)(t)=2/3+F_n(Tt)/3.                 (8)
```

Compose `(8)` in the order in `(6)`.  For every rational start `x`, not only
the root,

```text
F_(n+4)(alpha;x)=32/81-F_n(alpha;x)/81,
F_(n+4)(beta ;x)=50/81-F_n(beta ;x)/81.                 (9)
```

The negative return multiplier is the load-bearing orientation datum.

## 3. Exact counts, densities, and generating functions

At the root let `U_n,A_n,D_n` be the three angle-class counts at depth `n`.
Multiplying `(9)` by `3^(n+4)` gives

```text
U_(n+4)=32*3^n-U_n,       U_0..U_3=(0,1,4,10),
A_(n+4)=18*3^n-A_n,       A_0..A_3=(0,1,2,6),
D_(n+4)=32*3^n-D_n,       D_0..D_3=(1,1,3,11).          (10)
```

The fixed points of these affine recurrences are

```text
U_n/3^n ->16/41,       A_n/3^n ->9/41,
D_n/3^n ->16/41.                                        (11)
```

More precisely, put

```text
R_n=(41U_n-16*3^n,41A_n-9*3^n,41D_n-16*3^n).
```

Then `R_(n+4)=-R_n`, so the exact period-eight residuals are

```text
(-16,-9,25), (-7,14,-7), (20,1,-21), (-22,3,19),
(16,9,-25), (7,-14,7), (-20,-1,21), (22,-3,-19).       (12)
```

The ordinary generating functions are

```text
U(z)=(z+z^2-2z^3+2z^4)/((1-3z)(1+z^4)),
A(z)=(z-z^2)/((1-3z)(1+z^4)),
D(z)=(1-2z+2z^3-z^4)/((1-3z)(1+z^4)).                 (13)
```

The pole at `z=1/3` carries `(11)`; the roots of `1+z^4` carry the bounded
eight-periodic discrepancy.

## 4. Regular languages and harmonic subseries

Read a word from leaf to root.  For either wall, an undecided automaton state
records one of the four points in `(6)` and the current inequality
orientation.  The unique branch interval containing the wall advances to the
next phase; `B` reverses the orientation, while `A,C` preserve it.  Either
other branch enters an accepting or rejecting sink.  This gives at most
eight undecided states and two sinks for each wall.  Regular languages are
closed under reversal and Boolean operations, so all three root-to-leaf angle
languages are regular.

The uniformity in `(9)` is stronger than root counting.  For any prefix
state `x`, the centered suffix count `E_r(x)` satisfies

```text
E_(r+4)(x)=-E_r(x).                                     (14)
```

For `r<4`, the trivial bounds `0<=Q_r(x)<=3^r` are uniform in `x`; hence
every suffix cylinder of size `3^r` contains

```text
delta_T*3^r+O(1)                                       (15)
```

nodes of type `T`, uniformly over its prefix.

Enumerate the tree in shortlex order with any fixed ordering of `{A,B,C}`.
A lexicographic initial segment of one level is a disjoint union of at most
two prefix cylinders at each remaining depth.  Summing `(15)` over those
cylinders and the completed earlier levels gives

```text
#{m<=N:m in S_T}=delta_T N+O(log N),                    (16)
```

where `(delta_U,delta_A,delta_D)=(16/41,9/41,16/41)`.
Partial summation yields constants `C_T` with

```text
sum_(m<=N,m in S_T) 1/m
 =delta_T log N+C_T+O((log N)/N).                       (17)
```

The logarithmic coefficient is independent of the chosen letter order; the
constant need not be.  These are ancestry-shortlex densities and index
weights, not densities by hypotenuse or Farey denominator.

## 5. The lawful four-state object and its six pairs

The slope of the affine update in `(8)` is positive on `A,C` and negative on
`B`.  Both wall cycles therefore carry the edge-sign word

```text
(-,+,-,-),                                              (18)
```

whose product is `-1`.  Modulo two, the negative-edge indicators form

```text
(1,0,1,1),
```

the nonzero class in

```text
H^1(|C4|;F2)=F2.                                        (19)
```

Vertex sign changes add coboundaries and preserve the odd cycle sum.  This
holonomy is the conceptual source of the negative multiplier in `(9)`, the
factor `1+z^4`, and period eight.

The successor relation is a directed `C4`.  Among the six unordered pairs of
four phases, it sees the four cycle pairs and omits the two antipodal pairs.
There are four ways to orient those diagonals and complete the cycle to a
tournament, but none is invariant under phase rotation: the half-turn fixes
each diagonal setwise and swaps its endpoints.  Thus the intrinsic object is
a signed partial tournament, not a canonical `T4` or a tournament on six
vertices.

## 6. Corrected Fibonacci sidecar

Let `F_0=0,F_1=1` and begin with `k=2`.  If `k` is not `1 mod 3`, the pair
`(F_k,F_(k+1))` has opposite parity and its canonical Euclid slope is at
least `1/2>beta`; the node is D-obtuse.

If `k=1 mod 3`, both entries are odd.  THM-3339's primitive normalization is

```text
T(m,n)=((n-m)/2,(n+m)/2),       x'=(n-m)/(n+m).          (20)
```

Since `n<=2m`, one has `x'<=1/3<alpha`, so the normalized node is U-obtuse.
Therefore the canonical Fibonacci three-ray locus has exact chamber word

```text
D,D,U,D,D,U,... .                                       (21)
```

The raw Fibonacci slope is not a lawful substitute at the odd/odd indices;
MISTAKE-418 records the minimal witness `(3,5)->(1,4)`.

## 7. Evidence and stopping boundary

The exact companion uses rational arithmetic in `Q(sqrt(145))`, contains no
floating literals or `assert`, directly enumerates all `88,573` nodes through
depth ten, verifies `(10)--(13)` through depth forty, and checks `484`
four-step suffix recurrences from `121` prefix states.  It exhausts the four
tournament completions and all sixteen sign gauges.  Normal and optimized
runs match the stored transcript.

```text
python -B 04-computation/berggren_descendant_angle_language_density_20260816.py
python -B -O 04-computation/berggren_descendant_angle_language_density_20260816.py
```

This theorem proves no density by arithmetic height, no canonical map to the
separate Fibonacci `P^1(F3)` cycle, no RXTX tensor identity, no LRC word-
current or D5 flux class, no Jacobian statement, and no matrix-multiplication
bound.  LRC(14) remains open.
