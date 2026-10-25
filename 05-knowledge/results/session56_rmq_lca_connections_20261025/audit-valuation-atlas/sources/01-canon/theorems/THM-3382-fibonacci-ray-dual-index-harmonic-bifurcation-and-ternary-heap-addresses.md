---
id: THM-3382
title: "Fibonacci-ray dual-index harmonic bifurcation and ternary heap addresses"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  The three THM-3379
  Fibonacci/Berggren rays have
  explicit, strictly interlacing ternary-heap addresses satisfying affine
  base-nine and common order-two C-finite recurrences.  On the flipped ray,
  Fibonacci time gives the support {3r+4}, of density 1/3 and divergent
  harmonic mass, while canonical heap address gives
  {(33*9^r-9)/8}, of density zero and finite harmonic mass with an exact
  Lambert/Dirichlet profile.  Depth gives the odd integers but identifies the
  flipped and unflipped odd-depth rays.  The full heap T4 bit is a six-kernel
  3-automatic sequence of density 1/2 and divergent harmonic mass, so neither
  automaticity nor the abstract ray alone determines harmonic class: the
  injection into N is a necessary sidecar.  No encoding-invariant density,
  global tournament, LRC current, or JC flux follows.
source: codex-2026-08-14-dual-index-harmonic-bifurcation
audit: independent heap, language, recurrence, Dirichlet, kernel, harmonic, and scope audit
depends_on:
  - THM-3379-fibonacci-ray-local-t4-bit-is-mod3-boolean-projection
related:
  - THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar
  - THM-3364-cyclotomic-boolean-clocks-berggren-t4-xor-and-crt-phase
  - THM-2000-support-harmonic-abel-dini-figurate-surface
  - THM-2005-support-dirichlet-automatic-tournament-atlas
script: 04-computation/fibonacci_ray_dual_index_harmonic_bifurcation_thm3382.py
output: 05-knowledge/results/fibonacci_ray_dual_index_harmonic_bifurcation_thm3382.out
script_sha256: 1049d46e642e6d7f4640e625ba5fe1eb6330c34dc0a5ee34e8a96885fb1429ae
output_sha256: e3217b7857d26ba13f9b6d490cf8141fcc95f690b11b1e3064b08380d7c6d060
semantic_sha256: 67e8bf54060f0269d1d47e9b22c0928da7d0c1377ee0deba55f95d98be7a0b6a
hash_basis: LF-normalized bytes
---

# THM-3382 -- Fibonacci time and ternary address have different harmonic class

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Connection contract

THM-3379 identifies the flipped local-`T4` fibre on the three distinguished
Fibonacci/Berggren rays.  That statement uses the Fibonacci index as its copy
of the natural numbers.  A tree also has a canonical breadth-first copy of the
natural numbers, and the two weights are not equivalent.

| field | connection |
|---|---|
| source | the same ordered ray objects `(ray,r)` and their THM-3379 words |
| targets | Fibonacci time, word depth, and ternary heap address |
| maps | `3r+j`, word length, and the recursive address `(1)` below |
| preserved by time/heap | the ray object, chronological order, word, parameter, triple and `T4` bit |
| destroyed by depth | ancestry: `R0` and `R1` have the same depth image |
| changed | growth of the integer label, hence density and reciprocal weight |
| needed sidecar | which injection into `N` defines the harmonic series |
| hostile | the full heap bit is automatic but has density `1/2`, not zero |

This is an exact Fibonacci--Berggren realization of the carrier distinction
already guarded by THM-3359 and THM-2000/2005, not a claim that the general
principle is new.

## 2. Canonical ternary heap address

For a word `w` in `{A,B,C}*`, define

```text
h(empty)=0,
h(wA)=3h(w)+1,   h(wB)=3h(w)+2,   h(wC)=3h(w)+3.          (1)
```

Words of depth `d` map in lexicographic order onto the consecutive interval

```text
[(3^d-1)/2, (3^(d+1)-3)/2].                              (2)
```

This follows by induction: every positive address `n` has the unique parent
and branch selected by `n mod 3`, with residue zero assigned to `C`.  Thus `h` is the
usual breadth-first/shortlex heap bijection from finite ternary words to
`N_0`.  The initial `1,2,3` offsets matter; using bare base-three digits with
a leading `A=0` would lose injectivity at startup.

## 3. The three Fibonacci rays

THM-3379 gives, for `r>=0`, the root-to-child words

```text
R2: (BA)^r,          R0: A(BA)^r,          R1: C(BC)^r.   (3)
```

Appending `BA` and `BC` acts on every existing address by

```text
h(wBA)=9h(w)+7,                  h(wBC)=9h(w)+9.           (4)
```

Consequently:

| ray | Fibonacci time | depth | heap address |
|---|---:|---:|---:|
| `R2` | `3r+2` | `2r` | `h_2(r)=7(9^r-1)/8` |
| `R0` | `3r+3` | `2r+1` | `h_0(r)=(15*9^r-7)/8` |
| `R1` | `3r+4` | `2r+1` | `h_1(r)=(33*9^r-9)/8` |

The three affine laws are

```text
h_2(r+1)=9h_2(r)+7,   h_2(0)=0,
h_0(r+1)=9h_0(r)+7,   h_0(0)=1,
h_1(r+1)=9h_1(r)+9,   h_1(0)=3.                           (5)
```

Each is C-finite after homogenizing the constant:

```text
h_j(r+2)=10h_j(r+1)-9h_j(r).                             (6)
```

Direct subtraction gives the strict chronological interlacing

```text
h_2(r)<h_0(r)<h_1(r)<h_2(r+1).                           (7)
```

In particular, heap address is not an arbitrary permutation of the carrier.
It preserves the order of all three rays while exponentially dilating it.  The
ordinary generating functions are

```text
sum h_2(r)t^r = 7t/((1-t)(1-9t)),
sum h_0(r)t^r = (1+6t)/((1-t)(1-9t)),
sum h_1(r)t^r = (3+6t)/((1-t)(1-9t)).                    (8)
```

## 4. One fibre, three integer carriers

The flipped fibre is `R1`.  In Fibonacci time THM-3379 proves

```text
K={3r+4:r>=0}={4,7,10,...}.                              (9)
```

It has natural density `1/3` and

```text
sum_(k<=x,k in K) 1/k
 = (1/3)log x + gamma/3 + pi/(6sqrt(3))
   + log(3)/6 - 1 + O(1/x).                              (10)
```

Depth sends the same fibre to the odd positive integers `{2r+1}`.  This image
has density `1/2` and divergent harmonic mass.  However, `R0` has exactly the
same depth image and the opposite `T4` bit.  Depth is therefore a quotient,
not an injective realization of the three-ray carrier.

Heap address is injective and gives instead

```text
H={h_1(r):r>=0}
 ={(33*9^r-9)/8:r>=0}
 ={3,36,333,3006,...}.                                  (11)
```

Its counting function is `O(log x)`, hence its natural density is zero and it
is not ultimately periodic.  All three positive ray-address supports are
harmonically summable; for the flipped one the mass is treated exactly next.

This does not contradict THM-3359.  Although the sequence `h_1(r)` is
C-finite, `(11)` is its **value support**, not the modular index support of a
C-finite sequence.

## 5. Exact heap harmonic and Dirichlet profiles

Put `q=1/9` and `a=3/11`.  Factoring `(11)` gives the convergent Lambert sum

```text
M_H=sum_(r>=0) 1/h_1(r)
   =(8/33) sum_(r>=0) q^r/(1-aq^r).                      (12)
```

Equivalently, with `(a;q)_infinity=product_(r>=0)(1-aq^r)`,

```text
M_H=-(8/33) d/da log(a;q)_infinity at a=3/11, q=1/9.     (13)
```

If `P_R` is the sum of the first `R` terms and `T_R=M_H-P_R`, then monotonicity
of `1/(1-aq^r)` gives the exact rational enclosure

```text
(3/11)9^(-R) < T_R
 <= ((3/11)9^(-R))/(1-(3/11)9^(-R)).                    (14)
```

Thus finite mass is proved without a decimal approximation.  More generally,
for `Re(s)>0`, the binomial series and absolute convergence give

```text
D_H(s)=sum_(r>=0) h_1(r)^(-s)
 =(8/33)^s sum_(m>=0) (s)_m/m! (3/11)^m
                    /(1-9^(-(s+m))).                    (15)
```

The abscissa of convergence is zero.  In contrast,

```text
D_K(s)=3^(-s) zeta(s,4/3)                               (16)
```

has its harmonic pole at `s=1`.  The abstract flipped ray has not changed;
the injection controlling the reciprocal weight has.

## 6. Sparse automatic and Mahler forms

In ordinary base three, `(11)` has the regular address language

```text
{10} union 11(01)*00.                                   (17)
```

Hence its indicator is `3`-automatic.  For the sparse address generating
functions, with the root address omitted from `G_2`, equations `(5)` give

```text
G_2(z)=z^7+z^7 G_2(z^9),
G_0(z)=z  +z^7 G_0(z^9),
G_1(z)=z^3+z^9 G_1(z^9).                                (18)
```

So regularity and a finite automaton are compatible here with zero density
and finite harmonic mass.

## 7. The full-tree automatic hostile

Let `e(n)` be THM-3364's local labelled-`T4` bit at the unique word of heap
address `n`.  Its reset/XOR laws become

```text
e(0)=0,
e(3n+1)=0,       e(3n+2)=1-e(n),       e(3n+3)=1.        (19)
```

The full `3`-kernel has exactly six states.  Writing

```text
E=e, C=1-e, U=1-delta_0, D=delta_0, Z=0, O=1,
```

their residue transitions are

| state | digit `0` | digit `1` | digit `2` |
|---|---|---|---|
| `E` | `U` | `Z` | `C` |
| `C` | `D` | `O` | `E` |
| `U` | `U` | `O` | `O` |
| `D` | `D` | `Z` | `Z` |
| `Z` | `Z` | `Z` | `Z` |
| `O` | `O` | `O` | `O` |

Thus `e` is `3`-automatic.  Its level-`d` one-count is

```text
(3^d-(-1)^d)/2.                                         (20)
```

For a prefix statement, put

```text
t(n)=2e(n)-1,                 T(N)=sum_(1<=n<=N)t(n).
```

Grouping `(19)` into three-term blocks gives, with `T(0)=0` and `T(1)=-1`,

```text
T(3m)=1-T(m-1)       (m>=1),
T(3m+1)=-T(m-1)     (m>=1),
T(3m+2)=-T(m)       (m>=0).                             (21)
```

Induction yields `|T(N)|<=1+floor(log_3 N)`.  Therefore the full heap support
has natural density `1/2`.  Summation by parts also shows that
`sum t(n)/n` converges and

```text
sum_(n<=x,e(n)=1) 1/n
 = (1/2)log x + C_tree + O(log(x)/x),                    (22)
C_tree=gamma/2+(1/2)sum_(n>=1)t(n)/n.
```

Its ordinary generating function satisfies the exact Mahler equation

```text
E(z)=z^2(1+z)/(1-z^3)-z^2 E(z^3).                       (23)
```

Equations `(18)` and `(22)` are the decisive hostile pair: both supports are
automatic shadows of the same tournament-bit automaton, but one is lacunary
and summable while the other has density `1/2`.

## 8. Boundaries and verification

The word-direction hostile is immediate: root-to-child `BA` has heap address
`7`, while the reversed word `AB` has address `5`.  The first flipped object
has Fibonacci weight `1/4` but heap weight `1/3`, so no finite-prefix
monotonicity explains the eventual convergence.  The depth collision and the
full-tree density in Sections 4 and 7 block the two most tempting extensions.

The exact companion pins THM-3359 and THM-3379; checks the heap bijection on
`265,720` addresses; checks `2,403` ray formulas through `r=800`, including
interlacing, `(6)`, `(8)`, all `801` language words in `(17)`, and `(18)`;
gives an exact twelve-term bracket `(14)` and `80` factored Dirichlet terms;
checks `(19)--(21)` and `(23)` through depth eleven; and performs `1,594,320`
explicit six-kernel transitions.  It contains no floating-point literal or
`assert`, and ordinary and optimized runs LF-normalized-byte-match the stored
output.  An independent audit rederived `(15)`, `(17)`, `(21)` and the
harmonic asymptotics in addition to replaying the companion.

The theorem concerns three distinguished rays and one canonical heap gauge.
It does not make harmonic density invariant under relabelling, recover a
nonlocal tournament, orient equal-hypotenuse ties, or provide LRC phase/current
or JC flux.

**QED.**
