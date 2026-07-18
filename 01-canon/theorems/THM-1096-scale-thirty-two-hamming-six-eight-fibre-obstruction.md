---
id: THM-1096
title: Scale-thirty-two Hamming-six eight-fibre obstruction
status: PROVED STRUCTURAL + INDEPENDENT FINITE-EXACT — a deterministic algebraic-CRT Python primary and a separately developed standard-library C++ literal-CRT referee agree on all 12,281 hereditary words, 11,347,644 labelled support/order contexts representing 883,622,412,288 raw states, the complete 3,450-row scalar bank on 284 supports, all 20,700 Z/8 owner relaxations, and the terminal live-owner histogram 0:2802,1:456,2:192. The referee additionally checks every exact union bank and the literal-cover implication. Frozen outputs replay byte-for-byte.
source: codex-2026-07-18-S67 scale-thirty-two continuation
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860, THM-994, THM-1072, THM-1090]
related: [THM-983, HYP-6820]
verification:
  - 04-computation/lrc13_scale_thirty_two_hamming_six_eight_fibre_obstruction_codex_c32.py
  - 05-knowledge/results/lrc13_scale_thirty_two_hamming_six_eight_fibre_obstruction_codex_c32.out
  - 04-computation/lrc13_scale_thirty_two_hamming_six_eight_fibre_referee_codex_c32.cpp
  - 05-knowledge/results/lrc13_scale_thirty_two_hamming_six_eight_fibre_referee_codex_c32.out
---

# THM-1096 — scale thirty-two has a terminal eight-fibre deficit

The theorem is:

> **The primitive proper AP-centred common-scale-32 Hamming-six face is
> empty.**

The proof-facing carrier is the quotient

```text
Z/32 -> Z/8.
```

Every order dividing eight is retained exactly as a union of its eight thick
fibres.  Each order-sixteen or order-thirty-two provider is then allowed to
choose its unit independently outside that anchor union.  This is a sound
upper relaxation.  Across the complete scalar bank it can remain live at no
more than two of the six owners, so no shared unit word can cover every owner.

The Python primary enumerates every labelled row directly.  The C++ referee
starts from bounded literal CRT search, uses a different owner gauge and
enumeration order, and independently reconstructs the theorem-bearing counts.

## 1. Hereditary two-adic grammar

At `c=32`, the effective orders and unit counts are

```text
D                 1  2  4  8  16  32
# units            1  1  2  4   8  16.
```

For an order word `(D_1,...,D_6)`, hereditary leave-one-out lcm means

```text
lcm(D_i : i != j) = 32                 for every j.            (1)
```

Since only order `32` carries the full two-adic valuation five, (1) is
equivalent to

```text
at least two coordinates have order 32.                        (2)
```

Deletion preserves a top-order coordinate exactly under (2), proving both
directions without a search.  Consequently the exact unweighted grammar has

```text
6^6 - 5^6 - 6*5^5 = 12281                                  (3)
```

words.  The total unit weight of the five proper divisors is `16`, equal to
`phi(32)`.  Hence the literal-state count per support is

```text
32^6 - 16^6 - 6*16*16^5
  = 32^6 - 7*16^6
  = 956301312.                                                (4)
```

Across the `binom(12,6)=924` labelled supports, (3)-(4) represent

```text
11347644 labelled support/order contexts,
883622412288 raw labelled unit states.                         (5)
```

Both implementations check the literal leave-one-out predicate against (2)
on every order word and check (3)-(5).

## 2. Literal masks and scalar capacities

For provider label `a`, owner label `b`, effective order `D`, and unit `u`, let
`B` be the unique CRT representative

```text
B = D*a (mod 13),                 B = u (mod D).
```

The local sheet mask is

```text
M(a,D,u;b)
 = {t in Z/32 : <B*(b^(-1)+13t)>_(13D) lies in (-D,D]}.       (6)
```

Writing `r=a/b in F_13^*`, direct interval counting gives the cardinality
vectors below; columns are `D=1,2,4,8,16,32`:

```text
r= 1: 32 16  8  8  6  5      r= 7:  0 16  8  4  4  5
r= 2:  0  0  0  4  4  5      r= 8:  0  0  0  4  6  5
r= 3:  0  0  8  4  4  5      r= 9:  0  0  8  8  6  5
r= 4:  0  0  8  8  6  5      r=10:  0  0  8  4  4  5
r= 5:  0  0  0  4  6  5      r=11:  0  0  0  4  4  5
r= 6:  0 16  8  4  4  5      r=12:  0  0  0  4  4  4.       (7)
```

The general formula behind (7) is

```text
|M(a,D,u;b)|
 = (32/D) * #{x in (-D,D] : x = D*r (mod 13)}.                (8)
```

It is independent of `u`.  The primary constructs every base algebraically
and by bounded literal search.  The referee begins with literal search and
then proves its ratio/rotation normalization against every labelled
provider-owner mask before using the gauge.

Summing (7) over providers is a necessary scalar capacity test.  On all
contexts in (5), the number of owners reaching capacity 32 has histogram

```text
0:76548, 1:2800212, 2:4692582, 3:2946408,
4:743040, 5:85404, 6:3450.                                   (9)
```

Thus exactly `3450` rows, on `284` supports, can survive at all six owners.
None uses order one.  They have 23 order-multiplicity profiles, 1649 capacity
vectors, and represent `621084672` literal unit words.  The support histogram
by survivor count is

```text
0:640, 1:96, 2:48, 3:30, 17:24,
22:24, 33:12, 34:24, 38:24, 54:2.                            (10)
```

Every context omitted from this 3450-row bank is already impossible at some
owner by the scalar union bound.

## 3. The sound `Z/8` relaxation

Replacing `t` by `t+D` in (6) changes the CRT argument by a multiple of
`13D`.  Therefore an order-`D` mask is the pullback of a subset of `Z/DZ`.
When `D|8`, it is in particular periodic under `t -> t+8` and is an exact
union of thick fibres for `Z/32 -> Z/8`.

Fix an owner and retain every provider with `D|8`.  For one exact anchor union

```text
Q = union_(D_i|8) M_i(u_i),
```

arbitrary nonanchor choices satisfy

```text
|union_i M_i(u_i)|
 <= |Q| + sum_(D_i not|8) |M_i(u_i) minus Q|
 <= |Q| + sum_(D_i not|8) max_e |M_i(e) minus Q|.             (11)
```

Maximizing the right side over the exact anchor bank gives a sound upper bound
`U_8`.  Independent nonanchor maxima forget their shared unit word and mutual
overlaps; both losses only enlarge the attainable union.  Thus `U_8<32` is a
terminal owner obstruction.

Across all `3450*6=20700` owner obligations, the exact anchor banks have
sizes

```text
1:8064, 2:4008, 3:636, 4:3576, 6:2088, 7:240,
8:192, 9:36, 10:684, 12:984, 14:168, 26:24.                  (12)
```

The relaxed bounds are

```text
20:24, 21:396, 22:420, 23:564, 24:1416,
25:2352, 26:3252, 27:3780, 28:3708, 29:2208,
30:1260, 31:480, 32:624, 33:192, 34:24.                      (13)
```

Most importantly, live owners per scalar row have histogram

```text
0:2802, 1:456, 2:192.                                       (14)
```

No row remains live at three owners, much less all six.  A global unit word
would have to cover all 32 sheets at each owner, whereas every scalar row has
at least four owner projections that are impossible for **every** unit choice.
Equations (9), (11), and (14) prove the emptiness statement.

## 4. Exact-union sharpness and covariance

As a non-load-bearing sharpness audit, the C++ referee separately retains all
provider masks in an immutable-union DP.  It constructs `403733376` reachable
masks, with largest owner bank `227272`, and obtains

```text
exact maximum union:
20:24, 21:396, 22:432, 23:612, 24:1524, 25:3312,
26:3852, 27:4596, 28:3660, 29:1680, 30:276,
31:96, 32:240;

exact feasible owners/context:
0:3258, 1:144, 2:48.                                        (15)
```

Thus the exact answer independently has no all-owner row.  The relaxation
equals the exact maximum on 13272 owner obligations and exceeds it by at most
five on the others.  It has 600 threshold false positives, all harmless and
all already covered by (14).  The referee checks the literal implication

```text
exact cover => U_8 >= 32
```

on every one of the 20700 labelled owner obligations.

Multiplication by `F_13^*` gives scalar-context orbit sizes

```text
6:3, 12:286
```

and support orbit sizes

```text
2:1, 6:1, 12:23.
```

These are telemetry only; the primary visits every labelled row.

## 5. Tournament Analysis and the challenged vertex set

Runner tournaments are not the proof object.  The useful vertices are the six
owner-proof obligations.  The primary directs an edge by the lexicographic
pair `(U_8, scalar capacity)` and uses coordinate order as the tie Hamiltonian
path.  The referee independently uses the richer keys

```text
K8     = (U_8>=32, U_8, scalar capacity, anchor-bank size),
Kexact = (exact cover, exact maximum, scalar capacity, exact-bank size).
```

Both gauges are total orders, so every one of the 3450 tournaments is
transitive: score histogram `(0,1,2,3,4,5)`, zero directed triangles, six
singleton SCCs, and one Hamiltonian path.  Switching from `K8` to `Kexact`
changes between zero and five edges per row.

This quotient preserves proof-obligation hardness but destroys fibre owners,
nonanchor intersections, and shared-unit compatibility.  The actual proof
carrier is the eight-fibre anchor incidence, not the tournament.  This
challenges the assumption that tournament vertices should be runners, gaps,
or arcs: here even the owner tournament is only a transparent audit sidecar.

## 6. Independent replay and scope

The frozen primary source/output SHA-256 values are

```text
6c023306ac6878c2a0b1c3d9e2fa22bbfff2efd2becf8787ba34ebfd6aa87e28
2bdea8baf5cb502ae46a2bd77cddfa55877834582fd34da6c09359f8b0789be3
```

and the independent referee source/output values are

```text
7ce838b13bf773298e2770fc6dcc2199a2e25984c83b04e97fe7c207a773bad2
987ad51a551d4236ba484175cafe2947c1e43cff4ce8a61089495fb4b7581816
```

respectively.  Both outputs replay byte-for-byte.  The referee uses only the
C++ standard library and bounded literal CRT search; it does not import the
primary's tables, NumPy batching, SHA digests, or owner-loop implementation.

This theorem closes only the primitive proper AP-centred common-scale-32
Hamming-six face.  It does not prove uniform sporadic emptiness or LRC(14).
Scale 33 is the next untreated common-scale Hamming-six face.  ∎
