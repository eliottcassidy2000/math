---
id: THM-970
title: Scale-ten Hamming-six projective-prism owner-obligation obstruction
status: PROVED FINITE-EXACT — the complete 821,620,800-context primitive proper AP-centred common-scale-ten Hamming-six sheet bank is empty; exact reduction leaves only the 64 all-order-ten sign transversals, whose owner-obligation nerve is the triangular-prism graph K6 minus C6 and has no common face
source: codex-2026-07-17 scale-ten exact C++ certificate and independent Python referee
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860]
related: [THM-862, THM-957, THM-958, THM-960, THM-962, THM-963, THM-969, HYP-6820]
verification:
  - 04-computation/lrc13_scale_ten_hamming_six_frontier_scout_codex_c10.cpp
  - 05-knowledge/results/lrc13_scale_ten_hamming_six_frontier_scout_codex_c10.out
  - 04-computation/lrc13_scale_ten_hamming_six_referee_codex_c10.py
  - 05-knowledge/results/lrc13_scale_ten_hamming_six_referee_codex_c10.out
---

# THM-970 — scale ten dies on the projective owner-obligation prism

Let `R subset F_13^*` have six elements, put `P=[12] minus R`, and consider

```text
A=10P union {w_r:r in R},
w_r=10r (mod 13),              w_r>0,              w_r!=10r.      (1)
```

Assume that `A` is primitive and `M(A)<=1/13`.  For each replacement put

```text
D_r=10/gcd(10,w_r),
e_r=(D_r w_r/10) mod D_r.                                        (2)
```

Then no such packet exists.  Equivalently, the complete primitive proper
AP-centred common-scale-ten Hamming-six sheet bank is empty.  Every primitive
proper packet in (1) therefore has `M(A)>1/13`, and there is no surviving bank
on which to start a metric-height recursion.

## 1. The two-prime hereditary grammar

THM-860 supplies common-sheet coverage at each replacement owner and the
hereditary leave-one-out law

```text
lcm(D_s:s!=r)=10                         for every r in R.         (3)
```

The complete effective state alphabet is

```text
(D,e)=(1,0),(2,1),
      (5,1),(5,2),(5,3),(5,4),
      (10,1),(10,3),(10,7),(10,9).                              (4)
```

Because `10=2*5`, (3) says exactly that at least two coordinates carry the
factor two (`D in {2,10}`) and at least two carry the factor five
(`D in {5,10}`).  Inclusion-exclusion on the four divisor letters gives

```text
4^6 - 2(2^6+6*2*2^5) + 49 = 3,249                              (5)
```

admissible divisor words.  Here the simultaneous bad count is
`49=(1+6+6+30)+6`: either no `10` occurs and at most one `2` and one `5`
occur, or exactly one `10` occurs and all other letters are `1`.

For effective states, the four prime-membership cells have sizes

```text
neither:1,          two only:1,          five only:4,          both:4.
```

The same inclusion-exclusion is

```text
10^6 - (5^6+6*5*5^5) - (2^6+6*8*2^5) + 175 = 889,200,          (6)
```

where `175=(1+6+24+120)+24`.  Across the `binom(12,6)=924`
supports the exact labelled banks therefore have

```text
divisor contexts       924*3,249   =   3,002,076,
divisor/unit contexts  924*889,200 = 821,620,800.                (7)
```

These are complete finite counts, with no random, floating-point, asymptotic,
or height-cutoff step.

## 2. Literal masks and scalar capacity

For owner `o`, provider `r`, and state `(D,e)`, let `u` be the least CRT
representative satisfying

```text
u=Dr (mod 13),                    u=e (mod D).                   (8)
```

The provider covers sheet `ell in Z/10Z` exactly when the centred
representative of

```text
u(o^(-1)+13ell)                  modulo 13D                     (9)
```

lies in the half-open interval `(-D,D]`.  The certificate reconstructs every
mask from (8)--(9), and its output freezes the complete owner-one table.

Mask cardinality is independent of the effective unit.  In the normalized
ratio `q=r/o`, the four rows are

```text
D=1:   size 10 at q=1;                                      0 otherwise;
D=2:   size  5 at q in {1,6,7};                             0 otherwise;
D=5:   size  2 at q in {1,2,3,5,6,7,8,10,11};              0 otherwise;
D=10:  size  2 at q in {1,2,3,6,7,10,11};                  size 1 otherwise.
                                                                    (10)
```

Thus every owner first obeys the scalar necessary condition

```text
sum_(r in R) |E_(D_r,e_r)(r,o)| >= 10.                         (11)
```

Scanning all `3,002,076` labelled divisor contexts leaves exactly `1,200`
on `388` supports:

| divisor multiplicities | scalar-capacity contexts |
|---|---:|
| `1^0 2^0 5^0 10^6` | 64 |
| `1^0 2^0 5^3 10^3` | 120 |
| `1^0 2^2 5^0 10^4` | 36 |
| `1^0 2^2 5^2 10^2` | 48 |
| `1^0 2^2 5^3 10^1` | 48 |
| `1^0 2^2 5^4 10^0` | 12 |
| `1^0 2^3 5^0 10^3` | 344 |
| `1^0 2^3 5^1 10^2` | 336 |
| `1^0 2^3 5^2 10^1` | 144 |
| `1^0 2^3 5^3 10^0` | 48 |

Every other divisor context is impossible before its units are considered.

## 3. Owner-local reduction finds the projective object

For each scalar survivor, relax global compatibility by allowing a different
unit word at each owner.  Exact dynamic programming in the `1,024`-element
union-mask semilattice decides whether that owner can cover all ten sheets.
Requiring this owner-local condition at all six owners leaves precisely

```text
64 all-D=10 contexts.                                            (12)
```

Their supports are exactly the sign transversals

```text
R intersect (-R)=empty,             R union (-R)=F_13^*.         (13)
```

Equivalently, `R` is a section of the six-point quotient
`F_13^*/{+-1}`.  This explains the exact count `2^6=64`; the finite scan also
verifies the converse, namely that every such section is owner-locally
feasible.  The reduced literal fibre has

```text
64*4^6 = 262,144                                                (14)
```

global unit words.

The structural carrier is therefore not primarily a six-set of runners.  It
is a choice of signs over six projective residue classes.  Passing to runner
labels alone retains the section but obscures the quotient cycle exposed by
the final obstruction.

## 4. The faithful owner-obligation nerve is a triangular prism

Fix a sign transversal and define, for each owner `o`,

```text
O_o = {global unit words in {1,3,7,9}^R
       whose literal masks cover all ten sheets at o}.           (15)
```

Literal enumeration gives, uniformly over all 64 sections,

```text
|O_o|=32,                                                        (16)

|O_o intersect O_o'| = 0   if o'/o in {+2,-2,+6,-6},
                          4   otherwise,                         (17)

O_i intersect O_j intersect O_k = empty for all distinct i,j,k. (18)
```

On the projective quotient, the zero pairs in (17) form the multiplication-by-
two cycle

```text
[1]-[2]-[4]-[5]-[3]-[6]-[1].                                   (19)
```

Consequently the complete nonempty-intersection nerve has six vertices and
nine edges, no faces of dimension at least two, and one-skeleton

```text
K_6 minus C_6,                                                   (20)
```

the triangular-prism graph.  In particular, either edge of the complementary
zero-pair cycle already supplies two incompatible owner obligations, so the
six-fold intersection is empty.

The full `4^6=4,096` unit fibre has satisfaction profile

```text
unit words satisfying 0 owners: 3,940,
unit words satisfying 1 owner :   120,
unit words satisfying 2 owners:    36,
all higher counts             :     0.                           (21)
```

Equations (16)--(21) are checked for every one of the 64 sections, not inferred
from a single representative.

## 5. Tournament view, and why it is secondary

The pairwise observable is the zero-intersection predicate in (17).  Switch
each chosen residue to the positive representative of its projective class;
this gauges the sparse binary relation as the directed cycle

```text
1 -> 2 -> 4 -> 5 -> 3 -> 6 -> 1.                               (22)
```

All other pairs are ties.  Choosing the lexicographically first Hamiltonian
path in the tie graph `K_6 minus C_6` and orienting ties along that path gives a
completed tournament.  Across the 64 sections there are five joint
fingerprints; every tournament has

```text
sorted score sequence  (1,2,2,3,3,4),
directed triangles     6,
one SCC                 of size 6.                              (23)
```

The sparse-edge flip histogram is `{2:8,3:52,4:4}`, and the directed
Hamiltonian-path histogram is `{29:32,31:20,37:12}`.

This viewpoint is useful telemetry but not the faithful proof object.  The
completed tournament forgets that the six sparse pairs mean *disjoint owner
obligations*, forgets the intersection multiplicity four on the other nine
pairs, and forgets the absence of every triple face.  Owner obligations, or
equivalently projective classes decorated by their literal obligation sets,
are the vertices that preserve the LRC contradiction.

## 6. Complete certificate and independent replay

Owner-local feasibility is a necessary relaxation: if even an owner-dependent
unit choice cannot cover one owner, no globally shared unit word can work.
It is therefore exact to discard the other `3,002,076-64` divisor contexts.
The C++ certificate exhausts all `262,144` remaining literal words and finds

```text
global common-sheet survivors = 0.                              (24)
```

The Python referee independently rebuilds CRT masks, enumerates actual divisor
words by six leave-one-out `lcm` tests, uses Python sets for owner-local mask
reachability, replays the packed global fibres, reconstructs the projective
prism nerve, and reproduces the tournament telemetry.

Frozen SHA-256 hashes:

```text
C++ source     4170ddb73c48edaf09ae67f2e52d7f47c94d0f98255e5a1fe93387051a441c67
C++ output     69c20e00820e328b0ae14b1a19dafe1e1f5c7fb6d3bdea4524688b6702e81f4e
Python source  9c8cb86d9263761dfa601014194eaa07cc0d0028920f6714a22c00605fb6f537
Python output  262eedb2dda4441f52762071bd601ff0ceb75da428a8efc47b11512372795f5c
```

The referee additionally freezes the mask digest
`bef8724defec542a2d5da1fbfa93d621c65ff4bbbb36bf30cc345ebe2c6778b7`,
the owner-local bank digest
`db36a42aeda784871f22ac83ff9458d4a99e1173715c8b5ef794bdbe37058c11`,
and the concatenated obligation-profile digest
`e0acf11ba6826b6f9af49354c2a873679ec966e1abf51d6b9b925f730ea05538`.

Therefore the primitive proper AP-centred common-scale-ten Hamming-six bank is
empty.  ∎
