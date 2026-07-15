---
id: THM-850
title: The chi-seven face-charge, cyclotomic B3 recursion, and equitable operation-congruence census
status: PROVED (all-size face-charge intertwiner, augmented-character sign law, cyclotomic recurrence, residue census, finite Radon rank theorem, and address/carry odometer) + FINITE-EXACT (depth-399 replay, Paley fingerprints, 21-body flood incidence, three-needle obstruction, and depth-11 Radon alias)
source: codex-2026-07-15-S15
depends_on: [THM-830]
related: [THM-211, THM-741, THM-801, THM-841, HYP-3575, HYP-6895, HYP-6900]
verification:
  - 04-computation/thm850_chi7_face_charge_referee_codex_S15.py
  - 05-knowledge/results/thm850_chi7_face_charge_referee_codex_S15.out
  - 04-computation/lrc14_fano_chi7_flood_needle_obstruction_codex_S11.py
  - 05-knowledge/results/lrc14_fano_chi7_flood_needle_obstruction_codex_S11.out
  - 04-computation/chi7_radon_carry_alias_codex_S15.py
  - 05-knowledge/results/chi7_radon_carry_alias_codex_S15.out
---

# THM-850 - the chi-seven face-charge law

THM-830 turns the three staircase faces into a commuting tournament-operation
monoid.  Its ordered `A,B,C` mask has the distinguished binary/Fano gauge:

```text
q(A)=1,                 q(B)=2,                 q(C)=4       in F_7.       (1)
```

The three charges are the quadratic residues modulo seven.  This one choice
simultaneously explains the historical `++-+--+` word, produces an exact
cyclotomic specialization of the full `A+B+C-D-E-F+G` recurrence, and gives
an equitable seven-class partition of every staircase carrier.  It is an
operation-indexed congruence, not a blue/black or LRC codec.

## 1. The face monoid has an additive Fano charge

Use THM-830's positive tile coordinates and shifted address

```text
(u,g,w)=(b,a-b-1,n-a+1),
(x,y,z)=(w-1,u-1,g-1),             x+y+z=n-3.               (2)
```

The coordinates `x,y,z` count inverse lifts in the `A,C,B` directions,
respectively.  More precisely, THM-830's operation

```text
D_(alpha,beta,gamma)(G)_(p,q)
  = G_(p+beta+gamma,q+gamma)                               (3)
```

has inverse lift

```text
(u,g,w) -> (u+gamma,g+beta,w+alpha),
(x,y,z) -> (x+alpha,y+gamma,z+beta).                       (4)
```

Therefore

```text
q(D_(alpha,beta,gamma))=alpha+2 beta+4 gamma mod 7         (5)
```

is a monoid homomorphism.  If `zeta=exp(2 pi i/7)` and a tile address has
weight

```text
wt(x,y,z)=zeta^(x+4y+2z),                                  (6)
```

then an inverse `D_(alpha,beta,gamma)` lift multiplies its weight by
`zeta^q(D)`.  Equation (6) is the exact character intertwiner.

**Proof of (4).**  A lower tournament chord `(p,q)` is the lower tiling chord
`(p+1,q)`.  Its upper chord in (3) is
`(p+beta+gamma+1,q+gamma)`.  Substitution in (2), while increasing size by
`alpha+beta+gamma`, adds `gamma,beta,alpha` to `u,g,w`.  Composition adds the
three operation exponents, proving (5)-(6).  This also gives a direct
tournament proof rather than a numerical identification.  QED.

## 2. The seven Mobius signs are an augmented Legendre character

Let `chi_7(0)=0`, `chi_7(q)=+1` for `q in {1,2,4}`, and `chi_7(q)=-1` for
`q in {3,5,6}`.  The top Boolean face has charge zero but coefficient `+1`,
so the exact convention needed here is

```text
chiHat_7(0)=+1,                 chiHat_7(q)=chi_7(q) for q!=0.  (7)
```

`chiHat_7` is not claimed to be a multiplicative character at zero.  It is
the unique one-point augmentation used by the Boolean top cell.

For every nonempty *squarefree* subset `S` of `{A,B,C}`,

```text
(-1)^(|S|+1)=chiHat_7(sum_(s in S) q(s)).                   (8)
```

Indeed:

| mask | face subset | charge | `chi_7` | Mobius coefficient |
|---:|---|---:|---:|---:|
| 1 | `A` | 1 | +1 | +1 |
| 2 | `B` | 2 | +1 | +1 |
| 3 | `AB` | 3 | -1 | -1 |
| 4 | `C` | 4 | +1 | +1 |
| 5 | `AC` | 5 | -1 | -1 |
| 6 | `BC` | 6 | -1 | -1 |
| 7 | `ABC` | 0 | 0, augmented to +1 | +1 |

Thus mask order `1,...,7` gives exactly `++-+--+`.  Equivalently, if
`U_A=X,U_C=Y,U_B=Z`, THM-830's recurrence can be written

```text
T_n = sum_(empty!=S subset {A,B,C})
          chiHat_7(q(S)) U_S T_(n-|S|),                    (9)
```

for `n>=6`, where `U_S=product_(s in S) U_s`.  Singletons are the Fano/Paley difference
set, pairs are its complementary nonresidues, and the triple is neutral
because `1+2+4=0 mod 7`.  This proves the Fano lead in HYP-6895 at the level
of the actual THM-830 face operations.  It does not identify `chi_7` with a
tiling-parity or blue/black invariant, which HYP-6895 correctly refuted.

## 3. Cyclotomic specialization and exact period seven

Let

```text
F_n=T_n(zeta,zeta^4,zeta^2)
   =h_(n-3)(zeta,zeta^2,zeta^4),                           (10)
eta=zeta+zeta^2+zeta^4=(-1+i sqrt(7))/2,
etaBar=zeta^3+zeta^5+zeta^6=(-1-i sqrt(7))/2.              (11)
```

The variable order in (10) is `A,C,B`, as in (2).  The elementary symmetric
functions of the three quadratic-residue roots are `eta,etaBar,1`, so

```text
F_n=eta F_(n-1)-etaBar F_(n-2)+F_(n-3),                    (12)
F_3=1,                    F_4=eta,             F_5=-1.
```

This recurrence is exactly seven-periodic:

```text
F_(n+7)=F_n,
(F_3,...,F_9)=(1,eta,-1,etaBar,1,0,0).                    (13)
```

**Proof.**  With `Q={1,2,4}` and `N={3,5,6}`, the seventh-root factorization
gives

```text
(1-t) product_(d in Q)(1-zeta^d t) product_(d in N)(1-zeta^d t)=1-t^7.
```

Since

```text
product_(d in N)(1-zeta^d t)=1-etaBar t+eta t^2-t^3,
```

we obtain the stronger generating identity

```text
sum_(r>=0) h_r(zeta,zeta^2,zeta^4)t^r
 = (1+eta t-t^2+etaBar t^3+t^4)/(1-t^7).                  (14)
```

Equations (12)-(13) follow coefficientwise.  QED.

## 4. Every carrier has an equitable seven-charge partition

For `r=n-3`, define

```text
c_r(q)=#{(x,y,z) in Z_(>=0)^3:
         x+y+z=r and x+4y+2z=q mod 7},
M_r=binom(r+2,2)=|S_n|.                                   (15)
```

The complete all-size census depends only on `r mod 7`.  In the table,
`c_R` is the value at each one of `1,2,4` and `c_N` the value at each one of
`3,5,6`.

| `r mod 7` | `c_r(0)` | `c_R` | `c_N` |
|---:|---:|---:|---:|
| 0 | `(M_r+6)/7` | `(M_r-1)/7` | `(M_r-1)/7` |
| 1 | `(M_r-3)/7` | `(M_r+4)/7` | `(M_r-3)/7` |
| 2 | `(M_r-6)/7` | `(M_r+1)/7` | `(M_r+1)/7` |
| 3 | `(M_r-3)/7` | `(M_r-3)/7` | `(M_r+4)/7` |
| 4 | `(M_r+6)/7` | `(M_r-1)/7` | `(M_r-1)/7` |
| 5 | `M_r/7` | `M_r/7` | `M_r/7` |
| 6 | `M_r/7` | `M_r/7` | `M_r/7` |

In particular:

1. The seven fibre sizes always differ by at most one.
2. They are exactly uniform iff `n` is congruent to `1` or `2` modulo `7`.
3. The tiling cube has the literal gauge-dependent fixed-path factorization

```text
X_n ~= product_(q in F_7) F_2^(c_(n-3)(q)).                (16)
```

4. Its seven-block Hamming enumerator is

```text
sum_(t in X_n) product_q s_q^|t cap S_(n,q)|
  =product_q (1+s_q)^c_(n-3)(q).                          (17)
```

The small carrier rows are:

| `n` | `r` | `M_r` | `c(0)` | per QR | per NQR |
|---:|---:|---:|---:|---:|---:|
| 3 | 0 | 1 | 1 | 0 | 0 |
| 4 | 1 | 3 | 0 | 1 | 0 |
| 5 | 2 | 6 | 0 | 1 | 1 |
| 6 | 3 | 10 | 1 | 1 | 2 |
| 7 | 4 | 15 | 3 | 2 | 2 |
| 8 | 5 | 21 | 3 | 3 | 3 |
| 9 | 6 | 28 | 4 | 4 | 4 |
| 10 | 7 | 36 | 6 | 5 | 5 |
| 11 | 8 | 45 | 6 | 7 | 6 |
| 12 | 9 | 55 | 7 | 8 | 8 |
| 13 | 10 | 66 | 9 | 9 | 10 |
| **14** | **11** | **78** | **12** | **11** | **11** |

Thus the order-fourteen staircase has one neutral block of dimension 12 and
six nonzero blocks of dimension 11:

```text
X_14 ~= F_2^12 x (F_2^11)^6.                              (18)
```

**Proof of the census.**  The finite Fourier transform of (15) is

```text
sum_q c_r(q) zeta^(kq)
  =h_r(zeta^k,zeta^(2k),zeta^(4k)).                        (19)
```

At `k=0` this is `M_r`.  Multiplication by a quadratic residue permutes `Q`,
so (19) is the value in (13); multiplication by a nonresidue sends `Q` to
`N=-Q`, so it is the conjugate value.  Fourier inversion, using

```text
sum_(k in Q) zeta^k=eta,       eta+etaBar=-1,
eta etaBar=2,
```

gives the seven rows displayed above.  Equations (16)-(17) then follow from
the independent binary coordinate cube.  QED.

## 5. The first repeated-role neutral layer

Ignoring the depth-zero apex, there is no neutral operation address at depths
one or two.  At depth three the only neutral address is

```text
(alpha,beta,gamma)=(1,1,1),                 ABC.           (20)
```

At depth four, the neutral fibre first acquires three anisotropic addresses:

```text
(3,0,1), (1,3,0), (0,1,3),
A^3 C,    A B^3,   B C^3.                                  (21)
```

They are the three congruences `3*1+4=0`, `1+3*2=0`, and
`2+3*4=0 mod 7`.  The whole depth-four census is `3` at charge zero and `2`
at every nonzero charge.

Equation (21) is a proved statement about repeated face operations.  It is
only a target for the LRC `j=4` tree: no functor from four far-speed peel
histories to `A/B/C` face words has been defined, so the shared number four
does not itself transfer the statement.

## 6. Tournament Analysis and information loss

The seven nonempty Boolean face strata, not runners or arcs, form the natural
vertices here.  The map

```text
S -> q(S)=sum_(s in S)q(s)                                  (22)
```

is a bijection to `F_7`.  For two strata define the pairwise observable
`Delta(S,T)=q(T)-q(S)`.  The Legendre switch orients

```text
S -> T  iff  Delta(S,T) in {1,2,4}.                         (23)
```

This is exactly the Paley tournament `P_7`.  There are no ties; the declared
tie path `0->1->2->3->4->5->6` is itself Hamiltonian because every step is
`+1`.  Its fingerprints are

```text
score histogram {3:7},       directed 3-cycles 14,
SCC sizes [7],                Hamiltonian paths 189.        (24)
```

Multiplying every charge by a quadratic residue is a gauge automorphism and
flips no edges; multiplying by a nonresidue takes the converse and reverses
all 21 edges.

The assumption challenge considered runners, chords, gaps, fixed sections,
section boundaries, wall events, residues, Fourier modes, matroid circuits,
and proof obligations as vertices.  Face strata were chosen because (22)
preserves exactly the B3 Mobius coefficient and the operation character (6).
It destroys the exponent triple inside a charge fibre, the raw tiling, the
Hamiltonian-path presentation after node quotient, the seam word, and the
blue/black relation.  In particular, reflection swaps `x,y`, and

```text
q(sigma(x,y,z))=y+4x+2z                                    (25)
```

is not a function of `q(x,y,z)` alone.  Hence (18) is a coordinate grading,
not an iso-class node invariant or an edge-colour classifier.  It does not
preserve the LRC predicate.

## 7. Transfer audit and precise backlog

The following separates the theorem from three useful but incomplete
cross-repository perspectives.

### 7.1 Farey violation ladders and toothpick growth

THM-841 proves that a Farey-gap endpoint of denominator `i` contributes the
nested multiplier thresholds `lambda/(mi)`.  On the dyadic subchain the exact
congruence

```text
2^h mod 7 = 1,2,4,1,2,4,...                                (26)
```

cycles through the singleton face charges `A,B,C`.  This is the only current
theorem-level overlap.  It does not intertwine the objects: doubling `m`
halves a geometric threshold, while a face lift increments an address; the
Farey denominators and opposite endpoint are absent from (5).  THM-466's
2-adic digits of `H` use a third variable, the odd-cycle collection index,
and must not be identified with either multiplier.

A precise probe, rather than a claimed reduction, is the dyadic character
ladder

```text
W_(zeta,i)^(2)(s)
  =sum_(h:2^h i<=k) zeta^(2^h) 1[s<lambda/(2^h i)].        (27)
```

Its label is three-periodic while its spatial thresholds contract by eight
per full label cycle.  Compute its breakpoint tree and factorial correlations
with the right-end analogue before comparing it with any toothpick recursion.
The full THM-841 divisor ladder includes every multiplier, so a dyadic match
alone would still be a subladder theorem.

The chi-labelled dyadic subladder does have an exact restricted self-affinity.
Writing `a=lambda/i` and truncating after depth `H`, one has

```text
W_a^H(s)=zeta 1_(s<a)+zeta^2 1_(s<a/2)+zeta^4 1_(s<a/4)
         +W_a^(H-3)(8s).                                  (27b)
```

This follows from `2^3=1 mod 7`. THM-841's full multiplier ladder has the
additional odd-reciprocal birth term, so (27b) is a self-similar stalk plus a
non-self-similar source, not a restoration of full toothpick self-similarity.

### 7.2 The `j=4` flood tail

THM-741 remains `CLAIMED`; the shared log records `171/2002` clean bodies, not
a complete aggregate run, and its `Q(E)=emptyset` flood bodies are the
expensive tail.  A separate exact audit does classify that tail before the
recursive sweeps.  There are exactly 21 flood bodies:

```text
E_(a,b)={8,9,10,11,12,13,14} union {a,b},
1<=a<b<=7.                                                  (27a)
```

They are naturally the edges of `K_7`.  Each edge belongs to the unique Fano
line `{a,b,a xor b}`, and the 21 edges split into seven triples.  The three
negative Boolean masks `{3,5,6}` are themselves one of those lines.  This is
an exact organizational bridge between the flood bank and the `chi_7` signs,
but not a quotient symmetry: among all 168 elements of `GL(3,2)`, only the
identity preserves any one of the exact body observables `r(E)`, `m(E)`, or
the THM-741 threshold `V1(E)` on all 21 rows.

The exact edge-module structure is sharper than the stabilizer test. Write a
flood body as edge `{a,b}` of `K_7`. Then

```text
r_(a,b)=x_a+x_b,           x=(16,16,14,14,10,14,6).        (27c)
```

Consequently the `GL(3,2)` orbit span of `r` has rank 7 and centered rank 6.
The orbit spans of `m` and `V1` each have the full edge-space rank 21 and
centered rank 20. Concrete non-additive curls are

```text
m_12+m_67-m_16-m_27=1/21,
V1_13+V1_56-V1_15-V1_36=48.                               (27d)
```

The seven point-star rows have rank 7, the seven Fano-triangle rows have rank
7, and together they have rank 13, with only their common total shared. Thus
point and Fano marginals leave an eight-dimensional edge field invisible;
both `m` and `V1` have nonzero components there. The natural symmetry object
is the defect cocycle `d_g(e)=w(ge)-w(e)`, not a nonexistent Fano quotient.

There is also a local needle obstruction.  Exhausting all weak switch orders
of three one-sided Farey endpoint needles shows that such a path visits at
most four Boolean masks and at most two points of the negative line, including
coincident switches.  Three endpoint needles therefore cannot realize the
full seven-mask carrier or its negative Fano line.  A whole speed predicate
may contribute from both ends of a cyclic gap, so this is an obstruction for
the endpoint-needle model, not for arbitrary triples of speeds.

Equations (20)-(21) still say only that operation depth four is the first
depth with neutral repeated-role words beyond the balanced `ABC` word.  The
Fano edge partition organizes the 21 flood inputs but does not close their
exact interval sweeps.

The falsifiable next computation is to define, from the exact safe-interval
event complex, a canonical role map

```text
tau: (body, ordered peel history) -> {A,B,C}^*.             (28)
```

It must commute with subtree restriction and be independent of arbitrary
component ordering.  Only after those tests pass should the THM-741 rows be
stratified by `(alpha,beta,gamma)`, charge, `Q(E)=emptyset`, threshold, and sweep
count.  Charge alone must not be retained: (21) already shows three distinct
neutral words, and (25) shows that one charge loses the mirror channel.

### 7.3 Kakeya needles and a full operation Radon deck

At fixed depth `r`, the operation exponents reduce to the affine plane

```text
H_r={alpha+beta+gamma=r} subset F_7^3.                     (29)
```

The charge is one nonzero linear projection on this plane.  Its kernel
direction is exactly

```text
(2,4,1),       2+4+1=0,       1*2+2*4+4*1=0 mod 7.        (30)
```

Consequently the seven charge fibres are one pencil of parallel affine lines.
A finite-field Kakeya set contains a line in every projective direction, so
(29)-(30) are explicitly *not* a Kakeya reduction.  They explain the lesson
of the S659 finite-field scout: one scalar projection saturates early and
forgets position along its fibres, just as raw support forgets pinned and
concurrency data.

There is an exact rank theorem behind the proposed upgrade. Identify `H_r`
with `F_7^2`, first regard `f` as a function on that residue plane, and let
`P^1(F_7)` index the eight nonzero linear forms up to scale.  For a direction
`v` define its seven-offset Radon family

```text
R_v f(s)=sum_(x:v dot x=s) f(x),             s in F_7.      (31)
```

Over any characteristic-zero coefficient field, a deck of `d>=1` distinct
directions has

```text
rank {R_v:v in D}=1+6d,       kernel dimension=48-6d.      (32)
```

Indeed, the one-dimensional Fourier transform in `s` gives the Fourier slice
identity `hat(R_v f)(k)=hat f(kv)`.  Each direction reveals its six nonzero
frequency points, distinct projective directions intersect only at frequency
zero, and every deck shares the same total sum.  Formula (32) follows.  In
particular, the one `chi_7` pencil has rank 7 and a 42-dimensional kernel;
seven directions still leave a six-dimensional kernel; all eight directions
are necessary and sufficient to reconstruct an arbitrary function on the
residue operation plane.

There is a second, essential rank boundary before applying this theorem to
the integer address simplex. Let

```text
S_r={(alpha,beta,gamma) in Z_(>=0)^3: alpha+beta+gamma=r}.
```

For a canonical residue `xi=(a,b,c) in {0,...,6}^3` with
`a+b+c=r mod 7`, put

```text
L_xi=(r-a-b-c)/7.                                         (33)
```

Reduction modulo seven has the exact odometer decomposition

```text
S_r ~= disjoint-union_(xi in H_r) Comp_3(L_xi),
(alpha,beta,gamma)=xi+7(i,j,k),                            (34)
```

where the fibre is empty for `L_xi<0` and otherwise has size
`binom(L_xi+2,2)`. Addition of A, B, or C increments the corresponding
residue; a `6->0` wrap increments the corresponding carry coordinate. Thus
the carry fibre is itself a smaller three-part composition simplex with the
same B3 recursion.

At `n=14`, `r=11`, canonical residue sums are 4, 11, or 18. Consequently

```text
15 residue fibres have size 3,    33 have size 1,
(6,6,6) is the one empty residue.                         (35)
```

Let `P:K^S_11 -> K^H_4` be scalar residue pushforward, summing address values
within each residue fibre. Over characteristic zero,

```text
rank(P)=48,             dim ker(P)=30=15(3-1).             (36)
```

The all-direction Radon map `R` has rank 49 on the plane, so

```text
rank(RP)=48,            ker(RP)=ker(P).                    (37)
```

The complete scalar deck reconstructs the residue-aggregated value `Pf`, not
the original 78 address values. On each of the 15 triple fibres the missing
label is exactly which coordinate carries the extra seven: A, B, or C. The
30-dimensional alias kernel consists of the two independent contrasts in
each such triple. Explicitly, for every canonical `xi` with coordinate sum 4,
a basis is

```text
1_(xi+7e_A)-1_(xi+7e_C),
1_(xi+7e_B)-1_(xi+7e_C).                                  (37a)
```

Face triality gives a further exact organization of the eight directions:

```text
coordinate orbit: {(1,0),(0,1),(1,1)}                     size 3,
mixed orbit:      {(1,2),(1,4),(1,6)}                     size 3,
chi orbit:        {(1,3),(1,5)}                           size 2.             (38)
```

The coordinate-pencil line loads are, up to offset permutation,
`(17,15,13,11,9,7,6)`. Every mixed or chi pencil has loads
`(12,11,11,11,11,11,11)`. The charge triple `(1,2,4)` projects to direction
`(1,3)`, and its affine offset places the 12-address line at charge zero.

After splitting the 78 addresses into four channels, namely no carry and the
A/B/C carry, the same reflection-invariant six-direction set

```text
{(1,0),(0,1),(1,2),(1,4),(1,3),(1,5)}                    (38a)
```

has channel ranks `(33,15,15,15)` and total rank 78. Five directions cannot
work because their rank on the 33-point no-carry support is at most 31. Exact
enumeration shows (38a) is the unique reflection-invariant six-subset with
full channel rank. The formally triality-invariant six-set, coordinate orbit
union mixed orbit, has ranks `(32,15,15,15)` and misses one mode. Thus true
endpoint reflection permits a minimal six-pencil repair, while insisting on
formal face triality requires all eight directions.

The concrete metagraph target is a full finite Radon deck.  For every
direction in `P^1(F_7)`, retain all line sums of the address-decorated tiling
tensor on (29), together with upper/lower node classes, blue/black colour,
and pair-groupoid defect. At `r=11`, also retain the A/B/C carry role in every
triple residue fibre. The full scalar direction deck is information-complete
only after residue pushforward; equations (35)-(37) prove that dropping carry
loses 30 address dimensions. A structured metagraph tensor may occupy a
proper subspace, so the open test is its restricted rank: which directions
and carry contrasts match the `Omega+S2` classification through `n=8` and
remain closed recursively?

This Radon construction acts on Parikh/exponent vectors. It cannot recover an
ordered j=4 peel chronology unless the proposed role map (28) is first proved
to commute with the relevant history operations. Word-to-Parikh is sound for
THM-830's commuting face monoid; no corresponding commutation theorem exists
for THM-741 histories.

## Verification record

The independent referee uses exact arithmetic in
`Z[eta]/(eta^2+eta+2)`.  It checks:

- all seven signs with the zero augmentation stated explicitly;
- the recurrence and period through depth 140;
- the census formula through depth 399, including equitability;
- the neutral words at depths three and four;
- the fixed-depth kernel direction (30);
- all Paley fingerprints in (24), with Hamiltonian paths by subset DP.

The independent flood/needle audit additionally checks all 2,002 nine-body
subsets, the 21-body characterization (27a), the seven Fano edge triples, all
168 Fano automorphisms against `r,m,V1`, and all 104 distinct weak-order
three-needle mask paths. It also verifies the orbit ranks, point/Fano incidence
ranks, and curls (27c)-(27d).

The Radon/carry audit checks every subset of the eight directions, obtaining
plane rank `1+6d` and address rank `1+6d` for `d<=7`, followed by the strict
`49/48` plane/address split at `d=8`. It independently verifies (35)-(38) and
the exact A/B/C carry repair of all 78 addresses, including the minimal and
unique reflection-invariant six-pencil deck (38a).

Neither script closes a THM-741 sweep or tests an LRC transfer, and none is
asserted here.
