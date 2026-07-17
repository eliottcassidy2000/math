---
id: THM-969
title: Scale-nine Hamming-six two-orbit obligation-nerve obstruction
status: PROVED FINITE-EXACT — the complete 482,294,736-context primitive proper AP-centred common-scale-nine Hamming-six sheet bank is empty; exact reductions leave 64 all-order-nine signed cycles and one 12-context mixed orbit, and their owner-obligation nerves contain incompatible pairs
source: codex-2026-07-17 scale-nine exact C++ certificate and independent Python referee
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860]
related: [THM-862, THM-957, THM-958, THM-960, THM-962, THM-963, THM-970, THM-974, HYP-6820]
verification:
  - 04-computation/lrc13_scale_nine_hamming_six_frontier_scout_codex_c9.cpp
  - 05-knowledge/results/lrc13_scale_nine_hamming_six_frontier_scout_codex_c9.out
  - 04-computation/lrc13_scale_nine_hamming_six_referee_codex_c9.py
  - 05-knowledge/results/lrc13_scale_nine_hamming_six_referee_codex_c9.out
---

# THM-969 — scale nine is killed by two sparse obligation nerves

Let `R subset F_13^*` have six elements, put `P=[12] minus R`, and consider

```text
A=9P union {w_r:r in R},
w_r=9r (mod 13),              w_r>0,              w_r!=9r.       (1)
```

Assume that `A` is primitive and `M(A)<=1/13`.  For each replacement put

```text
D_r=9/gcd(9,w_r),
e_r=(D_r w_r/9) mod D_r.                                      (2)
```

Then no such packet exists.  Equivalently, the complete primitive proper
AP-centred common-scale-nine Hamming-six sheet bank is empty.  Every primitive
proper packet in (1) therefore has `M(A)>1/13`, and no metric-height recursion
is born on this face.

## 1. State grammar and exact size of the raw bank

THM-860 supplies common-sheet coverage at each replacement owner and the
hereditary leave-one-out law

```text
lcm(D_s:s!=r)=9                         for every r in R.         (3)
```

The complete effective state alphabet is

```text
(D,e)=(1,0),(3,1),(3,2),(9,1),(9,2),(9,4),(9,5),(9,7),(9,8).   (4)
```

Because `9=3^2`, (3) is equivalent to at least two of the six coordinates
having effective order nine.  There are therefore

```text
3^6-2^6-6*2^5 = 473                                             (5)
```

admissible order words, and

```text
9^6-3^6-6*6*3^5 = 521,964                                      (6)
```

admissible order/unit words per support.  Across the
`binom(12,6)=924` supports, the exact labelled banks have

```text
order contexts        924*473     =     437,052,
order/unit contexts   924*521,964 = 482,294,736.                (7)
```

No random, floating-point, asymptotic, or height-cutoff step occurs in these
counts.

## 2. Literal sheet masks and scalar capacity

For owner `o`, provider `r`, and state `(D,e)`, let `u` be the least CRT
representative satisfying

```text
u=Dr (mod 13),                    u=e (mod D).                   (8)
```

The provider covers sheet `ell in Z/9Z` exactly when the centred representative
of

```text
u(o^(-1)+13ell)                  modulo 13D                     (9)
```

lies in the half-open interval `(-D,D]`.  The certificate reconstructs every
mask directly from (8)--(9).  Its stored output freezes the full owner-one
table for all nine states and all twelve ratios.

Mask cardinality is independent of the unit `e`.  In owner-normalized ratio
coordinates the three order rows are

```text
D=1:  size 9 at ratio 1;                         size 0 otherwise;
D=3:  size 3 at ratios {1,4,5,8,9};              size 0 otherwise;
D=9:  size 2 at ratios {1,2,5,8,11};             size 1 otherwise. (10)
```

Thus every owner must first satisfy the scalar capacity inequality

```text
sum_(r in R) |E_(D_r,e_r)(r,o)| >= 9.                         (11)
```

Scanning all `437,052` labelled order contexts leaves exactly `1,186`:

| order multiplicities | scalar-capacity contexts |
|---|---:|
| `1^0 3^0 9^6` | 82 |
| `1^0 3^2 9^4` | 474 |
| `1^0 3^3 9^3` | 132 |
| `1^0 3^4 9^2` | 330 |
| `1^1 3^3 9^2` | 168 |

They lie on `316` supports.  Every other order context is impossible before
units are considered.

## 3. Owner-local reduction: exactly two structural banks

For each surviving order context, relax global compatibility by allowing a
different unit word at each owner.  Exact dynamic programming in the
`512`-element union-mask semilattice decides whether that owner can cover all
nine sheets.  Requiring this owner-local condition at all six owners leaves
only

```text
64 contexts of type D^6=9^6,
12 contexts of type D^6=3^2 9^4.                               (12)
```

The `76` contexts use `76` distinct supports.  They have exact structural
descriptions.

### All-order-nine bank

The 64 supports are precisely the signed-doubling six-cycles

```text
p -> o             iff             o/p in {2,-2} mod 13.       (13)
```

The relation is a directed Hamiltonian cycle on the support.  Writing

```text
v_(i+1)=(-1)^(a_i) 2v_i,
xor_i a_i=1,                                                   (14)
```

the usual rooted count gives `12*2^5/6=64`.  The finite reduction verifies
that these 64 cycles, and no other all-order-nine supports, are owner-locally
feasible.

### Mixed bank

The remaining 12 contexts form one free multiplicative orbit.  For a unique
`a in F_13^*`, their order-labelled supports are

```text
D=3 at {a,5a},
D=9 at {2a,-2a,3a,-3a}.                                      (15)
```

This is the first common-scale bank in the `c=4,...,9` run whose terminal
owner-local frontier is not solely the recurring 64-cycle family.  Keeping
the order colouring on the support is essential: the uncoloured six-set does
not preserve (15).

The equalities in (12)--(15) are finite-exact classifications from the literal
mask table.  They are not presented as an independent hand classification of
all 473 order words.

## 4. The faithful object is the owner-obligation family

Fix one of the 76 order-labelled contexts.  For each owner `o`, let

```text
O_o = {global unit words whose literal masks cover all nine sheets at o}. (16)
```

The unit word is shared across owners in (16).  Hence a common-sheet packet
would be a point of `intersection_(o in R) O_o`.

### The signed-cycle nerve is `3K_2`

For every all-order-nine signed cycle, exact literal evaluation gives

```text
|O_i| = 432,
|O_i intersect O_j| = 12   if j=i+3 (mod 6),
                         0   otherwise.                         (17)
```

No triple meets.  Thus the complete intersection nerve consists of the three
cycle-antipode edges:

```text
3K_2.                                                          (18)
```

The whole `6^6=46,656` unit fibre has satisfaction profile

```text
unit words satisfying 0 owners: 44,100,
unit words satisfying 1 owner :  2,520,
unit words satisfying 2 owners:     36,
all higher counts             :      0.                         (19)
```

In particular, adjacent owner obligations are disjoint, so the six-fold
intersection is empty.

### The mixed D=9 sub-nerve is `2K_2`

For the context (15), the two order-three owner obligations have size `1,152`
and the four order-nine owner obligations have size `144`.  Within the four
order-nine owners,

```text
O_(2a)=O_(-2a),             O_(3a)=O_(-3a),
O_(epsilon 2a) intersect O_(delta 3a)=empty
                         for epsilon,delta in {+1,-1}.          (20)
```

Indeed, each antipodal intersection has size `144`, while every cross-pair
intersection is zero.  The D=9 induced pair nerve is therefore `2K_2`, which
already kills the six-fold intersection.  For completeness, the D=3/D=3
intersection has size `192`; the eight D=3/D=9 intersections split four of
size `24` and four of size `48`.  The `2^2 6^4=5,184` unit fibre has profile

```text
unit words satisfying 0 owners: 2,928,
unit words satisfying 1 owner : 1,776,
unit words satisfying 2 owners:   336,
unit words satisfying 3 owners:   144,
all higher counts             :     0.                          (21)
```

Equations (17)--(21) are finite-exact identities checked for every structural
context.  Once those identities are known, the final nerve contradiction is
human-readable and does not require inspecting individual unit words.

## 5. Complete reduced literal certificate and independent referee

The owner-local step is a necessary relaxation: if even an owner-dependent
unit choice cannot cover one owner, no globally shared unit word can work.
It is therefore exact to discard the other `437,052-76` order contexts.

The C++ certificate then exhausts every state fibre above the 76 remaining
contexts:

```text
64*6^6                 = 2,985,984 all-order-nine words,
12*2^2*6^4             =    62,208 mixed words,
total                  = 3,048,192 literal words,
global common-sheet survivors = 0.                             (22)
```

This is an exact reduced literal certificate, not a claim that the executable
blindly loops over all 482,294,736 raw contexts.  Completeness comes from the
exact scalar and owner-local eliminations followed by (22).

An independently organized Python referee reconstructs the CRT masks, the
473 divisor words, scalar capacity, set-valued owner-local dynamic programs,
the 76-context classification, all 3,048,192 packed literal fibres, and the
nerve profiles.  It reproduces the zero verdict.  Two stable cross-language
payload invariants are

```text
literal mask digest       e9727a1660bc9b433ed1302673521a72ba12f3509e34dd8879372efd303438c2
owner-local bank digest   b0dbd8e63873c37f361704cfbc6c22cb05729f3e09cf5531f04510405d2fc695. (23)
```

## 6. Tournament analysis and challenged vertex choices

On a six-label support, take the pairwise observable `o/p`.  Declare the
sparse binary relation (13); the sign in `{+2,-2}` is the switch/edge colour,
and absent pairs are ties.  Complete ties along the lexicographically first
Hamiltonian path of the absence graph.

For the 64 signed-cycle supports this reproduces the familiar telemetry

```text
joint fingerprints             5,
score multiset                 (1,2,2,3,3,4),
directed triangles             6,
SCC profile                    (6),
sparse-cycle flip histogram    {2:8,3:52,4:4},
Hamiltonian-path histogram     {29:32,31:20,37:12}.             (24)
```

The 12 mixed supports give six joint fingerprints; their directed-triangle
histogram is `{0:1,1:4,3:2,4:2,6:3}` and their SCC histogram is
`{1+1+1+1+1+1:1,3+1+1+1:4,5+1:2,6:5}`.

These completions are telemetry, not the proof.  They preserve the signed
ratio scaffold but destroy the divisor colouring, literal sheets, obligation
sizes, and pair-intersection multiplicities.  In particular they cannot see
either `3K_2` or the mixed D=9 `2K_2`.

The vertex-choice audit is therefore:

```text
runner labels:       preserve ratios and signed C6; destroy units and coverage;
residue/mask rows:   preserve literal local coverage; remain larger than needed;
sheet vertices:      preserve one owner's cover; destroy shared-unit compatibility;
owner obligations:  preserve exactly the global compatibility predicate;
nerve edges:         destroy witness identities but preserve the zero obstruction. (25)
```

The challenged assumption is that tournament vertices should be runners.
Here the useful vertices are proof obligations.  Their intersection relation
is not naturally a tournament at all; forcing a completion discards the
decisive information.

## 7. Cross-scale structural sharpening

On the recurring signed-doubling support, the owner-obligation pair graph has
now evolved as follows:

```text
c=5: every owner pair meets, but eight owner triples fail;
c=6: K_6 minus the three antipodal edges;
c=7: the signed cycles die earlier by the row-product character;
c=8: K_3,3;
c=9: 3K_2;
c=10: K_6 minus C_6 (the triangular prism).                    (26)
```

Thus the sparse signed cycle is only the carrier.  The actual scale-dependent
object is the intersection nerve of the owner constraints.  Scale nine also
adds the mixed orbit (15), showing that a support-only recurrence would miss a
new divisor-coloured component.

## 8. Verification and formalization status

Reproduce the stored outputs with

```bash
c++ -O3 -std=c++20 -Wall -Wextra -pedantic \
  04-computation/lrc13_scale_nine_hamming_six_frontier_scout_codex_c9.cpp \
  -o /tmp/lrc13-scale-nine-c9
/tmp/lrc13-scale-nine-c9 | \
  cmp - 05-knowledge/results/lrc13_scale_nine_hamming_six_frontier_scout_codex_c9.out

c++ -O0 -std=c++20 -Wall -Wextra -pedantic \
  04-computation/lrc13_scale_nine_hamming_six_frontier_scout_codex_c9.cpp \
  -o /tmp/lrc13-scale-nine-c9-O0
/tmp/lrc13-scale-nine-c9-O0 | \
  cmp - 05-knowledge/results/lrc13_scale_nine_hamming_six_frontier_scout_codex_c9.out

python3 04-computation/lrc13_scale_nine_hamming_six_referee_codex_c9.py | \
  cmp - 05-knowledge/results/lrc13_scale_nine_hamming_six_referee_codex_c9.out
python3 -O 04-computation/lrc13_scale_nine_hamming_six_referee_codex_c9.py | \
  cmp - 05-knowledge/results/lrc13_scale_nine_hamming_six_referee_codex_c9.out
```

An AddressSanitizer/UndefinedBehaviorSanitizer build is also byte-identical.
The frozen SHA-256 values are

```text
7aa55fef8f99903b7fad531f02de47458189e5687f9ec5408f64fa5749b2c79b  C++ source
eedac03346bb2fd21f37510f41f0ac92b8aa9bceda37dda5014539cb3e1bdc6e  C++ output
9fbe8bb946560545e4c692b27afd68e0ae36b08b0dd4cb3eaabc51bb71334a31  Python source
83e7bc3efcf55646896b2604bf04546dd3329c97c6229a58c10fd8e01fd9fc55  Python output. (27)
```

`LRCScaleNineOwnerNerve.lean` now kernel-checks the terminal combinatorial
step.  It defines exact pair and induced-subfamily nerves; verifies by ordinary
`decide` that the all-order-nine table is `3K2` and the mixed order-nine
subtable is `2K2`; and proves abstractly that either reported nerve forces the
total owner intersection empty.  There is no `sorry` or `native_decide`; the
axiom audit is at most the standard foundational trio.

The module deliberately does not recheck the preceding nine-state mask table,
the `482,294,736`-to-`76` reduction, or the literal pair counts.  Those remain
the frozen C++/Python completeness certificate.  The next formalization layer
would be the mask-to-exact-nerve bridge, not another brute-force terminal
theorem.

## Scope guardrail

This theorem closes only the primitive proper AP-centred common-scale-nine
Hamming-six face under the `M(A)<=1/13` hypothesis that invokes THM-860's
common-sheet and hereditary-lcm conclusions.  This theorem alone does not
close common scales `10<=c<=840`, arbitrary non-AP Hamming-six packets, the remaining ramified
Hamming-five bank, deep-sheet languages outside this chart, the seven-wall
continuum-to-grid problem, the block-partition trichotomy, or global LRC(14).

THM-970 and THM-974 subsequently close `c=10,11`.  The next AP-centred
common-scale Hamming-six frontier is `c=12`.  ∎
