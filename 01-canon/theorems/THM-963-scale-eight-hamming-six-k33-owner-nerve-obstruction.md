---
id: THM-963
title: Scale-eight Hamming-six K3,3 owner-nerve obstruction
status: PROVED STRUCTURAL + FINITE-EXACT — all 215,728,128 hereditary-lcm-compatible AP-centred c=8 label/state contexts fail necessary common-sheet coverage; owner-local feasibility leaves the 64 signed-doubling supports, whose owner-obligation nerve is K3,3 and therefore has no six-fold intersection
source: codex-2026-07-17-S64 scale-eight divisor, literal-mask, and owner-nerve audit
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860]
related: [THM-861, THM-862, THM-957, THM-958, THM-960, THM-962, HYP-6820]
verification:
  - 04-computation/lrc13_scale_eight_hamming_six_frontier_scout_codex_c8.cpp
  - 05-knowledge/results/lrc13_scale_eight_hamming_six_frontier_scout_codex_c8.out
---

# THM-963 — scale eight is killed by a `K_{3,3}` owner nerve

Let `R subset F_13^*` have six elements, put `P=[12] minus R`, and consider

```text
A=8P union {w_r:r in R},
w_r=8r (mod 13),              w_r>0,              w_r!=8r.       (1)
```

Assume that `A` is primitive and `M(A)<=1/13`.  For each replacement put

```text
D_r=8/gcd(8,w_r),
e_r=(D_r w_r/8) mod D_r.                                      (2)
```

Then no such packet exists.  Equivalently, the complete primitive proper
AP-centred common-scale-eight Hamming-six sheet bank is empty.  Every primitive
proper packet in (1) therefore has `M(A)>1/13`, and no metric height recursion
is born.

## 1. Exact divisor grammar and hereditary lcm

THM-860 supplies common-sheet coverage at every replacement owner and the
hereditary leave-one-out law

```text
lcm(D_s:s!=r)=8                         for every r in R.         (3)
```

The complete effective state alphabet is

```text
(D,e)=(1,0),(2,1),(4,1),(4,3),(8,1),(8,3),(8,5),(8,7).          (4)
```

Since `8=2^3`, equation (3) is equivalent to saying that at least two of the
six coordinates have effective order eight.  Indeed, deleting any one
coordinate must leave another carrier of the maximal `2^3` valuation; the
converse is immediate.  Thus there are

```text
sum_(k=2)^6 binom(6,k) 3^(6-k) = 1,909                         (5)
```

order words.  There are four order-eight unit states and four states of
smaller order in (4), so the exact state-word count is

```text
sum_(k=2)^6 binom(6,k) 4^6 = 233,472.                           (6)
```

The two distributions are

| number `k` of order-eight coordinates | 2 | 3 | 4 | 5 | 6 |
|---:|---:|---:|---:|---:|---:|
| order words | 1,215 | 540 | 135 | 18 | 1 |
| order/unit words | 61,440 | 81,920 | 61,440 | 24,576 | 4,096 |

Across the `binom(12,6)=924` supports this gives

```text
labelled order contexts       924*1,909   =   1,763,916,
labelled order/unit contexts  924*233,472 = 215,728,128.        (7)
```

No asymptotic, floating-point, or random step occurs in these counts.

## 2. Literal local mask table

For an owner `o`, provider `r`, and state `(D,e)`, let `u` be the least CRT
representative satisfying

```text
u=Dr (mod 13),                    u=e (mod D).                    (8)
```

The provider covers sheet `ell in Z/8Z` exactly when the centred
representative of

```text
u(o^(-1)+13ell)                  modulo 13D                      (9)
```

lies in the half-open interval `(-D,D]`.  This is THM-860's common-sheet
mask specialized to `c=8`.

Normalize `o=1` and encode sheets `0,...,7` as an eight-bit hexadecimal
mask.  The complete table is

```text
(1,0):  1 -> ff; all other ratios -> 00.

(2,1):  1 -> aa; 6,7 -> 55; all other ratios -> 00.

(4,1):  1 -> 88; 3,4 -> 44; 6,7 -> 22; 9,10 -> 11;
        all other ratios -> 00.
(4,3):  1 -> 88; 3,4 -> 11; 6,7 -> 22; 9,10 -> 44;
        all other ratios -> 00.

(8,1):  1->0c, 2->04, 3->02, 4->03, 5->01, 6,7->80,
        8->40, 9->60, 10->20, 11,12->10.
(8,3):  1->09, 2->01, 3->20, 4->24, 5->04, 6,7->80,
        8->10, 9->12, 10->02, 11,12->40.
(8,5):  1->48, 2->40, 3->02, 4->12, 5->10, 6,7->80,
        8->04, 9->24, 10->20, 11,12->01.
(8,7):  1->18, 2->10, 3->20, 4->60, 5->40, 6,7->80,
        8->01, 9->03, 10->02, 11,12->04.                        (10)
```

Multiplying every label by `o^(-1)` transports (10) to an arbitrary owner.
For fixed `D,r,o`, mask cardinality is independent of the unit `e`.  This
makes the following capacity reduction exact.

## 3. Capacity and owner-local reductions

At every owner, common-sheet coverage first requires the scalar inequality

```text
sum_(r in R) |E_(D_r,e_r)(r,o)| >= 8.                           (11)
```

Scanning the `1,763,916` labelled order contexts using the unit-independent
cardinalities in (10) leaves exactly

```text
3,166 scalar-capacity contexts.                                (12)
```

Next, relax global compatibility by allowing a different unit word at each
owner.  Exact dynamic programming in the 256-element mask semilattice tests
whether each owner is locally coverable.  This stronger necessary condition
leaves exactly

```text
64 owner-locally feasible order contexts.                       (13)
```

All 64 have

```text
D_r=8                         for every r in R.                  (14)
```

Thus every mixed order word dies before any global unit compatibility or
metric height is tested.

The all-order-eight scalar-capacity bank itself has 66 supports.  It splits
exactly as

```text
64 signed-doubling six-cycles
 2 quadratic cosets Q={1,3,4,9,10,12}, NQ={2,5,6,7,8,11}.       (15)
```

The two quadratic cosets already fail owner-local unit feasibility.  The
other 64 supports are precisely those on which

```text
p -> o                iff                o/p in {2,-2} mod 13    (16)
```

is a directed Hamiltonian six-cycle.  This is the same 64-support scaffold
that organized scales two and four through six.  As usual, a rooted signed
cycle has an odd number of negative edges because `2^6=-1 mod 13`, and the
count is `12*2^5/6=64`.

The equality in (15), the exclusion of every mixed order context, and the
quadratic-coset owner-local failure are finite-exact statements checked from
the literal masks.  They are not being presented as a separate hand proof of
the full `1,909`-word classification.

## 4. Symbolic all-order-eight owner grammar

Fix one of the 64 cycles and write

```text
v_(i+1)=(-1)^(a_i) 2v_i,                  xor_i a_i=1,            (17)
```

with indices in `Z/6Z`.  Encode the odd order-eight unit by

```text
e_i=2x_i+1,                                x_i in Z/4Z,            (18)
```

and put

```text
A_(i,d)=a_i xor a_(i+1) xor ... xor a_(i+d-1).                   (19)
```

At owner `v_i`, the providers at cyclic distances `0,...,5` contribute:

```text
distance 0: fixed sheet 3 plus one even sheet,
distance 1: one even sheet,
distance 2: one even and one odd sheet,
distance 3: one even sheet,
distance 4: one odd sheet,
distance 5: fixed sheet 7.                                     (20)
```

The two odd sheets are `1` and `5`.  Their binary symbols are

```text
b_2=(x_(i+2) mod 2) xor A_(i,2),
b_4=(x_(i+4) mod 2) xor A_(i,4).                                (21)
```

Divide the four nonfixed even sheets by two and define their symbols in
`Z/4Z` by

```text
q_0 = 1-x_i,

q_1 = 1-x_(i+1)       if A_(i,1)=0,
      x_(i+1)+2       if A_(i,1)=1,

q_2 = x_(i+2)         if A_(i,2)=0,
      3-x_(i+2)       if A_(i,2)=1,

q_3 = 3-x_(i+3)       if A_(i,3)=0,
      x_(i+3)         if A_(i,3)=1.                              (22)
```

The literal owner mask covers all eight sheets if and only if

```text
b_2 != b_4
and
{q_0,q_1,q_2,q_3}=Z/4Z.                                       (23)
```

This is the scale-eight faithful local object: one parity-complement
constraint coupled to a four-symbol all-different constraint.

Each map from its `x` variable to the corresponding `q` in (22) is a
permutation of `Z/4Z`.  Hence the all-different part has `4!` assignments.
For each such assignment, (21) permits two values of `x_(i+4)`, while
`x_(i+5)` is free.  Every owner obligation therefore has exactly

```text
4!*2*4=192 unit-word witnesses.                                (24)
```

The verifier checks (23) against the literal masks for all

```text
64 supports * 4^6 unit words * 6 owners = 1,572,864 owner rows. (25)
```

Thus (23)--(24) are not an inferred pattern from a sample.

## 5. The owner-obligation nerve is `K_{3,3}`

Let `O_i subset (Z/4Z)^6` be the 192 unit words satisfying owner `v_i`.
Exact intersections of the symbolic/literal obligations depend only on
cyclic distance:

```text
|O_i intersect O_j| =  8,   distance(i,j)=1,
                         0,   distance(i,j)=2,
                        16,   distance(i,j)=3.                   (26)
```

Consequently two owner obligations are compatible exactly when their cycle
positions have opposite parity.  Their pair-intersection graph is

```text
C_6 plus its three antipodal edges = K_{3,3}.                    (27)
```

In particular every triple of owners contains two positions of the same
parity, hence a distance-two incompatible pair.  No unit word satisfies
three owners, let alone all six.  The six distance-two pairs are the complete
minimal obstruction family.

Equation (26) also gives the whole per-support satisfaction profile without
further search.  There are six distance-one pairs with eight witnesses and
three antipodal pairs with sixteen witnesses, so the number of unit words
satisfying exactly two owners is

```text
6*8+3*16=96.                                                   (28)
```

The six owner sets contribute `6*192=1,152` incidences.  Since no triple
meets, exactly `1,152-2*96=960` words satisfy one owner, and

```text
unit words satisfying 0 owners: 3,040,
unit words satisfying 1 owner :   960,
unit words satisfying 2 owners:    96,
all other satisfaction counts  :     0.                         (29)
```

The verifier obtains the same profile independently for every one of the 64
supports.  Equations (26)--(29) prove that no all-order-eight global unit word
exists.  Together with Section 3, this empties the complete common-sheet bank.

## 6. Direct literal census

Independently of the nerve deduction within the same executable, the verifier
packs the eight sheets at six owners into 48 bits and scans every context in
(7).  It finds

```text
literal contexts tested            215,728,128,
literal common-sheet survivors               0,
surviving supports                            0,
surviving state words                         0.                 (30)
```

This direct pass uses the CRT definition (8)--(9), not the symbolic
all-different predicate.  The reduction and the literal pass therefore check
the same zero verdict through different internal representations, although
they live in one source file and are not claimed as two independently written
implementations.

## 7. Tournament and alternate-carrier audit

On the 64 owner-local supports, the ordered-pair observable is the ratio
`o/p`.  Declare

```text
p -> o iff o/p in {+2,-2};                                     (31)
```

the sign is the binary edge colour/switch, and all absent unordered pairs are
ties.  Complete the ties forward along the lexicographically first
Hamiltonian path of the absence graph.  The resulting tournaments have

```text
joint fingerprints             5,
score multiset                 (1,2,2,3,3,4),
directed triangles             6,
SCC profile                    (6),
sparse-cycle flip histogram    {2:8,3:52,4:4},
Hamiltonian-path histogram     {29:32,31:20,37:12}.              (32)
```

These are the same bare completion fingerprints seen at scale six.  They
recognize the recurring signed-cycle scaffold but do not decide scale eight:
completion forgets which owner pairs have empty intersection and forgets the
intersection multiplicities `8` and `16` in (26).

The challenged vertex choice is therefore decisive.  Runner labels alone,
or a completed tournament on them, are lossy.  Take the six **owner
obligations** `O_i` as vertices and pair intersection as the observable.  The
resulting `K_{3,3}` nerve preserves the common-sheet obstruction and is already
sufficient to prove emptiness.  The information ladder at scale eight is

```text
runner labels / signed C6
  -> edge-coloured ratio cycle
  -> parity plus all-different owner obligations
  -> K3,3 pair-intersection nerve.                              (33)
```

The last object is smaller and more faithful than the tournament completion.
This explicitly refutes the assumption that the natural tournament vertices
must be runners: proof obligations are the useful vertices here.

## 8. Cross-scale interpretation

Scale eight reconnects the two structures separated at scale seven:

```text
c=4,5,6: the 64 signed-doubling supports organize the hard bank;
c=7:     a row-product character kills those 64, while the two quadratic
         cosets survive to the terminal square-sum contradiction;
c=8:     scalar capacity sees 64 signed cycles plus those same two quadratic
         cosets; owner-local units kill the cosets, and the K3,3 nerve kills
         every signed cycle.                                    (34)
```

Thus the quadratic-coset pair from scale seven is not an isolated accident,
and the recurring 64-cycle bank is not merely a low-scale artefact.  Capacity
at scale eight contains both, while the owner-local/global split separates
them cleanly.

## 9. Verification and formalization status

Reproduce the stored output with either optimization level:

```bash
c++ -O3 -std=c++20 -Wall -Wextra -pedantic \
  04-computation/lrc13_scale_eight_hamming_six_frontier_scout_codex_c8.cpp \
  -o /tmp/lrc13-scale-eight-c8
/tmp/lrc13-scale-eight-c8 | \
  cmp - 05-knowledge/results/lrc13_scale_eight_hamming_six_frontier_scout_codex_c8.out

c++ -O0 -std=c++20 -Wall -Wextra -pedantic \
  04-computation/lrc13_scale_eight_hamming_six_frontier_scout_codex_c8.cpp \
  -o /tmp/lrc13-scale-eight-c8-O0
/tmp/lrc13-scale-eight-c8-O0 | \
  cmp - 05-knowledge/results/lrc13_scale_eight_hamming_six_frontier_scout_codex_c8.out
```

An AddressSanitizer/UndefinedBehaviorSanitizer build is also byte-identical.
The frozen files have

```text
source SHA-256  af1d3d9d4b7b7537ca47ee3b54ec1c6f1494a61f733813d8ea9df78f49b6f5e2
output SHA-256  7c0f14abc913a8eb9776232790f82fd5cac959e591a0a621d094ee02bb6438a4. (35)
```

The bank-reduction certificate is exact native C++: integer CRT arithmetic, finite mask
dynamic programming, exhaustive 48-bit literal unions, and exhaustive finite
nerve checks.  It uses no floating point, randomness, heuristic pruning, or
`assert` whose removal under optimization changes the checks (`require`
terminates explicitly).  Normal and unoptimized builds are byte-identical,
and sanitizer replay is clean.

`LRCScaleEightOwnerNerve.lean` now kernel-checks the terminal symbolic
quotient.  Its reduced four-digit truth table has `2^6*4^4=16,384` rows and
uses ordinary `decide` to prove distance-two owner obligations disjoint after
forgetting the two private digits.  A separate six-cycle pigeonhole theorem
shows that any three distinct owners contain such a pair, hence the sixfold
intersection is empty.  The module has no `sorry` or `native_decide`; all four
public theorems audit to `propext`, `Classical.choice`, and `Quot.sound`.

The native checker still supplies the preceding completeness bridge: it
validates the symbolic quotient against every literal all-order-eight owner
row and separately proves that all mixed orders die owner-locally.  Thus the
terminal `K_{3,3}` contradiction is kernel-pure, while the 215-million-context
reduction remains a frozen finite-exact certificate.

## Scope guardrail

This theorem closes only the primitive proper AP-centred common-scale-eight
Hamming-six face under the `M(A)<=1/13` hypothesis that invokes THM-860's
common-sheet and hereditary-lcm conclusions.  This theorem alone does not
close common scales `9<=c<=840`, arbitrary non-AP Hamming-six packets, the remaining ramified
Hamming-five metric bank, deep-sheet languages outside this chart, the
seven-wall continuum-to-grid problem, or the block-partition trichotomy.

In particular, THM-963 makes a real finite frontier advance but is not a proof
of global LRC(14).  THM-969, THM-970, and THM-974 subsequently close
`c=9,10,11`; the next AP-centred common-scale Hamming-six frontier is `c=12`.
∎
