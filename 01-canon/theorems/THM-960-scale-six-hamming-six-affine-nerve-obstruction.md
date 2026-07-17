---
id: THM-960
title: Scale-six Hamming-six affine-nerve obstruction
status: PROVED STRUCTURAL + FINITE-EXACT — all 37,710,288 primitive proper AP-centred c=6 order/unit contexts fail necessary common-sheet coverage; an independent reduction leaves 64 all-order-six signed-doubling supports whose six affine owner four-flats have octahedral pair nerve and no triple intersection, so no metric recursion is born
source: codex-2026-07-17-S62 scale-six exact and structural audit
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860]
related: [THM-862, THM-957, THM-958, THM-962, THM-963, HYP-6820]
verification:
  - 04-computation/lrc13_scale_six_hamming_six_common_sheet_codex_S62.cpp
  - 04-computation/lrc13_scale_six_hamming_six_structure_codex_S62.py
  - 05-knowledge/results/lrc13_scale_six_hamming_six_common_sheet_codex_S62.out
  - 05-knowledge/results/lrc13_scale_six_hamming_six_structure_codex_S62.out
---

# THM-960 — scale six is killed by an affine-flat nerve

Let `R subset F_13^*` have six elements, put `P=[12]\R`, and consider

```text
A=6P union {w_r:r in R},
w_r=6r (mod 13),              w_r>0,              w_r!=6r.
```

If `M(A)<=1/13`, THM-860 supplies common-sheet coverage at every replacement
owner and the hereditary leave-one-out law for the effective orders

```text
D_r=6/gcd(6,w_r),             e_r=(D_r w_r/6) mod D_r,
lcm(D_s:s!=r)=6               for every r.
```

## Exact verdict

The primitive proper AP-centred common-scale-six Hamming-six sheet bank is
empty.  Therefore every packet above satisfies `M(A)>1/13`; no metric height
recursion is required.

The possible states are

```text
(D,e)=(1,0),(2,1),(3,1),(3,2),(6,1),(6,5).
```

The leave-one-out law is equivalent to at least two even orders and at least
two orders divisible by three.  Inclusion-exclusion gives

```text
6^6 - (3^6+6*3*3^5) - (2^6+6*4*2^5) + 91 = 40,812
```

admissible state words per label set.  Thus the literal finite bank has

```text
binom(12,6)*40,812 = 37,710,288
```

contexts.  Packing the six sheets at each of six owners into 36 bits and
testing the union of the six provider masks finds zero covers.

There is an independent structural reduction.  Forgetting units leaves 3,249
admissible order words and 3,002,076 labelled order contexts.  Of these, 5,838
pass the scalar sum-of-mask-cardinalities test at every owner.  Requiring only
that each owner be coverable by *some owner-dependent unit word* leaves exactly
64 contexts.  All 64 have

```text
D_r=6 for every r.
```

Thus all mixed order words fail before any global unit compatibility or metric
height is tested.

## Local mask grammar

Normalize the owner to one.  Write hexadecimal masks for sheets `0,...,5`.
The complete local table is

```text
(1,0):  1 -> 3f;                         all other ratios -> 00
(2,1):  1 -> 2a;  6,7 -> 15;            all other ratios -> 00
(3,1):  1 -> 24;  4,5 -> 12; 8,9 -> 09; all other ratios -> 00
(3,2):  1 -> 24;  4,5 -> 09; 8,9 -> 12; all other ratios -> 00
(6,1):  1 -> 20;  2,3 -> 10; 4,5 -> 08; 6,7 -> 04;
         8,9 -> 02; 10,11 -> 01; 12 -> 00
(6,5):  1 -> 20;  2,3 -> 01; 4,5 -> 02; 6,7 -> 04;
         8,9 -> 08; 10,11 -> 10; 12 -> 00.
```

The last two rows expose the all-order-six structure.  At every owner:

- the self provider fills sheet 5;
- exactly one provider with ratio in `{6,7}` fills sheet 2 independently of
  its unit;
- exactly two providers from `{2,3,10,11}` must fill sheets `{0,4}`;
- exactly two providers from `{4,5,8,9}` must fill sheets `{1,3}`;
- ratio `12` is forbidden.

Exactly 64 six-label sets have this pattern at every owner.  They are exactly
the recurring signed-doubling supports

```text
p -> o iff o/p in {2,-2},
```

on which this relation is a directed Hamiltonian six-cycle.  Their
multiplicative orbit sizes are

```text
4,12,12,12,12,12.
```

This identification is structural.  The unique `{6,7}` provider at an owner
is exactly its unique incoming signed-doubling predecessor, because
`{6,7}^{-1}={2,11}`.  A positive-indegree finite digraph contains a directed
cycle.  A cycle of length `m<=5` would force `(+-2)^m=1 (mod 13)`, which never
occurs, so the cycle uses all six vertices.  Conversely, from a signed
doubling cycle the five cyclic distances from an owner have ratio classes

```text
{+-2}, {+-4}, {+-8}, {+-(2^-2)}, {+-(2^-1)},
```

which distribute as `2+2+1` among the two reflected sheet pairs and the
central `{6,7}` bin.  Thus every signed cycle is owner-locally feasible.  The
usual rooted-cycle count `12*2^5/6=64` supplies the exact cardinality without a
`binom(12,6)` support scan.

## The faithful object is an affine-flat nerve

Encode `e=1,5` by a bit `x=0,1`.  At a fixed owner, the two providers serving
each reflected sheet pair impose one affine parity equation:

```text
x_p xor x_q = 1   when their ratios lie in the same ratio bin,
x_p xor x_q = 0   when their ratios lie in opposite reflected bins.
```

There are two such equations per owner.  They use disjoint provider pairs;
the owner's own bit and the unique `{6,7}` provider bit are free.  Hence each
owner obligation is an affine four-flat in `F_2^6`, with exactly 16 unit words.
The literal mask predicate and these two parity equations agree on all

```text
64 supports * 2^6 unit words * 6 owners.
```

The nerve can then be proved without enumerating those unit words.  Orient the
doubling cycle and write

```text
v_(i+1)=(-1)^(a_i) 2v_i,       xor_i a_i=1,            (7)
```

where the last identity is the cycle-closure law `2^6=-1 (mod 13)`.  If `x_i`
is the unit bit at `v_i`, the two equations for owner `v_i` reduce uniformly
to

```text
A_i: x_(i+1) xor x_(i-2)=a_(i-2) xor a_(i-1) xor a_i,
B_i: x_(i+2) xor x_(i+3)=a_(i+2).                     (8)
```

The equations `A_i` and `A_(i+3)` have the same left side, while their right
sides differ by `xor_j a_j=1`.  Thus opposite owner flats are disjoint.  For
any nonopposite pair the four rows in (8) are independent, so the intersection
is an affine two-flat with four unit words.  For each transversal owner triple,
the six rows have rank five and their unique row dependence has right side
`xor_j a_j=1`; hence no three owner flats meet.  This proves the nerve claims
below analytically.  Finally, the six flats contribute `6*16=96` owner-word
incidences and their twelve compatible pairs contribute `12*4=48` pair-word
incidences.  Since no triple meets, exactly 48 unit words satisfy two owners,
none satisfies one, and the remaining 16 satisfy zero.

Their intersection nerve is uniform over all 64 supports:

```text
unit words satisfying 0 owners: 16
unit words satisfying 2 owners: 48
all other satisfaction counts: 0.
```

Two owner flats are disjoint exactly when the owners are antipodal on the
doubling cycle, equivalently their label ratio lies in `{5,8}`.  These three
antipodal pairs form a perfect matching.  Every other pair intersects in an
affine two-flat of four unit words.  Thus the pairwise nerve is the octahedral
graph

```text
K_6 minus the three cycle-antipode edges.
```

No three owner flats meet.  The exact minimal inconsistent subfamilies are:

```text
the 3 antipodal pairs,
the 8 triangles of the octahedral graph.
```

Writing the doubling cycle as `(v_0,...,v_5)`, the eight triangles split as

```text
{v_i,v_(i+1),v_(i+2)} for i in Z/6Z,
{v_0,v_2,v_4}, {v_1,v_3,v_5}.
```

These are exactly the six consecutive and two alternating owner triples found
at scale five.  Scale six strengthens that same triple obstruction by making
the three cycle-antipode owner pairs themselves incompatible.

## Tournament and carrier audit

Completing the sparse signed-doubling cycle by the lexicographically first
Hamiltonian path of its triangular-prism tie graph reproduces the scale-four
tournament fingerprints:

```text
score multiset             (1,2,2,3,3,4)
directed triangles         6
SCC sizes                  (6)
sparse-cycle flips         {2:8,3:52,4:4}
Hamiltonian-path counts    {29:32,31:20,37:12}
joint fingerprints         5.
```

This completion recognizes the recurring 64-support scaffold but cannot see
the obstruction: it forgets the two affine equations attached to each owner
and therefore the octahedral intersection nerve.  Useful carriers form a
strict information ladder:

```text
runner labels / signed C6
  -> edge-coloured ratio bins
  -> six affine owner four-flats in F_2^6
  -> their octahedral nerve.
```

The last two are faithful for common-sheet compatibility.  A bare tournament
is not.

## Cross-scale interpretation

The same 64 signed-cycle supports now organize three consecutive scales:

```text
c=4: a rank-four two-triangle code leaves 4 full-cover unit words;
c=5: every pair of owner obligations is compatible, but 8 cycle triples fail;
c=6: only octahedral pairs survive; the same 8 triples fail and 3 antipodal
     pairs fail already.
```

This is a genuine ramification ladder on one underlying signed-cycle object,
not three unrelated enumerations.  The relevant invariant is the nerve of the
owner-obligation family as the sheet alphabet changes.

## Verification

```bash
c++ -O3 -std=c++20 -Wall -Wextra -pedantic \
  04-computation/lrc13_scale_six_hamming_six_common_sheet_codex_S62.cpp \
  -o /tmp/lrc13-scale-six-S62
/tmp/lrc13-scale-six-S62 | \
  cmp - 05-knowledge/results/lrc13_scale_six_hamming_six_common_sheet_codex_S62.out

python3 04-computation/lrc13_scale_six_hamming_six_structure_codex_S62.py | \
  cmp - 05-knowledge/results/lrc13_scale_six_hamming_six_structure_codex_S62.out
python3 -O 04-computation/lrc13_scale_six_hamming_six_structure_codex_S62.py | \
  cmp - 05-knowledge/results/lrc13_scale_six_hamming_six_structure_codex_S62.out
```

Frozen SHA-256 values:

```text
e6c197e66923bd10f34fd3359eaf6a521f0a09ced6ecd375e3a941de94bcd727  C++ source
188f3bf57c1b8dd3510e9728f962a718cd476a6d62edb15e012646b1c17b3a0a  C++ output
6c997963d4052da988e5bdd8f14cef6e0ee92bd93c2d1444fd83f9ade1eb9358  Python source
3f750387ec3fe71225406bc59539dccaefa8e7f5283972a5be6259bb6bb037d1  Python output
```
