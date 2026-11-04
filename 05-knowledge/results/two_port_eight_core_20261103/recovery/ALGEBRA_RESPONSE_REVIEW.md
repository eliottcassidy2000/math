# Independent audit of module inheritance, twins, and the full twelve-point response

**ACCEPTED with the stated finite dependencies.** Reviewed
[../algebra/PROOF.md](../algebra/PROOF.md) C1–C4, both algebra checkers, the
inherited complete order-ten modular sieve and its independent review,
and [../response/check.cpp](../response/check.cpp) with its wrapper and
normal/optimized receipts. Universal two-port reconstruction remains open.
This is not an independent complete order-twelve census.

## The all-order identities

C1(a) is the unconditional modular identity from the inherited B4 proof.
The uniform-run expansion has leading coefficient H(Q), an odd integer.
For every higher ordered k-path-cover polynomial, the internal directed
derivative skew is divisible by 2k!. Hence all k>=2 terms disappear modulo
four. Symmetry of C on a tournament therefore forces symmetry modulo four
on each of its uniform modules. This uses no hypothetical reversal of the
outside tournament and no descent of its entire response function.

C1(b) has a different target. For each exterior endpoint, the uniform-run
expansion has higher terms divisible by k! and one-run coefficient equal
to its quotient endpoint count times the odd H(M). Thus endpoint parity
outside the module agrees with the quotient. The missing quotient coordinate
follows from the exact equality of total start and end counts. Endpoint
parity passes to the quotient; it is not generally inherited by the module.

The C2 endpoint and ordered-cover identities are exact. If twins a->b are
nonconsecutive in a word or belong to different paths in a two-cover,
exchanging their names preserves every used arc: all their neighbors in
those positions are exterior and see the twins uniformly. In C_bw-C_aw,
the unmatched covers end their first path with the adjacent block a,b;
contracting that block gives C_zw(Q). In C_wa-C_wb, the unmatched block
starts the second path and gives C_wz(Q). These maps are bijections,
including singleton paths and covers with twins in different components.

Consequently the contrast

```
L(a,b;w) = C_bw-C_aw+C_wa-C_wb = C_zw(Q)+C_wz(Q)
```

is strictly positive when a->b and strictly negative for the same named
twins oriented b->a. The strict positivity theorem is correctly cited:
Session176's K(u,v) is exactly the present ordered C_uv. Its existence
proof also suffices for the sign without division: insert the freely
oriented adjacent block z,w into a Hamilton path on the other quotient
vertices, then remove its mutual edge. At least one of the two C entries
is positive. No strongness hypothesis is used.

The comparison allows arbitrary differences elsewhere between the two
quotients. The only requirements are that the same named pair consists
of twins in both tournaments and its internal orientation is opposite.
An exterior label exists for total order at least three. At order two,
the endpoint identities distinguish the orientations. Since ordered C
is recovered linearly from the complete two-port table, the sign also
rules out a positive proportionality of such tables. None of this claims
that a pair must be twins merely because its contrast has a particular sign.

C3 correctly identifies a special reverse direction for endpoint parity.
In M->a->b->M, a path starting in M has one M-run or two separated by the
entire segment a,b. The exact counts are S_M+B_2 and E_M+E_2, with no
additional factor. The fixed-slot squarefree identities give
B_2=E_2=S_M+E_M modulo two. Therefore the intrinsic start and end parities
are swapped in the larger tournament. This justifies descent to M in this
specific quotient; C1(b) alone would not justify it.

## The prime twelve-vertex conclusion

Suppose reversing a twelve-vertex module M in an arbitrary ambient
tournament preserves its complete two-port response. The inherited B3–B4
give intrinsic endpoint parity balance and C symmetry modulo four on M.
For a proper nontrivial submodule N:

- Odd orders at least three violate the all-order odd-current identity.
- Orders four and eight are excluded by the complete inherited class bank
  and the residue computation in this packet.
- Order six has just one C-symmetric-modulo-four class, code22693. Its
  pairs 5->0, 3->1, 4->2 are actual twins. A twin pair inside N remains a
  twin pair in M and in the entire ambient tournament. Reversing M flips
  that pair, contradicting C2's actual ambient sign.
- Order ten leaves a three-vertex quotient. C1(b) makes that quotient
  endpoint-parity balanced, so it is the directed triangle. C3 then makes
  N endpoint-parity balanced as well, contradicting the complete order-ten
  joint modular sieve. That sieve tests 9,733,056 classes, retains 12,670
  endpoint-parity candidates, and finds no joint survivor.
- Order two is already an opposite twin pair in the actual ambient pair,
  contradicting C2 directly.

This exhausts all sizes two through eleven and proves that M must be prime.
There is no ambient-order bound and no requirement that its outside arcs be
strong. The use of C2 is essential: the proof does not infer exact intrinsic
response equality through an unknown outside multiplier. The 336 rooted
six-module expansions independently corroborate the exclusion but are not
needed for this ambient prime-core conclusion. They support the separate
statement about congruence-only order-twelve survivors.

The same reasoning excludes every proper nontrivial submodule of size two
through nine in any hidden core at any order; odd larger submodules are
also excluded. Remaining proper nontrivial submodules must have even size
at least ten. These are necessary structural restrictions on a **single
uniform module reversal**, not a theorem about every pair of tournaments
with unrelated changed arcs.

## Exact finite controls and the full response audit

The scout's 311,872 inputs are specified labelled families and deterministic
streams, with duplicates allowed. The only bit-packed triangle code is
decoded at orders at most eight or small quotient orders. Actual twelve-point
objects use adjacency rows, preserving all 66 arcs. The fast parity DP is
the exact predecessor recurrence over F2; the separate integer and modulo-four
DPs check all 146 parity survivors and specified rejected controls. No sampled
stream is described as a complete class universe.

The independently implemented [algebra_review.py](algebra_review.py) uses
literal permutations and cut positions for small graphs, without importing
the scout. It tests exact twin endpoint/cover differences and signs, the
module and quotient directions in C1, the exact triangle identities and their
swapped parities, and the unique six-vertex hostile. It checks the finite
bank and the full order-ten source pins rather than claiming a second full
execution of those inherited computations.

Both saved twelve-vertex controls are prime. The producer checks all 4,082
proper nonsingleton subsets of each by comparing their interior vertices'
outside rows. The review separately checks every such subset by asking
whether each exterior vertex beats all or none of it; each subset fails.
The sparse and full-response controls therefore exercise the remaining
prime-core regime rather than a module family already ruled out by C4.

The response producer literally scans all 12! permutations for each of two
specified converse cores. Zero, one, and two deleted positions give the
Hamilton, two-path, and three-path tensors. The independent squarefree subset
recurrence computes the same S,E,M,C,U,V,D. All ordered endpoint indices are
retained; in particular a singleton middle component permits j=k in D.
The tensor formula separates adjacent ports, one intervening old path, and
two separated ports. Its terms are exactly the corresponding literal word
decompositions, with no extra ordering factorial.

Each port contributes at most two old-port incidence factors, so the response
has separate Boolean degree at most two in each port's mask. Its values at
all masks of weight at most two therefore determine it by ordinary Boolean
Möbius inversion. The complete battery has 79^2=6,241 entries per core;
the producer compares all entries to direct fourteen-vertex Hamilton DP and
also checks nine high-weight evaluations. The source-pinned normal and
optimized outputs agree.

Every integer range is safe. An order-twelve fixed ordered C entry is at
most 11!, and a three-path tensor is at most 12! times binomial(11,2).
The complete response is at most 14!. Each Boolean coefficient uses at most
16 battery values, and even a conservative sum of 6,241 such coefficient
bounds fits signed 64-bit arithmetic. Array indices and row masks fit their
declared ordinary unsigned types. There is no 64-bit encoding of 66 arcs.
The module checker casts every operand of its signed twin contrasts before
subtraction, avoiding unsigned promotion and implementation-defined conversion
of a wrapped negative value. The scout's unsigned subtraction is used only
for residues modulo four, where wrapping is harmless because four divides
the machine modulus.

For the specified parity-balanced core, 28 pairs have nonzero C skew modulo
four. Its converse differs at 6,240 of the 6,241 small-mask queries. The first
actual witness, with old-to-x mask0, old-to-y mask1 and x->y, is
627597 versus 1139039. The independent review recomputes this query by integer
Hamilton DP and checks all stored reversal identities and the complete
battery difference count. This is a hostile to parity-only faithfulness,
not a modular survivor or a collision.

Reproduce the independent review:

```powershell
python 05-knowledge/results/two_port_eight_core_20261103/recovery/algebra_review.py
python -O 05-knowledge/results/two_port_eight_core_20261103/recovery/algebra_review.py
```

Its receipts pin the reviewed sources and dependencies. They do not replace
the producer's much larger literal twelve-factorial audit or the inherited
complete-bank proofs.
