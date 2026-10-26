# Independent audit of the 1–2–4 tournament interface

2026-10-26. **PASS.** No proof or implementation defect found in the frozen
producer packet. Its counterexample to an inherited feedback-arc wording
error is independently confirmed. No producer file was edited.

Reviewed:

| File | SHA256 |
|---|---|
| `../tournaments/TOURNAMENTS_124.md` | `c6fb8d3696fb54822d2452a7fbd7a46b76f996dcc01dff84c6a7f96db5f6b2f9` |
| `../tournaments/check_tournaments_124.py` | `027d19f5eb64be7f4d83c9482c07131d863301e430eb02c2558b4670808d6f7c` |
| Both producer receipts | `889aee8d7513adc9d2c3b7049c2011bed837f12854a98b97720ce9f1bbda6067` |

The proof and computation scopes below are distinct. Inherited all-order
results were read at their current source; this audit does not claim to
rerun their larger historical censuses.

## 1. The literal multicomplex equivalence

The divisor exponent box minus its top is a downset. Evaluation at the
actual primes gives the divisor sum and evaluation at1 gives its cardinality.
The graph carrier retains exactly its degree-one monomials as vertices,
with absent coordinates treated as ghosts. This convention is essential
for a prime input: its proper-divisor set contains only1, so its graph has
zero vertices, rather than one unused vertex.

The three positive branches are exhaustive at total prime multiplicity at
most2: a prime, a prime square, or a product of two distinct primes.
They give respectively the empty graph, a single isolated vertex and K2.
For total multiplicity at least3, a repeated prime puts the nonsquarefree
monomial X_i² in the downset. It cannot be the deleted top unless the
original integer was precisely that prime square. In the remaining,
squarefree branch there are at least three coordinates; every pair is
allowed but the full set is forbidden. Any graph with all pairs allowed
has no edges and would allow the full set. Thus both directions of the
claimed equivalence are proved, with no numerical-evaluation substitution.

The n=1 void case is correctly separate because every graph independence
complex contains the empty face. A different encoding with new vertices
is not excluded by this theorem. The distinctions between a nonsquarefree
downset, a squarefree nonflag complex and an arbitrary graph are retained.

The coefficient controls are correct. Three independent vertices force
the entire Boolean three-face and already contribute27 at weight2. Hence
an independence polynomial with value7 or21 has degree at most2. Solving
the stated linear equation with the pair-count bound gives exactly the
listed possibilities. K3 and P4 attain these graph values, while the
separate tournament realization condition excludes them as whole cycle
conflict graphs. Their occurrence as induced subgraphs is not excluded.

## 2. The OCF source, path words and orbit counts

I read [Irving–Omar, Corollaries19–20](https://arxiv.org/html/2412.10572v3).
Its cycle convention is directed cycles modulo rotation, and its exponent
counts nontrivial odd cycles; fixed points account for vertices outside
those cycles. The complement has loops, but these do not alter paths with
distinct vertices. Its loopless part is the converse tournament, whose
Hamiltonian words are in bijection with the original words by reversal.
This gives exactly the producer's H=I(Omega,2), with literal directed
cycles retained even when two have identical supports. No loop vertices
are inserted into Omega.

I also read [Busch, Corollary1 and its proof](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf).
Its three residue-class formulas give the stated strong minima and
f(n+3)=5f(n), for the declared strong orders. The producer correctly cites
the inherited small-spectrum and strong-component argument for the
all-order H exclusions, rather than deriving them from the new small scan.

The current source files THM-1370 / h-spectrum-omits-7-21-all-n and
THM-4790 / presentation-fibres-7-and-21-absent-through-eleven-vertices were
read. They support exactly the two different scopes stated here:
ordinary H omits7,21 at every order, whereas the independently audited tc
exclusion covers orders through11. The all-order nonrigid tc question and
ordinary H spectrum completeness remain open. A group action on ordered
Hamiltonian words is free because fixing a word fixes every vertex; this
justifies H/|Aut| but does not identify it with H. The Paley7 example
(189,21,9) makes the loss explicit.

The field-size grouping formula retains every exponent residue modulo r.
It gives s(64)=63=3·21 without requiring the H spectrum to be closed under
division. The H63 complete-conflict witness and the arithmetic sum are
different realizations of63, as claimed.

## 3. Marked F8, its Fano line and its normal basis

The polynomial x³+x+1 is irreducible over F2 because neither0 nor1 is a
root. Its nonzero generator a has the displayed order-seven power table.
The marking r↦a^r is a multiplicative-group marking, not an additive map:
already the image of1+1 mod7 is a², whereas a+a=0.

With that marking, residue doubling is Frobenius and multiplication by4
is its square. The vectors (a,a²,a⁴) have coordinates(2,4,6) and sum0;
they span a two-dimensional subspace and its three nonzero vectors form
a Fano line. Replacing a by a³ produces coordinates(3,5,7). Their span
has dimension3, and squaring rotates them. They are consequently a normal
basis. Their sum is1. The same orbit length3 therefore supports two
different linear ranks; the packet does not confuse those predicates.

The ordinary Collatz cycle traversal1→4→2→1 equals multiplication by4
mod7 only on those three inputs. The Paley triangle has the reverse
traversal, and multiplication by4 preserves its arcs because4 is a square.
The inputs2 and9 are a valid hostile to a global residue-only Collatz
quotient. No dynamical conclusion outside the three-point cycle follows.

## 4. The feedback-five witness and the required inherited correction

For the row masks(106,60,89,48,97,68,10), I independently verified the
tournament axioms, H=127, automorphism order1, and the cycle counts

    length3:11,  length5:23,  length7:11.

There are45 literal odd cycles and9 disjoint pairs. No three can be
disjoint on seven vertices. Thus its profile is(1,45,9), giving
1+2·45+4·9=127.

The independent feedback oracle uses subset dynamic programming, appending
one last vertex at a time. Its result is5. A second enumeration over all
7!=5040 orders agrees. The zero-based order

    (1,2,4,0,5,6,3)

attains5. Exactly nine orders attain that minimum; none has fewer. The
full feedback histogram is recorded in `independent.normal.json`.

This refutes the old universal reading of THM-4748's finite-record sentence
that classes with H in{31,63,127,255} have feedback number3. It does not
refute the lower bound fa≥3 or the fan construction attaining3 at a
possibly larger order. A minimum over an H-fibre and a uniform value on
every class in that fibre are different quantifiers. The producer's
reported correction and its strongest survivor are exact.

## 5. Independent implementations and reproduction

The frozen producer passes8,494 explicit gates normally and with `-O`;
both reviewer replay files match its saved bytes exactly. The separate
`check_independent.py` imports no producer implementation and passes2,042
gates in both modes. Its independent receipt SHA256 is
`5317e13fa81d56b75273e0cd403e3d8d44aa614e1d067c00244a966391539ba5`.

Its distinct routes are: literal integer divisibility for graph faces
rather than an exponent box; subset minimization for the feedback number
rather than only permutation enumeration; all rooted cycles divided by
their rotation length rather than least-vertex-rooted cycles; and
polynomial convolution/reduction for F8 rather than the producer's
shift-and-reduce multiplier. Full binary span enumeration checks ranks2,3.

Reproduce from this audit directory:

```sh
python3 -B ../tournaments/check_tournaments_124.py --output replay.normal.json
python3 -B -O ../tournaments/check_tournaments_124.py --output replay.optimized.json
cmp replay.normal.json replay.optimized.json
cmp replay.normal.json ../tournaments/normal.json
python3 -B check_independent.py > independent.normal.json
python3 -B -O check_independent.py > independent.optimized.json
cmp independent.normal.json independent.optimized.json
```

The finite universes, primality weights, literal-cycle identities and
quotient boundaries agree with the written claims. No new all-order tc
result, H spectrum completeness, aliquot untouchable classification or
Collatz convergence result is certified by this audit.
