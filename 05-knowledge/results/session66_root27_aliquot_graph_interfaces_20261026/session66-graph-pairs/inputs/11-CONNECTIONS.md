# Session65: continuation interfaces, not one universal isomorphism

This is a proposed organizing framework supported by the exact interfaces
below, not a claim that the repository's open problems form one completeness
class. It continues Session63's integral-observer proposal and Session64's
continuation-exact summaries. The decisive extra question is now: **which
computation, clock, guard and witness survives the representation?**

## 1. Three carriers must remain distinct

**Finite tape carrier.** A binary tape window w is encoded by an ordered join
of transitive three-vertex tournaments for0 and cyclic triangles for1.
The word is recoverable from ordered strong components, even when successive
zero blocks merge into larger transitive clone modules. The inherited
order-join metric theorem makes the isomorphism edit distance exactly the
Hamming distance between equal-length words. A native write flips at most
one designated literal arc. Tape origin, head position, machine state and
transition table remain additional data.

This is an exact compositional code, not intrinsically a fractal, and not a
simulation by the repository's Hamiltonian K moves. The operation bank is
different. Its Hamiltonian-path count is3^(number of1bits), so that scalar
forgets symbol order. Words01 and10 have the same count but different tapes.
The Irregular Bell gives a stronger hostile: a guarded sweep changes one
bit, hence tournament distance1, while taking2k+1 native instructions.
No running-time complexity bound follows from endpoint edit distance.

**Arithmetic inverse-limit carrier.** The Marxen–Buntrock macro map has live
ternary residue branches0 and1 and halts at residue2. Its contracting inverse
branches give an actual binary Cantor full shift in the3-adics. Here
self-similarity is mathematical, not a visual analogy. A level-m safe word
is exactly a residue cylinder modulo3^m. The integer embedding, branch
address and macro/native clocks remain indispensable.

Under the signed ternary tournament from THM-4754, the induced survivor-set
tournament is transitive: alternate the complement of every other code bit
to obtain a lexicographic orientation. Thus a nontrivial fractal dynamical
system can have a very simple tournament shadow. That shadow does not decide
whether a nonnegative integer lies in the infinite survivor set.

**Uniform problem carrier.** A counterexample-search program is a computable
presentation of a quantified predicate. Its finite control graph is not the
space of possible configurations: tape contents and head locations make the
latter infinite. A25-state machine is not a25-vertex copy of Goldbach.

## 2. A hierarchy of claims, not interchangeable terminology

The evidence must say which of the following it establishes:

1. An injective encoding with an explicit decoder and required labels.
2. A labelled control-table isomorphism, including read/write/move/halting
   entries, not just a graph on target-state names.
3. A conjugacy or guarded macro simulation of full configurations, with a
   clock map and a converse describing which target states are realizable.
4. A halting or output equivalence for a specified source/input.
5. A uniform computable reduction between parameterized problem families.
6. A resource-bounded reduction, with representation size and running time
   analyzed, sufficient to transfer a particular complexity classification.

These are not all totally ordered in strength. Each names a different
preserved predicate. In particular, logical implications between mathematical
conjectures must be labelled separately from machine simulations or complexity
reductions. A single fixed conjecture is not by itself a uniform language
with a proved hardness classification.

The Goldbach case illustrates two genuinely different successful claims:
the wiki and pinned Lean25-state tables are isomorphic by state renaming;
the theorem-to-machine bridge instead says that the blank computation halts
exactly when strong Goldbach has a counterexample. The2016 construction
lineage had27 states; the25-state formalized compression is later. Reachable
two-state bounce padding preserves the full macroconfiguration and halting
while changing program size and native clock. Thus25 is not a semantic
invariant shared by every implementation of that mathematical problem.

## 3. Why the naive universal graph classification fails

A deterministic configuration graph has at most one outgoing edge per
configuration. On four or more vertices it cannot be a tournament. A single
visited run, after forgetting its state/tape labels, is merely a finite
halting path, a tail entering a cycle, or a nonrepeating infinite ray. All
nonrepeating nonhalting runs have the same unlabelled rooted ray shape.
Recognizing which case a finite program produces is still difficult; the
simple shape does not decide the case.

There is a second, independent symmetry obstruction. Finite tournament
automorphism groups have odd order: an involution would exchange a pair and
reverse its directed arc. Arbitrary control graphs can have involutions.
Hence there is no fully faithful equivariant encoding of all such graphs
into tournaments that injects their automorphism groups. This does NOT rule
out a general isomorphism-testing reduction, which need not preserve groups,
nor an explicitly ordered or symmetry-broken encoding with sidecars.

The useful repair is a typed interface category: objects contain the full
state carrier, guards, allowed operations, outputs and clock; morphisms must
state what they preserve, what they forget and how realizability is recovered.
The reduction atlas supplies concrete examples and boundaries rather than
declaring all objects equivalent.

## 4. Fifteen: an actual Frobenius obstruction, not a count coincidence

Let S_p be the set of geometric supersingular j-values, with Frobenius
permutation sigma. Its nontrivial orbits are b transpositions. Every
tournament orientation on S_p must fail Frobenius invariance on each swapped
pair, and all other unordered-pair orbits can be oriented consistently.
Thus the minimum orientation defect is exactly b.

The primary modular-curve genus formula identifies b with g(X_0^+(p)).
Consequently the fifteen Ogg primes are exactly those with defect zero.
At p=37, S_p={8,3+10t,3-10t}, t^2=2, and minimum defect is1. This is an
optimization over all orientations of the actual Frobenius set. It does not
construct a preferred arithmetic tournament or explain the cardinality15.

The general permutation theorem clarifies the shared mechanism: each orbit
of unordered pairs has a return sign. Negative-return orbits force odd
defect parity; positive ones force even parity. The product of their exact
local enumerators gives the full defect polynomial. Its minimum is the
number of even vertex cycles. This recovers the involution case and the
odd-order automorphism obstruction in one precise signed-action language.

Sloane's threshold15 is a different fact about zero-free ternary powers.
The historical15-state Erdős machine encodes a different digit conjecture.
Neither number is identified with the Ogg set by a proved isomorphism.

## 5. Finite digit information does not identify the intended integer orbit

For every k>=1, the exponent n_k=3^(k-1)-1 satisfies

    2^n_k = (3^k-1)/2 modulo3^k.

Its k trailing ternary digits are all1, so every finite tail depth admits
safe candidates for both digit filters. The compatible profinite exponent
has even parity and3-adic coordinate-1. It is not an ordinary nonnegative
integer exponent, and it is not the usual exponent-1 with its odd parity.
The parity sidecar cannot be discarded. Tail safety does not assert that
all digits of the full positive integer2^n are safe.

Likewise the Marxen–Buntrock survivor cylinders have measure(2/3)^m at
depthm, and the infinite survivor set has Haar measure0. Exact residue
counting proves density-one halting among nonnegative integer seeds and a
power-law upper bound on possible exceptional seeds. It does not prove
every-seed halting. A sparse prescribed orbit is not controlled by ambient
density. No equality between the two Cantor carriers is asserted.

## 6. Multiply-perfect numbers: an inherited map with a sharp scope

At order2^r with r>=2, the tournament arc cube has dimension
D_r=2^(r-1)*(2^r-1). Multiplicativity of sigma gives

    sigma(D_r)/D_r = sigma(2^r-1)/2^(r-1).

This is an exact multiperfect test on the specified family. For r=4 the
dimension120 is triply perfect; r=10 gives another triply perfect example,
523776. Prime Mersenne apex gives the familiar even-perfect case. The map
forgets every orientation and cannot represent hypothetical odd witnesses,
because all D_r here are even.

For odd n, the separate two-adic divisor-sum ledger gives necessary prime
exponent conditions. These do not supply missing odd divisibility or an
existence proof. A search for odd n dividing sigma(n) is a finite-witness
halting search, not an equivalence with this restricted tournament family.
`multiperfect/REPORT.md` supplies the proof, inherited status and hostiles.

## 7. What the recovered25+ session contributes

The old S62c25+ cross-field connection list is a genuine historical source,
not a certificate that25 transformations were proved. Its dihedral, Walsh,
Paley, tropical, cluster and matroid lenses are candidate methods. Later
atlases explicitly repaired missing coordinates and scope. Reuse the maps
and their counterexamples, not the cardinality of the list.

The organizing proposal for future work is therefore narrower and stronger
than “all open problems are isomorphic”: assemble a graph of **typed,
verified transfers**, retain exact witnesses and complexity overheads, and
leave absent converse reductions visibly open. The present session supplies
several such arrows and several rigorous obstructions to falsely adding them.
