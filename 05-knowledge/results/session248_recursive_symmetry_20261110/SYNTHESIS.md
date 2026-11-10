# Session248: how symmetry actually emerges as tournaments grow

**COMPLETE, 10 November 2026.** New proved-in-note mechanisms, inherited
inputs explicitly credited, and small new-formula probes. No independent
canon audit or mathematical-priority claim. Chan–Pak and LRC(14) remain OPEN.

The session changes the question from “does a larger tournament have more
symmetry?” to “which symmetry, on which carrier, survives which operation?”
Five answers now fit together: group replacement at births, supported tree
symmetry, arithmetic memory after symmetry loss, recursive invariant modes,
and exact thresholds for coarse observational symmetry.

## 1. Symmetry can grow by replacing the old group completely

Let U be T with newborn x, let r be the orbit size of x's incidence subset
under Aut(T), and let o be the orbit size of x under Aut(U). Then

    |Aut(U)| / |Aut(T)| = o/r.

The denominator is old symmetry broken by the new incidence; the numerator
is the new freedom to replace the newborn by an old vertex. The two
currencies are individually meaningful. Both are odd, but neither vertex
count nor automorphism order is monotone.

The sharpest example is explicit. Let R_q be the cyclic carousel on q odd
vertices, and let B=R_q minus one vertex. Then B is rigid and

    C3[B,B,B]  ->  C3[R_q,B,B]

is a one-vertex birth with r=3, o=q and group orders3 then q. At q=5 this
is a 12-to-13-vertex birth replacing C3 by C5. Every old nontrivial symmetry
is destroyed, yet the total symmetry increases. The construction works
for every odd q>=5.

Regular children explain a contrasting regime. Their deleted parent's
outdegrees recover the newborn incidence, so every parent automorphism
extends and r=1. Searching only regular symmetric children would therefore
miss simultaneous fracture and gain entirely.

To decide whether an old vertex can replace the newborn, an isomorphism of
deleted cards is not enough: it must also transport the full missing
incidence word. Compatible transporters form stabilizer bitorsors. Their
loop products may retain an actual stabilizer permutation even though all
scalar group-order ratios telescope around loops.

[Birth exchange, compatible transport, and constructions](birth/PROOF.md).

## 2. Recursive modules have a measurable depth of symmetry

The cyclic ternary tower R_d has vertices {0,1,2}^d; orientation is decided
at the first unequal digit. Its hierarchy-preserving group is the iterated
cyclic wreath group, of order

    3^((3^d-1)/2).

This is not merely a group-order formula. For any individually marked leaf
set S, let a(S) count the internal prefix nodes occupied by S. The surviving
marked group has order 3^((3^d-1)/2-a(S)). Thus the location of marks, not
just their number, determines the loss.

For a single reversed arc whose endpoints agree for l digits, precisely
2d-l-1 powers of3 disappear from the inherited group. Yet coherently
reversing the quotient edge at that scale requires 9^(d-l-1) physical arc
flips. A local defect can destroy a large amount of symmetry immediately,
while completing its intended coarse change is expensive. Undoing that
single defect still costs only one flip.

All hierarchical quotient choices form an exact weighted repetition code
inside the labelled arc cube. Its rectangles partition all arcs; majority
on each rectangle is the precise nearest-codeword map. This is a real
coding correspondence with a known metric and decoder, unlike a loose
identification of tournament score classes with linear codes. Its minimum
distance is1 if the bottom level is allowed to vary, an important boundary.

The same tower is genuinely 3-adic: first-difference orientation becomes
orientation by the first nonzero base3 digit of y-x. Actual carry and
carry-suppressed odometer actions both become tournament automorphisms.
At depth2 the old physical metacyclic and finite Heisenberg actions inhabit
the order81 hierarchy group, and their supported local rotations explain
the full join. This imports a niche digit-carrier result through an explicit
vertex map, not a shared vocabulary.

[Occupied-prefix laws, weighted code, and the 3-adic map](modules/PROOF.md).

## 3. One observer can erase inherited symmetry without erasing its charge

Coloring leaves by their adjacency to one added root is a dense binary
observation. The exact colored-type recursion shows something unusually
rigid: there are exactly TWO binary coloring types with trivial stabilizer
at EVERY depth of the ternary tower. They are built recursively from
[A,A,B] and [A,B,B].

Thus one carefully chosen root can erase the entire inherited internal
hierarchy group. A random binary root usually does not: the exact success
probability at n=3^d leaves is

    2*3^((n-1)/2)/2^n.

Meanwhile every actual automorphism subgroup K of a module M acts freely
on ordered path covers. Therefore

    |K| divides gcd_c(c! pc_M(c)).

Combining this with Session247's interface budget yields, for k depth-d
towers and q arbitrary exceptional vertices,

    3^[((3^d-1)/2) max(0,k-2q)] divides H.

The root can constrain every leaf, but a Hamilton path uses at most two
incidences at that root. These are different observation budgets. For
three nine-vertex towers plus one root, all inherited internal rotations
can be destroyed while81 still divides the Hamilton-path count.

This is an exact form of **arithmetic memory of symmetry**. It is not a
claim that the final tournament is fully rigid: new unmarked symmetries
might move the root or the chosen hierarchy. The surviving divisor does
not certify a surviving subgroup of the final tournament.

## 4. The metagraph has a recursive spectrum of newly born invariants

Now change carrier from one tournament's automorphism group to the space
of all isomorphism-invariant functions at size n. Let t_n count tournament
types. Define a concrete weighted metagraph move: keep a uniformly chosen
k-vertex induced tournament and resample every other arc fairly.

The old signed-even-graph Fourier basis diagonalizes this new operator.
A graph mode using v nonisolated vertices has eigenvalue

    lambda_v = binom(k,v)/binom(n,v),

with multiplicity t_v-t_(v-1). Hence the operator has rank exactly t_k;
its kernel has dimension t_n-t_k. Every new size-n invariant direction is
annihilated by one-vertex refresh, while older directions survive at known
rates. This describes invariant information, not growth of Aut(T).

The first nonconstant support is v=3. Consequently one-vertex refresh has
exact gap3/n on type functions, versus2/n on labelled functions. The
slowest type mode is centered cyclic-triangle count. Modes of the same
edge degree may have different vertex support, so edge-flip noise and
vertex-refresh noise give two commuting but distinct resolutions.

Birth and deletion have the same structure from another direction. Their
correct type measure is n!/(2^binom(n,2)|Aut(T)|), not uniform weight on
classes. They obey exact detailed balance and

    (n+1) B_n D_(n+1) - n D_n B_(n-1) = I.

The reason is elementary and physical: after birth, deleting the newborn
returns the original state, while deleting an old vertex is a fair refresh.

[Exact spectrum, inherited Fourier inputs and new kernel](bridges/PROOF.md).

### A genuine oscillator, not a metaphor

On the direct sum of invariant function spaces, multiply birth/deletion
by the appropriate square roots of the level. The resulting lowering and
raising operators satisfy exactly

    [a,a*]=I.

At size v there are t_v-t_(v-1) primitive “vacuum” modes. Each generates an
entire orthonormal oscillator ladder as size grows. If L is total size and
N=a*a counts growth steps above a vacuum, then S=L-N is its primitive
vertex support. Refreshing each vertex at rate1 has generator exactly -S.

This connects groupoid mass, random growth, signed graph modes, and
recursive information in one construction. The finite Heisenberg group
acting on ternary vertices in §2 is a DIFFERENT object from this real
oscillator algebra acting on observables across sizes. Both maps are
explicit, and neither is inferred from a name.

[Full oscillator proof and explicit triangle ladder](OSCILLATOR.md).

## 5. Coarse observations can acquire even symmetries that no tournament has

Return to the physical cube of invisible module switches. Session247 says
q exceptional vertices produce an H-response of degree at most2q in the
module-choice bits. An actual binary code now gives exact ensemble results:
if the minimum weight of a distinguishing dual character is d, the cosets
cannot be separated in mean whenever 2q<d.

For even versus odd parity on k bits, the threshold is exactly ceil(k/2).
Ordinal products of the old one-root two-module response attain it.

A richer bank uses k=2^m modules addressed by z in F_2^m. Record

    sigma(b)=(sum b_z, sum b_z z).

There are 2k syndrome sectors. Every fixed mean-H experiment with fewer
than k/4 roots gives the same answer on every sector. Exactly k/4 roots
suffice to separate any chosen pair, with a suitable quotient presentation
and observer template. A family of m+1 separate experiments distinguishes
all sectors. Their recursion is

    (p0,u0),(p1,u1) -> (p0+p1,u0+u1,p1).

Doubling the leaves doubles the necessary interface budget while adding
one hidden syndrome bit. These are explicit physical families and sharp
thresholds, not an appeal to coding distance by analogy.

Below the threshold, the response sectors have an elementary abelian
2-group translation symmetry. Each individual tournament still has an
odd-order vertex automorphism group. The even action changes module arcs
and permutes ensembles; it does not relabel one tournament while fixing
its edges. This is how symmetry can emerge under observation without
contradicting the odd-group theorem.

The result concerns means and the specified bounded joint moments, not
full laws or individual counts. Sharp pairwise constructions may choose
different quotient orders. A universal fixed-quotient realization at the
same budget remains open. The inherited six-vertex modules also retain
their charge3: visibility does not imply arithmetic freedom.

[Exact code transfer, sharp parity and affine-address banks](observers/PROOF.md).

## 6. What is new here, what is inherited, and what to do next

Inherited inputs include the equal-profile module switches and interface
degree bound, the orbit-splitting H barrier, the signed-even graph invariant
basis, and the actual digit actions. They are neither re-audited nor claimed
as new. Standard orbit-stabilizer, character cancellation and oscillator
algebra are tools, not priority claims.

The session's concrete new synthesis is the set of exact operation laws:
- symmetry exchange with an explicit completely replaced larger group;
- occupied-prefix defect loss and coherent block-flip distance;
- two rigid binary observer types and surviving arithmetic group charge;
- vertex-support metagraph spectrum, recursive multiplicities and the
  tournament growth oscillator;
- sharp physical observer thresholds for recursively coded ensembles.

High-value next steps:
1. Compare the repository's older tiling-line operator with the support
   projections. Its commutator with vertex refresh would identify exactly
   which recursive invariant scales its lines mix. No such commutation is
   assumed here.
2. Classify full automorphism groups of one-defect ternary towers. The
   exact common inherited subgroup and defect depth are now available.
3. Seek fixed-quotient, charge-free observer banks attaining the same
   thresholds. Session247 supplies candidate switches but not this theorem.
4. Track the joint birth data (r,o) rather than their ratio alone, and retain
   compatible deletion transporters when an actual symmetry is required.

The new proofs have scoped finite companion probes, not an independent
canon audit. No large old census was rerun. All intentional session files
are published; the session has ended.
