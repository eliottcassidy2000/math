# Wide coupled moves can have a three-state exact compression

23 October 2026. **PROVED scoped.** The incoming coupled-move theorem is
an exact criterion for a supplied bounded integer model. The new native
family below shows that no fixed local support bound suffices for its
one-hot frame encoding, even on a path; it also gives an exact compression.

## 1. A native tournament family

Take m named directed triangles, with owners(i,a), a in F3. Between consecutive
cells i,i+1 use the matching block B(a,b)=1[a=b], oriented from i to i+1.
Every other intercell block is fixed uniform. Internal automorphisms are
the three rotations a->a+t_i. The matching equations force t_i=t_(i+1),
so the full cell-preserving group consists of precisely three simultaneous
rotations, for every m>=1.

This has actual Hamiltonian ports. In cell i, visit the owners

    (t-i), (t-i+1), (t-i+2) modulo3,

and concatenate the cells in increasing i. The last owner of a cell equals
the first owner of the next, so every cross-cell step uses a matching arc.
The three simultaneous rotations carry these three actual paths into one
another. An indexed owner fixes their phase. Forgetting the phase preserves
the tournament class but not the literal ordered path.

## 2. The smallest move is necessarily global in the one-hot encoding

Use z_(i,a) in {0,1}, sum_a z_(i,a)=1, and retain the shared equations
z_(i,a)=z_(i+1,a). The integer fibre has just three points x^(a), each choosing
the same state at every cell. Every nonzero kernel vector of the equality
matrix has the same three coordinates at each cell, with their sum zero.
Its conformally primitive vectors are exactly

    1_m tensor (e_b-e_a), a != b.

Indeed any zero-sum integer triple decomposes conformally into unit transfers
from negative to positive coordinates. Hence each atom has support and l1
norm2m. At a feasible x^(a), the two moves to the other states are feasible.
No nonzero admissible change supported on fewer than m cells exists.

Thus a bank of individually feasible bounded-width cell switches misses
real changes for arbitrarily large m. This is an encoding-relative lower
bound, not a complexity lower bound for the problem: the fibre has size3.
Retaining one root phase t immediately solves it.

For a separable observation P(z)=sum_i sum_a w_(i,a)z_(i,a), form the exact
three-vector

    W(a)=sum_i w_(i,a).

P is constant iff W(0)=W(1)=W(2). A port-sensitive example is w_(0,a)=a
with all other weights zero, giving three different values. A cancellation
example uses row weights v_i(a)-v_(i+1)(a) with v_0=v_m, giving W=0 although
individual cell contributions vary. The shared root makes this cancellation
law valid. Forgetting shared phases gives false combinations.

More generally, on functional chart components use one authenticated root
state and the composed actual transports to evaluate each local factor.
Sum all contributions in that SAME state. If multiple root states remain,
compare their full projected values; if distinct components interact through
an observer or guard, retain that joint factor. This is a special efficient
consumer of the incoming owner-indexed compiler, not a replacement for it.

## 3. Fibonacci compression retains addresses, not only the final product

Suppose the local weights are prescribed by a fixed substitution word, such
as the incoming canonical C->CA,A->B,B->CA. For a block U retain its length,
its actual accumulated permutation rho_U of the root-state set, and its
observation vector W_U(t). For concatenation U V use

    rho_(UV)(t)=rho_V(rho_U(t)),
    W_(UV)(t)=W_U(t)+W_V(rho_U(t)).

These equations follow by evaluating the two consecutive blocks. They are
associative, give exact weighted finite-state composition, and work for
noncommuting permutations and vector-valued costs. Repeated substitution
blocks can be evaluated by a straight-line recursion without expanding
their whole words. Noncanonical occurrences with different choices cannot
share a cached block unless those choices are included in its key.

A queried interior port requires the prefix transformation at that actual
address. Retain the substitution decomposition and split recursively by
stored lengths. At each split carry the transported root state. The final
permutation and total W alone do not answer arbitrary such queries: for
example a nontrivial permutation followed by its inverse has the same final
permutation as two identity blocks but a different middle owner.

This connects the Fibonacci address structure to exact attachment and cost
composition. It does not identify the three letters with the three F3
owners: the incoming binary quotient explicitly identifies two letters.
The choice of a permutation and cost for each letter is declared data.
The checker's noncommuting S3 example includes a transposition. It tests a
generic finite-state adapter, not an internal automorphism of the native
C3-chain cells. Those native path witnesses are checked separately.

These prefix queries use ordered slot addresses. For positive-length tiles
stored lengths locate a slot; zero-length slots need their index or an
explicit before/after boundary convention as well. A zero-length phase
change can separate two states at the same geometric position. Exact
composition retains that event even though its length contribution is zero.

## 4. Retained constraints

Actual path positions, marks, converse and every required electrical guard
must be transported or checked in the resulting root state. A guard which
distinguishes root states is a retained predicate, not an average. The toy
linear observations here are not asserted to be the native mod8 invariant.
The old mod8 input remains inherited unchanged; no numerical coincidence
substitutes for it.

The [square-sum lane](square_sum/REPORT.md) supplies another actual graph
context: a pinned three-vertex simultaneous recolouring can permit a new
attachment although every allowed one-vertex recolouring fails. Its graph,
pins and colouring predicate differ from the tournament frame model, but
the preserved-boundary coupled-move mechanism is the same.
