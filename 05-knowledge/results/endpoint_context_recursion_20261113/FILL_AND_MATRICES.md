# Spatial multiplicities and matrices with vertex memory

**CITED inheritance / PROVED elementary controls / FINITE-EXACT.**
The spatial-filling idea becomes exact after specifying the space being
filled. It does not permit assuming that every allowed odd count occurs.

## Existing count fibres already have exact spatial geometry

The [Session24 multiplicity theorem](../actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md)
deletes singleton SCCs, leaving an ordered word W of r nontrivial strong
components with total size m and count product h. The r+1 singleton gaps
give the integer orthant N^(r+1). Its size-n slice contains
binom(n−m+r,r) ordinary types. Summing over all such W gives the exact
fixed-h multiplicity formula, eventually polynomial in n. This is
inherited, not rediscovered in this packet.

Incoming [fill-tree multiplicity](../fill_tree_matrix_20261113/multiplicity/PROOF.md)
adds a literal uniform-fill process: choose among all neutral named
attachments, giving transition weight(g_i+1)/(d+r+1). Every size-d
simplex cell is then equally likely. This explains precisely how existing
fibres fill as size grows. It leaves empty target fibres empty.

Let b_n be the least odd value other than7,21 absent at order n. Padding
makes b_n nondecreasing. The assertion b_n→infinity is exactly Chan–Pak,
not an independent premise usable in its proof. The current complete
small head b_3,...,b_8=5,9,17,35,63,611 is inherited from fill-tree.
Growing multiplicity of an already attained value says nothing by itself
about the next absent value.

For comparison with A000568, the neutral-context sequence J_n in
OBSTRUCTIONS is a rooted refinement; a_n counts all ordinary unrooted
types. Different observables need different quotient data. We add no new
term of A000568, A135244 or A002496 and make no unrestricted speed claim
for their enumeration. The practical improvement here is that once a
context is compiled, arbitrary long marked construction words can be
counted with four integers instead of a fresh2^n vertex-subset calculation.
That improves evaluation of a specified family, not enumeration of all
tournament types.

## Correcting the adjacency-power formulation

An n-vertex Hamilton path has n−1 edges; a Hamilton cycle has n edges.
Ordinary(A^k)_(i,j) counts walks of k edges, with repeated vertices allowed.
A transitive tournament has one Hamilton path and A^n=0. More generally,
for a finite digraph A^n≠0 iff the graph has a directed cycle: any n-edge
walk repeats a vertex, while a cycle supplies arbitrarily long walks.
Thus this power tests cyclic behavior, not Hamilton-path existence.

With X=diag(x_1,...,x_n) in a commutative algebra satisfying x_i²=0,

    H(D) = [x_1...x_n] 1^T X(AX)^(n−1) 1,
    HC(D) = (1/n)[x_1...x_n] tr((XA)^n).

Repeated visits now vanish. The cycle convention identifies rotations
but retains orientation. Equivalently use states(S,v), meaning the visited
vertex set and current endpoint; every transition adds a new vertex.
This lifted state graph is acyclic even when D is cyclic. These matrix
and squarefree formulations are classical; see
[Irving–Omar, section2](https://arxiv.org/html/2412.10572v3) and the incoming
[matrix/forest proof](../fill_tree_matrix_20261113/matrix/PROOF.md).

Independent controls here compare permutation enumeration, subset DP,
and squarefree coefficient propagation on all4,096 loopless labelled
digraphs on four vertices. A four-vertex tournament with rows(2,4,9,3)
has H5 but16 walks of length4; a transitive three-vertex tournament has
H1 but zero walks of length3. The oddness theorem is about H, not those
ordinary walk totals or the number of nonzero matrix entries.

## Why tournament completeness matters for7 and21

The forbidden values can occur in oriented graphs that are missing one
edge. Exact witnesses found and independently recounted are:

| H | Oriented graph rows | Complete by adding | Tournament H |
|---:|---|---|---:|
|7|(2,4,9,19,3)|2→4|9|
|21|(4,12,24,33,43,7)|0→1|37|

Thus oddness and the7/21 exclusions do not come merely from the lack of
opposite arcs. Having exactly one arc on EVERY vertex pair is essential.
This is a concrete map from tournaments to the larger digraph family:
delete one specified arc. It preserves vertices and acyclicity of each
individual path, destroys pairwise completeness, and changes the path
count. Its needed sidecar is which paths use that arc. Restoring the
arc restores the tournament restriction; it does not supply a way to
keep count7 or21 inside tournaments.

## Automata, grammar, and information that must survive

The visited-set machine is an exact finite-state encoding of one
Hamilton-path problem. The port compiler instead summarizes a fixed
composition boundary by four integer registers. Its rules form a weighted
grammar with literal graph constructors. Encoding a construction as a
Turing program or cellular automaton is possible, but does not prove
termination of its inverse search or coverage of all odd labels.

The useful additional invariant is the weighted branching sum in
LANGUAGE. It proves that the nine small renewable rules create too few
addresses relative to target growth: even perfect physical composition
has density-zero numerical support. This connects the computation model
to an actual obstruction, rather than treating universal programmability
as a theorem about tournament spectra.
