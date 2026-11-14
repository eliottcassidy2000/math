# Nonnegative cyclic run inverse beyond singleton joins

For nonempty tournaments A,B,C, let C₃[A,B,C] orient every cross edge
as A→B→C→A. This is strong even when a source block is not strong.
Let G_j(T) count ordered partitions of a tournament's vertices into j
nonempty directed paths; it is zero for j>|T|. If F_i allows empty slots,
then G_j = Σ_i (−1)^(j−i) binom(j,i) F_i.

Compress any full path in the cyclic join into maximal runs in the
same block. Its block labels advance cyclically: a run in A can only
be followed by one in B, then by one in C. Conversely any such label
sequence and three actual ordered path partitions give one full path.
For each k≥1 the triple of run counts is either (k,k,k), with three
choices of the initial label, or an unequal triple in {k,k+1}³,
with exactly one compatible initial label. Thus the target formula is

    H(C₃[A,B,C])
      = Σ_(k≥1) [ 3 G_k(A)G_k(B)G_k(C)
                  + Σ_ε G_(k+ε_A)(A)G_(k+ε_B)(B)G_(k+ε_C)(C) ],

where ε ranges over the six nonconstant binary triples. This is a
nonnegative expansion in physical supports, independent of signed
odd-cycle formulas. In particular the ordinary count is symmetric
under all permutations of the three block tournaments, despite the
directed orientation of the join. Converse can also be applied to
any individual block without changing its G_j or the total.

Each ordered partition of the n vertices into nonempty slot sets
gives at least one cover by Rédei, so G_j(T)≥j!S(n,j), with equality
for transitive T. In particular a module of order a contributes at
least 2^a−2 through its two-run term, in addition to the at least
three paths at (1,1,1). For prescribed H=N, every module therefore
has order at most floor(log₂(N−1)). This gives a finite exact inverse
for the entire uniform cyclic triple class and a direct size pruning
rule without assuming that any particular profile is physical.

## Two consecutive connector vertices

If C=TT₂, write h=G₁, p=G₂, q=G₃ for A and B. The run formula is

    H = 5h_A h_B + 3(p_A h_B + h_A p_B) + 7p_A p_B
        + 2(q_A p_B + p_A q_B + q_A q_B).

The rank-three kernel here is already recorded in [THM-5133's audit](../s92_prime_powers_quaternions_bases_20261111/audit/AUDIT_REPORT.md).
The use here is a direct run proof and exact target separation.

At A=TT₂ and B=TT_n, this reduces to H=4·3^n+5·2^n−11.
It gives 45 with n=2 and the prime 137 with n=3. At B=C₃
the same two-connector template gives 159. None of 45, 137,
159 is attainable with a singleton block in any order. Indeed
when C is a singleton the formula simplifies to

    N=(F₂(A)−h_A)(F₂(B)−h_B)+2h_A h_B,

which bounds each order by floor(log₂(N−1))≤7 for these
targets. The complete canonical class search at these orders
in replay_cyclic_runs.js checks the actual h,F₂ possibilities
and the literal witnesses. This is a genuine alternative at
these prescribed targets rather than another factorization
of the singleton equation.

The whole uniform cyclic triple class still misses 63.
A singleton block fails the same bounded exact search.
If all orders are at least two and at least one is at least
three, nonnegativity and the ordered partition lower bound
force H≥H(C₃[TT₂,TT₂,TT₃])=137. If all are two, the unique
output is 45. Hence any universal selector must also include
a route beyond uniform cyclic triples. This is not a scalar
nonexistence result for 63, which is independently attainable.

Run node replay_cyclic_runs.js to check arbitrary cyclic triple
formulas through total order seven, the literal 45/137/159
outputs, and all-order singleton exclusions for 45/63/137/159.

## One split seam reaches a hole of the whole uniform class

Take TT₂,TT₂,TT_n and reverse the cross arc from the first vertex
of the first TT₂ to the first vertex of the second TT₂. The TT_n
block is still uniform; the other four vertices are individual
marked sites. To count, compress the tail's appearances in a full
Hamiltonian path to j identical tail symbols. Repeated adjacent
symbols are not allowed by maximality; each occurrence must obey
the original entry and exit direction. The exact skeleton counts
for j=1,2,3 are 13,12,2. They can be checked by enumerating the
four marked vertices and tail symbols (companion replay). No path
can use four tail runs: each intermediate fragment between two
such runs starts in the first TT₂ and ends in the second, so
three intermediate fragments would require six marked vertices.

Multiply these coefficients by the respective counts of ordered
nonempty run assignments in TT_n:

    H = 13 + 12(2^n−2) + 2(3^n−3·2^n+3)
      = 2·3^n+6·2^n−5.

At n=3 this is the prime **97**, which C₃[A,B,C] with uniform
cross arcs cannot hit: a singleton is excluded by the finite
test, and all nonsingleton blocks give either 45 or at least
137. The split seam is still strong: the first TT₂ vertex
still beats the second vertex of the facing block, which
reaches TT_n and returns to each vertex of the first block,
and all other blocks remain reachable.

A smaller-count fallback within one apex is also explicit.
Start from C₃[TT₁,TT₁,TT₆] with H=65; flip the apex contacts
to the third and fourth vertices of TT₆. The normalized binary
contact word is 1001100, giving H=63 and rooted response R₊=77.
More generally the word 10011 followed by r≥1 zeroes gives
H=19·2^r−13 by the known one-apex recurrence (g(19)=5).
The 97 and 63 escapes replace the failed geometry, not just
the numerical factor pair. They are not a universal selector.
