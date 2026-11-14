# Least-missing obstructions: black-box module substitution and arithmetic symmetry

In the least-missing argument the numerical induction hypothesis is weak:
it supplies some tournament with every smaller admissible path count,
not any interface or special shape for such a tournament. Write
H(T) for the number of ordered directed Hamilton paths.

The existing [two-prime-lane reduction, THM-4094](
../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md)
is the starting point, not a new result: ordinal multiplication and
the known carry atoms 49,63,343 force a least missing admissible
count to be prime or seven times a prime. The new question is
what geometric restrictions those two arithmetic lanes impose,
and how much an inductive argument can do without a marked
interface.

## Fixed module substitutes cannot supply a nonmultiplicative H-only lift

Fix a quotient tournament and some ordinary fixed modules. In the
other r slots put arbitrary nonempty modules X_1,...,X_r. A
Hamilton path of the result collapses to a word of maximal runs in
their module labels. Adjacent distinct labels must follow quotient
arcs, and every module appears. For a fixed word, if slot i
appears j_i times, the number of internal realizations is the
product of ordered nonempty j_i-path-cover counts C_j(X_i),
including the corresponding counts in the fixed modules. No
cancellation is possible. The count is a finite sum for fixed
slot sizes; the set of potential run words need not be finite
before slot sizes are fixed.

Suppose this fixed template is H-only in an input (or jointly in
all its inputs): it is required to have the same output for all
witnesses of some fixed input count h, with the other inputs
held fixed. Replace any such witness T by its ordinal sum with
arbitrarily many singletons. This keeps H=h. A potential feasible
word revisiting that input j≥2 times grows unboundedly under
this padding (in all its slots if the same witness is used
more than once), contradicting constancy. Indeed assign the old
T path to the first run, and distribute the ordered padding
singletons surjectively to j named runs. All runs are directed
and there are j! S(m,j) such distributions at height m.
All other positive factors can be maintained by padding as
needed.

In particular if a template using r copies of one arbitrary
tournament T is required to depend only on H(T) on all sizes,
there can be no run pattern with any of these r copies revisited.
Thus every copy is traversed exactly once and the only allowable
lift is c H(T)^r. With separate inputs it is c ∏ H(X_i).
A coefficient c=7 cannot give an all-input one-copy lift:
putting transitive singletons in the variable slots would realize
H=7. Likewise no universal affine lift with a nonzero intercept
can be obtained this way. This obstruction concerns a fixed
unmarked whole-module replacement recipe, not a marked or
target-selected construction.

## Quantitative padding sparsity for fixed input counts

Fix attainable h_1,...,h_r. At a fixed H, finite nontrivial SCC
cores suffice: the path counts of strongly connected components
multiply, each nontrivial factor is at least three, and
vertex-pancyclicity together with the positive odd-cycle packing
formula gives H(C)≥|C|−1 for a strong component C.
After its singleton SCCs are stripped the word therefore comes
from a finite catalogue. If B is the stripped word (possibly
empty), inserting m singleton SCCs anywhere multiplies its
ordered q-cover count (empty runs permitted) by q^m, since each
run traverses its SCCs in order. Hence every ordered nonempty
j-cover count satisfies the exact formula
$$
 C_j(B_m)=\sum_{q=1}^j(-1)^{j-q}inom jq P_B(q)q^m
         \asymp j^m\quad(j\ge2).
$$
The asymptotic constant is positive.

For a fixed module recipe and fixed stripped choices, classify
each slot i by whether some feasible word revisits it. Fixed
extra modules impose fixed run-capacity bounds. For every slot i,
the finite-capacity query reduces to a finite antichain of
componentwise-minimal variable run-count vectors among words with
j_i≥2. This is Dickson's lemma: every set in N^r has finitely many
coordinatewise minima (an infinite sequence contains an infinite
subsequence nondecreasing in coordinate one, then in coordinate
two, and so on). Thus which repeats are feasible depends only
on sizes truncated at a fixed constant. In each such regime use
one fixed witness word for each repeated slot; its factor grows
exponentially in that slot's padding length, and all other factors
are positive integers. Hence an output at most X has O(log X)
possible padding heights at each repeated slot. Height at any
never-repeated slot is irrelevant, since its H is fixed and
feasibility cannot depend on its size beyond one. The fixed bank
has O((1+log X)^r) distinct output counts up to X. This does not
bound a bank whose input counts or quotient recipes themselves
range with X.

## Arithmetic prime factors and actual geometric symmetry

Automorphisms act freely on ordered Hamilton paths, so
|Aut(T)| divides H(T). Conversely arithmetic divisors need not
produce geometric symmetry. For an automorphism g, each vertex
orbit of length d induces a vertex-transitive tournament,
therefore every vertex there has outdegree (d−1)/2. Its number
of cyclic triangles is exactly
$$
 inom d3-dinom{(d-1)/2}{2}=d(d^2-1)/24.
$$
All the empty packings and single-triangle packings on disjoint
vertex orbits can coexist independently (whether or not the
orbits are homogeneous modules in T). The odd-cycle packing
formula gives the useful budget
$$
 H(T)\ge\prod_{	ext{orbits}}(1+d(d^2-1)/12).
$$
Thus any witness with H=p prime p≥5 is rigid: an order-p
automorphism, which would exist if |Aut|=p, has a p-orbit and
that orbit's triangle bound alone exceeds p. For H=7p with
p≥11 prime one similarly rules out p dividing the group order,
since 1+p(p²−1)/12>7p. Its automorphism group therefore has
order 1 or 7. These are restrictions on possible witnesses,
not a representation theorem for the two prime lanes.

The distinction is not just negative. The eight possible
circulant tournaments obtained on a 7-orbit have H equal to
175 or 189 (literal dynamic-programming check), sharpening
the per-orbit budget there. The circulant with forward steps
1,2,3 modulo 7 has (H,P_2)=(175,1988), where P_2 is the
ordered two-path cover count allowing empty slots. Ordinally
append a transitive k-tape to form B. Now take the cyclic
three-module substitution C3[1,1,B]. Its path count is
H(B)+P_2(B): packings avoiding the two singleton portals
contribute H(B), and all remaining odd cycles use both
portals plus a nonempty odd-order path in B. Their complement
weights sum to half of P_2(B) because the signed two-cover
convolution vanishes on nonempty B (each odd packing block
has two opposite signs). Thus
$$
 H=175+1988\cdot2^k=7(25+284\cdot2^k).
$$
At k=1 and k=3 the cofactors 593 and 2297 are prime (trial
division). Each resulting tournament carries an actual
7-rotation of the circulant, and the budget proves that
the entire automorphism group has order exactly seven.
No assertion about infinitely many prime cofactors is needed.

Run python replay_symmetry.py.
