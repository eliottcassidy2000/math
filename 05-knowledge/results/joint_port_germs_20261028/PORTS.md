# Port registers, strong tags, and exact behavioral minima

**PROVED scoped elementary deductions; same-agent audit.** Local labels J1-J6
belong to this file. A tournament here is an actual binary relation, not a
tournament imposed on unrelated data.

## J1. A recognizable ordered port register

Let T be a nonempty tournament on n vertices and let p_1,...,p_r be an ordered
tuple of its vertices, with r>=1 and repetitions allowed. Construct R(T;p):

1. Append a universal sink a.
2. For j=1,...,r append z_j, with z_j -> p_j and v -> z_j for every other
   vertex v already present.

Thus an old vertex v has outdegree

    d_T(v) + 1 + r - multiplicity(v),

a has outdegree r, and z_j has outdegree 1+r-j.

**Theorem.** The unlabelled isomorphism class of R, together with r, determines
the isomorphism class of (T,p_1,...,p_r). Conversely a port-preserving
isomorphism of the inputs extends to their registers. The encoding costs
exactly r+1 added vertices.

**Decoder proof.** For r>=2, the only vertex of outdegree one is z_r, except
possibly an old vertex v which is a universal sink of T and is selected
at all r ports. In that exceptional case the two degree-one vertices are
z_r and v, with z_r -> v; their induced source identifies z_r uniquely.
Its unique out-neighbor recovers p_r. Delete z_r and repeat.

At r=1, a and z_1 have degree one. If p_1 is not a universal sink of T,
they are the only such vertices; a -> z_1 identifies z_1 as their induced
sink. Record its out-neighbor, delete z_1, then delete the unique universal
sink a. This recovers T and p_1.

If p_1 is a universal sink, the degree-one vertices are exactly
p_1,a,z_1. They form a terminal directed triangle. Contract that entire
triangle to the recovered vertex p_1; every outside vertex dominates it.
Any previously recorded port in the triangle maps to that same recovered
vertex. No equivariant choice of one triangle vertex is needed.
The construction and all these steps commute with isomorphism. QED.

The one-dent construction without a is not a universal decoder of an
unmarked pointed core: denting the source of the transitive three-tournament
and denting a vertex of C3 both give the same strong four-tournament.
The temporary sink repairs this particular ambiguity. We do not claim the
r+1 bound is optimal among every possible register design.

**Exact symmetry loss.** Let G=Aut(T;p_1,...,p_r). The register group is G,
except when p_1 is the universal sink of T and no later p_j equals p_1.
In that case it is G x C3. Indeed the decoder fixes every appended tail
vertex; the only remaining kernel consists of the rotations of the terminal
triangle. A later port at p_1 singles out its actual triangle vertex and
kills that kernel. Every element of G extends by fixing the added vertices,
so the quotient splits and the surviving C3 commutes with G.
This is a class decoder, not an unqualified automorphism-exact codec.

## J2. The exceptional-sink menu has a discrete minimum

A state in X_n is a tournament on positions 0,...,n-1 containing the forward
Hamilton path. There are 2^((n-1)(n-2)/2) such literal presentations.

Start with the previous singleton/native-complement/fresh-re-marking language.
Add D_i: append a vertex beating only the selected current position i and
losing to every other current vertex. It is legal as a Hamilton extension
when i is not the current final position. The observations include the
ordinary unlabelled tournament class; the old blue flag may also be retained.

**Theorem.** For this menu, with arbitrary choice of the index i at each
step, behavioral equivalence on X_n is equality of the literal relation.
Its smallest deterministic carrier is therefore the discrete X_n itself.
Any two unequal states are separated after at most three added vertices.

**Proof.** Choose an old pair i,j on which the states differ. Append one
universal sink, then D_i, then D_j. Both selected indices are old positions,
so both additions are legal Hamilton extensions. J1 decodes the final
ordinary class into (T,i,j), including the direction of the old arc.
Thus the final classes differ. Conversely the literal relation computes
every permitted operation and every fresh-reset law exactly. QED.

This is an all-size lower bound on required information, not just a sufficient
representation. At n=7 the minimum rises from the inherited32272 cells to
32768 literal states. For uniformly distributed X_n this requires
((n-1)(n-2)/2) bits in any lossless fixed-length state code. Particular
presentations or restricted families may still have shorter descriptions.

The index choice is an explicit operation input, costing up to log2(n)
bits to specify. A fixed successor selector or a fixed boundary-width menu
does not satisfy this hypothesis. The old minimum is not contradicted.

## J3. Genuine C3 substitutions already force the same minimum

Permit substitution of C3 for any selected current vertex, with all external
arcs uniform as inherited from the old vertex. Choose its internal order
0->1->2->0, so a Hamilton presentation remains a Hamilton presentation.
The selected position is movable, and substitutions may be iterated inside
previously inserted modules.

**Lemma.** Every strongly connected module M of a tournament is a strong
module: another module cannot overlap M without containing it or being
contained in it. To prove this, suppose N overlaps M and contains an
outside vertex v. That vertex has one uniform direction to M. The module
property for N forces every arc between M intersect N and M minus N to
have that same direction. This contradicts strong connectivity of M.

Build an odd-order strong tournament A_m, m>=3, by repeatedly replacing
the first vertex by C3, starting from a singleton. Its identity order is
a Hamilton path and it is strongly connected. It has automorphism group
C3: at each nonbottom C3 quotient, the larger child is fixed by size,
forcing the two singleton children to be fixed; only the bottom triangle
can rotate.

For n>=2 choose

    p=n+1+(n mod2),    q=2n+1.

Both are odd, p>n, and q>p+n-2. For two distinct old vertices i,j, substitute
A_p at i and A_q at j. This is a word in the single allowed C3 operation.
The final order is

    n-2+p+q = 4n+(n mod2).

**Decoder.** The q-vertex bag is a strongly connected module and hence
cannot overlap any other module. It is the unique module of size q:
a disjoint candidate would lie in the remaining p+n-2<q vertices, a
properly contained candidate would be smaller, and a containing candidate
larger. Contract it, retaining that contracted vertex as the second port.

Now the p-bag is likewise the unique p-vertex module, because only n-1<p
vertices lie outside it. Contract it as the first port. The original
tournament and the ordered pair i,j are recovered.

Consequently two unequal Hamilton presentations can be separated by a
C3-substitution context of order at most4n+(n mod2). This proves that the
C3-only arbitrary-position substitution language already has the discrete
minimum X_n; adding native complement and fresh re-marking cannot merge it.
There is no claim that this size bound is the least possible context.

The two recognized modules give

    Aut(output) = Aut(T;i,j) x Aut(A_p) x Aut(A_q),
    |Aut(output)| = 9 |Aut(T;i,j)|.

The quotient recovers the rooted core group; the two internal C3 factors
are a known sidecar. Actual module searches test the decoder through order4
and on the stated selected order5 controls. The all-size assertion uses
the strong-module proof, not a finite extrapolation.

## J4. The seven-vertex selector, with costs stated separately

Use the inherited presentations A=1072 and B=4176, with rows

    A=(66,12,41,17,39,75,30),
    B=(66,20,9,51,37,71,30).

Their old singleton carriers agree forever. This is inherited from
[minimal-blue M1-M5](../minimal_blue_selectors_20261028/THEORY.md), which
already proves that one extra interior marked port is minimum in port arity.

Here are actual attachment realizations with ordinary class observations:

| Operation at position1, after native complement | Added vertices | H on the two outputs |
|---|---:|---:|
| Exceptional sink D_1 | 1 | 307,293 |
| Substitute C3 | 2 | 777,981 |

Different Hamilton-path counts certify nonisomorphic outputs. A Boolean
readout H=307 suffices for the first context. Within this additional-attachment
cost, one added vertex is the least positive cost. For a nontrivial strong
tournament atom, order3 is minimum, so C3's two net added vertices are sharp.
This does not declare a universal weakest Boolean predicate or compare
port arity to arbitrary information-bit costs.

At position3 the C3 outputs are nonisomorphic even though their H counts
agree, both before and after native complement. Thus H is a convenient
separator at the chosen successor port, not a complete attachment observer.

There is an exact path identity for a dent at p. Write h(S) for the number
of Hamilton paths of T[S], h(empty)=1, and h_p^start / h_p^end for endpoint
conditioned counts. Then

    H(D_p T)
      = H(T) - h_p^end(T)
        + sum_(S subset V minus {p}) h(S) h_p^start(T[V minus S]).

If the new vertex is last, its predecessor cannot be p, accounting for
the subtraction. Otherwise split the path before the new vertex; its next
vertex is necessarily p, and every possible previous endpoint dominates it.
This proves the formula. Omitting the endpoint subtraction fails already
on the one-vertex core; that failed first probe was repaired before use.

## J5. Spectral blindness and observation-algebra saturation

For a fixed core T with acting group G=Aut(T), attach r named vertices with
their internal tournament fixed and let all nr cross arcs vary. Mark old
and new roles. The action on cross coordinates is r copies of the vertex
permutation action of G, with no orientation inversions. Hence

    F_T,r(z) = (1/|G|) sum_g product_(vertex cycles C of g)
                                      (1+z^|C|)^r.

It depends on the permutation action of G, not on all arcs of T.
Naming old ports replaces G by their joint stabilizer and changes nothing
else in this formula. The previous signed-spectrum theorem applies to
this particular unsigned cross-coordinate action.

Both seven-vertex cores and their native complements are rigid. Therefore

    F_A,r(z)=F_B,r(z)=(1+z)^(7r)

for every r, even after naming ports. Actual chosen attachment outputs
still differ by J4. This is a theorem about these marked attachment spaces;
it does not count distinct unmarked output tournaments after forgetting
the embedded core.

The81/27 hostile gives the complementary test: its automorphism actions
differ, and F(1) gives176 versus512 one-vertex attachment orbits. It still
has67 identical labelled cyclic triples and462915 Hamilton paths each.
Thus one hostile tests spectral sensitivity to symmetry, the other tests
spectral blindness to the underlying relation.

There is a positive operator conclusion. For each free old pair i,j, J2
or J3 provides a context whose output-class readout is its orientation bit.
These bits generate every function on X_n under sums and products.
Over characteristic different from2, their signs generate all Walsh
characters, so the resulting observation algebra has dimension |X_n|.
All Fourier modes are observable in this richer language. A trace or an
eigenvalue-multiplicity polynomial alone cannot replace those readout maps.

## J6. Comparison with the old K and what remains open

The old minimal carrier
(E(T),E(native complement T),B(T)) is still exact for singletons, native
complement and fresh uniform/counted re-marking. Its blue interface has
4n-8 states for n>=5, by the inherited theorem. Fixed boundary ports fail
on the padded pair at order2w+5; arbitrary position selection is stronger.

J2 and J3 settle the minimum for the two explicit enriched menus above.
They do not classify every intermediate menu: a bounded set of allowed
positions, bounded substitution depth, position-oblivious attachment, or
a weaker output statistic can preserve nontrivial compression. Finding the
least observation-preserving state there remains a meaningful question.
The finite horizon bound n+3 in J2 and4n+(n mod2) in J3 concerns distinguishability,
not a polynomial-time canonical-isomorphism or enumeration algorithm.
