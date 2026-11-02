# A complement-reciprocal carrier for square-sum path defects

**Status: PROVED elementary identities; FINITE-EXACT small controls.** This
is an explicit representation, not a new proof of Gerbicz's existence
theorem or an equivalence with the raw tournament omissions7,21.

For a simple undirected graph G on n>=1 vertices, define

    P_G(z)=sum_(vertex permutations pi) z^b_G(pi),

where b_G(pi) counts consecutive pairs of pi that are NOT edges of G.
This retains the missing-edge mask as a formal weight. Unlike an ordinary
tournament completion, it records whether its displayed spanning path
actually uses allowed graph edges.

## 1. Exact maps and the information they preserve

Let A be the fully fair weighted tournament, Aij=1/2 for i!=j. Give each
arc a separate multiplier1 on a G-edge and z on a nonedge. Its weighted
Hamiltonian sum is exactly2^(1-n) P_G(z). Equivalently this is the expected
guarded Hamiltonian sum in a uniformly independently oriented tournament.

At z=1 every graph collapses to the SAME fair tournament and
P_G(1)=n!. At z=0 the constant coefficient counts ordered Hamiltonian
paths of G. The matrix at z=0 is generally NOT a tournament: a missing
edge has total weight zero across its two directions. The guard cannot be
discarded while claiming to preserve Hamiltonicity.

For n>=2, [z^0]P_G is twice the number of unoriented Hamiltonian paths;
for n=1 it is1. The exact complement identity is

    P_(complement G)(z)=z^(n-1) P_G(1/z).                (1)

Each consecutive pair changes from allowed to forbidden or vice versa,
so its exponent changes from b to n-1-b. This is a polynomial identity,
including at z=0 after expansion. It is not a Riemann-zeta functional
equation. All coefficient identities survive reduction modulo every
integer modulus; positivity/nonvanishing as an INTEGER does not.

## 2. The lowest power is a minimum path-cover defect

A spanning linear forest is a spanning subgraph whose components are
paths, including isolated vertices. For such a forest F let k(F) be its
number of components and s(F) its number of isolated vertices. Define

    C_k(G)=sum_(F with k components) 2^(k-s(F)).

Then

    P_G(z)=sum_(k=1)^n k! C_k(G) z^(k-1)(1-z)^(n-k).    (2)

**Proof.** In the product along a permutation replace each adjacency
factor by z+(1-z)1_(edge in G) and expand. The selected allowed edges
form a spanning linear forest. Its k components can be ordered in k!
ways and each nontrivial path oriented in two ways. It has n-k selected
edges and k-1 unselected connecting positions, giving its stated term.
Conversely every such oriented ordered forest gives exactly one expansion
term. This proves(2), without a disjointness assumption between different
forests.

Let pc(G) be the minimum number of paths in a spanning path cover.
The first nonzero term in(2) is at k=pc(G). Therefore

    ord_(z=0) P_G = pc(G)-1,
    [z^(pc(G)-1)] P_G = pc(G)! C_(pc(G))(G).             (3)

The same valuation is the minimum number of missing edges that must be
used by a spanning vertex ordering. Cutting at missing edges proves one
inequality; ordering a minimum path cover proves the other. If two ends
of different minimum-cover paths were adjacent they could be joined,
contradicting minimality. This also explains the leading coefficient.

Thus Hamiltonicity is the zero-defect stratum of a universally defined
carrier. It is useful before a Hamiltonian path exists, rather than only
after choosing one as observer order. Formula(2) is a nonnegative
Bernstein-type expansion, but computing all its coefficients is not
claimed efficient.

## 3. The three four-vertex three-edge shapes

The star K_(1,3) and triangle-plus-isolate are complements. Both have

    P(z)=12z+12z^2,

whereas the four-vertex path has

    P(z)=2+10z+10z^2+2z^3.

The first pair has the same unmarked path-cover profile and defect1;
the latter is Hamiltonian and self-complementary. This gives an exact
meaning to the user's complementary-shape intuition while also displaying
a limit: equality of this polynomial is not graph isomorphism or equality
of all marked attachment responses.

For an explicit port hostile label the star center0 and its leaves1,2,3;
in its complement,0 is isolated and1,2,3 form the triangle. Attach a new
leaf4 to port0 in each graph. The first becomes K_(1,4), with pc=3;
the second becomes K3 disjoint-union K2, with pc=2. Thus equal unmarked
polynomials can separate under the SAME port operation. Exact gluing must
retain marked linear-forest boundary states, as already required by the
overnight packet l's correction to the cyclic product rule.

## 4. What the square-sum exceptions mean in this carrier

Apply this to Q_n, whose vertices are1..n and whose allowed edges satisfy
x+y=square. Inherited [THM-4851](../../../01-canon/theorems/THM-4851-the-prehistory-of-the-square-sum-numbers-unique-cycles-at-32-and-33-terminal-spectrum-with-the-hole-at-39.md)
records Hamiltonian paths exactly at n in{1,15,16,17,23} union[25,infinity)
(the infinite existence claim is CITED/reproduced there). Gerbicz also states the result directly in his January21,2018 comment on [OEIS A090461](https://oeis.org/A090461), checked again November2,2026; the linked original forum endpoint was inaccessible in this session. In this carrier,
those are exactly the orders where ord_0 P_(Q_n)=0.

The older [THM-4510](../../../01-canon/theorems/THM-4510-sum-graphs-as-reflection-orbits-pythagorean-zigzags.md)
already explains15,16,17 through the tight three-reflection system with
targets9,16,25; THM-4505 gives the exact degree obstruction and unique
15-path. The exceptional order23 and failure24 require further path/port
constraints, not a missing tournament completion. A tournament completion
exists at EVERY order, and its unguarded Hamiltonian existence is automatic.

The sets are differently typed:15,16,17,23 are graph VERTEX ORDERS with
allowed paths;7,21 are forbidden RAW PATH COUNTS across all tournaments.
No map transporting one forbiddenness predicate to the other is supplied
by their small size. The common mechanism that survives is exact
incidence-constrained path assembly, retaining the mask and ports.

## 5. Relation to the other lanes

FRACTRAN likewise has a total priority order before any computation is
performed. What chooses the next rule is the divisibility guard, including
the failure of earlier guards. The incoming
[least-clock theorem](../stability_priority_spectrum_20261102/FRACTRAN_CLOCKS.md)
and Session135's word acceleration already give exact arithmetic versions
of that distinction. A run can revisit rules; it is not a Hamiltonian
path through its instruction list.

Weighted stability uses a different retained coordinate: the common
rounding baseline B and assigned light mass U. The314 proof succeeds
because each survives long enough to pay a later cost. The analogy is the
operation of controlled forgetting, not a reduction among these problems.

Two concrete next questions: can the square-sum lift propagate a specified
low-defect marked linear-forest state, not just a path-existence bit; and
can the stability proof use its remaining local square terms to exclude
the m8 scalar cell at312? These are separate, testable extension problems.
