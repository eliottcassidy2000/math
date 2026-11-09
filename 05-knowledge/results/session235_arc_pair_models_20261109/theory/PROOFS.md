# Pair-state bridges, parity, and the missing-pair boundary

Status: **PROVED**, with separate **FINITE-EXACT** checks. Local labels PB1–PB7
belong only to this file. No novelty claim is made. All digraphs are finite,
loopless, labelled unless stated otherwise. H counts ordered directed
Hamiltonian paths; H(empty)=1. Diagonals may be added freely without affecting H.

Inheritance: Session234 PROOFS.md DG3/DG6/DG7 and GENERAL_GRAMMAR.md give the
semicomplete SCC product, weighted tournament resolutions, and bounded strong
cores. The hostile triangle with an unused reverse arc forbids assuming strict
monotonicity. The least-used sidecar here is the generating polynomial recording
which ambiguous pairs a path actually uses.

## PB1. The three meanings of “11 means no edge”

For i<j the raw code 10 means i→j and 01 means j→i.

1. If 00 is forbidden and 11 means empty, the three codes describe precisely
   oriented graphs (at most one arc per pair), bijectively. Raw semicomplete
   matrices are being **decoded**, not used as adjacency matrices.
2. If both 00 and 11 mean empty, the decoded adjacency is
   B[i,j]=A[i,j](1-A[j,i]). This again describes oriented graphs, but is not
   injective: an oriented graph with m empty pairs has exactly 2^m raw encodings,
   before independent diagonal choices. Uniform raw codes do not give the
   uniform three-state oriented-graph model.
3. If 00 means reciprocal and 11 means empty, while 10/01 keep their directions,
   the decoding is the involution τ(A)=J-I-A^T off the diagonal. It exchanges
   empty and reciprocal pairs, fixes single arcs, and describes **all** digraphs
   bijectively. It agrees with case 1 on raw semicomplete matrices, but not on
   arbitrary raw matrices.

In cases 1/2 a tournament is exactly a graph with no decoded empty pairs.
Reinterpreting symbols does not preserve raw H in general: the two-vertex raw
code11 has H=2 as adjacency and H=0 as an empty-pair code. A one-pair 00 code
has H=0 under case2 and H=2 under case3. The diagonal remains irrelevant only
because H is defined on distinct-vertex paths, not because 11 has a diagonal
interpretation.

## PB2. One polynomial controls deletion, reciprocal lift and completion

Let O be oriented, with m missing unordered pairs M. Let L add **both** arcs
at every pair in M. Thus L is semicomplete. Define the multiaffine polynomial

P_O(z)=Σ_{Hamiltonian paths p in L} ∏_{e∈M used consecutively by p} z_e.

A vertex-simple path uses an unordered pair at most once. Consequently

H(O)=P_O(0),             H(L)=P_O(1),
Σ_{tournament completions T of O} H(T)=2^m P_O(1/2).

Proof: a path using k missing pairs dictates k completion choices and is
compatible with exactly 2^(m-k) of the 2^m completions. Unused pairs remain
free. This is the same double count as Session234 DG6, now viewed from below.
Every completion contains O and lies in L, hence

H(O) ≤ min_T H(T) ≤ average_T H(T) ≤ max_T H(T) ≤ H(L).

Moreover the intersection of their Hamiltonian-path sets is exactly the path
set of O, and their union is exactly that of L. A path using a missing pair is
excluded by a completion choosing its reverse; every L-path can be extended
to a completion because it cannot demand both directions of one pair.

Conversely, deleting **both directions** of all reciprocal pairs of any
semicomplete L produces O. Thus these are two endpoints of the same bridge,
not unrelated recipes. A tournament resolution of L is a completion of O.
Scalar resolution counts alone lose which paths overlap. Neither their sum
nor their average generally equals H(O) or H(L).

## PB3. Inclusion-exclusion and exact response

Let N(S) count paths in L that use every pair in S⊆M consecutively, in either
direction. Then

H(O)=Σ_{S⊆M} (-1)^|S| N(S),
N(S)=(∂_S P_O)(1).

This follows by inclusion-exclusion on forbidden pair-use events, or by
expanding P_O at 1 and evaluating at 0. More generally deleting just the pairs
R⊆M replaces z_e by0 for e∈R and leaves other coordinates1. Coefficients of
P_O give nonnegative exact counts by the set of used ambiguous pairs. For a
single deleted reciprocal pair e the loss is precisely N({e}); for several
pairs one must retain overlap terms. This is why separately measured losses
cannot simply be added.

The ordinary directed path polynomial
F(x)=Σ_{permutations p} ∏_i x_(p_i,p_(i+1))
likewise specializes to every binary adjacency matrix; under case2 decoding
substitute x_ij=a_ij(1-a_ji), and under case3 substitute x_ij=1-a_ji.
These are exact polynomial evaluations, not a statement that raw and decoded
path counts coincide.

## PB4. Complement–reverse parity, for every digraph

For every loopless digraph A,

H(A) ≡ H(τ(A)) (mod 2),          τ(A)=J-I-A^T.

In fact this holds on every induced vertex set. An elementary stronger identity
is available. Let f_c(A) count spanning directed linear forests contained in A
with c unordered nonempty path components, including isolated vertices. Then,
for n>=1,

H(τ(A)) = H(Ā) = Σ_(c=1)^n (-1)^(n-c) c! f_c(A),   Ā=J-I-A.

Proof: expand Σ_p ∏_(consecutive uv in p)(1-A_uv) by inclusion-exclusion.
A selected subset of consecutive arcs of a permutation is a spanning directed
linear forest F with n-c arcs. Conversely F occurs in exactly c! permutations,
one for every ordering of its path components. Group the expansion by F.
Reversing all vertex orders proves H(τ(A))=H(Ā). Modulo2 the terms c>=2 vanish,
leaving f_1(A)=H(A). The empty case is1=1 by convention. This proves parity
on any induced vertex set as well.

For an independent proof retaining all induced counts simultaneously, work
algebraically as follows; neither proof requires the tournament conjecture.

Work in R=F₂[x_1,...,x_n]/(x_1²,...,x_n²), and put X=diag(x_i),
Ā=J-I-A. In characteristic two signs are immaterial. Since every product of
more than n vertex variables vanishes, I+XA is invertible by a finite geometric
series. The expression

G_A=1+1^T(I+XA)^(-1)X1

has coefficient H(A[S]) mod2 at x_S=∏_(i∈S)x_i. Indeed, expansion enumerates
all directed walks weighted by their visited vertices; walks repeating a
vertex vanish, leaving exactly directed paths on each nonempty S. Its constant
term1 represents the empty path.

The matrix determinant lemma over the commutative ring R gives

G_A=det(I+XA+XJ)/det(I+XA).

Entrywise, I+XA+XJ=diag(1+x_i)(I+XĀ): off-diagonal entries agree since
(1+x_i)x_i=x_i, and diagonal entries are1+x_i. Therefore

G_A= [∏_i(1+x_i)] det(I+XĀ)/det(I+XA).

Every element q of R with constant term1 has q²=1: in characteristic two
cross terms vanish, and every nonconstant squarefree monomial squares to0.
So each determinant in this formula is its own inverse. The right side is

[∏_i(1+x_i)] det(I+XĀ) det(I+XA),

symmetric in A and Ā. Thus G_A=G_Ā. Extract x_S, then use the bijection reversing
all vertex orders to obtain H(Ā[S])=H(Ā^T[S])=H(τ(A)[S]).

Consequences:

- For an oriented O and its reciprocal lift L=τ(O), H(L)-H(O) is even and
  nonnegative. For semicomplete raw A, its 11-empty decoding is τ(A), so its
  decoded count is at most its raw count and differs by an even integer.
- For a tournament τ(A)=A. The congruence is then tautological: **this proof
  alone does not prove Rédei's oddness theorem**.
- On unrestricted raw codes, case2 deletion is not τ, so no such parity
  assertion follows. The graph with sole reciprocal pair12 and sole arc23
  has raw H=1, whereas deleting12 leaves H=0. Explicitly adjacency
  010/101/000 has the
  unique Hamiltonian path012 and decoded H=0. This is a hostile control.
- No order inequality is asserted between A and τ(A) in the unrestricted
  case, since τ both adds and removes arcs.

The congruence is **Berge's complement-parity theorem**, followed by the
path-reversal bijection. The two displayed proofs are sufficient for its use
here; no source priority is claimed. Current canon
[THM-5071, statement and proof (3)](../../../../01-canon/theorems/THM-5071-weighted-cycle-formula-for-digraph-path-totals-and-one-reciprocal-pair-exclusions.md)
records a weighted-cycle refinement derived from the cited Grinberg-Stanley
formula. That refinement and the present linear-forest proof retain different
sidecars; neither should be presented as a new parity theorem.

## PB5. Every nonnegative integer occurs for oriented graphs

H=0 occurs on an edgeless pair; H=1 occurs on a transitive tournament.
For h>=3 the directed h-cycle has exactly h Hamiltonian paths, obtained by
cutting its cycle at any one of its h vertices. To get H=2, take two directed
cycles of lengths at least3 sharing exactly one vertex v, with no other arcs.
Removing v leaves two disjoint directed path arms. A Hamiltonian path must
cross between arms exactly once, through v, and must traverse each arm in its
forced direction. There are exactly two choices for which arm goes first.
This graph is strongly connected and oriented.

Therefore the all-order spectrum of oriented graphs, and hence of unrestricted
digraphs, is exactly the nonnegative integers. This is **not** a closure of
the tournament spectrum: missing pairs are essential in these witnesses.
For h=7 or21 a directed h-cycle already supplies a reciprocal-free example.
Ordered sums with transitive tournaments preserve H and give examples at every
larger order. These are existence families, not minimum-order assertions.

## PB6. Arbitrarily large strong oriented cores at every fixed h>=2

A general cut-vertex identity explains the failure of the semicomplete grammar.
Let A and B meet only at v, with no arcs between A-v and B-v. For a lobe A define
p_A as the number of Hamiltonian paths of A with last vertex v, and q_A as the
number with first vertex v; define p_B,q_B similarly. If both lobes have a
vertex other than v, then

H(A∪B)=p_A q_B + p_B q_A.

Proof: after deleting v, any Hamiltonian path can have vertices in at most two
consecutive lobe blocks, with v between them. Since both lobes are nonempty,
v cannot be an endpoint. The two possible orders give the two products.
Conversely each pair of endpoint paths concatenates at v uniquely.

Fix h>=3. Let C be a directed h-cycle and partition its vertices into two
nonempty sets I and O. Add a hub v and arcs i→v for i∈I and v→o for o∈O.
This lobe B is oriented and strong: the cycle connects all its vertices, and
the hub has both an incoming and outgoing arc. Its endpoint counts are
p_B=|I| and q_B=|O|, since every directed-cycle Hamiltonian path has one
specified start and end, each running once through all h vertices.

As other lobe A choose any directed m-cycle, m>=3, sharing only v. It has
p_A=q_A=1. Their union is strong, oriented, has h+m vertices, and has

H(A∪B)=|I|+|O|=h.

For h=2 use two directed cycles sharing v as above, allowing either cycle
length to grow. Thus every fixed h>=2 admits strong oriented witnesses of
unbounded order. In particular **the finite strong-core alphabet fails for
H=7 and H=21 themselves**. This directly locates the semicomplete hypothesis
used in Session234: strong oriented graphs need not have a spanning directed
cycle, nor satisfy H>=|V|. A structural classification based on other,
possibly compressed, objects is not ruled out by this counterexample.

## PB7. What the Chan–Pak question does and does not say

Inherited source: Session234 RESEARCH_CONTEXT.md and PROOFS.md cite Chan–Pak,
*Computational complexity of counting coincidences*, TCS1015 (2024),114776,
§6.5 Conjecture6.2, arXiv:2308.10214v2. The tournament spectrum assertion that
every positive odd integer except7 and21 is attained remains **OPEN**.
The established tournament exclusions7/21 and parity imply that an oriented
graph with either count must have a missing pair. They do not prevent those
counts for oriented graphs. PB5/PB6 settle only the enlarged oriented/digraph
spectra and show why transferring tournament or semicomplete fibre grammars
would be invalid. Tournament completion bounds and parity do not supply a
completion with the same H.

Reproduction: run `python3 theory/check_bridges.py` from the Session235 packet.
The script exhausts orders1–4 for parity and all oriented graphs through4 for
polynomial identities; it checks the sparse constructions independently by
subset dynamic programming and permutation enumeration at small orders.
