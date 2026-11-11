# Session260: arithmetic sources, physical seed moves, and domination defects

**11 November 2026. PROVED / FINITE-EXACT / CONDITIONAL / OPEN as individually
labelled below. Chan–Pak and quasiperfect existence remain OPEN.**

This session joined the requested historical sequence with the strongest
current tournament constructions, then pursued the operations that the
connection actually permits. It moved through divisor sums, modular
substitution, marked seed arithmetic, local surgery, graph domination, and
quadratic norms. The work is not a claim that many suggestive formulas add
up to a complete proof.

## 1. The strongest new structural bridge

For any tournament T and marked vertex z, put A=N⁺(z), B=N⁻(z). Make the
bipartite graph G_z with shores A,B and edges a—b exactly when a→b in T.
Write H for Hamilton-path count and d for dominating-set count.

**PROVED:**

    H(T) ≥ d(G_z), with equality if and only if T−z is transitive.

This is a literal injection. A dominating set S determines
L=(A\S)∪(B∩S) and R=(A∩S)∪(B\S). Domination forces the terminal strong
component of T[L] to meet B and the initial strong component of T[R] to
meet A. Choose endpoint paths and concatenate L,z,R. The path recovers S.

Splitting all Hamilton paths at z gives the stronger bookkeeping identity

    H(T)−d(G_z) = Σ_L (e_B(L)s_A(V\L)−ε(L)),

whose **individual terms are nonnegative**. Equality first forces both
shores transitive through the all-selected dominating fibre. A remaining
mixed triangle forces a fibre with at least two paths. Thus equality is
exactly the transitive-deletion locus, not a finite observation.

For a fixed bipartite graph, even arbitrary choices of internal shore
orientations can give zero defect **if and only if the graph is a chain
graph**. An induced 2K2 is the smallest structural obstruction. Its count9
lifts to13 or15 at this interface, never9.

This supplies a real map from the solved domination side toward Chan–Pak,
and identifies precisely what it fails to preserve. A future transfer must
control the positive defect, replace the source graph by another graph of
the same count, or change the interface. Shore sorting alone cannot do it.

Proof: [BIPARTITE_DEFECT.md](bridges/BIPARTITE_DEFECT.md).
Focused review: [review/BIPARTITE_DEFECT.md](review/BIPARTITE_DEFECT.md).

### A second payoff: reverse path capacities form a circulation

Let U=T−z. For each original arc u→v, let L_uv count old Hamilton paths
using that adjacency and R_uv count paths gained by reversing it. Then

    H(T)−H(U)=Σ_(u∈A,v∈B,u→v) (L_uv+R_uv).

This **PROVED positive directed-cut law** is supported only on original
arcs. Compare it with Session258's raw even-path endpoint capacities c_uv.
The exact edge weight is2(c_uv+c_vu), and the reverse capacity c_vu placed
on each original arc u→v is a nonnegative circulation. This explains how
the old nonarc terms redistribute without changing any newborn count.
On a transitive host the weights are2 at adjacent vertices and2^(d−1)
at distance d≥2, recovering the chain and apex laws directly.

The stronger conjecture H(T)−d(G_z)≥H(U)−1 remains **OPEN**. The positive
law turns it into a precise injection problem for marked path covers; its
finite survival is not promoted to proof. A selected-fibre lower bound
and spectator-product transport are proved independently.
See [DEFECT_QUANTITATIVE.md](bridges/DEFECT_QUANTITATIVE.md).

## 2. What20,104,464,650,1952 actually mean

The recovered repository sources place these in the **excess-two** equation
σ(n)=2n+2. Quasiperfect numbers instead satisfy σ(n)=2n+1. The Chowla
function c(n)=σ(n)−n−1 therefore sends these five terms to

    21,105,465,651,1953.

The distinction matters:21 is an inherited forbidden tournament count;
105,465,651,1953 now have explicit constructors in this packet. The old
classification through divisor count17 and its next term130304 are retained
as inherited results, not re-audited or silently broadened.

For the formal family n_a=2^a(2^(a+1)−3),

    n_a+1=(2^a−1)(2^(a+1)−1).

The primality guard on2^(a+1)−3 is needed for the divisor-sum source, not for
the graph constructor. The whole successor family is decided for a≥2:
a=2 gives the genuine21 exception, a=3 has a literal105 tournament, and
every a≥4 follows by the ordinal product of two inherited Mersenne fans.
The separate650 branch has a concrete651 construction.

The cyclic chamber T(a,b)=C3(TT_a,TT_b,1) has

    H(T(a,b))=(2^a−1)(2^b−1)+2.

Reversing its last-A/first-B edge yields the physical shifted-product law

    4H−11=(2^(a+1)−3)(2^(b+1)−3).

This is a specified graph operation and identity, not merely arithmetic
factorization. The general assertion that every even excess-two source
other than20 has a realizable successor remains **OPEN**.

Details and historical paths: [arithmetic/RESEARCH.md](arithmetic/RESEARCH.md).

## 3. Prime and Carmichael arithmetic really acts on tournaments

Replace a vertex v of a rooted tournament by a transitive m-vertex module,
m≥2.
Its actual path-cover expansion is Σ_r a_r r!S(m,r), with a₁=H(T).

**PROVED:** the resulting H is congruent to H(T) modulo m for every rooted
tournament exactly when m is prime or Carmichael. It suffices to test a
specified family of strong rooted tournaments. Their triangular skeleton
coefficients recover the full prime-divisor criterion and exclude repeated
prime factors. The weak C3 test alone would incorrectly admit341; the
third strong test rejects it. Carmichael561 survives the universal test.

This is a genuine arithmetic classification of a physical substitution.
It is a congruence theorem, not an exact-count coverage theorem. Cyclic
roots also obey an exponential growth lower bound; a congruence-preserving
large clone cannot be treated as a small numerical correction.

Proof: [bridges/RESEARCH.md](bridges/RESEARCH.md).

## 4. Physical multiplier arithmetic and its missing return operation

The current seed mechanism retains H−W=1 and U=q, not only H. Six new
all-length strong seed families are certified. Two have the particularly
useful physical transitions

    (H,U) → (2H+1,3U+2),     (H,U) → (2H+1,3U+6).

Selected-clone placement matters: marking the first or second clone is
legal in the classified families, whereas an arbitrary interior mark is
not. An exact finite atlas has627 presentations and212 odd slopes, with
three explicit rays whose below-start heads are retained.

The larger two-marked-clone grammar has state(q,δ), δ∈{0,1,3}. Its physical
left and right insertions are

    L(q,δ)=(2q+2δ,0),     R(q,δ)=(3q+2δ,δ).

The selected gap is exactly0 or1. The twelve resulting odd-slope sequences
are classified and sparse. There is no free binary or ternary alphabet:
L erases the state needed for a corrected next digit.

**CONDITIONAL bridge to coverage:** a finite seed bank plus a specified
same-q rearming operation to δ=1 would give all sufficiently large even
multipliers by binary induction. A further physical deflation would then
give the missing exact multipliers. These operations are not supplied;
small universal versions are false. The statement identifies an actual
remaining constructor rather than replacing it with an arithmetic wish.

Proofs: [constructors/RESEARCH.md](constructors/RESEARCH.md) and
[coverage/PROOF.md](coverage/PROOF.md).

The final incoming-work integration coarsens actual adjacent twin blocks,
reducing the enlarged grammar to roots q=29,47. Besides L and R it has a
one-use gap birth G:q→(13q+8δ)/5. All three directions commute on their
legal domains. The whole carrier has δ>0 exactly on its odd-q boundary;
its even-q interior has δ=0. Therefore same-q rearming is impossible
**within this carrier**, even after the enlarged presentation. Reversing
the last L gives60→29, not60→30. This is an exact boundary obstruction,
not an all-tournament exclusion. The eight unified slope sequences absorb
the earlier twelve-sequence presentation and its additional boundary heads.
See [INCOMING_CONNECTIONS.md](coverage/INCOMING_CONNECTIONS.md).

## 5. When local repair fails, reconstruction can still succeed

For T(a,a+1), the unmodified count is n_a+3. The desired successor is two
lower. **PROVED:** every decreasing single-edge reversal drops at least
3·2^(a−1)−2, for every a≥2. The crossing-edge theorem and exact apex-edge
formulas establish this globally, not only at small a.

The stronger **all-apex theorem** permits every orientation of every edge
incident to the apex. The maximum is n_a+3 and its next distinct value is
lower by the same gap. Hence none of these2^(2a+1) choices makes n_a+1 at
any a≥2. A binary-word sorting proof supplies the all-size statement.

The multi-edge exploration is deliberately separate:

- **FINITE-EXACT:** complete two-edge shells for a=2,...,6 hit the target
  only at a=3, in exactly two ways.
- A specified three-edge exploration and selected fourth-edge continuations
  find a four-edge465 repair at a=4.
- **PROVED:** the natural extrapolated four-edge rank pattern has an exact
  exponential formula and misses the target for every a≥5. Its success at
  a=4 is isolated in that pattern.

These obstructions do not refute the successor family: section2 already
constructs every allowed member by changing the architecture. They explain
why continuing to flip the same few edges is the wrong uniform proof.

Proofs: [LOCAL_RIGIDITY.md](LOCAL_RIGIDITY.md),
[CROSSING_FLIP.md](bridges/CROSSING_FLIP.md), and
[TWO_EDGE_REPAIR.md](constructors/TWO_EDGE_REPAIR.md).

## 6. The domination route has its own arithmetic and exceptions

For a binary word in A,B, connect each A to later B vertices. The chain
graph and its apex tournament have exactly the same count. If x is the
binary prefix with A=1,B=0 and h is the count, append operations are

    A: (h,x)→(h,2x+1),
    B: (h,x)→(h+x+(x mod2),2x).

This gives a count-preserving compiler with a necessary second register.
It does not cover every allowed count:29 is missing from all chain graphs
and their disjoint unions, with an all-order size-bound proof.

Moving beyond chains produces another exact result:

**PROVED WITH FINITE EXACT BASE:** over bipartite graphs of arbitrary
order, the domination spectrum below225 is every positive odd number
except7,13,21,43. Every attained value has an explicit witness on at most
nine vertices. The finite base is complete because a rooted-tree state
compiler proves every ten-vertex connected graph has at least225
dominating sets; larger connected graphs inherit that lower bound.
Disconnected components are handled by products. In particular43 is a
true all-order exception, not a search hole.

Lawful moves retain boundary data. A leaf at a nonisolated root can never
implement d→2d+1. False twins can implement it on the genuine family
3·2^m−1, but keep its fixed coefficient. Shore hubs have a closed
six-coordinate interface and an explicit two-shore affine-bilinear family.
The forced-all clone cases are essential; dropping them falsely creates
forbidden counts.

These are constructive positive families, not just negative tests. But
the bipartite support problem and the positive tournament-defect problem
are separate, and neither scalar recurrence solves the other.

Proofs: [DOMINATION_TRANSFER.md](bridges/DOMINATION_TRANSFER.md) and
[BIPARTITE_COVERAGE.md](constructors/BIPARTITE_COVERAGE.md).

## 7. Quasiperfect numbers contribute a real norm interface

For a nonempty tournament A, let Q₂(A) count ordered spanning two-path
covers, allowing empty slots. Put x=Q₂−H and y=H. Literal cyclic
substitution gives

    H(C3[A,B,1])=x(A)x(B)+2y(A)y(B),
    H(C3[A,A,1])=x²+2y².

The inherited odd-square necessity for a hypothetical quasiperfect n=m²
puts σ(n)=1+2m² in the same quadratic norm. Equality of norm values does
not identify their coordinates or prove the existence of either source.

**PROVED:** cutting Hamilton paths at every gap gives Q₂≥(|A|+1)H, hence
x≥|A|y. This supplies a complete finite decision procedure for every
target in the diagonal cyclic grammar: enumerate norm representations,
then only graphs of order at most⌊x/y⌋. A prime target1+2m² cannot arise
there for m>1. This prime barrier does not obstruct quasiperfect existence,
whose divisor sum would already be composite.

Switching a conjugate factor of norm3 requires Q₂=5H±3. The cut bound
reduces its entire physical classification to a justified finite base:
only the singleton, TT3, and C3 occur, giving outputs3,51,99. Pell
continuation first fails at the exact missing profile H=1,Q₂=42.
Automorphism content and inert-prime valuation parity retain further
structure, but no quasiperfect nonexistence conclusion is claimed.

Proof: [QUASIPERFECT_NORMS.md](arithmetic/QUASIPERFECT_NORMS.md).

## 8. Integration with the moving frontier

The session started from Session258's endpoint seed synthesis and retained
the inherited odd-cycle/path-cover machinery rather than re-auditing it.
Concurrent main commits were preserved, read, and used as mathematical
signal:

- 48ffb162415d completed the residual d=3 and d=21 grids and classified
 stable-port families. Those completed heads are credited there, not
 counted here as newly discovered numerical coverage.
- 8cc23df02699 proved arbitrary finite strict-seed closure and eventual
 upper-frontier stabilization in the summand/multiplicand graph. That
 solves the arithmetic closure question once its operations exist; it
 does not install our missing tournament rearming operation.
- d0383cbc2486 added the d=35 residual family at exponent1 and exponents≥6,
 leaving its four uniform heads open, and interpreted our selected-clone
 moves as activation of two seed birth roots. The source/mark information
 remains essential; scalar coefficient agreement alone is insufficient.

The later incoming e4173071bc62 strengthens the arithmetic continuation to
eventual additive saturation and a two-seed product saturation classification.
It is preserved at close-out; it still does not implement a marked-tournament
return.

No numbered canon claim or external novelty claim is made in this packet.
New proofs, inherited inputs, bounded experiments, conditional reductions,
and open conclusions are kept distinct. The remaining targets are stated
concretely in [OPEN_QUESTIONS.md](OPEN_QUESTIONS.md).
