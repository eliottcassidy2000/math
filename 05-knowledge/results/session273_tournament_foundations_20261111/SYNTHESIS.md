# Four ways tournament structure carries, transports, or loses mathematics

**Session273, 11 November 2026. Research packet complete. Chan–Pak remains OPEN.**

The productive interpretation of “tournaments form a basis” is not that every
open problem is secretly Chan–Pak. It is that a precise binary relation,
its admissible histories, and its boundary responses can provide a common
computational carrier. This session produces four actual bridges and identifies
where the carrier stops being faithful. No new large numerical census was run.

## 1. Constructive advance: move a resource without changing the observable

Old lonely-runner telescoping suggested allowing local changes rather than
insisting that every building block preserve H. For a strongly connected finite
transition graph, positive edge-weight replacements preserve the product on
**every closed word if and only if** their ratios have the form

    w'_uv / w_uv = λ_v / λ_u.

This is standard graph cohomology used as a finite all-word conservation test.
Prime valuations turn it into integer resource transport. Open words remember
the boundary ratio; closed words cancel it.

The new physical implementation uses Session268's four strong rigid ten-vertex
cells with H-values215,179,245,201. Add equal-order triangle/singleton fillers
so that the off-diagonal weights become

    179·3^k, 245·3^(K−k),          0≤k≤K.

Each cell has10+3K vertices. Compatible circular binary words of lengthN give
exactly

    2+(K+1)(2^N−2) ordinary tournament classes,
    floor(N²/4)+2 distinct Hamiltonian counts.

The complete multiplicities and labelled counts are proved, not inferred from
finite trials. The class quotient retains the linear SCC order: rotations of
the binary word are not silently identified.

More strongly, the gauge preserves the **entire path-cover response**, not just
H. Every uniform-substitution context therefore remains blind to its position.
Yet a marked terminal seam exposes it. Two26-vertex equal-profile inputs with
H=131565 give27-vertex outputs131565 and219275. This uses the known3→5 move;
the advance is the all-order resource transport and its exact invisibility and
exposure laws, not a new claim about that small move.

**Why useful for Chan–Pak:** one can rearrange a successful construction while
keeping its count and all unrooted substitution data, then expose a useful
rooted resource. This is a concrete repair architecture. The still-missing
lemma is renewal: after exposure, replenish the endpoint statistics needed for
the next independently chosen target step. The displayed finite-prime-support
family cannot itself cover arbitrary odd values. No such implication is made.

[Proof, class formulas, generating functions and controls](interfaces/PROOF.md).

## 2. A000568 becomes an exact growth process, but only with a boundary state

Let Q_N mean: choose an ordinary tournament class uniformly, then label it
uniformly. Let P_N be independent fair arcs. Orbit–stabilizer gives

    dQ_N/dP_N = |Aut(T)|/(1+δ_N),
    δ_N=A000568(N) N!/2^binom(N,2)−1.

The Q_N do **not** form consistent marginals. Deleting a random vertex from a
uniform four-vertex class gives cyclic-triangle probability1/4; the uniform
three-vertex class law gives1/2.

Nevertheless a precise finite-horizon growth process exists. Retain

    h_n^(N)(t)=E[|Aut(T_N)| | T_n=t].

Each new incident row has probability2^(−n)h_(n+1)/h_n. The product telescopes
to the desired terminal law. A conditioned signed-pair-orbit Burnside formula
computes h exactly; no claim that this is always fast is made.

Session271's census error envelope now has a new interpretation: uniformly
for every prefix size n≤N, including growing n,

    TV(Q_N restricted to n vertices, P_n) ≤ δ_N/(1+δ_N)→0.

Thus the entire sequence, not merely selected subsequences, has the fair
infinite tournament as its local limit. Its cylinder measure is unique. Almost
every infinite sample is the homogeneous random tournament, with abundant
automorphisms, although finite tournaments are rigid with probability tending
to one. Infinite automorphisms need not preserve any named finite prefix.

This is a complete finite-to-infinite demonstration of the boundary,
consistency, uniqueness and symmetry issues raised by the previous Yang–Mills
session. It is a proof-mechanism bridge, not a map from gauge fields to arcs or
a solution of a Yang–Mills objective.

The infinite tournament contains every finite tournament. Therefore its finite
induced H-spectrum is exactly the unknown global spectrum. But its extension
property manufactures prescribed **arc patterns**, not a pattern with a
prescribed **count**. The note supplies a quantitative witness-amplification
bound once a witness exists, without using that bound to assume existence.

[Conditional growth, exact controls and age/bulk distinction](orbits/PROOF.md).

## 3. Why continuum-looking structure does not determine exact arithmetic

For every real λ≥0 there are n-vertex tournaments T_n such that

    normalized arc-edit distance to TT_n = O_λ(1/n),
    log H(T_n)/n → λ.

Even divergent entropy can occur with distance tending to zero. Thus the same
transitive macroscopic limit supports every prescribed nonnegative finite
Hamiltonian entropy, or infinite entropy. The proof uses actual ordered
components and the elementary random-tournament mean; it does not appeal to a
major conjecture. Triangle chains give the explicit exampleH=3^(n/3) with
onlyn/3 arc reversals from transitive.

For a fixed nontransitive core S hidden by neutral padding, a fixed-size
unlabelled k-sample notices it with probability

    c3(S)(k)_3/(n)_3 + O_(S,k)(n^-4).

Hence independent unlabelled k-sample access requires ordern³ samples for
constant distinguishing advantage. Marked/adaptive/full-matrix access is
outside this statement. Renormalized strong-motif densities and the measure of
strong components restore parts of the lost information, with an explicit
uniform-integrability requirement for varying cores.

This turns the previous contact-data analogy into an exact discrete model:
macroscopic convergence may be perfect while the target information lives in
vanishing-scale defects. It explains why A000568 prediction, graph limits and
Chan–Pak integer realization demand different observers.

[Arbitrary entropy theorem and sharp sampling scale](microscopic/PROOF.md).

## 4. A real bridge to lonely runners and planar sign geometry

At generic time t, integer speeds give the intrinsic tournament

    i→j iff sin(2π(v_j−v_i)t)>0.

Collision and antipodal walls remain boundary symbols, not arbitrary tie
orientations. For a triple with primitive adjacent gaps p,q, exactly one
quadratic Walsh correlation survives. The exact cyclic probability is

- (1+1/(pq))/4 if p,q are odd;
- (1−1/[q(p+q)])/4 if p is even;
- (1−1/[p(p+q)])/4 if q is even.

It lies sharply between1/6 and1/2. The lower extremes are gap ratios1:2 and2:1;
the upper is1:1. This transfers the old arithmetic square-wave decoder into
an exact tournament distribution rather than a generic random model.

For any fixed speed row, the sequence at t,2t,4t,... has an exact finite hidden
Markov presentation using rational phase intervals. Every accepted state path
is physically realizable. The number of lengthK tournament words lies between
2^K andc(V)2^K, so entropy is exactly one bit per scale. The number of states
depends on the speeds, not just the number of runners.

The labelled root-edge signs reconstruct all phases digit by digit. At depthK
they bracket the actual loneliness margin with error≤2^(−K). Every strict
witness is eventually certified; equality walls require a separate argument.
Two cells for speeds(0,1,3,7) have the same instantaneous strong4 class/H5 but
next doubled-time valuesH1 versusH5: class identity does not determine the
continuation.

Real planar rank-two sign tournaments have a binary cut-switch compiler, giving
at most2^(n−1) ordinary classes, a vanishing fraction of A000568. This is a
capacity bound, not a reason a sparse family cannot cover H-values. It also
blocks conflating real positivity with finite-field Paley constructions.

[Exact distributions, hidden-state compiler and phase decoder](wildcard/PROOF.md).

## 5. What connects the four lanes?

Three precise principles describe the common direction without pretending to
solve all four subjects:

1. **Boundary-conditioned growth.** The correct next-state law may depend on a
   terminal horizon, a phase interval, or rooted endpoints, even when the
   current unmarked scalar is identical.
2. **Transport before observation.** A coboundary can move a real resource while
   leaving a whole closed-observer category invariant. A deliberately marked
   seam can expose it. This gives a way to search beyond scalar-preserving
   gadgets without abandoning control.
3. **Resolution-specific completeness.** The class census, its infinite local
   limit, a microscopic component measure and an exact target count are
   different levels of information. Completeness at one does not promote itself
   to completeness at another.

These are named working principles supported here by explicit maps, controls
and counterexamples. They are not new universal mathematical laws. The old
repo sources and lost-coordinate ledgers are recorded in each proof and the
[research board](RESEARCH_BOARD.md).

**Bottom line:** the session improves constructive interfaces, exact recursive
prediction and the understanding of emergent symmetry. It does not prove
Chan–Pak, a new LRC case, planar JC, or Yang–Mills. The most actionable next
Chan–Pak target is the rooted-resource renewal lemma, stated concretely in
[open questions](OPEN_QUESTIONS.md), rather than another extension of numerical
coverage.
