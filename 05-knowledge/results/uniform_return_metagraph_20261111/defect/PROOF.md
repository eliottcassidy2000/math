# Fixed defect: a finite grammar for every branching stratum

**Proof candidate under independent review; small automata FINITE-EXACT.**
No all-height source-return or Chan–Pak theorem is claimed.
The new scope is all fixed positive ready defects with at least two inports,
including nontransitive residuals; the older atlas treats two inports with
transitive residual and neutral start/delete response.

## 1. Inheritance and the positive quantity

For an actual strong rooted tournament (T,r), put h=H(T), q=H(T-r),
s=s_T(r), F=E_r(2), and b=F-h. Let I=N^-(r) and d=|I|.
For U=T-r, W_ac(U) counts two nonempty ordered spanning paths, the first
ending at a and the second starting at c. The inherited identities in
[ADAPTIVE section3](../../adaptive_metagraph_reformulation_20261103/ADAPTIVE.md#3-the-defect-detects-indegree-and-the-legal-inverse)
and [READY section1](../../two_entry_funnels_20261103/READY.md) give

    b = sum_(a in I) s_U(a) + sum_(a!=c in I) W_ac(U).             (1)

The first sum is q-s: prefixing r accounts for exactly those deletion
paths whose first vertex lies outside I. In particular (1) is a sum of
counts of specified objects. A capped algorithm must compute this sum
directly, not subtract two separately clipped large path counts.

The inherited pair-capacity bound is

    W_ac(U)+W_ca(U) >= 2 H(U-{a,c}).                             (2)

It is a bound on the PAIR, not each individual W term. An individual
ordered term may be zero. This distinction is essential for weak deletions.

## 2. A bounded residual kernel for every branching defect

Suppose b is fixed and d>=2. Write R=T-({r} union I). Induced insertion
monotonicity gives H(U-{a,c})>=H(R), including H(empty)=1. Summing(2),

    d(d-1) H(R) <= b,       d(d-1)<=b.                          (3)

Thus d is bounded and H(R)<=floor(b/[d(d-1)]). In the unique ordered SCC
decomposition of R, every nontrivial strong component has at least three
Hamilton paths; their product is H(R). Each such component has order at
most its path count, and there are at most floor(log_3 H(R)) components.
Consequently only finitely many ordered nontrivial strong kernel lists
can occur. An elementary bound on their total vertices is

    K floor(log_3 K),     K=floor(b/[d(d-1)]),                  (4)

with zero nontrivial vertices when K<3. Sharper inherited Busch bounds
can replace this deliberately crude bound. Arbitrarily many singleton
components remain in each gap of the ordered kernel list.

This is the same finite-kernel/variable-gap decomposition used in the
[capped metagraph lane](../metagraph/README.md); the important new map is
from a bounded ROOT DEFECT to a bounded PATH COUNT after deleting r and
all its inports. The original h itself can grow without bound.

## 3. Effective regular languages, with exact isomorphism quotient

Fix d, an ordered residual kernel list, literal representatives of those
cores, all their incidences with the d named inports, and the tournament
on the inports. Treat those finitely many vertices as feedback vertices.
The remaining vertices are the ordered singleton gaps. They form a
transitive backbone, with a finite letter recording their incidences
with the inports and their gap position. Incidences with the residual
cores are forced by SCC order. The root beats the backbone and all
residual cores and loses to every inport.

For each path object in(1), list the feedback vertices in their path
order and, for two paths, choose the split between their lists. Every
backbone vertex belongs to exactly one increasing run between consecutive
feedback vertices or at a path boundary. The inherited
[finite-feedback circuit](../../session176_tournament_information_20261104/spectrum/PROOF.md#4-a-growing-nonuniform-interface-circuit-that-counts-all-orphans)
extends directly to these specified path boundaries: each run is empty,
nonempty with its last vertex pointing to its right boundary, or nonempty
with the wrong last incidence. A first run vertex must satisfy its left
boundary. A final empty internal run requires the feedback-feedback arc.
Forbidden endpoint gaps are kept empty. This is a finite nonnegative
integer weighted automaton counting precisely (1), with no orphan loss.

Cap every register at b+1. Truncation commutes with nonnegative sums and
products, so the exact equality test b is retained. The capped vector
therefore gives a finite deterministic automaton. Strongness is regular
as well. A tournament fails strongness exactly when a nonempty proper
vertex set dominates its complement. In this encoding, its backbone
vertices must form an initial segment. A finite automaton can guess its
subset of the fixed feedback vertices (including r), guess the backbone
cut position, and check every required cross-incidence before and after
that position. It also tracks that both sides are nonempty. Complement
this regular failure language to recognize strongness. This is an effective finite
construction; no uniform polynomial state bound in b is asserted.

The residual SCC list is intrinsic after removing the intrinsic set I.
Any rooted isomorphism preserves its order and every singleton position.
For each fixed list the finite group

    S_d times product_i Aut(C_i)

acts on the inport template AND contact letters. There is no permutation
of equal core occurrences at different SCC positions. This gives exactly
the rooted-isomorphism quotient of the representations, not just a
presentation count. For each group element, fixed templates and fixed
letters impose regular restrictions. Burnside's finite average therefore
shows that the ordinary rooted-class generating function, by vertex
order, is rational for every fixed b and d>=2. The finite union over d
remains rational. This is a classification algorithm for each fixed
branching defect, not a fixed finite bank as b increases.

The d=1 case has no bound(3). It retains the inherited endpoint-budget
funnel grammar and is not included in this rationality conclusion.
Nor does the theorem say that a ternary output returns to a low-defect
stratum: the inherited H35 obstruction already refutes that general step.

## 4. Complete first branching strata, including the missing nonneutral roots

For b<=5, (3) forces d=2 and H(R)<=2, hence R is transitive. Order the
two inports a->c. Encode a backbone vertex by the pair of bits saying
whether it beats a,c. The code in `two_inport.py` uses integers
0=00,1=10,2=01,3=11 (low bit is the first coordinate).
The source is strong exactly when the word is nonempty, its last letter
is nonzero, and at least one letter has first bit1. These conditions
ensure that the root reaches a, then c, and that every backbone vertex
can reach an inport; every inport returns directly to the root.

Four disjoint two-run skeletons count(1): the two port orders in a path
starting at an inport, and the two ordered path-bin endpoint choices.
Each has nine raw empty/good/bad run states. Exhaust all reachable
vectors capped at b+1, then retain exactly states that can reach an
accepting strong word. Removing the listed self-loops leaves a DAG.
Every remaining path is retained as an exact regular expression in the
receipt; these disjoint paths certify ALL word lengths.

The resulting ordinary rooted-class series are

    b=2: (x^4+2x^5+2x^6)/(1-x),
    b=3: 0,
    b=4: x^4+x^5+(6x^6+10x^7)/(1-x),
    b=5: 2x^5.                                                (5)

The b2 series is the inherited ready pentagon. The infinite part at b4
is exactly the inherited sixteen neutral tails. The two extra b4 words
are `10` and `01,10`, of orders4 and5. At b5 the complete list is
`10,01` and `01,11`, both of order5. These are nonneutral branching roots;
they were outside the old q=s atlas, not counterexamples to it.
Rooted isomorphisms fix the ordered two inports and the unique residual
chain, so different words are different rooted classes and all are rigid.

At b6 the transitive-residual d2 chart has four additional isolated words
and the inherited31 neutral tails. It is NOT the complete branching b6
stratum: d3 and a C3 residual first become possible at that budget.

## Reproduction and scope

Run `python -B defect/two_inport.py` from the packet directory.
All341 literal words of lengths0..4 (including weak hostiles) are checked
by independent subset endpoint counting against the positive automaton.
The full capped state spaces at budgets2,3,4,5,6 have respectively
736,1217,1972,2765,3681 reachable states. These are not claimed minimal.
The computed finite automata certify their infinite languages; the
general rationality theorem rests on the construction in section3.
