# Realizing an unsigned curl field: positive parity phases and negative NAE constraints

Status: PROVED elementary exact encoding; FINITE-EXACT exhaustive input census through order5. The underlying Boolean encoding is standard constraint compilation, not a claim of a new general SAT method. This is a separate bounded third-round addendum; the previous frozen proofs and receipts were not changed.

## 1. Inheritance and target

The recovered C3-structure theorem of Boussaïri–Ille–Lopez–Thomassé identifies realization fibres: two labelled tournaments have the same cyclic triples exactly when they are related by homogeneous-interval reversals; for an indecomposable realization the fibre is just the tournament and its converse. That is a CITED classification, not a recognition algorithm supplied by this note. Primary source: https://perso.ens-lyon.fr/stephan.thomasse/liste/bilt.pdf , Theorem2 and Corollary1. Its Hamiltonian-response consequence was proved in `CURL-AND-CLAN-COMPLETION.md`.

A targeted repository search found the parity union-find method already used for a different object in `math/05-knowledge/results/tiling_line_metagraph_self_lines_twins_depth_20261020.md`: permutation-edge constraints for self-lines. Signed NAE appears elsewhere, including `THM-862-scale-three-hamming-six-decision.md`, in an unrelated arithmetic provider test. Neither result supplies the present C3-realization theorem. No repository or literature priority claim is inferred from that search.

Here the input is an arbitrary labelled3-uniform hypergraph `H` on `[n]`: precisely its triples are demanded to be cyclic. The question is whether it is the unsigned triangle-curl field of ANY tournament, and, if so, to recover its complete labelled realization fibre.

## 2. The exact Boolean system

For every unordered pair `i<j`, introduce

    x_ij=1 iff i->j.

For an increasing triple `i<j<k`, the arc bits in the cyclic direction `i->j->k->i` are

    (x_ij, x_jk, 1-x_ik).                                (1)

The triple is cyclic if and only if the three bits in (1) are all equal. Consequently:

- for each positive triple in `H`, impose

      x_ij XOR x_jk = 0,
      x_ij XOR x_ik = 1;                                (2)

- for each negative triple outside `H`, impose

      NAE(x_ij, x_jk, 1-x_ik),                           (3)

where `NAE` means not all three arguments are equal.

The sign on the long edge `ik` is load-bearing. Replacing (1) by three unordered edge bits all with the same sign gives a false system.

### Parity-DSU reduction

Apply parity union-find to the equations (2). If they are inconsistent, there is no realization. Otherwise, each of the `p` connected phase components has one free bit `z_r`, and every edge has a recorded expression

    x_e=z_(r(e)) XOR epsilon_e.                          (4)

Substitute (4) into all clauses (3). This yields a signed NAE formula on exactly `p` phase bits, allowing repeated variables within a clause.

**Theorem1 (exact realization criterion).** The input is tournament-realizable if and only if its positive parity system is consistent and its reduced signed NAE formula is satisfiable. More strongly, satisfying phase assignments are in bijection with all labelled tournament realizations, via (4).

Proof. Every edge assignment uniquely specifies a tournament. The positive equalities are exactly the cyclic conditions on positive triples. A consistent parity component has one independently selectable root bit, so (4) is a bijection between its `2^p` assignments and all edge assignments satisfying those positive conditions. Under that bijection, (3) is exactly the requirement that every remaining triple be transitive. This proves soundness, completeness and uniqueness of decoding. QED.

Thus positive information is a genuine affine phase reduction. Negative information is not merely a check of parity consistency: the residual clauses can couple the phases globally.

## 3. Gauge, repeated literals and complement pairing

Changing the chosen root orientation of a DSU component replaces one phase variable by its complement and toggles the recorded signs. It changes the coordinates, not the decoded tournaments or satisfiability.

Repeated-literal clauses can be simplified exactly:

- if a variable and its complement both occur, the clause is a tautology;
- three identical literals give immediate contradiction;
- `NAE(a,a,b)` is exactly the binary XOR condition `a != b`;
- otherwise the clause genuinely involves three distinct phase variables.

Binary conditions may be propagated into a second parity reduction, but this is optional and does not eliminate the general residual obstruction below.

For `n>=2`, simultaneously complementing all phase bits complements every edge bit, hence takes a tournament to its converse. Every clause is invariant under this move and no edge assignment is fixed. Therefore every nonempty realization fibre has even size and splits into exact converse pairs. One chosen phase may be fixed to0 to enumerate one representative per pair. At `n=1`, there is one realization and no edge variable, so this exceptional boundary must not be folded by two.

This is a gauge reduction for search or representation, not a theorem that one converse representative is canonically the desired tournament for another labelled consumer.

## 4. Two distinct impossibility mechanisms

### 4.1 Positive parity itself may be inconsistent

On four vertices, demand that `012,013,023` all be cyclic. From the first two triples,

    x_02=1-x_01,   x_03=1-x_01.

The third demands `x_02=1-x_03=x_01`, contradiction. Negative triples play no role in this failure.

### 4.2 Consistent positive phases can leave an unsatisfiable NAE system

On five vertices, demand exactly the two cyclic triples

    H={014,023}.                                         (5)

They share one vertex but no edge. Their four positive equations are independent, so the ten edge variables reduce to `p=6` freely assignable phase bits. Every one of the eight negative triples becomes a genuine ternary clause on three distinct phase variables. Thus neither a positive parity check nor repeated-literal simplification rejects the input.

Nevertheless no realization exists. Here is an elementary certificate.

By global converse and, if necessary, exchanging labels2 and3 (an automorphism of the prescribed hypergraph), normalize the two demanded cycles to

    0->1->4->0,    0->2->3->0.

Avoiding a cycle on `013` forces `3->1`. Avoiding one on `024` forces `4->2`. Now `1->4->2`, so avoiding a cycle on `124` forces `1->2`. But then

    1->2->3->1

is a forbidden cycle on `123`, contradiction.

**Finite-exact minimality in order.** Every cyclic-triple hypergraph on at most four vertices whose positive XOR system is consistent is realizable. Thus (5) is a smallest-order counterexample to “positive parity consistency suffices.” It also uses only two positive triples; zero or one positive triple is always realizable: use a transitive tournament or one cyclic triangle as an ordinal block among singleton blocks.

This is the first failed implication, not evidence of arbitrary3-SAT hardness. The derived clause families have special tournament incidence structure.

## 5. Indecomposable rigidity still needs negative information

If the input has an indecomposable tournament realization and `n>=2`, the cited C3-structure theorem forces its complete solution set to contain exactly two assignments, paired by converse. But the positive DSU need not already have `p=1`.

A concrete five-vertex input is

    H={012,014,023}.                                     (6)

The positive equations leave four independent phase bits `a,b,c,d`. In ascending edge order `01,02,03,04,12,13,14,23,24,34`, write

    (a, 1-a, a, 1-a, a, b, a, 1-a, c, d).                (7)

Of the seven negative clauses, four are tautologies. The remaining clauses successively force

    d=1-a,    b=a,    c=1-a.                              (8)

Hence there are exactly two realizations. Their ascending-edge bit masks, with bit1 meaning the ascending arc, are117 and906, complements in ten bits. The checker verifies directly that they have no nontrivial homogeneous interval.

Thus an indecomposable response class can have four positive phase components but only one global orientation bit after the negative constraints. The saved classification is not “one component implies prime”; it is that the full phase system has one converse pair whenever an indecomposable realization exists.

## 6. Exhaustive census and honest reduction size

The standalone checker enumerates every candidate hypergraph, not only those already known to arise from tournaments. Its results are:

| order n | candidate hypergraphs | parity-inconsistent | parity-consistent but NAE-unsatisfiable | realizable |
|---:|---:|---:|---:|---:|
|1|1|0|0|1|
|2|1|0|0|1|
|3|2|0|0|2|
|4|16|5|0|11|
|5|1024|771|45|208|

For the208 realizable order-five hypergraphs, the positive-DSU phase count histogram is

    p=1:12, p=2:70, p=3:15, p=4:70,
    p=6:30, p=8:10, p=10:1.

The original edge space has ten bits. Its median phase count over these208 input classes is four, but the empty cyclic-triple hypergraph retains all ten bits. The latter's NAE system cuts the1024 edge assignments down to the120 transitive tournaments, one for each linear order. Thus there is no uniform nontrivial phase reduction.

Among order-five input classes having indecomposable realizations, the histogram is

    p=1:12, p=2:60, p=4:60,

and every class has exactly two actual realizations. These are class-weighted finite counts, not random-tournament statistics or an asymptotic complexity claim.

The reduced system has `binom(n,3)-|H|` raw NAE clauses. Positive processing uses `2|H|` parity equations on `binom(n,2)` edge variables; standard parity union-find implements this stage in near-linear time in those sizes. Exhausting the residual phase assignments costs up to `2^p` times the clause count. No polynomial-time recognition theorem, general speedup, or NP-hardness result follows from this presentation.

## 7. What this changes and what remains outside scope

- An unsigned curl/complete Hamiltonian-response state now has an explicit exact realization certificate: consistent phase identifications plus a satisfying residual orientation assignment.
- It separates impossible local phase cycles from globally incompatible transitive requirements, with a smallest hostile exhibiting the latter.
- Fixing global converse is safe for enumerating response-state realizations. It is not safe if the target consumer observes signed curl, marked paths, port orientations or a chosen metagraph wiring.
- The theorem recognizes and enumerates labelled realizations of a supplied cyclic-triple hypergraph. It does not reconstruct the tournament metagraph or solve the root's path-transport constraints.
- A future practical target is structured SAT or modular-decomposition propagation on these special signed NAE systems. No such algorithmic improvement is claimed here.

## 8. Reproduction

Run

    python3 session36-context/verify_realization_phases.py

The verifier uses no repository imports. It generates all tournaments directly, computes their cyclic-triple hypergraphs without DSU, and compares each resulting exact fibre against the independent parity-DSU/NAE solver. It also checks complement pairing, every reported histogram, the empty-hypergraph `n!` boundary, all homogeneous subsets for indecomposability, and the two-hostile distinction.

Outputs: `realization_phase_receipt.out`, `realization_phase_receipt.json`. This addendum and its companion are the only new third-round files; no previous frozen artifact was rewritten.
