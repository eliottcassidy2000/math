# Six edges, positive witnesses, and integer descent

**Session 104, 2026-10-29.** This packet develops exact connections requested after the all-order tournament-word classification. Proofs and independently coded audits distinguish a real shared mechanism from numerical coincidences. It does **not** prove Collatz, LRC(14), or abc.

## 1. Six edges are a bivector space, but the six-state machine is not its wedge action

K4 has six edge coordinates, naturally the basis of Λ²Q⁴. Its incidence boundary has a three-dimensional kernel (cycles) and a three-dimensional image (zero-sum vertex charges). This is a genuine 6→3 structure.

The actual tournament language has a different hierarchy: three linear counting dimensions, four acceptance states, and six states when legality must also be reported. Its six-state `ee` transition has rank **two**. Exterior squares of four-dimensional linear maps have ranks only **0,1,3,6**. Thus no change of basis can identify this legality machine with an exterior-square representation of four-dimensional letter operators. The coincidence 3/4/6 is not a transition-preserving isomorphism.

There is nevertheless a genuine tournament wedge carrier: at four vertices the native two-run observer has twelve ordered-pair ports. Its antisymmetric run-swap sector is exactly Λ²Q⁴, with six dimensions, and vertex re-markings act by exterior squares. This is distinct from the six legality states. The rational symmetric/antisymmetric split also retains a nontrivial integral gluing issue. As vertices are inserted the wedge space grows, and preserving a kernel vector requires the separate guarded transport proof.

There is an honest alternative: lift the three-state counting matrix C to 1⊕C and take its exterior square. This gives C⊕Λ²C, an explicitly chosen six-dimensional linear model. It is not the legality automaton. Moreover an explicit determinant-one integer change of basis identifies Λ²C with C−I. The dual cubic is the real seventh-cyclotomic polynomial: on its roots, the transformation x→−1/(x+1) equals x→x²−2. This is angle doubling on seventh-root traces, returning after three doublings because 2³≡1 modulo7. It is an exact squaring/doubling bridge, with the counting growth connection to seven inherited from the earlier recurrence work.

The [exterior note](exterior/K4_EXTERIOR_BOUNDARY.md) states the exact boundary, distinguishes cross products from wedges and Hodge stars, and provides reproducible matrix checks.

## 2. The useful common mechanism is positive evidence that cannot cancel

A signed rank-three projection loses half of K4's edge information. Yet it can still detect every nonzero **nonnegative** edge flow, provided the orientation is acyclic. A directed triangle is the hostile control: its three positive edge charges have zero incidence boundary.

More generally, for a finite matrix L the following are equivalent:

- no nonzero nonnegative vector lies in ker L;
- zero is outside the convex hull of its columns;
- a linear functional y makes every entry of yL strictly positive.

This gives a scalar positivity certificate, not a faithful representation for arbitrary later operations. For an oriented graph it becomes the classical directed-cycle/topological-potential alternative. It explains the successful Session103 construction: the auxiliary DAG supplies strictly positive excess for every off-cycle path edge, so no cancellation can hide it. Equality then identifies the canonical cycle-cut cover and produces unit-triangular observers.

The [positive-interface note](positive/POSITIVE_INCIDENCE_AND_ATOMS.md) also implements owner-tagged endpoint insertion for LRC. Adding a runner filters old accepted tags **and births new endpoint tags**. Forgetting the birth term is a false descent argument: accepted tag count can increase from 2 to 6 even though the safe set shrinks. For speeds 1,...,13 at threshold 1/14, the inherited endpoint certificate has six safe locations and twelve accepted owner tags among 182 candidates, while ordinary measure is zero.

### Concurrent work: common ranks survive only compatible port operations

The [concurrent fixed-port packet](https://github.com/eliottcassidy2000/math/blob/ba6481180a0c988e8ed8a11e08dd42155eb0dbf9/05-knowledge/results/wedge_descent_168_20261029/THEORY.md) derives a concrete recurrence from Session103. For the canonically port-matched ordered-k observer B=I+N on n vertices,

    N^(k(n-k)+1)=0,
    B^d=Σ_(j=0)^(k(n-k)) binom(d,j) N^j.

Thus repeated multiplication of this specified table has bounded polynomial depth. This is not an assertion that arbitrary tournament substitution equals table multiplication.

For a finite menu of nonnegative matrices with diagonal one on the same labelled carrier, every product is unipotent exactly when their off-diagonal support union is acyclic. This is the same positive-cycle/common-potential alternative as above. It gives either a common integer grading or an explicit short failing word.

Our [independent narrow audit](audit_incoming/INCOMING_COMMON_RANK_AUDIT.md) verifies these two incoming proofs and reconstructs the native K4 control. Re-marking vertices1 and3 makes a product of two individually unipotent twelve-port tables have trace22. Even trace12 is insufficient: seven re-markings with that trace still have cyclic support unions, and the checker constructs failing longer words. Re-marking is therefore a real loss of common rank, analogous to an unpaid reset, not an innocuous change of names at an unrecorded seam. Other claims in the concurrent packet are outside this audit's scope.

## 3. A literal decreasing sequence of natural numbers

The elementary termination principle is sound: a natural-number rank that strictly decreases on every nonterminal transition proves termination. On an infinite arithmetic graph, acyclicity alone is weaker; n→n+1 has no cycle but runs forever. A decreasing real potential can also run forever if its decrements approach zero.

For the repository's auxiliary (3,3) return controller, the inherited potential is

    P(n,Q)=(n+1)(3/2)^Q,   rho=27/28.

Here n is the actual positive odd state, Q is an active queue, and the native step uses the exact valuation of 3n+1. Branches increase Q by two and contract P by at least rho; leaves decrease Q by one and do not increase P. The following exact integer rank specializes the earlier October11 energy/credit mechanism:

    B(P)=min{k≥0 : rho^k P<3},
    R(n,Q)=3B(P(n,Q))+Q,   with terminal rank R(n,0)=0.

Every active step decreases R. The actual run

    (13,1) → (5,3) → (1,2) → (1,1) → (1,0)

has ranks **163 → 162 → 38 → 4 → 0**. For a queue starting at one, the proof also gives an explicit bound of **60 times the binary length of n, minus56**, on its native return time. This proves each auxiliary return finishes. It does not prove Collatz: a completed return can increase the original number, as **3→5** does, and restarting the queue incurs a reset cost.

### A precise compiler for paid resets

Suppose a block retains immutable source n, has local rank 0≤r≤B(n), and lowers r internally. Suppose every reset is authenticated to have smaller new source n'<n and refreshed rank r'≤B(n'). Then

    G(n,r)=Σ_(m<n)(B(m)+1)+r

is one ordinary natural-number rank for the entire process. It decreases both internally and across resets. For this controller the budget envelope even permits the closed formula **60n·bit_length(n)+r** instead of the prefix sum. This converts local termination plus genuine source descent into global termination without hiding a lexicographic or clock reset. What remains open for the Collatz application is exactly the source-payment obligation, not the impossibility of infinite integer descent.

The [descent note](descent/POSITIVE_ATOMS_AND_INTEGER_DESCENT.md) retains exact carries and guards. Its source debt

    L_J=(2^D−b^J)n0−C = 2^D(n0−n_J)

detects whether a block really ends below its original source. A slope below one alone is insufficient: the 5x+1 word (1,1,5) has slope125/128 and a positive fixed point13. Its source debt is zero, and its native cycle defeats any strict well-founded rank.

## 4. Positive atoms can turn quantitative averages into a rank, but not for free

With fixed positive weights w(x) on every integer source, survival mass below w(x) certifies that source has already terminated. More quantitatively, an independently proved uniform bound

    Σ_x w(x) min(tau(x),K) ≤ M  for all K

implies tau(x)≤floor(M/w(x)). With effective rational data this gives a bounded simulation that computes the exact hitting-time rank. Positive atoms make exceptions visible; they do not establish the required tail or budget estimate. Universal termination need not imply finite mean for a chosen prior. These distinctions preserve the actual target quantifier rather than replacing it by a generic or almost-everywhere statement.

This is a sharpened interface to inherited atomic coverage and October11 rank work, not a claim that the required Collatz estimate has been obtained.

## 5. The old three pieces are recovered, and they explain the even starts

The [recovered September22 discussion](summands/THREE_CHAINS_RECOVERED.md) separates two objects that had been conflated:

- Strict distinct-parent closure from {2,3} misses exactly **{1,4,6}**. Those are holes, not connected components.
- The square-sum graph Q_12 really has three path components, founded at **{1,2,4}**. Vertices 13 and 14 merge them; connectivity persists thereafter.

The complete formulas behind the displayed balances are

    D1(a;k)=a−k²,
    D2(a;k)=(a+k)(a−k(2k+1)),

where Dp is the sum of k+1 consecutive p-th powers minus the following k. Thus even-number balances begin at **2k²**, and consecutive-square balances begin at base **k(2k+1)**. The next square example is

    21²+22²+23²+24² = 25²+26²+27².

A useful exact connection emerges: the even starts 2,8,18,32,... are precisely square-sum graph vertices whose self-pair would sum to a square. The graph excludes these diagonals. This explains the common formula, not recurring isolation: 32 already has two earlier neighbours, 4 and17.

Doubles and squares are also the swap-fixed diagonals of addition and multiplication. The [summand note](summands/SUMMAND_SQUARES.md) retains distinct-parent versus shadow-graph semantics, gives all-k proofs and finite-difference identities, and explains why Gaussian squaring moves prime information into other coordinates rather than preserving all arithmetic information.

## 6. K4 contains the whole abc problem, not an automatically easier subfamily

Four integer potentials give edge differences d_ij=x_j−x_i. Their decomposable bivector satisfies

    d01 d23 + d03 d12 = d02 d13.

This is a real addition/multiplication connection, using the three opposite-edge matchings of K4. In fact it carries **every primitive positive abc triple** after gcd reduction. Given coprime A<B with C=A+B, put d=B−A and choose potentials

    (0,d,2d,2B).

The three products are (2dA,2dB,2dC), with gcd2d. The remaining primitive equal case (1,1,2) comes from potentials(0,2,3,6), with products(6,6,12).

Consequently a universal abc bound proved for this unrestricted reduced K4 family would already prove abc in full. The geometry is a faithful encoding of the additive relation; it supplies no radical or height improvement by itself. This is a constructive scope test, stronger than a vague warning that similar numbers need not be related.

## 7. A unifying research contract

The productive common theory is not that all these objects are isomorphic. It is a set of independently testable obligations:

1. **Represent the relation exactly.** Exterior products, summand hyperedges, and graph shadows retain different data.
2. **Retain what future operations can observe.** A scalar count, an acceptance state, and a legality state have different minimal costs.
3. **Give actual targets positive visibility.** A signed projection or diffuse average can erase precisely the exception being studied.
4. **Pay every transition and reset.** A natural-number rank must follow native clocks and immutable sources, including source payment at block boundaries.

The exterior rank obstruction, positive triangle circulation, isolated LRC endpoints, unpaid3→5 return, 5x+1 cycle, and nonleaf32 are hostile controls for different failed implications. The surviving mechanisms are proved in their stated scopes, rather than transferred by analogy alone.

## Verification and provenance

The [algebra audit](audit_algebra/REPORT.md) and [descent/positivity audit](audit_descent/REPORT.md) use independently coded checks without producer imports. Producers and audits run normally and under Python optimization with explicit guards. Source manifests record inspected local documents and live GitHub blobs; the three-chain recovery distinguishes a preserved condensed framing from a raw user transcript.

Classical ingredients and prior repository results are explicitly credited. No scarce theorem identifier has been reserved. Exact computational checks corroborate the proofs; they do not replace all-order arguments or close the named open conjectures. Run `python3 verify_all.py` from this directory (or invoke it by absolute path). The aggregate replay script and manifest document the final verification scope.
