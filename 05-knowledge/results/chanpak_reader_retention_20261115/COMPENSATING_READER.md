# A compensating two-closer window

**Proved:** two consecutive odd counts at each height of a common clock, from changing the host and its genuine cover response. First pair: **29,311 and 29,313**. This is a width-one window, not overlapping placement or Chan–Pak.

A deterministic sample of old notes found the [Fibonacci-box matching identities](https://github.com/eliottcassidy2000/math/blob/main/07-reflections/fibonacci-boxes-farey-and-the-k4-matching-triple-opus-20260728.md) and the [four-point turnpike note](https://github.com/eliottcassidy2000/math/blob/main/07-reflections/four-point-turnpike-xor-gram-branch-language-and-harmonic-sidecars-codex-20260815.md). The former exhibits a readout retained by balanced sums and differences rather than equality of its primitive products; the latter separates a projected invariant from missing sidecars. These suggested matching a *selected linear combination* of the two physical responses, not matching both independently.

## Exact colored-bin calculation (inherited transfer; written in scoring coordinates)

A marked host is a tournament B with marked set O. Its closing C(B,O;A) puts B before module A, then a vertex v; the arcs v→B go exactly to O, and the other inter-module arcs are forward. For a subset S of the host let h(S) count Hamiltonian paths (h(empty)=1), E(S) count those ending in O, and define M(S) as the number of ordered two-path covers of S whose second path starts in O (first may be empty). Let W(S) count such covers with first path nonempty and ending in O, and r(S)=h(S)−W(S). Then for nonempty A,

    H(C(B,O;A)) = r H(A)+M Q(A),

where r,M are evaluated on all B and Q counts ordered two-path covers with empty slots permitted. This is the inherited THM-5130 closing law. For **empty** A the correct expression is r+M−E, because the join at v no longer allows a host path ending at a marked vertex immediately before v. This exceptional term matters in the next color.

Write C_3 for ordered three-path covers with empty slots allowed, and define the folded host statistics

    R = sum_{S⊆B} r(S) h(B−S),
    U = sum_{S⊆B} M(S) h(B−S),
    F = sum_{S⊆B} E(S) h(B−S).

Now choose one of two path bins as the bin containing v. Write S for its host vertices and T for its module vertices. The second bin, without v, has h(B−S)h(A−T) paths. The v-bin uses the ordinary closing expression except at T empty, where subtract E(S). Summing all choices gives

    Q(C(B,O;A)) = 2R Q(A)+2U C_3(A)−2F H(A).

This is the two-bin specialization of the fuller marked-bin transfer of [Session193](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session193_triangulated_proof_20261105/profile_substitution/PROOF.md), in converse orientation. The correction term must be retained, even for a nonempty physical A, because its contribution to the bin containing v can be empty.

## Compensation rather than interface equality

Use the outer closer with host TT2 and both vertices marked. Its known state is r=−1,M=3, hence it reads the inner tournament as 3Q−H. It recognizes the host only through

    H(outer(inner(A))) = 6U C_3(A)+(6R−M) Q(A)−(r+6F) H(A).

Two hosts can preserve U and 6R−M while neither R nor M is preserved. The exact finite host certificates are:

| order | orientation code | mark mask | M | R | U | r | F |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 8388864 | 65 | 81 | 96 | 1674 | −27 | 306 |
| 7 | 1025 | 73 | 177 | 112 | 1674 | −25 | 306 |

Code bits correspond in ascending lexicographic order to pairs i<j, least-significant bit first; bit 1 says i→j. Mask bit i marks vertex i. The host with order 8 has corrected offset 1809 and the order 7 host has 1811. Both have 6R−M=495. Thus for every t≥1, set A=TT_t to get the consecutive odd pair

    10044·3^t+495·2^t−1811,
    10044·3^t+495·2^t−1809.

The smaller count uses the smaller host. At t=1,2,3 the pairs are (29311,29313), (90565,90567), and (273337,273339). In the **inner** tournaments, h8−h7=−96·2^t−2 and Q8−Q7=−32·2^t; their changes cancel except for the 2 in the actual outer reader. In particular the mechanism is not equal Q and cannot be justified using only a multiplicative iso-profile. The motion remains exactly 2H(A) for a general nonempty A; it is exactly 2 only for a unit-count module. Consecutive t-windows are disjoint by direct subtraction and grow geometrically apart.

For comparison, the order-7 pairs (1382658,38),(67680,76) have genuine same-(M,R,U,F)=(397,464,2705,433) and r=−17,−11. These have same Q and h-gap 6 for a transitive module, not 2. This is consistent with the modulo 3 cover law Q≡(−1)^n H: a same-order equal-Q pair of odd counts can only have gap divisible by 6. Compensation removes the requirement to match Q; the explicit witness also changes order parity.

**Control:** using instead the outer 6Q−H, codes (71836161,45),(143663362,45) look like an adjacent pair if F is omitted. They match U=7750 and 12R−M=−19423, with r=−313,−311, but F=1310,1322. The true constant spacing is 146. This is a literal counterexample to the tempting empty-bin shortcut.

## Transitive marked hosts: exact 2-adic choice, sparse absolute image

Take a fully transitive host B=TT_n with ordered vertices 1..n and nonempty O, clock A=TT_t, t≥1, and fixed outer anti-reader αQ−H from any known anti-seed of fixed slope α≥1. The linear folded parameter is

    U(O) = sum_{i∈O} 2^{i−1}3^{n−i}.

For a specified marked start i in the first v-bin's two-path cover, the place i has to be in the marked path, earlier vertices have two other bins available, later ones have any of the three. Thus each bit i−1 has a single coefficient of that 2-adic valuation. Binary peeling proves that O↦U gives a bijection mod 2^n, although not a real interval. The empty mask supplies the residue zero in that statement.

There are only Θ(X^{log_3 2}) **presentations with nonempty O** in this family up to X. This is also an upper bound on distinct counts; no distinct-count lower bound is claimed. Indeed W(S)≤M(S), since W counts the subset of M with the additional marked end, so r(S)≥−M(S). For a nonempty clock subset T counted in the bin containing v, Q(T)≥2h(T), giving ordinary closing count at least M(S)Q(T)/2. Empty-T contributions to the full counting argument are nonnegative. Consequently the inner Q is at least U(C_3(A)−H(A))=U(3^t−1). Also outer H≥(α−1/2)Q and hence at least a fixed positive constant times U3^t.

If the first mark is i, let l=i−1 and d=n−i. There are 2^d choices for its suffix mask and U≥2^l3^d. For counts ≤X this gives 2^l3^(d+t)≲X. Write J=floor(log_3(CX)) with a fixed suitable C. For a given d there are O((J−d+1)^2) possible pairs (l,t) under this inequality, so summing 2^d over d≤J with this polynomial factor gives O(2^J). Conversely, choose first mark at i=1 and t=1. We have R≤2^n, U≤3^n, F≥0, so all 2^{n−1} mask presentations at these parameters have inner Q≤4·2^n+6·3^n; their outer count is at most α times this. Choose n=floor(log_3(X/(10α))) to get the matching presentation lower bound. Empty O contributes only O(log X) **distinct counts** (its separate presentation count is polynomial in log X). This is a clock-and-transitive-reservoir limitation, not a bound on the general two-hub spectrum or on nontransitive compensation hosts.

## Replay

Run `python verify_compensating.py`. It independently counts end/start-marked induced subpath partitions, reconstructs literal inner and outer tournaments and checks three heights. Optional exact host search: `reader_census.cpp` takes class-list rows `order code` in this orientation convention and an optional outer α argument; it folds endpoints from a subset DP and prints collisions in the exact key `(U,2αR−M)`, with corrected offset `r+2αF`. The mathematical statements need only the displayed literal witnesses, not census completeness.
