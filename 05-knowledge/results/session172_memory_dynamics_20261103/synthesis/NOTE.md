# Contextual memory growth: what a complete observation still cannot prove

**Status: proved synthesis and explicit open targets.** S172 constructs stronger physical memory and sharper observation interfaces. Its Collatz parity and finite-stopping-density ingredients are classical; it does not prove convergence, a new Hamilton-count spectrum classification, or a new square-sum existence range. The useful synthesis is a separation between finite-horizon information, continuation-compatible memory, and arithmetic control.

## 1. A small theorem organizing the session

Let X be a nonempty source set with allowed maps f_1,...,f_m:X→X and a finite-valued observer o. For k≥0 let x~_k y mean that o(f_w(x))=o(f_w(y)) for every allowed word w of length at most k, including the empty word. Here all words in the declared alphabet are allowed. Let q_k be the number of equivalence classes, assumed finite at the horizons under discussion. Define the **contextual memory requirement**

    m_k=ceil(log_2 q_k).

This is the worst-case fixed-length binary length needed to name a response class. For any source probability law, H(class)≤log_2 q_k≤m_k. Equality throughout requires uniformly distributed classes and q_k a power of two. Thus the memory requirement is not generally a Shannon entropy.

**Proposition (memory growth and closure).**

1. Any deterministic record supporting every horizon-k response needs at least q_k states, and class labels attain this bound. The sequence q_k is nondecreasing.
2. Exact positive stochastic encodings with deterministic output channels also require at least q_k hidden states. Randomizing the encoding does not reduce this cardinality.
3. If ~_k=~_(k+1), then ~_k is preserved by every f_i and equals full future-response equivalence. In particular, q_k=q_(k+1) implies this conclusion when q_k is finite.
4. If q_k is unbounded, no finite autonomous record can preserve every permitted future response. Conversely, if all q_k are finite and bounded by one common integer, the chain stabilizes and supplies a finite exact autonomous quotient.

**Proof.** Distinct response classes must have distinct deterministic records; storing the class suffices. Adding contexts refines classes. For positive encodings, a hidden state in the support of two sources must give both sources' deterministic response under every output channel: extremality of probability-simplex vertices forces those responses to agree. Hence different classes have disjoint nonempty supports, proving (2).

If x~_k y and ~_k=~_(k+1), every word of length at most k following f_i gives equal outputs, so f_i(x)~_k f_i(y). Induction then preserves all future outputs. Finite equal class counts in a refinement imply equal partitions. An unbounded family of finite-horizon class counts cannot fit into a finite autonomous record. A bounded nondecreasing integer sequence eventually has consecutive equal terms, proving the converse. ∎

The statement combines the nonlinear quotient in [observers](../observers/NOTE.md) with the positive-state argument in [carry](../carry/NOTE.md). It is elementary, not a priority claim. Horizon indexing differs from the observers packet: its k-symbol single-map prefix includes times 0,...,k−1, whereas the proposition includes all words through length k.

This criterion requires contexts closed under the specified continuation alphabet. A convenient finite experiment bank, or a scalar that presently separates some states, need not satisfy it. Graph extensions also change the labelled vertex set; their interfaces are naturally graded by order. They must not be silently treated as endomorphisms of one fixed set of n-vertex graphs.

## 2. Collatz: maximal information growth is not a descent theorem

The [b-adic observer theorem](../observers/NOTE.md) covers affine residue branches

    T(x)=(u_a x+v_a)/b,  x≡a mod b,

where u_a is a b-adic unit and u_a a+v_a≡0 mod b. Every length-k residue word occurs on exactly one class modulo b^k. Therefore an exact k-symbol record has b^k states and needs ceil(k log_2 b) binary bits. Its response metric is exactly the b-adic metric.

A dynamical step consumes one digit of available precision: the class modulo b^k determines the next class modulo b^(k−1), but not modulo b^k. Finite snapshots are complete at their declared horizon yet do not form same-horizon autonomous quotients.

The [Collatz tape](../collatz/NOTE.md) realizes every k-parity word using 5k tournament vertices, with the same full path-cover profile for all 2^k states. A specified one-vertex Hamilton context reads a requested bit through the exact 17/19 seed response. This is a literal finite storage theorem, not a physical autonomous implementation of Collatz.

On the cylinder with word b_0,...,b_(k−1), the shortcut map satisfies

    T^k(n)=(3^s n+a_k)/2^k,
    a_(j+1)=3^(b_j)a_j+b_j2^j,  a_0=0.

The pair (s,a_k) distinguishes all words of fixed length. At fixed weight s, all binom(k,s) additive offsets are distinct. A simultaneous unweighted Hamilton read gives only s and leaves log_2 binom(k,s) bits hidden under the uniform conditional law. Thus the additive cocycle is exactly where the omitted order information resides; two integer registers are not two bits.

Arithmetic descent needs an additional inequality:

    T^j(n)<n iff (2^j−3^(s_j))n>a_j.

The packet gives the exact threshold for every finite cylinder, including its small exceptions. Its survivor-count dynamic program reproduces classical density-one finite stopping, not universal convergence. The classes n≡−1 mod2^K have K strictly increasing steps, so no uniform bounded descent horizon exists.

The observer theorem itself cannot resolve this gap. The 1x+1, 3x+1 and 5x+1 shortcut maps have identical full-shift observer geometry, yet 1x+1 decreases trivially and 5x+1 has the positive cycle 13→33→83→208→104→52→26→13. Its parity word corresponds under 3x+1 to the noninteger 2-adic point 19/101. Symbolic conjugacy does not preserve the positive-integer subset or its height.

## 3. Strong physical memory: a genuine advance over S171's construction

The [strong wrapper](../spectrum/NOTE.md) resolves S171's stated strongly connected storage target at the level of strong, rather than 2-vertex-strong, tournaments. Begin with M≥2 ordinal five-vertex seed blocks. Add a vertex a beating the first block and beaten by all remaining blocks. The three groups form a directed cycle.

This single feedback vertex preserves full-profile blindness and labelled Hamming distance, while making every stored graph strong. Its deletion is the unique single-vertex deletion destroying strong connectivity, so the ordered seed sequence remains canonically recoverable. All 2^M states have

    order 5M+1,
    H=75·84^(M−1)−57·9^(M−1).

One-vertex local read heads still distinguish every addressed bit. The polynomial wrapper identity proves constancy of the entire path-cover profile, not merely H. Its arithmetic interface uses both H and the two-slot count F(2), exposing why a Hamilton scalar alone is inadequate for this composition.

This is an actual improvement to the physical construction, not merely another observer analogy. Its limits are also structural: the unique articulation remains a vulnerability, and the new additive-exponential count family is not asserted to contain previously unrealized integers. Count-measurement separation and adjacency-error protection are different noise guarantees.

The strong tape can encode arbitrary finite parity words as well as codewords. A Reed–Muller restriction may protect encoded information, but cannot be imposed directly on all raw Collatz words, since every binary word occurs.

## 4. Two read heads, with different permission to continue

The [carry interface](../carry/NOTE.md) considers marked upper-triangular matrices

    G_b=[[a,b],[0,d]].

Every trace formed using only upper-triangular contexts is independent of b, even with repeated unknown factors and defined inverses. An external probe P reads b exactly when its lower-left entry is nonzero. In particular,

    J=[[0,−1],[1,0]],  tr(JG_b)=b.

This enlarges the experiment alphabet and reads retained information. It does not recover b from the old trace alone. More importantly, J is not supplied as a legal Collatz step.

The exact accelerated-odd words (1,3) and (3,1) have identical diagonal data but carries 5 and 11 and authentic source guards 19 mod32 and 29 mod32. Reading the carry distinguishes the formal branch record; checking its source and intermediate states establishes whether the word is executable. Neither operation supplies descent.

The tournament read head is an actual graph extension; J is an external matrix context; parity continuation is an actual arithmetic iteration. They share a response-fiber theorem, not an interchangeable operation algebra.

## 5. One-step graph memory can be small and still fail at two steps

The [graph attachment theorem](../squares/NOTE.md) characterizes every undirected one-new-vertex Hamilton response. If the new vertex has old neighborhood S, then

    R_G(S)=2Σ_(i∈S)E_i+2Σ_({i,j}⊆S)K_ij,

where E counts Hamilton endpoints and K counts ordered spanning two-path pieces joining at specified ports. The exact linear response rank over all labelled n-vertex graphs is n(n+1)/2 for n≥3. This is a signed-linear dimension, not a count of positive hidden states or a running-time bound.

The packet gives nonisomorphic five-vertex graphs with identical complete path-cover profiles but a fixed leaf readout of 0 versus 4 Hamilton paths. Thus contextual recovery can expose existence, not just change one positive count to another.

It also proves a continuation failure. The edgeless five-vertex graph and a graph of two disjoint edges plus an isolated vertex both have E=K=0. Every one-vertex response agrees. Adding two independent vertices each adjacent to every old vertex gives counts 0 and 48 respectively. The horizon-one record is complete for its declared bank but is not closed under repeated extension.

For square-sum graphs the actual new neighborhood supplies an arithmetic specialization. The inherited 22→23 and 24→25 arrivals read K=3 and K=10, giving 6 and 20 oriented Hamilton paths even when the previous Hamilton counts are zero. These rechecked examples explain the missing boundary memory; they are not new arithmetic realizations of the abstract equal-profile seeds.

## 6. The research agenda now has explicit proof obligations

- **Collatz arithmetic control.** Build the exact cylinder-descent tree on n>1, with 1 as induction base. Supply a well-founded reduction of every unresolved cylinder or a new restriction on the positive-integer parity image. Better reading or smaller exceptional density alone does not do this.
- **Stronger tournament storage.** Remove the canonical articulation while proving a full-profile identity, nonisomorphism recovery and bounded-cost readout. A second feedback vertex is a candidate, not a proved solution.
- **Arithmetic-legal graph interfaces.** Determine response ranks for specified square-sum source and arrival families, then prove a continuation law. The rank for all labelled graphs cannot simply be transferred to an arithmetic subset.
- **Legal carry experiments.** Identify which nontriangular probes, if any, have an authenticated arithmetic implementation; distinguish that from externally measuring a matrix coordinate.
- **Memory-growth measurement.** For a declared continuation family, measure response-class growth or linear rank growth as the horizon increases. Keep positive state cardinality, encoded bits, signed-linear dimension, physical graph size and count precision separate.

The session establishes exact response growth for full-branch dynamics, a stronger physical tournament carrier, a sharp one-step graph interface, and explicit failures of continuation closure. A complete observer can preserve every relevant finite detail while leaving the central arithmetic inequality entirely unproved. That separation is the productive connection across these problems.

## Source and verification boundary

This synthesis reads and uses the five sibling source notes linked above. Historical references and exact finite-check certificates remain with their owning packets. All numerical values here are proved identities or those packets' verified finite examples. The proposition in Section 1 is proved directly here; no new empirical extrapolation is used.
