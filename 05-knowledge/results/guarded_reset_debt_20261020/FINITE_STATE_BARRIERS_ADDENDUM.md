# Finite-state debt barriers: exact cycle criterion and a coverage obstruction

Date: 2026-10-20. Status: **PROVED elementary finite-graph results and precisely scoped Collatz obstructions**. This is a separate extension; RESET_DEBT_COMPOSITION.md and its original verification artifacts remain unchanged. No global ROOT or Collatz theorem is claimed.

## 1. Exact criterion for positive state barriers

Let G be a finite directed multigraph with a nonempty finite vertex set and finitely many edges. Loops are allowed. An edge e:u→v has real parameters `lambda_e>0`, `delta_e>=0`, and carries the forward resource estimate

    x_out <= (x_in+delta_e)/lambda_e.

A positive barrier vector is a choice `b_v>0` satisfying every edge inequality

    delta_e <= lambda_e b_v-b_u.                         (1)

**Theorem 1 (iff).** Such a vector exists if and only if every directed cycle has multiplier product at least one, and every cycle whose product is exactly one has zero debt on every edge. It suffices to check simple directed cycles, including loops.

This is a finite-state feasibility theorem. It says nothing by itself about which graph edges are physically available to a given source.

### Necessity and the exact path potential

For a chronological path e_1,...,e_m from v_0 to v_m, define

    Lambda_0=1, Lambda_j=product_(i<=j) lambda_i,
    Delta_m=sum_(j=1)^m Lambda_(j-1) delta_j,
    sigma_j=lambda_j b_(v_j)-b_(v_(j-1))-delta_j.

Then exact telescoping gives

    Lambda_m b_(v_m)-b_(v_0)-Delta_m
       =sum_(j=1)^m Lambda_(j-1) sigma_j >=0.            (2)

In particular,

    x_m-b_(v_m) <= (x_0-b_(v_0))/Lambda_m.               (3)

For a cycle based at v, (2) gives `Delta <= (Lambda-1)b_v`. Since `b_v>0` and `Delta>=0`, its product cannot be below one. If the product is one then Delta=0, and every coefficient in the debt sum is positive, so every edge debt is zero.

The same argument applies to closed walks. Alternatively, decompose a closed walk into simple cycles: their positive multiplier products multiply, so products at least one survive; a product-one walk consists entirely of product-one simple cycles and therefore has zero debt throughout.

### Sufficiency: difference constraints followed by scaling

For each simple cycle C containing at least one positive-debt edge, let

    p_C = number of its edges with delta_e>0.

The hypothesis gives `log Lambda_C>0` on each such cycle. There are finitely many simple cycles. If any contain a positive-debt edge, choose

    0<epsilon<=min_C (log Lambda_C)/p_C;

if none do, choose any epsilon>0. Define modified edge weights

    w_e=log lambda_e-epsilon*1_(delta_e>0).

Every simple directed cycle has nonnegative total w-weight, and hence so does every closed walk. Add a new source with weight-zero edges to all vertices and take shortest-path distances d_v. These distances are finite: deleting a nonnegative cycle never increases a path's weight, so a shortest path may be taken simple. They satisfy `d_v<=d_u+w_e`.

Set `r_v=-d_v` and `b_v^0=exp(r_v)`. Then

    r_u-r_v <= w_e,
    lambda_e b_v^0 >= exp(epsilon*1_(delta_e>0)) b_u^0.

Thus `s_e=lambda_e b_v^0-b_u^0` is nonnegative on every edge and strictly positive wherever delta_e is positive. A common finite scaling factor

    K >= max(1, max_(delta_e>0) delta_e/s_e)

makes `b_v=K b_v^0` satisfy (1). The maximum over an empty set is omitted. This proves sufficiency, including zero-debt product-one cycles and positive-debt edges that lie on no cycle.

### A finite exact linear-programming criterion

When lambda and delta are rational, feasibility is equivalent to the rational linear system

    b_v>=1 for every v;
    lambda_e b_v-b_u>=delta_e for every e.               (4)

Indeed any positive feasible vector can be scaled upward until its smallest coordinate is at least one; scaling preserves all inequalities because delta is nonnegative. Conversely (4) already gives positive barriers. A nonempty rational polyhedron has a rational feasible point, so a finite rational certificate exists. Exact linear elimination, exact LP, or exact cycle-product comparisons can therefore validate the criterion without treating floating-point logarithms as equality tests. The logarithms above prove existence, not a numerical certification shortcut.

### Fully rational constructive barriers

For rational data one can remove logarithms from the construction entirely. For each positive-debt simple cycle C let p_C be its positive-debt edge count. Choose an integer

    N>=max(1, max_C (2^p_C-1)/(Lambda_C-1)),
    kappa=1+1/N.

The maximum over no such cycles is omitted. Since N>=1, binomial expansion gives

    kappa^p <= 1+(2^p-1)/N,

so `kappa^p_C<=Lambda_C`. Put `mu_e=lambda_e/kappa` on positive-debt edges and `mu_e=lambda_e` on zero-debt edges. All cycle products of mu are at least one.

Add a source with multiplier-one edges to all vertices. Let d_v be the minimum product along its paths to v. Deleting cycles whose products are at least one can only lower a path product, so the minimum is attained among finitely many simple paths. Each d_v is a strictly positive rational and satisfies `d_v<=d_u mu_e`. Thus `b_v^0=1/d_v` has nonnegative original edge slack, strictly positive on positive-debt edges. The same rational scaling K used above produces an explicit rational barrier certificate. This supplies a direct constructive proof of rational feasibility without a floating-point logarithm or a generic polyhedron theorem.

## 2. What strict cycle products actually buy

If every simple directed cycle has product strictly greater than one, every infinite graph path has multiplier product tending to infinity, in fact exponentially in path length.

To see this, let V be the number of vertices and let mu be the minimum positive logarithmic product over simple cycles. An m-edge walk decomposes into simple cycles of at most V edges and a simple residual path of at most V-1 edges. The number of extracted cycles is at least `(m-(V-1))/V`. The finitely many residual paths have a uniform lower bound on their log weights. Therefore

    log Lambda_m >= mu*(m-(V-1))/V-C.                   (5)

If the graph is acyclic there is no infinite path, so the assertion is vacuous. Formula (3) then gives asymptotic control by the terminal-state barrier. If x is the actual positive integer height and `B=max_v b_v`, an infinite lawful path would eventually lie in the finite head `x<=floor(B)`. This still requires a separate analysis of that head, its legal outgoing edges, and ROOT. If x is another resource, an additional height or target implication is needed.

The inequalities do not force a graph path to exist. They also do not exclude bounded non-target dynamics in an arbitrary affine system: the one-state map `x→(x+3)/2`, with lambda=2, delta=3 and barrier b=3, fixes the positive integer 3 forever. This is an abstract affine hostile, not a claimed Collatz branch or Collatz counterexample.

## 3. Necessary hypotheses and small hostiles

- One loop `(lambda,delta)=(1/2,0)` has no positive barrier: `0<=-b/2` is impossible.
- One loop `(1,1)` has no barrier despite product one; zero debt on neutral cycles is essential.
- One loop `(1,0)` has every positive barrier and supplies no strict progress.
- Two edges `u→v:(1/2,0)` and `v→u:(2,0)` have a neutral cycle and valid barriers `(b_u,b_v)=(1,2)`.
- Replacing either of those zero debts by a positive number makes the neutral cycle infeasible.
- Two edges `u→v:(1/2,1)` and `v→u:(3,2)` have cycle product 3/2 and valid barriers `(4,10)`. An individual edge may have lambda below one when source and target barriers differ.
- Positive-debt acyclic edges cause no obstruction: a neutral zero-debt loop at u together with `u→v:(1,5)` admits `(1,6)`.
- Allowing negative debt invalidates the criterion: the loop `(1/2,-1)` admits b=1 despite product below one.
- Allowing a zero barrier invalidates the necessity of product at least one: `(1/2,0)` admits b=0.
- Finitely many edges are essential. Infinitely many parallel edges `u→v:(1,n)`, n=1,2,..., have no cycles but admit no finite barriers, since they require `b_v-b_u>=n` for every n.

## 4. No finite bank of reset-containing words covers all odd sources

Here and only here, a Collatz word means a nonempty finite list of exact odd valuations `(a_1,...,a_k)`, each a_i>=1, for

    U(n)=(3n+1)/2^v2(3n+1).

A reset is a valuation at least two. The complete exact residue guard is the one in section 4 of RESET_DEBT_COMPOSITION.md: if `(S n+C)/M` is the affine endpoint, require `S n+C=M mod 2M`, including the odd-terminal bit, and respect ROOT at proper prefixes.

**Theorem 2.** Let a finite nonempty bank consist of nonempty valuation words, each containing at least one reset. Let L be their maximum odd-step length. No word in the bank is lawful from any source in the infinite progression

    n=2^(L+1)t-1,   t=1,2,... .                         (6)

In particular, no such finite bank covers all positive odd inputs.

**Proof.** For `0<=j<=L`, direct induction gives

    U^j(n)=3^j 2^(L+1-j)t-1.                            (7)

For j<L,

    3U^j(n)+1=2[3^(j+1)2^(L-j)t-1],

whose bracket is odd. Thus each of the first L actual valuations is exactly one. Every bank word demands a reset by position L, so its guard fails. These actual endpoints strictly increase and stay above ROOT. QED.

**Corollary.** A finite bank of nonempty unit-budgeted physical height words cannot cover all odd sources. For an all-ones word of odd length k, the inverse multiplier is `lambda=(2/3)^k<1`, while its additive debt is positive. Such a word cannot satisfy `delta<=lambda-1`, so every unit-budgeted word contains a reset and Theorem 2 applies.

The same obstruction applies to any finite bank of genuinely descending odd-valuation words, whether or not it has a common unit budget: an all-ones word strictly increases every positive input. This is a bounded-horizon obstruction, not a refutation of adaptive unbounded-depth certificates.

## 5. Stronger: finite-state height barriers cannot give universal word coverage

**Theorem 3.** Suppose a finite graph has V>=1 states and finitely many edges, each labelled by a nonempty actual odd-valuation word of length at most an integer L>=1. Give each edge its **physical height** parameters

    lambda_e=2^A_e/3^k_e, delta_e=C_e/3^k_e>0,

so the exact branch is `(x+delta_e)/lambda_e`. Suppose positive barriers satisfying (1) exist. Then, from every source

    n=2^(VL+1)t-1,  t>=1,                               (8)

and from **any** initial graph state, there is no lawful path of V edges. ROOT is not reached during the relevant actual prefix. Consequently the graph cannot provide a next legal block at every unresolved state along every positive odd orbit.

**Proof.** The first VL actual valuations are all one, by (7) with L replaced by VL. If a path of V legal edges existed, their total odd length would be at most VL, so every letter of every edge word would be one. The V-edge state walk visits V+1 states and therefore contains a directed closed subwalk. Every edge on it has multiplier `(2/3)^k_e`; the whole nonempty closed subwalk has product below one. This contradicts Theorem 1. The actual all-ones trajectory strictly increases from n>=3, so an early ROOT stop cannot rescue coverage. QED.

This is stronger than saying that closure does not *prove* coverage: for this finite-state, finitely-many-fixed-words, physical-height model, universal continued coverage is impossible.

The restrictions are load-bearing. The theorem does not forbid finitely many **parameterized** rule schemes with unbounded valuation-one run lengths, infinite-state models, additional unbounded state coordinates, or a different signed/matrix resource. Nor does it forbid a finite bank from certifying large or useful subclasses. In particular, the infinite first-reset bank of the base note avoids the finite maximum-length hypothesis, but its arbitrary-source coverage remains unproved for independent reasons.

## 6. Validation and the surviving route

The separate `verify_finite_state_barriers.py` compares the exact cycle criterion with independent rational Fourier–Motzkin feasibility, tests the listed hostiles, and replays the all-ones guarded obstructions. It does not modify the base checker or its receipts. The finite experiments are controls; Theorems 1–3 are proved above.

The constructive use is to certify chosen reachable subgraphs and state-dependent debts. The universal route must retain an unbounded entry mechanism or a different resource, and it must prove actual guards and the final ROOT/finite-head implication. Merely replacing one scalar budget by finitely many height barriers cannot erase the long valuation-one obstruction.
