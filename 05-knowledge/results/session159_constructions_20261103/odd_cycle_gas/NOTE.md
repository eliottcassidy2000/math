# Odd-cycle gas on transitive tournament blow-ups

**Status.** Exact cycle-packing polynomial for a standard strong tournament family, obtained by a labeled hard-core gas expansion. The Hamilton-path specialization recovers the already known S158 exponential generating function; the general determinant formula, refined activity parameter and the two-C-vertex formulas below are the contribution here. No novelty relative to the literature is claimed. The finite checks are controls, not proofs.

## General determinant identity

Let S be any fixed finite tournament with adjacency matrix B on m vertices. Replace each vertex i by an internally transitive block of size n_i, orienting edges between blocks according to S. Let J_n(t) be the odd-cycle packing polynomial of this blow-up. Set

    U(x) = diag(e^{x_1}−1,...,e^{x_m}−1),
    D(x) = det(I−U(x)B),  s = x_1+...+x_m.

Then in Q[t][[x_1,...,x_m]],

    sum_{n_i≥0} J_n(t) product_i x_i^{n_i}/n_i!
       = e^s [D(−x)/D(x)]^{t/2}.                         (G)

Indeed, the maximal block-runs in every directed cycle follow a closed directed walk in S. Each run is a nonempty set of labels in its block, in the uniquely determined increasing order, and contributes e^{x_i}−1. Therefore

    C_all(x) = sum_{r≥1} tr((U(x)B)^r)/r = −log D(x)

is the multivariate labeled EGF for directed cycles. The factor 1/r correctly removes the choice of starting run even for a periodic color sequence: all nonempty run-sets are pairwise disjoint labeled sets, so no nontrivial rotation fixes the labeled structure. The odd-total-degree projection is (C_all(x)−C_all(−x))/2 = (1/2)log[D(−x)/D(x)]. Finally exp(s+t C_odd) counts vertex-disjoint cycle collections and unused vertices by the labeled exponential formula. As induced blow-ups depend only on subset cardinalities, this application is valid for the internally transitive blocks.

## Directed triangle base

Let T(a,b,c) have three linearly ordered, internally transitive blocks A,B,C of sizes a,b,c and all edges A→B, B→C, C→A. If a,b,c≥1 then T is strong. Set

    J_{a,b,c}(t) = sum_F t^{|F|},

where F ranges over vertex-disjoint collections of directed odd cycles, each cycle counted up to rotation. In particular, J_{0,0,0}=1. Write s=x+y+z and Q=(e^x−1)(e^y−1)(e^z−1). Then the full trivariate exponential generating function is the formal power-series identity

    sum_{a,b,c≥0} J_{a,b,c}(t) x^a y^b z^c/(a!b!c!)
        = e^s ((1+e^(−s)Q)/(1−Q))^(t/2).                (1)

The rational power means its formal binomial expansion over Q[t][[x,y,z]], where the parenthesized series has constant coefficient 1. Each exponentially normalized coefficient a!b!c! [x^a y^b z^c] on the right is in Z_{≥0}[t] by the combinatorial definition.

## Proof of (1)

A directed cycle necessarily meets all three blocks: any two blocks form a transitive tournament. As one traverses a cycle, consecutive maximal block-runs follow the cyclic order A,B,C,A,...; suppose there are k runs of each type. Vertices inside a run must appear in increasing order. For any selected sets of a,b,c labels, an ordered partition into k nonempty runs of each color therefore gives all cycles with this run count; simultaneous cyclic rotation of the k triples of runs is the only overcount. Thus the number of cycles on the *specified* selected labels is

    D(a,b,c) = sum_{k=1}^{min(a,b,c)} (k!)^3 S(a,k) S(b,k) S(c,k)/k,  (2)

where S is a Stirling number of the second kind. The formula has value 0 when any block is empty. Its exponential generating series is

    C_all(x,y,z) = sum_{k≥1} Q^k/k = −log(1−Q).

The number of vertices on a cycle is odd precisely in the odd-total-degree part, so

    C_odd = [C_all(x,y,z)−C_all(−x,−y,−z)]/2
          = 1/2 log((1+e^(−s)Q)/(1−Q)).

For each block, the induced tournament on any chosen subset is canonically isomorphic to the same three-block construction with smaller block sizes. Hence the labeled exponential formula applies to collections of disjoint cycles, with activity t per cycle, while isolated unused vertices have generating function e^s. This gives exp(s+t C_odd), namely (1). The identity counts disjoint cycles, not all cycles independently: the exponential formula partitions the labels among components.

## Hamilton-path specialization

The [odd-cycle-packing identity](https://github.com/eliottcassidy2000/math/blob/3d6628b1804cf80457da16aa52489f016dfc99c4/05-knowledge/results/session157_gaps_20261103/cycles/README.md) H(T)=J_T(2) yields

    sum H(T(a,b,c)) x^a y^b z^c/(a!b!c!)
      = (e^s+Q)/(1−Q).                                  (3)

This recovers the triangle-block Hamilton-path generating function in [S158](https://github.com/eliottcassidy2000/math/blob/6bb391f902df4a6fd34b947a7923f64a0d5e6057/05-knowledge/results/session158_creative_20261103/cyclic_blocks/NOTE.md): writing U=e^x−1, V=e^y−1, W=e^z−1 and subtracting the empty path from (3) gives exactly (U+V+W+UV+VW+WU+3UVW)/(1−UVW). For explicit coefficients, expanding by powers of Q gives

    H(T(a,b,c)) = sum_{k≥0} (k!)^3 S(a+1,k+1)S(b+1,k+1)S(c+1,k+1)
                 + sum_{k≥1} (k!)^3 S(a,k)S(b,k)S(c,k).  (4)

Both sums terminate. The +1 in the first arguments accounts for each e^x(e^x−1)^k coefficient. Formulas include empty blocks with the convention H(empty)=1. A naive ordinary generating function that assigns one path to each block-word would be wrong: the same block can occur in multiple runs, and its labels can be partitioned among them in many ways.

## Complete polynomial with two vertices in C

Fix a,b≥1 and c=2. Each directed cycle meets C, so J_{a,b,2} has degree at most 2. Put

    A_n = 3^n−2^(n+1)+1,   E_n = 1+(−1)^n.

Then J=1+p_1 t+p_2 t², where

    p_2 = [A_a A_b + 2(2−2^a)(2−2^b) + E_a E_b]/4,        (5)
    p_1 = [1+3(2^a−1)(2^b−1)+A_a A_b−E_a E_b]/2.         (6)

For (5), each of two disjoint cycles uses one distinguished C vertex and a nonempty subset of each of A and B. Because its C vertex is unique, each such cycle has one run in each block, so it is uniquely determined by those subsets. Its length is odd exactly when the two subset sizes have the same parity. For a block of size n the signed assignment polynomial onto unused, first cycle, second cycle, requiring both cycle subsets nonempty, is

    f_n(u,v) = (1+u+v)^n−(1+u)^n−(1+v)^n+1.

The oddness conditions are imposed by (1/4) sum_{u,v=±1} f_a(u,v)f_b(u,v). Since the two C vertices distinguish the two cycles there is no extra factor of 2. The four values of f_n are A_n, 2−2^n, 2−2^n, and E_n; this proves (5). Formula (6) follows from (4) and H=1+2p_1+4p_2. These expressions are integers by the counting proof.

Examples: J_{2,2,2}=1+14t+4t² and H=45; J_{2,3,2}=1+44t+12t² and H=137. These are strong tournaments. Formula (1) provides all coefficients for arbitrary block sizes, but does not classify the Hamilton-path spectrum.

## Independent verification

Run `python3 check.py` and `python3 -O check.py`. A direct simple-cycle permutation enumerator plus an independent vertex-disjoint-packing recursion computes J for all 165 nonnegative triples a+b+c≤8, including the empty tournament. A separate endpoint-subset dynamic program computes H. All agree with the Stirling-cycle recursion derived from (1), equation (4), H=J(2), and (5) whenever a,b≥1 and c=2. The outputs agree under both interpreter modes. The computations do not prove the general identities.

For a separate general-base control, `python3 check_general.py` (also with `-O`) checks (G) after substituting (x_1,x_2,x_3,x_4)=(q,2q,3q,4q) through q^6 at t=0,1,2,3,4, for a four-vertex strong base. The left side comes from direct cycle and packing enumeration; the right side comes from exact rational polynomial determinant and series arithmetic. This weighted univariate specialization is a useful control, not a coefficientwise proof.
