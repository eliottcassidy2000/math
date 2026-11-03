# Matched sources for transported edge jumps

**Status.** Exact literal constructions and profile identities. These match explicit sources, not an arbitrary prescribed Hamilton-path count. No global optimality or new spectrum inclusion is claimed. The old source vertices are retained when the indicated edge is reversed.

Write H(X) for the directed Hamilton-path count and F_X(r) for ordered r-slot path covers, empty slots allowed. Let T_k be the k-vertex transitive tournament, T_0 empty, and let G(X) be the triangle substitution X→c→d→X. We use the earlier identities

    F_{T_k⊕X}(r)=r^k F_X(r),
    H(G(X))=H(X)+F_X(2).

For every nonempty X, G(X) is strongly connected. Both displayed formulas apply to the two physical orientations of the same edge. These facts are given in S159's [positive profile operator](https://github.com/eliottcassidy2000/math/blob/7589463e54e20585aff6a2e68019ee4a4af582cd/05-knowledge/results/session159_constructions_20261103/profile_operator/NOTE.md), and the compact pair below can also be checked directly.

## 1. Four vertices yield exact strong +18 and +130 pairs

Start with 0<1<2<3, all edges pointing upward. In X− reverse 0–2. In X+ reverse both 0–2 and 0–3; the only edge changed between X− and X+ is 0–3. Their full cover profiles are

    F_{X−}(r)=r^4+2r^2,
    F_{X+}(r)=r^4+4r^2.

To prove these, use the odd-cycle packing expansion from [S159](https://github.com/eliottcassidy2000/math/blob/7589463e54e20585aff6a2e68019ee4a4af582cd/05-knowledge/results/session159_constructions_20261103/cycle_profiles/NOTE.md): the only possible odd cycles on four vertices are triangles, contributing 2r² each. X− has one triangle and X+ has two. Hence (H,F₂)=(3,24) and (5,32), respectively, and their difference at every slot count is 2r².

For any k≥0, G(T_k⊕X−) and G(T_k⊕X+) are strong, differ in the same old edge, have order k+6, and realize the exact transition

    3 + 24·2^k  →  5 + 32·2^k,
    jump = 2 + 8·2^k.

Thus k=1 gives the literal strong seven-vertex transition 51→69, a +18 jump. At k=4, 387→517 is a +130 jump on ten vertices. At k=0 the pair is 27→37 (+10). These are source-and-target statements, not a free +18 or +130 instruction at an arbitrary source.

## 2. A compact full-profile match for −126

Use the two-apex convention of [THEORY](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/apex_wedge_constructor_20261103/THEORY.md): the base vertices 0,...,6 are transitive; apex a beats base i exactly when w_i=1, apex b beats base i exactly when u_i=1, and Y+ has a→b while Y− has b→a. Take

    w=0000111,  u=1000010.

The complete ordered path-cover polynomials are

    F_{Y−}(r)=r^9+16r^7+46r^5+32r^3+2r,
    F_{Y+}(r)=r^9+14r^7+50r^5+32r^3+2r.

A finite, reproducible proof of these explicit identities is the odd-cycle packing enumeration in `check.py`: for Y− the numbers of directed cycles of lengths 3,5,7,9 are 8,23,16,1, and for Y+ they are 7,25,16,1; there are no disjoint pairs of odd cycles. The checker enumerates all directed odd cycles (anchored at their least vertex) and all disjoint packings, and separately computes H and F₂ by subset-and-endpoint path dynamic programming. Thus

    (H−,F₂−)=(97,4292),  (H+,F₂+)=(99,4164),
    F_{Y+}(r)−F_{Y−}(r)=2r^5(2−r²).

For every k≥0, actual strong tournaments G(T_k⊕Y−), G(T_k⊕Y+) have order k+11 and counts

    97+4292·2^k  →  99+4164·2^k,
    jump = 2−128·2^k.

In particular k=0 gives the strong eleven-vertex transition 4389→4263, a −126 jump. The Y± themselves do not need to be strong: G makes each wrapper strong. The proof works for every k; enumerations are controls for the finite polynomials and literal small graphs.

The complete difference 2r⁵(2−r²) coincides with the incoming balanced P=6,R=8 response in [PROFILE P3](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/apex_wedge_constructor_20261103/PROFILE.md). The incoming pair instead has (H−,F₂−)=(2625,131004) and (H+,F₂+)=(2627,130876), yielding wrapped transition 133629→133503. Equality of the response polynomial persists in any common one-hole context built from fixed ordinal-sum factors and G because these act linearly on the profile difference. It does **not** assert equality of the individual tournaments or of their absolute counts.

## 3. Why unrestricted matching is impossible

S157's [gap proof](https://github.com/eliottcassidy2000/math/blob/main/05-knowledge/results/session157_gaps_20261103/gaps/gaps-note.md) establishes that 7 and 21 are not tournament Hamilton-path counts. Consequently no construction can support a +2 jump at every realizable source (5→7 fails), a +18 jump at every source (3→21 fails), or a −126 jump at every source (133→7 fails). Counts 3 and 5 occur above, and 133 is realized by S160's W(2,2). The cited excluded outputs explain the scope issue; they do not rule out a +130 all-source rule, nor solve the remaining spectrum problem.

## Controls

Run `python3 check.py` and `python3 -O check.py`. It independently enumerates explicit odd cycles and disjoint odd-cycle packings, obtains H and F₂ from a separate Hamilton-path subset/endpoint DP, and directly counts literal wrapped examples, including both ±126 sources and the 387→517 pair. Assertions remain active under `-O`. These finite computations do not establish global minimality or universal source matching.
