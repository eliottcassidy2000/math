# A common infinite fan with both physical ternary children

Session 197, 2026-11-05. **All-order paired construction from an explicit common parent family.** This is not a universal ternary grammar or a claim about eventual spectrum coverage. The two child branches use the stated different initial port pairs; neither a shared returning marked interface nor closure of the child family is asserted.

## 1. Statement

For every positive integer m there is an explicit strong tournament P_m with a distinguished vertex b such that

    H_m=H(P_m)=18*2^m+13,
    (H_m,e_b,R_b,H(P_m-b))=(H_m,13,H_m,13).            (1)

The same actual two-sided two-port gates, addressing the same named old vertices at every m, construct strong tournaments with counts

    3H_m+2=54*2^m+41,
    3H_m+4=54*2^m+43.                                (2)

At m=1 these are the named H49 input and the physical H149/H151 extensions. The formulas hold for arbitrary m because the leftward gates commute with a proved one-seam chain construction, not because a few cloned examples fit a recurrence.

For contrast, the [S196 equal-response theorem](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/global-induction/PROOF.md) excludes both values in (2) from every P_m when it starts marked only at b and is restricted to its arbitrary-gadget one-root language. The present construction changes which old vertices are exposed. It does not contradict that theorem.

The proof combines an inherited one-seam transfer with two different physical core extensions. It also gives a general prefix-triangle compiler for the digit `2*s_a(S)`, for every actual rooted source S. No literature-priority or new numerical-record claim is made.

## 2. Explicit cores and fixed leftward gates

All adjacency rows are outgoing bitmasks on zero-based named vertices. Let

    C=[26,16,3,6,12], distinguished a=2.

Its exact rooted profile is

    (H,e_a,R_a,H(C-a),s_a,L_a)=(13,2,18,3,3,15),      (3)

where s_a counts Hamilton paths starting at a and L_a counts ordered two-covers whose first path ends at a. In R_a the second path starts at a and the first can be empty; in L_a the first is nonempty and the second can be empty.

A leftward birth with exception set S adds a vertex losing to precisely the old ports in S and beating every other old vertex. Forgetting a port changes only its exposed status, not the graph. The two old ports are ordered, and a mask selects their positions; the forgotten position is in `(old0,old1,new)`.

**Gate A, for the +2 branch.** Expose old vertices `(0,1)`.

1. Add a universal source u and forget old port 1: `(L,0,1)`.
2. Add v losing only to the retained old vertex 0, and retain both current ports: `(L,1,2)`.

On C, the resulting core is

    A=[90,16,3,6,12,31,62], a=2,
    (H,e_a,R_a,H(A-a),s_a,L_a)=(41,6,54,9,3,99).      (4)

Here u=5,v=6, and the final exposed pair is `(0,5)`.

**Gate B, for the +4 branch.** Expose old vertices `(1,4)`. Add one vertex losing to both and beating every other old vertex, then retain the old pair: `(L,3,2)`. The resulting core is

    B=[26,48,3,6,44,13], a=2,
    (H,e_a,R_a,H(B-a),s_a,L_a)=(43,8,54,13,6,72).      (5)

Thus both core gates triple R_a, while their H outputs are respectively `3*13+2` and `3*13+4`. The supplied verifier independently computes all six entries in (3)–(5) both by subset endpoint DP and by literal permutations and cuts. These are finite exact core facts; their all-order use is justified next.

C is strong, for example along the directed cycle `0,1,4,3,2,0`. Gate A's final graph is strong: old vertices reach0, then v, then u, which reaches the old core. Gate B adds a vertex with both in- and out-neighbors in the strong core, so B is strong as well.

## 3. The exact one-seam fan

Given a nonempty rooted core `(X,a)` and m>=1, build F_m(X,a) as follows:

- Add a distinguished vertex b.
- Add a transitive chain `z1 -> ... -> zm`.
- Orient every core-to-chain edge forward, every chain-to-b edge forward, and every core-to-b edge forward except for the single reverse seam `b->a`.

All internal core arcs remain unchanged. Equivalently the cross layout is

    X -> Z -> b, with X -> b except b -> a.           (6)

**Lemma 1 (one-seam fan identity).** If h=H(X) and R=R_a(X), then

    H(F_m)=h+2^m R,
    ends_b(F_m)=H(F_m-b)=h,
    R_b(F_m)=H(F_m).                                 (7)

The first identity is the [S195 rooted one-seam rule](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session195/rooted-interfaces/PROOF.md) with a source-first transitive `(m+1)`-vertex gadget marked at its sink. Here is a direct proof of all entries.

A Hamilton path not using `b->a` must end at b, since b has no other out-neighbor. Its earlier word is a core Hamilton path followed by the whole chain, giving h choices. A path using the seam splits the core into an ordered R_a-cover `(P,Q)`, Q starting at a. Each chain vertex can independently belong to the increasing run before b or the increasing run after Q; empty runs are allowed. These 2^m assignments give `2^m R` paths, and restriction recovers each choice uniquely. This proves H. The endpoint-b statement follows from the first case. Deleting b leaves the ordered join `X -> Z`, of count h.

For an R_b-cover, either its second path is the singleton b, contributing h, or it starts `b,a`. In the latter case the old core portions form an R_a-cover and the chain vertices again choose their two increasing runs independently. Thus `R_b=h+2^m R=H(F_m)`. QED.

If a reaches every vertex of X, the fan is strong: every vertex reaches the chain and b, b reaches a, and a reaches the whole core. In particular this holds for all three strong cores above.

## 4. Same named gates at every chain length

Set `P_m=F_m(C,2)`. For literal labels, keep core vertices 0 through 4, set `b=5`, and use chain vertices `6,...,m+5` in increasing order. The whole construction has m+6 vertices. At m=1 its rows are exactly

    [122,112,67,102,108,4,32].                         (8)

Applying Gate A to P_m always starts with the same old pair `(0,1)`. Its new source has label `m+6`, its second new vertex has label `m+7`, and the final ports are `(0,m+6)`. Applying Gate B always starts and ends with the same old pair `(1,4)`, adding its new vertex at label `m+6`. Neither gate refers to an internal chain position or changes with m.

Every newly added gate vertex beats b and every chain vertex, because the exceptions belong to the unchanged old core. Consequently the following is an equality of actual graph constructions up to the displayed natural relabeling:

    Gate A(F_m(C,2)) = F_m(A,2),
    Gate B(F_m(C,2)) = F_m(B,2).                     (9)

For A, the new core vertices5,6 of its canonical fan correspond to the two gate births, its fan root7 corresponds to the old b=5, and its chain corresponds in order to the old chain6 through m+5. For B, its new core vertex5 corresponds to its single gate birth, its fan root6 corresponds to old b=5, and the chains likewise match. The verifier checks these graph equalities at adjacency level, not only their counts.

Substitute the three exact core pairs

    (h,R)=(13,18), (41,54), (43,54)

in (7). This proves (1)–(2) for every m. The root b remains the same **old physical vertex5** after either gate. Its resulting four-coordinate states are

    Gate A: (3H_m+2,41,3H_m+2,41),
    Gate B: (3H_m+4,43,3H_m+4,43).                   (10)

In particular the deletion count is not silently kept at 13. Both final graphs are strong by (9) and the strongness of A and B.

For m=1 the gates give exactly the previously checked rows

    H149: [378,112,67,102,108,4,32,127,254],
    H151: [122,240,67,102,236,4,32,109].               (11)

The first six parents are 49,85,157,301,589,1165. The construction covers their two literal children and every later member of the same exponential parent family.

### Exact comparison with the older interface

Equation (1) is an equal-response state with d=13, h>=49, and

    h-3d=18*2^m-26>=10.

The S196 arbitrary-gadget equal-response theorem therefore gives, from that input marked only at b, complete one-root descendant support through `3h+4` equal to

    {h,2h-13,3h}.                                   (12)

Neither child in (2) belongs to (12). The state is primitive as well: `gcd(h,13)=1` for every m. Thus the new construction escapes an infinite family of primitive marked-state obstructions, not only the initial H49 example.

This is a comparison of different future interfaces. The two leftward branches expose different old vertex pairs, fixed within each branch for all m. They are not both operations on a single unchanged two-port state. The output cores A and B are not shown to carry the same two gates again. Hence (2) supplies no universal ternary induction, and no claim that all larger odd counts are covered follows from it.

## 5. A general source-side triangle compiler

The +2 branch has a broader physical explanation. Let `(S,a)` be any nonempty rooted tournament, and put

    g=H(S), s=s_a(S), E=e_a(S), Q=R_a(S), L=L_a(S),
    B=P_(S-a)(2), d=H(S-a), delta=[|S|=1],

where P(2) counts ordered spanning two-path covers with either slot allowed empty.

First add x with `x->S-a` and `a->x`, retaining the distinguished a; call this X0. Then add a universal source u, followed by v losing only to x and beating every other old vertex; call the result X1. The vertices x,v,u form a directed triangle. Their edges to S all point toward S except `a->x`.

The one-seam transfer for a rooted prefix G feeding into S is

    H = g H(G)-s e_x(G)+L R_x(G),
    e_a = E H(G)-delta e_x(G),
    R_a = Q H(G)+B R_x(G),
    H(delete a)=d H(G).                              (13)

The first row separates ordinary concatenations from the unique reverse seam. In the R_a row, covers not using the seam put all of G in the first path and give QH(G); covers using it begin their second path `a,x`, giving an R_x-cover of G and an arbitrary ordered two-bin cover of S-a. Deletion removes the only reverse edge. This proves the displayed laws directly; they are the inherited S195 returning interface.

The rooted singleton prefix has `(H,e,R)=(1,1,1)`, while the rooted triangle has `(3,1,3)`. Thus, writing `(h,e,R,d0)` for X0,

    X1: (3h+2s, 3e+2delta, 3R, 3d0).                (14)

After adding the same positive fan chain to X0 and X1, (7) gives the actual scalar transformation

    H(F_m(X1,a)) = 3 H(F_m(X0,a)) + 2s.              (15)

When S has a unique Hamilton start at a, s=1 supplies +2. When s=2 it supplies +4. Each is an all-order construction from the stated source class, not a promise that an arbitrary desired parent has such a presentation. Even if S is weak, a reaches the whole of X0 through x and the whole of X1 through x,v,u, so both fans in (15) are strong.

For the common family above, use the actual source

    S=[8,1,3,6], a=1,
    (g,s,L,E,Q,B,delta,d)=(5,1,9,2,6,12,0,3).

Adding x and relabeling it as old vertex 0 gives C; replacing the prefix by the triangle gives A. Hence the +2 branch is explained uniformly by (14), not only by the core profile table.

### Why the +4 core is a genuinely different synthesis in this presentation

Keep that exact S and its mark. Allow **any** nonempty rooted prefix G in the one-seam construction (13), of arbitrary order. An output with `R_a=54` must satisfy

    6 H(G)+12 R_x(G)=54,
    H(G)+2R_x(G)=9.

Since H(G) is positive odd and `R_x(G)>=H(G)`, the only numerical candidates are `(H,R)=(1,4)` and `(3,3)`. The scalar formula then gives respectively `41-e_x(G)` and `42-e_x(G)`, both less than 43. Therefore no prefix in this entire factory supplies the H43,R54 core. Gate B does so by a different leftward two-port incidence pattern.

This is a restricted synthesis distinction, not a claim that the core B has no other representation or that the general s=2 compiler cannot attain 43 elsewhere.

## 6. Verification and scope

    python3 verify.py > receipt.json
    python3 -O verify.py > receipt.optimized.json
    cmp receipt.json receipt.optimized.json

The standard-library verifier uses explicit exception checks. It independently counts the three full core profiles by subset DP and literal permutation/cut enumeration; checks both literal gate constructions; constructs the parent and both children at m=1 through 6, verifies the exact physical graph relabelings, marked root states, and strongness; and checks the general triangle compiler on every labelled source of orders 1 through 4 and every root, including the singleton endpoint correction.

Those 285 rooted source checks control the proof of (14), not replace it. The finite core facts and universal one-seam identities prove the infinite fan. No finite-to-infinite support extrapolation is made. All files were created in the session-specific scratch; no dirty checkout or worktree was changed.
