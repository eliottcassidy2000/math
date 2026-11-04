# Terminal-port strongness and an exact ordinal-component grammar

2026-11-04. **PROVED:** an exact seven-coordinate transfer under ordered sums, a necessary and sufficient strongness test for arbitrary suffixes, and an infinite clock-uniform +4 construction with nonzero defect. **NOT PROVED:** two-branch closure, unrestricted realization, native-operation legality, or the global Hamilton-path spectrum classification.

## 1. Sources and precise increment

This is a specialization and constructive extension of the inherited ordered path-cover machinery, not a new general substitution theory. Targeted sources used were S181 `multiport/PROOF.md`, `grammar/PROOF.md`, and `older/PROOF.md`. The latter identifies THM-3121 (complete path-cover substitution) and THM-3134 (endpoint jets). The newer `zero_ready_strong_switch_20261103/structure/README.md` already has a general decorated maximal-run interface and an affine sufficient criterion; no priority is claimed for those general principles.

The useful incremental observation is that the S181 assumption that the suffix is strong can be removed exactly. Nonstrong suffixes with a terminal port produce strong full tournaments. Their ordered component assembly has a closed finite-coordinate transfer that retains the joint correction, and gives the explicit uniform +4 family below.

A targeted pinned-main survey found no occurrence of this particular TT3→TT5 family or seven-coordinate law in the surveyed sources; this is not an exhaustive novelty claim. No startup helper or forbidden mistakes file was used.

## 2. Exact strongness, with no strong-suffix assumption

Let A and B be nonempty tournaments, W a subset of B, and r fresh. As in S181, let M(A;B,W) have A→B, W→r, and r→A and r→B\W. Then

**M is strong if and only if W meets the terminal strongly connected component of B.**

Proof: the strongly connected components of a tournament are linearly ordered, with every earlier component beating every later component. If the terminal component Z misses W, Z beats neither r nor any vertex outside Z, so remains a terminal cut in M. Conversely, if Z meets W, every vertex of B reaches that port inside B, then reaches r, then reaches A. From A every vertex of B is reachable. These routes connect all vertices, including different vertices inside A. This also handles singleton components. In particular W cannot be empty.

Thus requiring strong B is sufficient but unnecessary. This does not rescue an earlier sink-containing candidate at its *unchanged* port: a global sink must actually belong to W, which changes its evaluator and must be checked.

## 3. A seven-coordinate realizable signature

For an actual named tournament B with actual port subset W, define

- h = H(B).
- p = P_B(2), ordered covers by two directed-path slots, either allowed empty.
- f = the number of those covers whose first slot is nonempty and ends in W.
- g = the number whose second slot is nonempty and starts in W.
- j = the number whose two slots are nonempty, with first ending in W and second starting in W.
- e = the number of Hamilton paths ending in W.
- s = the number of Hamilton paths starting in W.

Write sigma(B,W)=(h,p,f,g,j,e,s). Empty W is allowed for component assembly and makes f,g,j,e,s zero. This signature is sufficient for the particular operation below, not for arbitrary births, reversals, physical re-markings, or native guards.

Take two nonempty named tournaments C,D on disjoint names, and their actual port subsets U,V. Let B=C→D, W=U union V. If sigma(C,U)=(h,p,f,g,j,e,s) and sigma(D,V)=(H,P,F,G,J,E,S), then

    sigma(B,W) = (
      hH,
      pP,
      fH+pF,
      gP+hG,
      jH+hJ+eS+gF,
      hE,
      sH
    ).                                                   (1)

All coefficients are nonnegative integers. This is a constructive transfer: given representatives, the ordered union and port union literally realize the result. Hence a finite or infinite bank of actual components generates an unbounded class of actual suffixes, not just formal signature vectors.

Proof: every path slot in an ordered sum is its C segment followed by its D segment. Two unrestricted slots factor, proving pP. A first-slot endpoint in C forces all of D into slot two, giving fH; an endpoint in D allows C in either slot, giving pF. The start formula is symmetric. For j, split the two named boundary ports into the four locations (C,C), (D,D), (C,D), and (D,C). Their contributions are respectively jH, hJ, eS, and gF. The last case allows both slots to cross from C to D and must not be omitted. A full Hamilton path traverses all of C then all of D, giving the remaining formulas.

Associativity on realizable signatures follows either from these bijections or associativity of ordered union. Physical names, rows, and port sets are retained alongside the signature, so iterating (1) has actual representatives at every step. To certify strong full attachments one tracks whether the last component's terminal SCC meets its ports. If components are themselves strong, this condition is simply that the last component has a nonempty port set.

## 4. Full count, transported defect, and a constructive sufficient criterion

The S181 identities, whose proofs do not require strong B, give

    H(M(A;B,W)) = H(A)(h-j)+P_A(2)f,
    E_r(M,2)-H(M) = H(A)j.                              (2)

For transitive A=TT_t, t>=1, put c=h-j. Then H_t=c+2^t f and the physical root defect is j. Nonstrong suffixes and negative c are permitted; only the actual full strongness condition in Section 2 is required.

Consequently a literal component-word rewrite from sigma to sigma' realizes the same digit d for every common transitive prefix if and only if

    f'=3f,          h'-j'=3(h-j)+d.                     (3)

The sufficiency includes strongness only if both actual terminal-port conditions hold. Counts at two different t already imply the two equations. A proof of both-digit induction would additionally need a recursively specified class of named component words, an actual rewrite for each d=2,4 from every word, preservation of the terminal-port condition, and return to that same class. Equation (1) makes checking a proposed rewrite exact; it does not produce missing rewrites. Native Hamilton-complement implementations would need their separate legal guards.

This is narrower and more concrete than the inherited general affine-interface criterion: its operation, realization, and boundary strongness are explicitly fixed.

## 5. Explicit infinite +4 family using nonstrong suffixes

Let B=TT_3, numbered 0→1→2, with W={0,1,2}. Let B'=TT_5, numbered 0→1→2→3→4, with W'={0,2,4}. Repeated (1), or direct two-slot enumeration, gives

    sigma(B,W)   = (1,8,7,7,6,1,1),
    sigma(B',W') = (1,32,21,21,12,1,1).

Both terminal vertices are ports. Thus for every t>=1 both full attachments are strong, and

    H_t = -5+7*2^t,
    H'_t = -11+21*2^t = 3H_t+4.

The physical root defects are 6 and 12, respectively. This is not a zero-ready construction, and the defect is not inferred from an observer phase. The first four pairs are 9→31, 23→73, 51→157, 107→325.

The transformation preserves all old names and all old suffix edges: add two successive terminal suffix vertices 3,4, make 4 an in-port, and remove old vertex1 from the port set. In the full tournament, that final change reverses the one old root edge 1→r to r→1. This is a direct tournament construction, **not** a certified legal native Hamilton-complement program. All old prefix and root names survive. An all-strong direct adjacency schedule is available: first add vertex4 as a terminal suffix vertex and port, then insert nonport3 after every old suffix vertex and before4, then reverse 1→r. Each completed step has an actual terminal port, namely4 after its birth. This verifies strongness of these adjacency edits, not native Hamilton-complement guard legality. The independent audit checks this schedule at eight prefix lengths.

The example illustrates why restricting to strong suffixes discards useful valid attachment families. It neither continues its own output under both digits nor says anything about the actual 170807 lineage.

## 6. Binary form for arbitrarily many physical ports

For B=TT_m, encode W by f=sum_(i in W)2^i, since a first slot ending at i can place each earlier vertex in either slot and must put every later vertex in slot two. Scan positions from 0 upward, retaining (f,g,j,e). When a new terminal vertex is born with port bit b in {0,1}, update

    f_new=f+b*2^m,
    g_new=2g+b,
    j_new=j+b*(g+e),
    e_new=b.                                           (4)

Start with the empty accumulator f=g=j=e=0; this is only an initialization convention, not use of (1) with an empty tournament. After the final bit, a strong full attachment exists exactly when b=1. The current g is the bit-reversal value and e is the latest bit; hence g+e is even, so j remains even. Formula (4) transports a genuinely changing joint defect while allowing arbitrary port sets and unbounded suffix order.

The example in Section 5 is the transition f=7 to f'=21 with j=6 to j'=12. Shifting both port patterns k places by prepending k nonport singleton components gives f=7*2^k and f'=21*2^k with the same j values, yielding the same clock-uniform +4 relation. No general digit-closure claim follows from this bit encoding.

## 7. Independent finite controls and limits

`verify.py` uses ordered permutation enumeration to compute all seven observables. It checks (1) for all labelled factors through order3 with combined order at most5, including empty component port sets. It independently checks Section 2 and both equations (2) for every labelled suffix through order4 and every port set, using a one-vertex prefix. Reachability tests are independent of component decomposition. There are 1,380 ordered-sum transfer checks and 1,098 attachment/strongness/defect checks. A further 62 transitive-suffix/port checks verify (4), and the two displayed full signatures are checked explicitly. Finally an independent subset-path DP checks the four displayed +4 pairs and full strongness. Normal and optimized outputs byte-match. All verification guards are explicit exceptions and survive Python optimization.

This packet supplies an exact, realized unbounded *ordered-component* grammar and one nontrivial infinite digit map. A global classification still needs coverage and both-digit closure, or some different global induction. Neither seven-coordinate evaluability nor permitting arbitrary component words proves those missing existence statements.
