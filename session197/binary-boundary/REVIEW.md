# Independent review of the binary deletion boundary

**PASS, no correction.** This review binds the frozen `session197-orbit-wildcard/binary-boundary/PROOF.md` SHA256 `a6c358fc78204c266d05ad04d458877abe2aa980e25403740511d1af2ab0abf7`, verifier `ca4a3c6d2d6556c5e4846733617731979e939a38b4effa3a098e6fe77a7769ff`, and receipt `875b1a2a1dba4bd7b4fef22c871f67e04758f0e0e4b0fa76609b7a13c80137db`. The source files were left unchanged.

## All-order proof review

The physical single-edge flip is exact. In the original ordered join I->O, every lost Hamilton path uses I_a->O_a as its only crossing and has e squared choices. In the reversed graph every gained path uses O_a->I_a exactly once. The two side restrictions form the stated q arrangements. Their remaining cross joins cannot require the missing I_a->O_a arc because both marked endpoints have already been used at the central edge. Empty complementary paths cause no extra multiplicity or missing join, including the singleton-side case. Thus the deletion count h squared plus q squared minus e squared is established physically rather than extrapolated from a positive clone length.

The endpoint reciprocity substitution d=-e is consistent with the independently reviewed squarefree proof and has the correct sign. Splitting q by path-length parity gives h squared plus four times q_even q_odd. Since the singleton endpoint path has a positive complement count, q_odd is positive; q_even is zero exactly when the marked vertex has no in-neighbour. This proves the stated equality criterion, including weak side tournaments. The mod 4 assertion also uses the established oddness of tournament Hamilton-path counts.

The Hamilton-cycle count at zero clone length is f squared for side order at least two. Every spanning cycle must use the only O-to-I edge, and the other join is valid because neither endpoint can then be the same marked vertex. At side order one the two-vertex graph has no cycle, exactly as separated in the proof. At positive clone length, two clone runs would require two disjoint uses of the single exceptional crossing edge. The unique run and the two complete side paths give f squared with no additional cycle-rotation factor.

For the single-backward-arc side S_k, the literal subset allocation around its backward edge gives h=q=2^k+1 and e=f=1. The displayed spanning cycle uses every vertex once and every claimed arc has the correct orientation. Thus the boundary family is strong, self-converse and has exactly one Hamilton cycle. The count formula, the first six values and the congruence k=0 or 2 modulo 3 for divisibility by 7 are correct. Divisibility is not equated with prime quotient: k=6 is explicitly composite. Every positive clone length has a nontrivial odd square factor; the exceptional possible value 49 for a 7p square is ruled out by the additional factor at least 3. No returning grammar, infinitude or full-spectrum conclusion is claimed.

## Independent computation

`check.cpp` was written independently and imports no source checker. Its graph construction places clones last rather than first and builds each unordered pair directly from its block membership.

1. Literal directed-word DFS checks every 285 rooted labelled side tournament through order 4. It counts the new graph, paths using the newly reversed edge, the original ordered join, and paths lost through the old edge separately.
2. The same literal DFS checks all 570 positive-length contrast graphs, including their anchored Hamilton-cycle counts.
3. The six quadratic-family graphs, of orders 6 through 16, are counted in two independent ways: literal Hamilton-word DFS and inclusion-exclusion of fixed-length directed walks over induced vertex subsets. All path counts match, and a separately anchored cycle enumeration gives exactly one cycle each.
4. Every arc of the paired anti-involution is checked on each family graph. The complete 14-vertex 2177 row list agrees with the frozen source receipt.
5. Complete trial division confirms prime quotients 7,23,311, and direct multiplication checks 8449=7*17*71.

The independently written checker agrees under optimized and unoptimized assertion-enabled compilation. The supplied source verifier was also compiled and replayed under both modes; both outputs exactly match its frozen receipt. No scratch scan or finite test is used to replace the all-order arguments above.

This review covers only the frozen deletion-boundary packet. A later fixed-word or automorphism refinement is not silently included. No literature-priority claim is made or audited here.
