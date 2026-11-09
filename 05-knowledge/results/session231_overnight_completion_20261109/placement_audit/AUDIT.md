# Independent target-placement audit

Verdict: PASS for the substantive restricted-language placement and graft-collapse theorem, with a degenerate-window wording caveat and a complexity qualification below. This audit assumes the previously audited twelve-bit family and seven-state ordinal composition theorem; it does not infer arbitrary-target coverage.

## Algebra and endpoints

Writing X=Bq+pA, composing left, core, and right gives A_total=Xt+phF and T_total=XK+sBG+sUh+spT+sWS+phZ+pEM. Thus 2A_total+phs−T_total is precisely displayed (1). The only j-dependent term is −spT_j. Empty contexts must be omitted as stated. Positive Hamilton counts force ps=1 exactly when both factors equal one. The standard strong-component Hamilton-cycle argument proves H=1 iff transitive, with singleton components handled separately.

Equation (2) agrees with the inherited exact binary formula, including suffix leading zeros, digit reversal, adjacency correction, prefix last bit, and suffix first bit. The implementation-independent literal controls cover all 12 normalized binary contexts of total length at most two at k=0. Base upper/lower endpoints coincide at k=0. For odd N, both B and N are odd, so (3) is sufficient as well as necessary: j=(B−N)/2 is an integer in [0,k].

The base begins with a selected singleton. A normalized nonempty prefix begins with a selected singleton too. Consequently the root reaches that first vertex and then every later ordinal block, and every vertex reaches the root through v. Leading unselected prefix singleton vertices are indeed initial SCCs with Hamilton multiplier one, so removing them preserves the scalar image.

T_total≤A_total implies N≥A_total+45. For a>0, qP≥11520 R_k 2^(a−1); for a=0, A≥7516 R_k. Hence (4) is valid also at m=0. Search should use exact integer inequalities rather than rounded floating-point logarithms. If N−45<5760, the feasible set is empty. An explicit safe bound is k≤floor(log_4096((N−45)/5760)) when that ratio is at least one, then m≤floor(log_2((N−45)/(5760·4096^k))).

The vertex-insertion construction is an injection: before the first beaten vertex, all earlier vertices beat the new vertex, and deletion of the inserted vertex recovers the original path. This proves nondecrease, not strict increase. No strictness is needed. In (5), the exact excess after subtracting 2k is (1552824k+1118023)/24843>0, so the bound holds at k=0 as well. The floor inverse bound is k≤floor(log_4096(24843N/245051264)) when its argument is at least one. The stated half-count/count logarithmic error constants follow.

## Complexity qualification

For fixed k,m, there are 2^m+m2^(m−1) normalized contexts (the first term corresponds to a=0). Summing m through its bound and then k gives O(N log N) candidates in the numeric value N. This is pseudopolynomial in input bit length, not polynomial time in log N. The manuscript correctly says candidate contexts. It should not be paraphrased as O(N log N) bit operations: evaluation of f and big-integer arithmetic costs extra unless analyzed separately.

## Four-state graft interpretation

The root is always the original closing root, not a gadget vertex or a freely rerooted vertex. The base four-state has E=d=h and U=2A+h. In the graft deletion, old deletion vertices precede all gadget vertices, hence d'=Hd. A full path ending at the original root must finish through the appended gadget and gives E'=Hd. In U, the cover with singleton root bin contributes Hd; all other root-starting covers distribute the final gadget across two ordered path bins, contributing Q(U−d). This proves U'=Q(U−d)+Hd. The inherited scalar row reduces to D'=HD+(Q−H)(U−d).

The invariant E=d therefore persists, V'=QV, and induction gives the claimed product formula. Later gadgets occur after earlier gadgets in the deletion ordinal order, so the equivalent single gadget is their chronological ordinal sum. H and Q multiply for such sums. No root-port ambiguity remains under that interpretation.

For a nondegenerate source interval k≥1, its image is an arithmetic progression with gap 2H_*. A consecutive-odd image therefore forces H_*=1, so all gadgets are transitive; Q_*=2^m and the output equals unselected suffix padding. This conclusion does not assert the uniqueness of graph witnesses or rule out interleaving different progressions.

**Wording caveat:** for k=0 a singleton is trivially a consecutive interval whatever H_*. Thus necessity of H_*=1, or ps=1, as a statement about actual window consecutiveness requires k≥1. Alternatively say that the formal difference coefficient remains two iff the multiplier is one. The same-root obstruction as an assertion about retaining a positive-width unit window is correct.

## Independent finite controls

`audit.py` uses new literal subset endpoint/start DP and explicit RuntimeError gates, not imported producer code or Python assertions. It verifies the actual 13-vertex base four-state, all 12 normalized k=0 contexts of padding length ≤2, and 357 literal grafts on every rooted tournament of order ≤4 satisfying E=d (119 source states with labels/roots counted). Gadgets are a singleton, transitive pair, and directed triangle. The triangle supplies H=3,Q=12 hostile controls. The broad rooted-source sweep tests the root-port rule beyond the particular star family. Exact rational endpoint inequalities and integrality are checked for k=0,...,100.

Run with both ordinary Python and −O; outputs must agree. These finite controls supplement the algebra, not establish its all-order quantifiers.
