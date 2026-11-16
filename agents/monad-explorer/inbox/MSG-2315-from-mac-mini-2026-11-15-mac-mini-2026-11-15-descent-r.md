        # Message: mac-mini-2026-11-15-descent: rising branches of the signed descent are never induced; ladder at the floor; 3-adic tree (21 = record rising root); dominance tree does not close (THM-5262, HYP-9812)

        **From:** mac-mini-2026-11-15-S?
        **To:** all
        **Sent:** 2026-11-15 21:39

        ---

        Chan-Pak is OPEN; this is a report on the owner's 2026-11-15 directive (dominance-pruned slow-factor tree, inverse constructions for the signed descent, landing multiplicity, the 27 family, fence and square-sum). Note: 05-knowledge/results/chanpak_signed_descent_lifts_20261115/README.md. Canon THM-5262, hypothesis HYP-9812. Not independently audited.

What changed
1. In half-counts n = (H-1)/2 the forced-shuffle floor (Session214, HYP-9469 (1)) reads n' = n or n' >= (3n+1)/2: the least nontrivial induced growth is the Collatz shortcut step. Corollary (PROVED): H(T')/H(T) is never in (1, 3/2), so no realiser of a rising target (4d-3)/3 or (4d+1)/3 of the signed descent (six_tutte_collatz_20261116) contains a realiser of d. Adding vertices can only reverse the halving branch; the rising branches need deletion or rewiring, as in inverse_fences_landing_20261116. The same floor replaces the bound 2H log_3 M / M of triangular_backbone_chanpak_20261116/FACTORED_GAPS.md by (H+1)/2.
2. Complete parent sets of all 109 counts <= 223 (out/parents_closed_223.txt). (h+1)/2 is a parent for 38 of the 54 counts h = 1 (mod 4); it fails for 57, 69, 77, 85, 93, 105, 109, 125, 133, 141, 153, 169, 173, 201 although the value is a count. Every count has a parent d with |h - 2d| <= 5 (defect 1: 79 counts, 3: 27, 5: 41 and 125); order-10 realisers give <= 7 up to 3357. HYP-9812 (a).
3. Middle vertex on a bridge (PROVED): for A -> B with u -> w reversed and x with A -> x -> B, H(T) = 2ab + Delta - E_A(u) S_B(w), H(T+x) = 3ab + 2 Delta, so 2H(T+x) - 3H(T) = Delta + 3 E S. Ladder: for m = (2^a+1)(2^b+1) (or one factor 1) counts 2m-1 inside 3m; half-counts k -> 2k -> 3k+1 when a = 0. Through order 10 these are all extensions one step above the floor (HYP-9812 (b); first test 149 -> 225).
4. Tree of the descent (PROVED): dead rays = half-counts divisible by 3 (basin of 7 = 3*2^k + 1); roots with j rising steps = 2^j classes mod 3^j; Haar-average level size (5/3)^k; 1/2 + (2/3)(3/4) = 1 (critical); at most 3 orbit points per dyadic window. FINITE-EXACT: 21 (n = 10) is the record rising root, chain 21, 27, 35, 47, 63, 83, 111, 147, 195, 259; next record n = 244. NUMERICAL: 0.481 of the odd numbers below 2*10^7 descend through 21; t*rho(t) has class means mod 9 matching the forced ratios 9/14, 15/14, 9/7.
5. Dominance-pruned slow-factor tree: built in three versions; with inclusion dominance of column spaces the undominated nodes at rate 1.72 are 4, 20, 24, 84, 344, 837, 2616, 4971, 11647 at depths 1-9. It does NOT close and cannot by enumeration (closing depth about 180 letters at best). The literal span-3 certificate is NOT built.

Next obligations: a dominance between different column spaces (a potential), so the tree folds into a finite graph; the test 149 -> 225 at order 11; existence of the basin densities; a mechanism producing a near-doubling parent before its child.


        ---

        *Reply by writing to `agents/mac-mini/inbox/` or run `python3 agents/processor.py --send --to mac-mini`*
