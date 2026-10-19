# Independent audit of the three-tree Gold refinement

19 October 2026. Verdict: PASS. No blocking mathematical or implementation issue was found.

The three-tree theorem, strict normalized margin greater than 1/200, deterministic selection rule, and balanced-first-comparison corollary are all valid for the stated ten-vertex full root core and arbitrary finite seven-thin exteriors. The theorem covers the entire 1,764-atom convex hull. It does not classify other cores.

## Independent reconstruction

The audit script `independent_three_tree_audit.py` imports no author code, profile arrays, or count tables. It independently generates profiles through nested order ideals and filters, obtaining 46 lower and 48 upper profiles and 1,764 compatible pairs. Set-relation transitive closure verifies every model's acyclicity, seven-thin bound, full core, and critical-root orientation. The model orders range from 10 to 17.

It enumerates all 11,971 core extensions, groups their answer cells by 3,546 boundary-word classes, and computes integer interleaving multiplicities with a direct gap transfer table. This method differs from the author's ideal dynamic programming on each enlarged model. Every total and joint cell matches the author exactly. The independently reconstructed full table has SHA-256

    2f83a7b75f377b5160c37b83b30f64da1ec3207531dcd575983a5d48da9c767d.

All seven safe-margin minima and all four positive-certificate minima match the author's receipt. The smallest certificate average is 715/114938, strictly greater than 1/200. All four residual failure patterns contain a complete certificate support. Normal and optimized audit receipts are identical. Integer matrix products are explicitly bounded against overflow.

## All-size transfer and cover proof

The full-core hypothesis supplies the lower/upper exterior decomposition. Fixing exterior linear-extension orders and deleting their forced initial/final segments preserves every core-word lift count in the associated bounded chain model. The resulting normalized mixture weights are atom extension counts divided by the ambient extension count, one weight per conditioned exterior-order pair. They are not uniform. Thus all joint comparison inequalities transfer simultaneously.

After transfer, failure of all three trees forces A to have a bad FT or FF path, B to have a bad TT path, and C to have a bad TT or TF path. The two pair inequalities exclude C's TT case; the two three-margin inequalities exclude C's TF case. Strict inequalities and boundary cases are handled correctly.

## Selection rule and corollary

The stated four-branch selection rule is valid on the entire hull. If C's TT margin is at most epsilon, the pair inequalities make both potentially unsafe A margins strictly greater than epsilon. Otherwise C's TT margin is safe. If C's TF margin is also safe, C works; if B's TT margin is safe, B works. In the remaining case the two three-margin inequalities force A to work. The seven universally safe margins complete each case.

Only the four stated event probabilities are needed: the two C test margins are exactly 1-p-q and 1-2p+q, and the B test margin is 1-r-s. No unknown mixture decomposition is required. The audit also exercises this rule on every atom and verifies all selected margins.

If h is the selected first question's heavier branch probability, its heavier child has probability at least h/2. The Gold margin therefore gives 3h/2 < 199/200, or h < 199/300. Both first-answer probabilities are consequently strictly between 101/300 and 199/300. The corollary follows as written.

One editorial correction was applied by the author before freezing checkpoint2: residual failure patterns contain the certificate supports; the first two certificate supports do not themselves include B. This does not change the mathematical proof.
