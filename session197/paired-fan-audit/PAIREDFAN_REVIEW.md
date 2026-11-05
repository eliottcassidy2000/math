# Independent review: common infinite paired fan

Verdict: **PASS**, 2026-11-05. No mathematical correction required.

Reviewed exact frozen files under `/workspace/scratch/07d5532f2fce/session197-hybrid-constructions/paired-fan/`:

- PROOF.md SHA256 `2615ac3e894719c2aec9fb9289b5b4b99233ea128e71b94427f1eae5c1bf9134`
- core.py SHA256 `9c104d778f7c38eee81629a81f9f9e88c2f1197a06c15c6b2d8afd424071d8aa`
- verify.py SHA256 `bc3439a6911ca3ea376d3ed63971945074a6c65e14050c41424ccaf8b961ab83`

## Proof audit

1. The six-coordinate core profiles and all literal gate masks agree. Gate A adds a source, then a vertex beaten only by retained old vertex0; Gate B adds a vertex beaten exactly by1 and4. The final root in the fan remains the original vertex5, not either gate's new vertex. The two initial port pairs differ, as explicitly stated.
2. The positive chain condition is essential and correctly present. For m>=1, no-seam Hamilton paths end at b and contain the entire chain immediately after the core, so there are H(X) of them. Seam paths correspond bijectively to rooted R_a covers and independent allocations of the chain vertices into the two increasing runs. The corresponding R_b cover dichotomy gives R_b=H. At m=0 the claimed h contribution would need an endpoint correction; the theorem does not apply there.
3. Both gates' new vertices beat every chain vertex and b, and their exceptional arcs address only unchanged named core vertices. Thus the adjacency-level gate/fan commutation is exact under the displayed relabelings for every m, not only at the six checked values. The root/deletion outputs are correctly changed to41 or43 rather than held at13.
4. All three cores are strong; the fan's return through b->a preserves strongness. For the general source compiler even a weak S works: a reaches x and then S-a, and in the enlarged prefix it reaches x,v,u and then S. The positive fan supplies return paths.
5. The S196 equal-response theorem was read at the named published source. Here h=18*2^m+13, d=13, h-3d>=10, so its exact low-window support is indeed {h,2h-13,3h}. Neither literal child can occur in that older marked-only-at-b language. The present two-pair comparison does not contradict that theorem.
6. The general prefix identities were rederived independently. Without the reverse seam, H is g*H(G) minus the invalid boundary s*e_x(G). A seam path is an L_a cover of S and an R_x cover of G. For endpoint a, the only invalid concatenation occurs when S is the singleton, giving the delta correction. In an R_a cover the seam can only occur in the second path, immediately a,x; the remaining S-a vertices occupy two unrestricted path bins, giving B*R_x(G). Deletion leaves an ordinary ordered join. Substituting prefix profiles (1,1,1) and (3,1,3) yields exactly (3h+2s,3e+2delta,3R,3d).
7. The fixed-source factory exclusion is valid for arbitrary prefix order. R=54 forces h_G+2R_G=9 with positive odd h_G and R_G>=h_G, leaving only(1,4),(3,3). The H-values have upper bounds41 and42 before subtracting a nonnegative endpoint count, and cannot be43. This is correctly limited to that exact source and mark.
8. The infinite family proves two children from each common parent, but not a common returning marked interface or closure of the child family. The scope statements accurately preserve that distinction.

## Independent execution

Ran the supplied verifier normally and with `python3 -O`, storing outputs separately in this reviewer's scratch. The JSON outputs are byte-identical. Complete core profiles also agree with literal permutations and cuts inside the verifier, independently of the subset DP.

The six parent/children triples are:

- 49 ->149,151
- 85 ->257,259
- 157 ->473,475
- 301 ->905,907
- 589 ->1769,1771
- 1165 ->3497,3499

All285 labelled rooted source cases through order4 pass the universal prefix-triangle compiler, including its singleton endpoint correction. Both physical graph relabelings, all marked root profiles, exposed port pairs and strongness controls pass.

Outputs: `pairedfan-review-normal.json` and `pairedfan-review-optimized.json` beside this report. Hashes were rechecked after the runs and remained equal to the frozen values. No repository or worktree was changed.
