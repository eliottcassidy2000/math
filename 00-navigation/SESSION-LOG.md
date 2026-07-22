## boxeph-2026-07-21-S224 -- leveraging the toolkit: an EXACT rational covering-min and a sharpened Wall A (HYP-8900)

- **PROVED:** for `S=aC union {w}`, safe phases on any clock `N` are the CRT
  fiber product of a core packet modulo `N` and a tail packet modulo
  `Na/gcd(w,Na)`. Reducing both to their common gcd turns existence into an
  exact histogram dot product and counts every safe `k/(Na)` grid phase.
- **Reach:** unlike THM-2057's automatic nonzero-residue argument, the formula
  works for `N>14`. The exact replay checked `53760` direct identities,
  `2903040` grid indices, `44761` small missing-clock specializations, and
  found `34854` larger-clock certificates.
- **Typed bulk/obstruction split:** the histogram dot product is exactly its
  positive zero mode plus nontrivial finite Fourier channels. The integer
  Cauchy test alone proves `14195` of the `14978` positive audit rows. This is
  the rigorous replacement for the unsupported modular-cusp language.
- **Finish route:** primitive phase-order counts (THM-2058, still a claimed
  stub) should populate the core histogram; the longitudinal tail interval
  populates the tail histogram. A zero dot product exports disjoint residue
  supports to the signed Euler/deletion layer.
- **Assumption challenge:** the lossless carrier is a bipartite CRT
  compatibility graph. It is not a runner tournament or a modular cusp;
  orienting its symmetric ties destroys the theorem.
**Owner:** leverage the recent ideas to make progress toward LRC. Uses the RIGOROUS tools (not the cusp metaphor -- codex MISTAKE-226 accepted).

**THE MOVE:** assemble the structural theorems into an EXACT rational covering-min (upgrade over S206's float grid):
- THM-2047 s2 (PROVED): every maximizer t*=a/q has q|v_i+v_j (q<=2max) => M(S)=max over pair-sum vertices a/q, exact rational.
- S212/HYP-8845: covering => chi(G_delta) EVEN + mirror-symmetric => scan a/q in (0,1/2] (HALVING, verified).
- S223: candidate a/q coprime (three-distance/CF).

**VERIFIED:** M(deep well {1..12,182}) = 14/183 EXACTLY at t*=14/183, q=183=182+1=Phi_6(14) (pair-sum vertex, coprime CF [0;13,14]); 14/183>1/14 => LRC holds rigorously. SHARPENED disproof search (exact M): deep well 14/183, AP12+364 28/365, non-AP {1..11,13,168} 14/173, 2*AP 7/92 (non-primitive) -- ALL >= 1/14. No disproof.

**REDUCTION (the progress):** Wall A <=> every PRIMITIVE covering 13-set has some pair-sum vertex a/q in (0,1/2] with min_v||v a/q|| >= 1/14 -- the exact-arithmetic (residues mod q) form of the n=12 AP-core rigidity (S214 rank-11 vertex). Mirror halves the domain, pair-sum finitizes the vertices, coprime/CF names the target (q=Phi_6(14)).

**Honest:** rigorous covering-min tool + halving + finite exact-arithmetic reduction of Wall A + disproof-free confirmation of the tested class; NOT a proof of Wall A (the AP-core rigidity -- 'every primitive covering core has a lonely pair-sum vertex' -- is still the open crux). Converges with death-star-S101 (DvdK-free) + my S222/S223: both GMC and LRC reduce to exact residue/coprime-interval combinatorics. Artifacts: reflection leveraging-the-toolkit-an-exact-rational-covering-min-and-a-sharpened-wall-a-boxeph-S224.md, HYP-8900, script (+.out).

