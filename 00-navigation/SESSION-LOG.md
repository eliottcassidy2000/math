## boxeph-2026-07-21-S223 -- one-dimensional coprime intervals complete the DvdK bypass (HYP-8895)
## codex-2026-07-21 -- THM-2059 exact CRT packet join

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

## boxeph-2026-07-21-S222 -- bypass GMC(2)'s DvdK dependency via the saddle-point/Watson method (HYP-8890)
## death-star-2026-07-21-S101 -- SHARP DvdK-free criterion (refines S100 a lot): a UNIQUE minimal balanced channel (unique tournament-zeta primitive cycle) => GMC(2) DvdK-free, coefficient-independently. 84% of supports. HYP-8878.

**Owner:** keep going on the DvdK bypass; think one-dimensional coprime intervals.

**KEY:** DvdK (THM-1630) is a ONE-VARIABLE theorem, and GMC(2) uses it only in 1 variable. In 1D the constant-term non-vanishing is elementary COPRIME-INTERVAL / NUMERICAL-SEMIGROUP combinatorics (no residues/Liouville). For f=sum c_k z^k, support S (0 not in S): Newton polytope = interval [min,max] (two-sided iff 0 interior); period d=gcd(gaps); the RETURN set R={m:0 is an m-fold sum} is a numerical semigroup.

**TWO POLES (verified, one_dimensional_coprime_intervals_return_semigroup_boxeph_S223.py):**
- bare coprime PAIR {-q,p}: PERIODIC, R=(p+q)Z (= THM-1840 seed m0=p+q, ONLY multiples, Frobenius=inf).
- FILLED coprime interval (endpoints+interior, gcd gaps=1): COFINITE, R = all m > Frobenius# ({-2,1,3}: all m>=3; {-2,2,3}: all m>=4 with a gap). Adding one interior exponent collapses period p+q -> 1.
- POSITIVE coeffs: CT(f^m)!=0 <=> m in R (verified); mixed signs: R minus sporadic cancellations, cofinite (saddle S222).

**RESULT:** DvdK's 1D content ('two-sided => CT!=0 for some m') is SHARPENED to 'all m > Frobenius#' -- exact, EFFECTIVE (Frobenius# = the open effective-DvdK m0), self-contained. Combined with S222 (mixed-sign saddle) + S208 (confluent cusp), the whole 1-variable GMC(2) non-vanishing is elementary; DvdK demoted from imported premise to the periodic special case.

**BONUS:** coprime intervals also drive LRC's three-distance -- the gaps of {kt mod 1} take 3 values by the CF of t (= coprime-interval structure); t*=14/183=[0;13,14] coprime partial quotients; danger arcs = coprime intervals. 1D coprime intervals = the shared engine of GMC constant-term returns AND LRC three-distance arc geometry.

**Honest:** coprime-pair periodicity + filled-interval cofiniteness + CT<=>R (positive) verified; mixed-sign uses the S222 saddle. Completes the DvdK-bypass 1D core (elementary+effective), not a full THM-2022 rewrite. Artifacts: reflection one-dimensional-coprime-intervals-complete-the-dvdk-bypass-boxeph-S223.md, HYP-8895, script (+.out).

