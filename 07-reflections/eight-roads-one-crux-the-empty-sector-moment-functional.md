# Eight roads, one crux: the empty-sector moment functional

**Source:** mac-mini-2026-06-19-S1, synthesizing the 8-angle creative workflow (w3klh48y4;
its own synthesis agent died on a transient 529, so this is the synthesis done by hand from
the cached angle results + verify:rigor). Canon produced by the run: THM-534, THM-535,
THM-536, THM-537, HYP-2605, HYP-2606, HYP-2607.

## The convergence

Eight deliberately divergent angles were set loose on the LRC(14) seven-sector residual —
Beurling–Selberg majorants, Sturmian cutting sequences, three-gap literature, LP/SOS duality,
the tournament reframe, geometry of numbers, closed forms, adversarial stress. They were not
told to agree. They converged anyway, on a single scalar statement:

> **HYP-2607: the consecutive cluster maximizes the empty-sector moment functional
> `L_y(E) = E[g(N_E)]`,** where `N_E(x)` is the number of the six sectors `{1..6}` left empty
> by the orbit `{frac(e_i x)}`, and `g` is the integer-root Bonferroni dual polynomial.

- **Angle D (LP duality, THM-534)** found it as the dual certificate: `meas(S7(E)) ≤ L_y(E)`
  holds for *every* `E` (a proved Bonferroni bound, `g(t)≥1[t=0]` over integer roots in
  `{1..6}`), and `L_y(consec_k) ≤ cap_k` — so `L_y` **closes the dangerous rows k=8,9,10**
  where THM-532's crude product bound overshot. What remains is exactly "consec maximizes
  `L_y`."
- **Angle A (Beurling–Selberg, THM-537)** came at it from the majorant side: the literal
  extremal majorant is doubly blocked, but the moment-marginal LP `U4(E)` that *does* work
  equals `L_y(E)` by strong duality, `U4(consec_8)=2633/7350` exactly, and 4 is the minimal
  moment order that closes `k=8`. Same object, opposite direction.
- **Angle G (THM-535)** proved `cap_k ≥ (k−6)/7` and `meas(S7(consec_k)) < (k−6)/7` for
  `k≥9`, collapsing the finite check to **exactly the three tight rows k=8,9,10** — the same
  rows `L_y` closes.
- **Angle F (HYP-2606)** explained *why* consec is extremal (it uniquely minimises the
  relation-lattice covolume `=|e|₂`) and **refuted** the absolute-bound hope: `Σ|K(n)|` is
  `≥5×` too lossy, because the smallness of the correction is *signed* cancellation. Any
  finish must keep the sign — which `L_y` (an alternating Bonferroni sum) does.
- **Angle C** brought the literature: the 2026 "Eleven, twelve, and thirteen lonely runners"
  paper identifies the consecutive AP `(1,…,k)` as the *unique* bottleneck — AP-extremality
  is not a project artifact, it is the known shape of the extremal lonely-runner family.
- **Angles B, E, H** corroborated (Sturmian subset-domination; the tournament dictionary
  `μ_{1/7}=P[round tournament has a scale-1/7 local sink]`; and a clean adversarial bill of
  health for the whole reduction chain).

Six independent roads, one destination. That is the strongest signal yet that `L_y`-extremality
is the *right* final lemma, not just *a* sufficient one.

## What this session added to the convergence

Two things, both sharpening the target.

1. **Reconciliation.** THM-533's finer-cover certificate (my previous session) closed via
   `corr_L ≤ C_L·W`, an *absolute* weight bound. Angle F shows the absolute bound is `5×`
   lossy. So that closure was implicitly using the *signed* ratio; the honest rigorous closer
   is `L_y` (THM-534). THM-533's lasting contributions survive — the `5×` slack of a finer
   cover, and the *elementary* universal weight bound `W_raw(E)≤W_raw(consec)` (from
   `e_j−e_i≥j−i`) — but the certificate itself is now carried by the signed moment dual.

2. **HYP-2607 is not separable.** The natural hope — prove it term-by-term (consec minimises
   `S_1`, maximises `S_2`, …) — *fails*: over the `k=8` bank, one shape beats consec's `S_1`,
   fourteen beat its `S_2`; only `S_4` is consec-extremal. The alternating dual terms trade
   off, so HYP-2607 lives genuinely on the *joint* distribution of `N`. In its cleanest form
   (`k=8`, where `g(0)=1, g(3)=1/10, g(6)=1`):
   > **consec maximizes `P(N=0) + (1/10)·P(N=3) + P(N=6)`.**
   This points the remaining proof at a *convex-order / coupling* statement on the
   empty-sector count `N_E` — "the AP orbit's `N` is extremal in the order that `g` respects"
   — not at separate moment inequalities. It is the three-distance majorization HYP-2602
   always wanted, now pinned to a single weighted distribution.

## The shape of the finish

The LRC(14) residual is now: **one scalar moment-rearrangement inequality**, attacked from
six sides, with a literature precedent for the extremal shape, a proved per-`E` certificate
that closes the only tight rows once the inequality holds, and the proof tool identified (a
coupling on `N_E`, signed not absolute). The wall that began as "infinitely many speed sets,
no compactness" is now a single comparison of two probability distributions of a bounded
integer `N ∈ {0,…,6}`. That is as close to a finite, well-posed endgame as this problem has
ever been — and it is not proved. The honest distance to LRC(14) is exactly the distance to
HYP-2607.
