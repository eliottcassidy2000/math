# The covering-bound cap is a pair-normalized Pascal mass — and a web of connections to push on

*mac-mini-2026-06-27-S63. Owner: merge "pair-normalized Pascal mass" if applicable, then generate my own
web of out-of-the-box connections as I push the frontier. The merge lands exactly on the live cap thread
(codex HYP-3090) and my gK8 binding (HYP-3085). Then a curated web — each spoke with a real hook, marked
proof-relevant vs exploratory.*

## The merge (VERIFIED, `lrc_pair_pascal_cap_margin_macmini_S63.py`)
The covering-bound cap is a **pair-normalized Pascal mass**:
```
cap_k = C(k+1,2)/C(14,2) = C(k+1,2)/91   EXACTLY for k = 10,11,12,13        (codex HYP-3090, verified)
      = P( a uniformly random pair {i,j} of the 14-clock both lie in a (k+1)-block )
      = the SECOND factorial moment of the block-occupancy indicator.
```
and the **covering-bound margin is the pair-complement**:
```
1 - cap_k = (C(14,2) - C(k+1,2))/91 = (# pairs OUTSIDE the (k+1)-block)/91     [exact for k>=10]
          = 36/91, 25/91, 13/91, 0   for k = 10,11,12,13.
```
The two binding rows sit a small **dip** below the pure pair-Pascal mass: `dip_9 = 1/4004 = 1/(44·91)`,
`dip_8 = 1081/76440 ≈ 0.0141` (the largest — the hardest row). **The dip is the higher-Pascal (Krawtchouk
degree ≥ 3) tightening** — exactly my gK8 `L_yK8 = ... + 10 S2 − 9 S3 + 6 S4` correction (S60), where k=8 was
the "S2/S3/S4-balance" row. So:

> **`C(k+1,2)/91` = the degree-2 (pairwise) value of the cap; it is EXACT once the config is dense (k≥10),
> and the binding obligation at the sparse rows k=8,9 is precisely the finite higher-Pascal dip.**

This re-reads the whole gK8 thread: the covering bound is "pairwise occupancy ≤ pairwise capacity" (a
pair-Pascal inequality), and the only non-pairwise content is the finite k=8,9 dip. The Krawtchouk
decomposition (HYP-2716: `W_j = Σ_{|R|=j} d_R`, binary Krawtchouk on the 6-cube) is the exact organ: j=2 is
the pair-Pascal mass, j≥3 is the dip.

## The web (my out-of-the-box connections; ⊕ = proof-relevant, ○ = exploratory)

### A. Combinatorial-geometric (the 91 pairs as a space)
- **⊕ Johnson scheme J(14,2) / Eberlein polynomials.** The 91 pairs carry the Johnson association scheme;
  `cap_k` is a **degree-2 Delsarte/Eberlein value** and the dip is higher-degree. HOOK: re-run the gK8 LP in
  the Eberlein (not Krawtchouk-on-7) basis — the pairs are the natural ground set, so the LP may *diagonalize*
  and the dip may get a closed Eberlein form. (Ties HYP-2716, HYP-3085.)
- **○ Hypersimplex Δ(2,14) / Gr(2,14).** The 91 pairs = vertices of the 2nd hypersimplex; `cap_k` = the
  fraction in a `(k+1)`-face. HOOK: cap as a normalized face-volume; positivity/cluster structure of the
  pair-space; the dip = a triangulation defect of Δ(2,14).

### B. Probabilistic / moment (Pascal as a distribution)
- **⊕ Truncated moment problem / Hankel positivity.** The factorial moments `S0..S4` are a moment sequence;
  the cap is the **degree-2 truncation** of the moment problem, exact when the degree-2 relaxation is tight
  (k≥10) and dipping when degree-4 tightens (k=8,9). HOOK: prove the dip = the degree-2→degree-4 Hankel gap,
  giving the binding rows a **finite Hankel-determinant certificate** (the scissors form of CRUX 1).
- **○ de Moivre–Laplace / Edgeworth (Pascal → Gaussian).** The miss-count `N` is asymptotically Gaussian; the
  dip is the **skew/kurtosis (S3/S4) non-Gaussian correction**. HOOK: k=8,9 = the most non-Gaussian (sparsest)
  rows; an Edgeworth bound on the dip.
- **○ Pair correlation / Montgomery–GUE.** `S2` = the orbit's **pair correlation**; "pair-normalized mass" =
  pair-correlation statistic. HOOK: does the binding sector pair-correlation match a universal (Poisson/GUE)
  form, explaining why k≥10 is exactly pairwise?

### C. Number-theoretic / arithmetic (the denominators)
- **⊕ The mod-41 Farey / Dehn scale (S62).** Last session's equidecomposability invariant `D=41` is the
  *reassembly* (Dehn) face; the pair-Pascal cap is the *occupancy* (volume) face. HOOK: both are pairwise —
  unify them as the two pair-invariants (volume `meas`/cap vs Dehn `D`) of the lonely set's fiber (HYP-3091).
- **○ Bloch group / dilogarithm.** Scissors congruence → the Bloch group → the dilogarithm `L2`, which is the
  *pair* of the log (`Li_2`, the weight-2 polylog). HOOK: is the pair-Pascal cap a "dilogarithm volume" and
  the Dehn `D` its Bloch companion? (Dupont–Sah scissors–K-theory; weight 2 = pairs.)
- **○ Apéry / ζ(2), ζ(3).** The witness floor `m_P = 14249/252252` and the cap denominators (`91 = C(14,2)`,
  `4004 = 44·91`) are central-binomial-flavored; Apéry's `ζ(3)` proof runs on `C(n,k)^2` (pair-Pascal) sums.
  HOOK: do the dip rationals have an Apéry/`ζ` denominator structure?

### D. Topological / spectral (the lonely set as an object)
- **○ Persistent homology of the lonely arcs.** The arcs as a barcode over the level `1/14`; the witness floor
  = the longest `H0` bar. HOOK: the cap = the death-time distribution; the dip = a higher-persistence feature.
- **○ Beurling–Selberg extremal functions.** The Delsarte/gK8 dual IS a one-sided Fourier (Beurling–Selberg)
  extremal of the sector indicator. HOOK: the cap = the degree-2 Beurling–Selberg `L1` value; the extremal
  majorant gives the dip directly.
- **○ Quasicrystal / cut-and-project diffraction.** The lonely set is a model set; its **diffraction
  (autocorrelation) is intrinsically pairwise** = the `S2`/pair-Pascal mass. HOOK: cap = the central
  diffraction-peak ratio of the sector point set.

## What to push first (the frontier move)
The single highest-value spoke is **B-Hankel × A-Eberlein**: prove the **dip = the degree-2→degree-4 gap of
the pairs' (Johnson/Eberlein) moment problem**, i.e. give the k=8,9 binding rows a finite Hankel/Eberlein
certificate. That would close the only non-pairwise content of the covering bound: for k≥10 the cap is the
pure pair-Pascal mass `C(k+1,2)/91` (clean, exact), and for k=8,9 the dip is a finite degree-4 determinant.
Combined with the closed-form margin `(91 − C(k+1,2))/91`, the covering bound's analytic core becomes
**one pair-Pascal inequality + two finite higher-moment dips.**

Honest status: the merge (cap = pair-Pascal mass, margin = pair-complement, dip = higher-Pascal) is VERIFIED
arithmetic and a clean reframing; the Hankel/Eberlein closure of the dip is the proposed (open) push. The web
is a map of leads, marked by proof-relevance — logged to the backlog.

Related: HYP-3092 (the verified merge + margin), HYP-3090 (codex triangular caps), HYP-3085 (gK8 = the S2
pairwise binding), HYP-2716 (Krawtchouk shadow), HYP-3091 (the lonely set's fiber; cap = the volume face),
[[the-four-faces-of-14-why-the-exceptional-structures-crowd-into-lrc]], [[three-notions-of-sameness-are-the-lonely-sets-fiber]].
