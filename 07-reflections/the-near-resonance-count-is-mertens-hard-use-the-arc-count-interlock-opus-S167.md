---
source: opus-2026-07-09-S167
status: (1) the near-resonance COUNT is roughly UNIFORM (~150) across longest-AP L -- so "dissociated
  => few resonances" is NOT a count phenomenon; the discriminator is the SIGNED cancellation (r_N),
  the MERTENS situation (absolute bound hopeless, cancellation essential, kps-S92 20x). (2) The clean
  closure is mac-mini's ARC-COUNT route, whose rho* input IS opus-S158's D3_inf^{(L)} -- DECREASING in
  L, hence LARGEST for dissociated (low L) with a +0.10..+0.36 a-priori margin, dropping below #arcs/V
  exactly at L~10 where klein's Dirichlet/LEM-012 near-AP branch takes over. The two branches INTERLOCK
  via D3_inf^{(L)} monotonicity. (3) Mertens = genuine structural analogy; Hadwiger = a covering-number
  analogy (+ a speculative theta=1/7 <-> Hadwiger-Nelson 7-coloring thread).
tags:
  - lrc14
  - good-period
  - dissociated
  - mertens
  - hadwiger
  - density-floor
  - arc-count
---

# The near-resonance count is Mertens-hard; use the arc-count interlock

**opus-2026-07-09-S167.** Owner: work the near-resonance count and the single dissociated-branch
inequality; consider Mertens and Hadwiger. Findings redirect the last mile onto a clean route.

## (1) The near-resonance COUNT does not discriminate — it is the SIGNED cancellation (Mertens)

The dissociated inequality is `r_N = |Corr_N|/(N(6/7)^k) < 1` at `N = ceil(7(k-1)/6)`,
`Corr_N = sum_{m!=0} Whatfreq(m) G_N(m/V)` (opus-S165), near-resonances the `m` with `||m/V||` small.
Computing the near-resonance count (small balanced relations `n`, `n_i` coprime-to-7, `||n.e/V|| <
1/(2N)`): it is **roughly UNIFORM (~150) across all `L`** — it does NOT track longest-AP. So the
kps-S92 framing "dissociated => few near-resonances" is not a COUNT phenomenon. What IS controlled is
the **signed** `Corr_N` (verified `r_N <= 0.85` uniformly, S165), via cancellation: `|Corr_N|` is
`~4-5x` below the sum of `|terms|` (e.g. AP: `|Corr|=0.49` vs `2.26`; dissociated: `0.15` vs `0.63`),
and kps-S92's a-priori absolute bound is `~20x` the target.

**This is exactly the MERTENS situation.** `sum_{n<=x}|mu(n)| = x` (trivial, uniform), but the SIGNED
`M(x) = sum mu(n)` has deep cancellation (`M(x) = o(x)`, conj. `O(x^{1/2+eps})`), and PROVING the
cancellation is the hard part. Here `sum |What|` diverges (opus-S154) / is `20x` too big, yet the
signed `Corr_N` cancels to `r_N < 1`. So the r_N route is a MERTENS-type problem: **no absolute bound
works; the cancellation is the theorem.** (The arc-Fourier `b(m)=(1-e(m/7))/(2 pi i m)` vanishes at
`7|m`, so the resonance product is supported on `n_i` coprime-to-7 — an inclusion-exclusion over
residues mod 7 = a MOBIUS structure; klein-S197's "x7 collapse" exploits exactly this.) The
no-cancellation extremal is the STRUCTURED case (the complete-residue AP = the tight `M=1/k` instance,
opus-S164) — the analog of Mertens' worst case being a conspiracy of zeta zeros.

## (2) The clean route: mac-mini's ARC-COUNT, fed by opus-S158, interlocking with LEM-012

Because the r_N route is Mertens-hard, the right closure is mac-mini-S61's route (c): a good period
EXISTS once `#arcs < rho* V` (`rho*` = good-set density, `#arcs` = blocked-period count), which needs
NO cancellation. Its two inputs are a-priori, and — the key point — **`rho* >= D3_inf^{(L)}`, opus-S158's
decorrelated density floor, which is DECREASING in `L`:**

| L | 2 | 4 | 6 | 8 | 9 | 10 |
|---|---|---|---|---|---|----|
| `rho* >= D3_inf^{(L)}` | .855 | .839 | .761 | .601 | .524 | .465 |
| margin over `#arcs/V <= 0.5` (S58) | +.36 | +.34 | +.26 | +.10 | +.02 | **−.04** |

So `rho*` is **LARGEST exactly for the dissociated (low-`L`) clusters** the arc-count handles — a big
`+0.1..+0.36` a-priori margin — and it drops below `#arcs/V` right at `L ~ 10`, **exactly where the
other branch (near-AP, klein-S196/197 LEM-012 Dirichlet clustering) takes over.** The two branches
**INTERLOCK via the monotonicity of `D3_inf^{(L)}`** (opus-S158): the density-floor result I proved
for the k=11 tail is precisely the a-priori `rho*` that closes the good-period dissociated branch. One
theorem, both legs.

> **Net: the dissociated branch closes a-priori via arc-count (`rho* >= D3_inf^{(L)} > #arcs/V` for low
> `L`), NOT via the Mertens-hard signed `r_N`; the crossover `L ~ 10` is exactly the LEM-012 boundary.**

## (3) The conjecture connections

- **MERTENS — genuine structural analogy.** The r_N route is a Mertens-type problem (absolute bound
  hopeless, signed cancellation essential, mod-7 Mobius, extremal = structured). This is *why* the
  fleet should not chase the r_N a-priori and should use arc-count. Useful: it names the difficulty and
  redirects. (A literal Mertens theorem is not invoked; the analogy is the shape of the hard direction.)
- **HADWIGER — a covering-number analogy.** The arc-count route IS a covering statement: "`#arcs`
  same-length arcs cannot cover the period grid when `#arcs < rho* V`" — a packing/covering count, the
  kind Hadwiger's covering (illumination) conjecture concerns. The load-bearing content is elementary
  (measure/pigeonhole), so Hadwiger is an analogy, not a tool. A speculative deeper thread: the
  `theta = 1/7` / 7-arc covering threshold echoes the Hadwiger–Nelson plane 7-coloring (project
  THM-418/419) — both pin "7" as a covering constant of a homogeneous space; worth a look but not
  load-bearing here.

## Ledger / next

- FINDING: the near-resonance count is uniform (~150), so the dissociated inequality is a Mertens-type
  SIGNED-cancellation problem (r_N route), NOT a count — redirect to arc-count.
- SYNTHESIS: mac-mini's arc-count route closes dissociated a-priori via `rho* >= D3_inf^{(L)}`
  (opus-S158, decreasing in L => large for low L), interlocking with LEM-012 at `L ~ 10`. The k=11
  density-floor theorem IS the good-period `rho*` — one shared result.
- NEXT: make `#arcs <= c(L)·spread` (mac-mini S58) and `rho* >= D3_inf^{(L)}` fully a-priori/Lean; the
  r_N route is deprioritized (Mertens-hard). Files: `lrc14_near_resonance_count_mertens_opus_S167`
  (+out). -> kps-S92 (r_N/near-res), mac-mini-S61 (arc-count c), klein-S196/197 (LEM-012), opus-S158
  (D3_inf^{(L)}), opus-S154 (L^1 divergence).
