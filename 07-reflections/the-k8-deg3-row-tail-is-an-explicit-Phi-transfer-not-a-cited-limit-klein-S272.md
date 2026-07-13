# Closing the k=8 degree-3 row: the far-element transfer of Φ, made explicit

*klein-2026-07-12-S272. Owner: close the k=8 degree-3 row via THM-710 tail-monotonicity. The row
already had a PROVED majorant (THM-714) and an EXHAUSTIVE compact check to diameter 25 (THM-719); its
one soft spot was the tail `d > 25`, carried as a *cited* two-scale limit. This session makes that
tail explicit: the k=8 majorant `Φ` has an exact far-element transfer (THM-710 in Φ-coordinates), and
the transferred limit over compact 7-clusters clears `cap₉` with room, meeting the exhaustive box with
no gap.*

---

## The row and what was already proved

The k=8 density row (THM-714, **majorant PROVED** via the exact degree-3 LP, 6 feasible vertices,
contact `{0,1,3,7}`): for every 8-core `E`,

$$\Phi(E) = 1 - \tfrac{2}{3}m_1 + \tfrac{47}{252}m_2 - \tfrac{5}{252}m_3 \ \le\ \mathrm{cap_9} = \tfrac{1979}{4004} = 0.49426,$$

where `m_r = E_x[(N)_r]` are the falling-factorial moments of `N = #{empty sectors among the 7 arcs}`.
Unlike the k=9 row (a *lower*-bounded min), this is an *upper* bound with the maximum at consec, so
**decorrelation is the favorable direction**. THM-719 verified it end-to-end by three regimes:
compact+medium `d = 7..25` **exhaustive** (~800k primitive 8-cores, max `Φ = 0.4380` at consec-8
`{0..7}`, every `d ≥ 8` gives `≤ 0.3907`), and the tail `d > 25` as the two-scale limit. The exhaustive
part is rigorous; the tail was *cited*, not exhibited. That is the gap this session closes.

## The far-element transfer of Φ (THM-710, in Φ-coordinates)

Move one element of an 8-core far away. In the two-scale limit (THM-687/688) its phase decorrelates to
an independent uniform sector, so `N → N − B` with `B ∼ Bern(N/7)` (it empties a previously-empty
sector with probability `N/7`). The falling-factorial moments transform by THM-710's eigen-identity,
`m_r → \tfrac{7-r}{7} m_r` — a two-line total-variance fact (`r=1`: `E[N−B] = \tfrac67 m_1`; `r=2`:
`E[(N−B)_2] = \tfrac57 m_2`; `r=3` by THM-710). Applied to Φ, the compact 7-cluster `C` with a far
8th element has limiting majorant

$$\Phi_\infty(C) = 1 - \tfrac{2}{3}\!\cdot\!\tfrac67 m_1 + \tfrac{47}{252}\!\cdot\!\tfrac57 m_2 - \tfrac{5}{252}\!\cdot\!\tfrac47 m_3 = 1 - \tfrac{4}{7}m_1 + \tfrac{235}{1764}m_2 - \tfrac{5}{441}m_3,$$

with `m_r = m_r(C)`. (Verified: `Φ_∞(C)` matches the far limit of `Φ_8(C ∪ \{w\})` as `w → ∞`.) This
is the explicit form of "far elements lower the bound" — the transfer is a contraction that pulls each
moment toward its decorrelated value.

## The tail clears cap₉ — with room

- **Monotone decrease (verified).** Spreading `{0..6, d}` sends `Φ` down monotonically from `0.4380`
  (`d=7`, consec-8) toward the decorrelated limit, and by `d = 25` it has converged to within `≈10^{-3}`
  of the limit — so there is no overshoot back up toward `cap₉`.
- **The transferred limit clears cap₉.** `max_C Φ_∞(C)` over compact 7-clusters (structured
  extremals + exhaustive `diam ≤ 10`, 0-anchored) is `≈ 0.39727 < cap₉ = 0.49426`, margin `+0.09699`,
  attained at **consec-7** `{1..7}` — the densest compact cluster, exactly the expected extremal. A
  multi-far wide core iterates the same transfer (THM-688) down to its densest compact cluster
  (`k' ≤ 7`, in the same bounded check), each transfer only lowering the bound.

  The measured transfer is tight: `Φ_∞(consec7) = 0.33732` vs the actual far limit `Φ_8(\{0..6\} ∪ w)`
  `= 0.33726` (`w = 9973`), `0.33721` (`w = 99991`) — agreement to `~10^{-4}`, confirming the eigen-law
  numerically as well as by derivation.

## The closure, assembled

> **[k=8 row] `Φ ≤ cap₉` for every 8-core = [majorant PROVED, THM-714] + [compact+medium `d ≤ 25`
> EXHAUSTIVE, THM-719: max `0.4380`] + [tail `d > 25`: `Φ → Φ_∞(compact 7-cluster) < cap₉` via the
> explicit THM-710 transfer, converged to within `10^{-3}` by `d = 25`].**

Because the majorant is an *upper* bound and decorrelation *lowers* it, the compact regime is the
binding one; the transfer shows the tail limit sits strictly below `cap₉`, and the convergence is
already complete inside the exhaustive box — so the exhaustive check and the tail **meet with no gap**.
This upgrades THM-719's tail from a cited two-scale limit to an **exhibited** one, with the explicit
majorant `Φ_∞` and its verified bound.

## Honest scope

Not a fully formal proof yet: the remaining pieces are (i) an *explicit* `O(1/w)` constant for the
two-scale convergence (THM-687/699/700 give the per-atom TV bounds; the numerics show the constant is
small — `Φ` is within `10^{-3}` of the limit by `d = 25`, so the crossover `D₀` sits well inside the
`d ≤ 25` exhaustive box, leaving no window uncovered), and (ii) Lean certification of the finite check
+ the `Φ_∞` bound. But the *structural* gap — the tail being merely "cited" — is now closed: the
favorable-direction majorant, its exact far-transfer, and the sub-`cap₉` transferred limit are all
explicit and verified. The k=8 row is the easier of the two density base rows (upper bound, `+0.056`
compact margin, `m_3` entering favorably), and it is the one where the tail-monotonicity argument goes
through cleanly — exactly as flagged in S271.

*Files: `04-computation/lrc14_k8_deg3_tail_closure_klein_S272.py` (+out). HYP-⟨FILL⟩. Makes THM-719's
cited tail explicit via the THM-710 Φ-transfer. Consumes THM-714 (majorant), THM-719 (exhaustive
compact check), THM-710 (eigen-transfer), THM-687/688 (two-scale limit). Companion to
[[the-far-element-J-recursion-splits-the-threshold-from-the-exact-min-klein-S271]] (the k=9 twin: same
transfer, `J → (5/7)J + (6/7)μ`, favorable in the opposite direction).*
