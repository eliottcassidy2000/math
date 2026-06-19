# The LRC residual is a round-tournament event

*mac-mini-2026-06-18-S7 (Angle E — tournament/parity home turf)*

This project is, at root, about tournaments and Rédei parity. The LRC(14) work drifted
far into circle geometry and singular series. This reflection records the moment the two
came back into contact: **the open LRC residual is, literally, an event about a round
tournament that lives on the cluster orbit.**

## The bridge

Fix the cluster offset set `E` (`0∈E`, `|E|=k`). For each phase `x∈[0,1)` the `k` points
`p_i(x)=frac(e_i x)` sit on the circle `R/Z`. The **difference-winding map**
`i→j ⟺ frac((e_i−e_j)x)∈(0,1/2)` is a tournament `T(x)` (HYP-2576). When tie-free it is a
**round (local) tournament**: order the points clockwise, and each beats exactly the points
in its open clockwise semicircle. Score `s_i = #points in the clockwise half-circle after
`p_i`. This is the canonical "points on a circle" tournament — the most studied family of
non-transitive tournaments.

## The dictionary (all exact, 0 failures, k up to 10)

- **`maxgap(x) > 1/2 ⟺ T(x) has a Condorcet winner ⟺ T(x) is not strongly connected.`**
  An empty arc longer than a semicircle means one point sees everyone in its semicircle.
  The circle-gap statistic and the tournament's strong-connectivity are the *same* fact.

- **`maxgap(x) > 1/7 ⟺ some point's clockwise 1/7-arc is empty` — a "scale-1/7 local sink".**
  This is the LRC-relevant one. The teeth-half-width-1/14 safe-point condition, the
  1/7-gap, the seven-sector miss — all are the existence of a **local dominator at scale
  1/7** in the winding tournament. Therefore the crux

  > `μ_{1/7}(E) = P_x[ T(x) has a scale-1/7 local sink ]`

  and the target `μ_{1/7}(E) ≥ thr_k` reads: *a positive density of phases must produce a
  winding tournament with a scale-1/7 local dominator.* The lonely runner is, at the
  optimal phase, a **Condorcet-like winner at scale 1/7** in the difference-winding order.

- **`c3(T(x)) = C(k,3) − Σ_i C(s_i,2)`** — the standard round-tournament 3-cycle identity,
  so cycle content is a pure function of the score (gap) vector.

## What the parity side says about extremality

The folklore of this project is "AP = transitive = extremal". The tournament side
**corrects** this. Over the winding ensemble, the AP `{0,…,k−1}` is not the *most
transitive* — it is the *most cyclic*: it **maximizes the average Hamiltonian-path count**
`E_x[H(T(x))]` (exact global max over the primitive box at k=7). The dangerous LRC cluster
is the one whose winding tournament carries the **most directed odd cycles on average**, not
the fewest. Sector-fill `meas(S7)` co-moves with this cyclicity — though, honestly, only as
a strong trend (36% discordant pairs), not an identity, which is exactly why AP-extremality
of `meas(S7)` is delicate (and fails at k=12,13, per HYP-2604). The parity invariant `H`
(Rédei: always odd) and the sector measure are *correlated by the same circle geometry* but
are not the same functional.

## Why this is the right home, even though it didn't close LRC(14)

Two independent frameworks converged on the same object. The singular-series people found
that "relation-rich / low-height / AP-like" clusters are the dangerous ones. The tournament
side says the dangerous clusters are the **maximally cyclic winders**. The 7-vanishing
(THM-503, `hat a_T(7m)=0`) is the statement that sector-occupancy is blind to the
7-divisible winding modes — the modes that would build a 7-fold-symmetric (odd) rotation.
The "7" in LRC(14) and the "7" in the sector cover are the *same* 7, and on the tournament
side it is the period of the rotation that the round tournament cannot see.

When the circle geometry of LRC and the parity geometry of tournaments name the same
extremal object from opposite directions, that is not a coincidence — it is the project's
recurring signal that the residual has a single underlying shape. The shape here is a
**round tournament with a scale-1/7 local dominator**. LRC(14) is not proved, but the
crux now sits squarely on home turf: prove that *enough phases force a 1/7-scale Condorcet
winner in the difference-winding tournament of any 13-speed cluster.*
