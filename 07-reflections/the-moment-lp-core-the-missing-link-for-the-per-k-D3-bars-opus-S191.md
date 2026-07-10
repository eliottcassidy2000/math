---
source: opus-2026-07-09-S191
status: FORMALIZED the one-sided MOMENT-LP CORE (THM-661 analytic link), the missing piece between the
  rational D3 certs and the good-set measure floor. LRCMomentLP.lean (kernel-pure): integral_le_measure_pos
  -- on a probability space, measurable W, and any p with p(w)<=1_{w>0}, integral p(W) <= (mu{W>0}).toReal;
  measure_pos_ge_of_moment_ge -- b <= integral p(W) => b <= mu(GOOD). This REDUCES the per-k bars
  rhoGlobFloorRat(k) <= mu(GOOD) to (a) the rational moment bound bar <= Sum c_i m_i (LRCD3FloorCert,
  native_decide, done k=11), (b) the moment identity E[W^i] = integral W^i (Farey-cell integration), and
  (c) the min over ALL k-clusters (compact check + decorrelation tail = the coupled-region crux). I did
  NOT fully prove the bars -- I formalized the analytic INEQUALITY that connects the (certified rational)
  D3 to the measure, leaving (b) and (c) as the shared open analytic content.
tags:
  - lrc14
  - moment-lp
  - density-floor
  - d3-bars
  - formalization
  - markov-krein
---

# The moment-LP core: the missing link for the per-k D3 bars

**opus-2026-07-09-S191.** Owner: prove the per-k D3 bars for k=8..12. These are the density-floor bars
`rhoGlobFloorRat(k) <= mu(GOOD) = P(W>0)` (`W = Sum (gap_i - 1/7)_+`), the genuine analytic crux. I could
not prove them fully -- they are the coupled-region content the whole fleet is working -- but I formalized
the KEY MISSING LINK and reduced the bars to it plus two identified pieces.

## The state I found

- `GoodSetBridge` (death-star, AP certs): `mu(GOOD) >= m_P` for BOUNDED DIAMETER (<=75), kernel-pure, via
  the Farey-roof AP20/30/44/76 certificates.
- `LRCD3FloorCert` (kps-S89): computes the EXACT rational moments `m1,m2,m3` and `D3 = m1/M +
  (m1-m2/M)^2/(m2-m3/M)`, and `native_decide`s `bar <= D3(E)` for the k=11 anchor shapes -- the RATIONAL
  moment bound.
- The MISSING piece: the analytic inequality `mu(GOOD) >= D3(E)` connecting the rational `D3` to the actual
  Lebesgue measure. This is THM-661's proof core (`E[p(W)] <= P(W>0)` for a feasible `p`), and it was NOT
  formalized.

## Formalized (LRCMomentLP.lean, kernel-pure)

`integral_le_measure_pos` : on a probability space, for measurable `W` and any `p : R -> R` with `p(w) <=
1_{w>0}` pointwise, `integral p(W) dmu <= (mu {W>0}).toReal`. Proof: `p(W x) <= indicator{W>0} 1 x`
pointwise, then `integral_mono` and `integral_indicator_one`. Plus `measure_pos_ge_of_moment_ge` (chains a
rational bar `b <= integral p(W)` to `b <= mu(GOOD)`).

Instantiated with a feasible degree-`d` polynomial `p` (`p <= 1_{w>0}` on `[0,6/7]`) and `integral p(W) =
Sum c_i E[W^i]`, this is EXACTLY `D_d(E) <= mu(GOOD)` -- the rigorous reduction of the density floor to a
moment bound.

## What this reduces the per-k bars to

`rhoGlobFloorRat(k) <= mu(GOOD)` now follows from three pieces, of which the first two are done/routine and
the third is the shared crux:

1. **(a) rational moment bound** `bar <= Sum c_i m_i` (feasible `p`): LRCD3FloorCert `native_decide`, done
   for k=11; routine to extend to k=8..12 (same Farey-cell engine).
2. **(b) moment identity** `E[W^i] = integral W^i dmu` (`m_i` = the Farey-cell integral): the analytic
   moment computation -- connecting the exact rational `m_i` to the actual integral. Involved but standard
   (piecewise-polynomial integration over Farey cells).
3. **(c) min over ALL k-clusters** `min_E D3(E) >= rhoGlobFloorRat(k)`: compact check (bounded diameter,
   native_decide) + the DECORRELATION TAIL (diam > 75), whose coupled band diam in [18,35] is the fleet's
   known open analytic step (LEM-005 / the aliasing route THM-665/666).

So the per-k bars are NOT closed; the moment-LP core (this session) supplies the missing measure-from-moment
inequality, reducing them to (a) [done/routine], (b) [Farey-cell integration], (c) [the coupled region].

## Ledger

- FORMALIZED the one-sided moment-LP core `integral_le_measure_pos` + `measure_pos_ge_of_moment_ge`
  (LRCMomentLP.lean, kernel-pure, root-wired) -- the analytic link `Sum c_i E[W^i] <= mu(GOOD)` between the
  rational D3 certs (LRCD3FloorCert) and the good-set measure floor.
- HONEST: this does NOT prove the per-k bars; it reduces them to (a) rational moment bound [native_decide,
  done k=11], (b) the Farey-cell moment identity, (c) the min over clusters (compact + coupled-region tail,
  the shared open crux). -> THM-661, LRCD3FloorCert (kps-S89), GoodSetBridge (death-star), LEM-005, THM-665/666.
- File: LRCMomentLP.lean.
