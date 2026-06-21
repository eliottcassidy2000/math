# L7 closure: the sole open LRC(14) lemma reduces to a classical 2D torus-line discrepancy bound

**kind-pasteur-2026-06-21** (HYP-2730; converges with codex HYP-2729, mac-mini HYP-2726).
Builds on the kps-S23 ledger (`lrc_q108_full_ledger_kps-S23.md`): L1-L6 PROVED/feasible, **L7**
(balanced multi-cluster, the joint r>=2 Erdos-Turan-Koksma constant) the SOLE gap, narrowed to the
bounded ratio window `f2/f1 in (1,2.15)`.

## The reduction (r=2; r>=3 splits pairwise into r-1 such windows)

For `E = B u {f1,f2}`, `B subset {0..14}` the bounded base, balanced `gamma=f2/f1 in (1,2.15)`:

**Step 1 (limit).** As `f1 -> inf`, `p0(E) -> p0_inf(B,gamma)`, and the far pair
`(frac(f1 x),frac(f2 x))` converges in law to `(frac(q v),frac(p v))` for `gamma=p/q` (coprime),
`v` uniform and independent of the base phases (standard equidistribution; the base elements are
bounded). For irrational/large-`q` `gamma` the `(q,p)` curve fills the torus -> the far pair is
independent -> `p0_inf = P2(B)` = the decorrelated plateau (mac-mini's Delsarte-LP side).

**Step 2 (exact cell form).** With `sector(y)=floor(7 frac(y))`,
```
   p0_inf(B,p,q) = sum_{i,j in Z/7} mu_{p,q}(i,j) * g_B(i,j),     P2(B)=sum_{i,j}(1/49) g_B(i,j),
   mu_{p,q}(i,j) = v-measure{ (sector(qv),sector(pv))=(i,j) },    g_B(i,j)=x-measure{ base(x) u {i,j} = Z/7 },  0<=g_B<=1.
```
So the resonance correction is `R(p/q) = p0_inf - P2 = sum_{i,j}[mu_{p,q}(i,j)-1/49] g_B(i,j)`, hence
```
   |R(p/q)| <= max_{i,j} g_B(i,j) * D_{p,q}  <=  D_{p,q},     D_{p,q} := sum_{i,j}|mu_{p,q}(i,j)-1/49|.
```
`D_{p,q}` is the **L1 cell-discrepancy of the `(q,p)` torus geodesic** -- a classical object. Its 1D
marginals are EXACTLY uniform (`qv`,`pv` sweep full periods), so `D` is pure 2D correlation.

**Step 3 (finite atlas + tail).** Computed exactly (`lrc_q108_L7_{resonance_atlas,tail_discrepancy}_kps.py`):
- `sup_q (D_{p,q} * q) = 12/7 = 1.714` over the whole window (attained at `3/2`); `D_{p,q}=0` exactly
  whenever the apex prime 7 aligns the grid (e.g. all `q in {7,14}` ratios).
- The largest `q` with `D_{p,q} >= margin (0.21)` is **`q=4`**: the FINITE ATLAS is the 6 ratios
  `{2/1,3/2,4/3,5/3,5/4,7/5...}` with `q<=4` (and `q<=8` to be safe), each exact-checked:
  `p0_inf(B,p/q) <= 0.247 (k=9) / 0.089 (k=8) << cap_9=0.494 / cap_8=0.381`, margin `>= 0.247`.
- TAIL: **`D_{p,q} <= 14/p`, PROVED ELEMENTARILY** (no three-gap citation needed) -- see the proof
  below; verified 0 violations over 1248 ratios, true constant `sup(D*p)=20/7`. So `|R| <= 14/p`,
  tail `p >= 67` safe (`14/67<0.21`); atlas `p <= 66` exact-checked (practically `p<=13` via `20/7`).

### Elementary proof that `D_{p,q} <= 14/p`  (the last analytic gap, now closed)

`mu(i,j)=meas{v: {qv} in I_i, {pv} in I_j}`, `I_t=[t/7,(t+1)/7)`. Fix `i`: `{qv} in I_i` is `q`
intervals of length `1/(7q)`; on the `m`-th (`m=0..q-1`), `v=(i+7m+s)/(7q)`, `s in[0,1)`, and `{pv}`
sweeps the arc `[a_m, a_m+L)`, `L=p/(7q)`, `a_m={pi/(7q)+pm/q}`. **Since `gcd(p,q)=1`, `{pm/q mod 1:
m=0..q-1}={0,1/q,...,(q-1)/q}` -- the `a_m` are EXACTLY equally spaced, gap `1/q`.** With
`h_j(a)=overlap([a,a+L),I_j)` (trapezoid, `Var(h_j)=2/7` since `L>1/7`), `mu(i,j)=(q/p)(1/q)sum_m
h_j(a_m)`. Koksma on equally-spaced points (`D*<=1/q`): `|(1/q)sum h_j(a_m)-int h_j|<=Var(h_j)/q=2/(7q)`,
`int h_j=L/7=p/(49q)`. Hence `|mu(i,j)-1/49|=(q/p)|err|<=2/(7p)`, and summing the 49 cells
`D_{p,q}<=49*2/(7p)=14/p`. QED. Script: `lrc_q108_L7_discrepancy_proof_kps.py`.

So the tail is RIGOROUS and ELEMENTARY; only the finite atlas (`p<=66`, exact `p0_inf<cap`) remains a
finite computation. (Earlier `D_{p,q}<=(12/7)/q` was the empirical `D*q` constant; the proven per-`p`
form `14/p` is what closes the tail without citing classical discrepancy theory.)

**Step 4 (finite f1).** `p0(B u {f1,f2}) - p0_inf = O(1/f1)` exactly the single-far comb rate
`(6/49)V/f1` (THM-546 PROVED), verified `(p0-p0_inf)*f1` bounded (`<=0.72` at `gamma=2`). So for
`f1 > W*` the finite value is within margin of `p0_inf`; `f1 <= W*` is the existing finite window.

## Exhaustive atlas verification (k=8..12)

Per-`k` tail threshold `p* = 14/(cap_k - P2(B))` (proven tail `D<=14/p` makes `p>p*` safe):

| k | base | P2 (plateau) | cap_k | margin cap-P2 | p* | atlas p<=p* |
|---|---|---|---|---|---|---|
| 8 | [0,2..10] | 0.0839 | 2243/5880 | 0.298 | 47 | 0 violations (sup p0_inf 0.089) |
| 9 | [0,2..12] | 0.2465 | 1979/4004 | 0.248 | 57 | 0 violations (sup 0.247) |
| 10| [0,2..14] | 0.3989 | 55/91     | **0.205** | 68 | 0 violations (sup 0.401, checked p<=72) |
| 11| [0,2..16] | 0.4807 | 66/91     | 0.245 | 57 | 0 violations (sup 0.482) |
| 12| [0,2..18] | 0.5586 | 6/7       | 0.298 | 47 | 0 violations (sup 0.559) |

**Binding case is k=10 (margin 0.205); the margin DIPS at k=10 then RECOVERS** (not monotone -- no
blow-up at large k). All `p<=p*` exact-checked (0 violations); `p>p*` covered by the proven `D<=14/p`.

**CAP CORRECTION (flag for the team):** correct caps are `cap_11=66/91`, `cap_12=6/7`. Some repo files
carry `cap_11=25/91`, `cap_12=2243/5880` (=cap_8, a copy-paste) -- those are WRONG (below the plateau
`P2`: `P2(k=11)=0.481 > 25/91=0.275`). With the correct caps L7 is comfortable at every k.

## Status

L7 = **[finite atlas, exact] + [tail `p>=67`: elementary torus-line discrepancy
`D_{p,q}<=14/p`, proved above] + [finite-f1 O(1/f1) = THM-546] + [base-size domination,
kps-S23]**.  The sharper empirical constants `sup(D*q)=12/7` and `sup(D*p)=20/7` make the
practical atlas tiny, but the proved `14/p` bound already closes the analytic tail.

This converts the kps-S23 ledger's "the joint r>=2 ET-Koksma constant is the ONE unproven input"
into an elementary discrepancy estimate plus finite exact checks.  The tail is no longer the open
part; the remaining work is proof packaging of the finite atlas and finite-`f1` window. **This is
the most reduced LRC(14) has been: one elementary discrepancy lemma + finite exact checks.**

## FINAL: L7 computationally CLOSED (all binding k)

- **Atlas exhaustive, k=8..12, p<=p* (per-k), 0 violations** (the resonance peaks all comfortable;
  binding k=10, margin 0.205). Tail `p>p*`: `R<=14/p` PROVED elementarily.
- **Worst base validated**: for k=10 the dense even AP `[0,2..14]` is the max over 400+ bounded bases
  (span<=14) at the binding ratio `2/1` (`p0_inf=0.40061`). r=2 dominates r>=3 (base-size domination,
  verified). So the worst balanced multi-cluster config is pinned and comfortable.
- **Finite-f1**: `p0(finite) -> p0_inf` at `O(1/f1)` = THM-546 (PROVED), verified.

L7 (the sole open lemma of the kps-S23 LRC(14)-S3 ledger) is therefore CLOSED modulo (i) the standard
far-element equidistribution limit (the basis of THM-546) and (ii) the previously-established L1-L6
and the S1/S2 cases / end-to-end audit. The analytic content of L7 -- the joint discrepancy constant
-- is now an elementary `D<=14/p`. Not a claim that LRC(14) is fully proved (the chain needs an
end-to-end audit), but the last *analytic mystery* of the sector route is resolved.
Scripts: `lrc_q108_L7_{resonance_atlas,finite_convergence,tail_discrepancy,discrepancy_proof,r3_domination,atlas_exhaustive}_kps.py`.

Convergence: codex HYP-2729 (same finite-atlas + 2D-ETK reduction, + Delsarte/tail45 packet
classification), mac-mini HYP-2726 (the Delsarte LP IS the `P2` decorrelated side), codex HYP-2708
(two-far = the r=2 instance). -> OPEN-Q-108, L7, THM-546/547, THM-557, HYP-2726/2729.
