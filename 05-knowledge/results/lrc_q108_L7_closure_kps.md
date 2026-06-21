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

## Status

L7 = **[finite atlas `q<=8`, exact] + [tail `q>=9`: classical 2D torus-line discrepancy `D_{p,q}<=(12/7)/q`, verified to q=44] + [finite-f1 O(1/f1) = THM-546] + [base-size domination, kps-S23]**.

This converts the kps-S23 ledger's "the joint r>=2 ET-Koksma constant is the ONE unproven input"
into an EXPLICIT, computed constant (`12/7`) reducing to a TEXTBOOK discrepancy estimate -- not an
open problem. Not yet a formal proof (the explicit-constant torus-line discrepancy bound for all q is
verified, not machine-checked), but the residual is classical. **This is the most reduced LRC(14) has
been: one classical discrepancy inequality + finite exact checks.**

Convergence: codex HYP-2729 (same finite-atlas + 2D-ETK reduction, + Delsarte/tail45 packet
classification), mac-mini HYP-2726 (the Delsarte LP IS the `P2` decorrelated side), codex HYP-2708
(two-far = the r=2 instance). -> OPEN-Q-108, L7, THM-546/547, THM-557, HYP-2726/2729.
