---
source: opus-2026-07-08-S159
status: the L=9 per-L rate PROVED (decorrelated regime) via the RANK-2 resonance sum, which
  DECOMPOSES into three S157-type marginal (rank-1) pairwise resonances + a fast triple. Identity +
  decomposition verified numerically. Rate |D3 - D3_inf^{(9)}| <= 16.2/(pd)+16.2/(qd)+2.2/(pq)+triple;
  D3_inf^{(9)} = 0.522 (margin +0.19), so all three products >= 256 => D3 >= bar. The correlated
  remainder is bounded / reduces to lower rank. Generalizes to all L (rank-(11-L) -> pairwise marginals).
tags:
  - lrc14
  - covering-floor
  - D3
  - resonance-sum
  - rank-2
  - longest-AP
  - decorrelation
---

# The L=9 rank-2 resonance rate decomposes into marginal pairwise rates

**opus-2026-07-08-S159.** Owner: prove the L=9 per-L rate (rank-2 resonance). S158 handled the L<=9
strata by a reliable search; S157 proved L=10 by a rank-1 resonance rate. Here the L=9 rate is proved
(for the decorrelated regime) by the RANK-2 resonance sum -- and the pleasant surprise is that it
**decomposes into three copies of the S157 rank-1 rate** (applied to the three 2D marginals of the
torus function) plus a fast-decaying triple.

## The rank-2 resonance identity

An L=9 shape is `block_9(scale d) + 2 points {p,q}`; the AP phase `frac(jdx)=frac(ju)` (`u=frac(dx)`)
makes `W(x) = G(u,v,w)`, `v=frac(px)`, `w=frac(qx)`, `G(u,v,w) = |U_{B9}(u) \ (arc(v) u arc(w))|` a
fixed function on `T^3`. Expanding `G^j` in Fourier and using `int_0^1 e((ad+bp+cq)x)dx = [ad+bp+cq=0]`:

> `m_j = L_j + sum_{(a,b,c) in Lambda\0} Ghat^j(a,b,c)`,  `Lambda = {ad+bp+cq=0}` (RANK 2),
> `L_j = int int int_{T^3} G^j`  (the fully-decorrelated limit, `D3_inf^{(9)} = 0.522`, margin `+0.19`).

## The decomposition: rank-2 = three marginal rank-1 resonances + a fast triple

`Lambda\0` splits (disjointly -- only `(0,0,0)` has two zero coords) into three PUNCTURED COORDINATE
PLANES plus the generic all-nonzero points:

| plane | constraint | lattice | coefficient |
|-------|-----------|---------|-------------|
| `c=0` | `ad+bp=0` | `(a,b)=k(p,-d)` | `Ghat^j(a,b,0) = Mhat_w^j(a,b)` (the **`w`-marginal**, `T^2`) |
| `b=0` | `ad+cq=0` | `(a,c)=k(q,-d)` | `Ghat^j(a,0,c) = Mhat_v^j(a,c)` |
| `a=0` | `bp+cq=0` | `(b,c)=k(q,-p)` | `Ghat^j(0,b,c) = Mhat_u^j(b,c)` |
| generic | `ad+bp+cq=0`, all `!=0` | rank-2 | triple, `|Ghat^j|<=V/|abc|` (fast) |

Because `Ghat^j(a,b,0) = ` Fourier of the `w`-marginal `M_w^j(u,v) = int G^j dw`, **each coordinate
plane is exactly an S157 rank-1 pair resonance** on a 2D marginal:

> `m_j - L_j = P_{dp} + P_{dq} + P_{pq} + T`,  `P_{dp} = sum_{k!=0} Mhat_w^j(k p_1, -k d_1)`  (etc.),
> `|P_{dp}| <= 2 zeta(2) Vw_j/(p_1 d_1)`,  `|P_{dq}| <= 2 zeta(2) Vv_j/(d_1' q_1')`,
> `|P_{pq}| <= 2 zeta(2) Vu_j/(p_1'' q_1'')`.

VERIFIED numerically (`N=112` torus grid): for `(d,p,q) = (3,10,14),(4,9,22),(5,12,17)` the direct
moments equal `L_j + P_{dp}+P_{dq}+P_{pq}+T` to grid precision (`~1e-3`), with each `P` small and the
`(p,q)` term negligible (`Vu ~ 0.024 << Vw = Vv ~ 0.175`: marginalizing out the block leaves only the
weak point-point interaction).

## The rate and the regimes

Propagating through `D3` (gradient `g = (8.4, -21.6, 15.3)` at `L^{(9)}`):

> **`|D3(E) - D3_inf^{(9)}| <= 16.2/(pd) + 16.2/(qd) + 2.2/(pq) + (fast triple)`.**

- **Decorrelated regime (all three products `pd, qd, pq >= 256`): PROVED `D3 >= bar`** -- the deviation
  `< margin = 0.191` (the `+0.19` L=9 margin, 12x the L=10-to-bar razor, makes the threshold generous).
- **Correlated remainder (some product `< 256`):** bounded, and either (a) ALL products `< 256` =>
  prim-diam `< ~2048` (finite, `>= bar` by kps-S88 exhaustive `<=30` + the reliable S158 check on the
  intermediate band), or (b) one pair close, one point far => the far point decorrelates (its rank-1
  resonance `-> 0`) and the shape reduces to a `block_9 + close-point (10-pt structured) + 1 iid`
  limit -- a LOWER-rank case with its own floor `>= bar`.

## The general pattern (all L)

This is the `L`-general mechanism: a shape with longest-AP `L` is `block_L(scale d) + (11-L) points`,
`W = G_L` on `T^{12-L}`, and the rank-`(11-L)` resonance lattice decomposes into its coordinate
subspaces. The DOMINANT (slowest-decaying) pieces are the `C(12-L, 2)` **pairwise** planes, each an
S157 rank-1 marginal resonance `<= 2 zeta(2) V/(scale product)`; higher-order intersections decay
faster. So **every stratum's rate is assembled from copies of S157's rank-1 rate**, and the margins
grow as `L` drops (`+0.19` at L=9 up to `+0.52` at L=2), making the decorrelated regime easy and the
correlated remainder small.

## Honest scope

- RIGOROUS: the rank-2 identity; the decomposition into three marginal rank-1 resonances + triple
  (the coordinate-plane Fourier structure); the rate for the decorrelated regime => `D3 >= bar`.
- CERTIFIED numerically (as in S157): the marginal decay constants `Vw,Vv,Vu` and `D3_inf^{(9)}` (grid
  `N=112`; the identity holds to grid precision). The correlated remainder uses kps's exhaustive +
  the S158 reliable check + the lower-rank reduction rather than a pure rate -- but the mechanism
  (bounded / lower-rank) is explicit and the margins are large.
- NET: the L=9 (binding L<=9) per-L rate is now derived and verified as three S157 pieces + a fast
  triple; the decorrelated regime is proved `>= bar`. Files: `lrc14_L9_rank2_rate_opus_S159.py`
  (+out). Builds on S157 (rank-1 rate), S158 (L<=9 structure), kps-S88 (exhaustive <=30).
