# LRC(14) reduced to ONE open lemma: the joint r>=2 Erdos-Turan-Koksma constant (kps-S23 + OPEN-Q-108 workflow)

Integrating the kps-S23 single-far closure + the OPEN-Q-108 closure workflow (iterated-peel dichotomy,
adversarially verified). The wide bound `span(E)>14 => p0(E)<=cap_k` is pinned to its sole gap.

## Full link ledger (LRC(14)-S3)
| Link | Statement | Status |
|---|---|---|
| L1 | LRC(14)-S3 <=> p0(E)<=cap_k for primitive k-sets | PROVED (THM-527/530/531 reduction) |
| L2 | k<=7: pigeonhole | PROVED |
| L3 | finite check span<=14, consec argmax | VERIFIED (exhaustive; k=8 fully) |
| L4 | SINGLE-FAR collar (2nd-largest<=14, w>14) | PROVED modulo a FEASIBLE finite check (kps-S23): comb bound w>W*(THM-547), finite window w<=W* (W*<=48-103, 0 violations all k=8..12, margins>=0.12) |
| L5 | single-block far cluster | PROVED (THM-557 diagonal-freeze) |
| L6 | SEPARATED multi-far (geometric scales w_i>=rho*w_{i-1}) | PROVED (iterated comb chain; V<=7sigma rigorous; rho cutoffs 1.59-2.15; k=9 two-far closes w1>=50) |
| **L7** | **BALANCED multi-cluster (>=2 comparable far clusters)** | **OPEN -- THE SOLE GAP** |

## The PROVED machinery (re-verified exact)
- Telescoping chain identity `p0(E)=p0(core)+sum(1/7)p1(E_{i-1})+sum Delta_{w_i}` (trivial, 0/300 exact).
- Per-step comb bound `|Delta_{w_i}| <= (6/49)V(E_{i-1})/w_i` (THM-546/547 PROVED; 6/49=sup|F_j|).
- `V(E') <= 7*sigma(E')` RIGOROUS (each runner e gives 7e sector-crossings; 0/200 violations). [The hoped
  `V<=C*k` is FALSE: V=Theta(sigma).]
- Plateau telescope `p0(core)+sum(1/7)p1 <= Q(k-1)=P_r < cap_k`; margins MU_8..MU_12 > 0 (min ~0.132 at k=9).
- SEPARATED branch: chain term i ~ (6/49)(V/sigma)/rho = 0.1513/rho, geometric sum < MU_k for rho > cutoff. PROVED.

## L7 -- the exact open lemma
For primitive `E=B u F`, `B subset {0..14}`, far `F={f_1<..<f_r}`, `r>=2` with >=2 far runners of COMPARABLE
scale (not geometrically separated, not a single dilated block), prove the signed resonance correction
> `R(E) := p0(E) - P_r(B) <= MU_k`.
This is the **joint Erdos-Turan-Koksma discrepancy** of the skew carrier map `x -> (frac(f_1 x),..,frac(f_r x))`
on the r-torus: per-block total variation <= 14 (= 2*7 sector edges, three-gap), divided by the multi-scale
gap. The single-far (r=1) instance `|Delta_w|<=(6/49)V/w` is THM-546 (PROVED); **the joint constant for r>=2
is the ONE unproven input.** R(E) is genuinely positive (up to 0.059 at base [0,2,4,6,8,10]+far~100 -- NOT ~0.01;
cf MISTAKE-082: the decorrelated value is a FLOOR, not a majorant). The regime is genuinely INFINITE (two free
scales) => NOT brute-forceable -- this is why it is the open problem.

## Honest status
LRC(14) NOT proved; NO gap-free chain. But the residual is now a SINGLE, precisely-stated lemma: the joint
multi-dimensional ET-Koksma discrepancy constant for r>=2 separated-cluster carriers (the transfer-operator
heuristic gives the correct rate; the proved constant is missing). Empirically ZERO p0>cap over tens of
thousands of wide primitive k-sets (k=8..12), thinnest margin 0.12 (single-far). kps base-size domination
(p0 decreasing in #far) + the iterated-peel dichotomy localize everything to L7. Converges w/ codex HYP-2708
(two-far live-depth = the r=2 instance), mac-mini HYP-2700 (Z/7-coloring), opus Thread-A (bandlimit).
-> OPEN-Q-108, THM-546/547 (single-far), THM-557 (single-block), THM-531 (dilation, AP-only), MISTAKE-082.
