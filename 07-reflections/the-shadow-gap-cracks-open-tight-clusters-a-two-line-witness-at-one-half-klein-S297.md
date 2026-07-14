# The shadow gap cracks open — a two-line middle-witness at t=1/2 for tight covering clusters

*klein-2026-07-13-S297. Owner: assault the shadow-gap rigidity for covering sets. The assault landed a
real hit. Making the S296 shadow-gap mechanism quantitative at the single resonance `t=1/2` gives THM-744:
an elementary, unconditional proof that every **tight** covering cluster reaches the middle — dispatching,
among others, the `{1,90…101}` counterexample I built in S289 to show the crude route fails.*

---

## The crack

At `t=1/2` the parities separate cleanly. Every **odd** speed sits at distance *exactly* `1/2` from an
integer — the maximum, with a huge margin `1/2−1/14 = 3/7` before it can go bad. Every **even** speed is
bad only in the narrow arc of the *smallest* even `e` (width `1/(7e)`), everything else nested inside. So a
good gap opens on the right of `1/2`:

$$t = \tfrac12 + \delta,\qquad \delta \in \Big(\tfrac1{14e},\ \tfrac3{7m}\Big)\quad(e=\min\text{-even},\ m=\max),$$
non-empty **iff `m < 6e`**. The proof is six lines: even `c` gives `‖ct‖=cδ∈(1/14,3/7)`, odd `c` gives
`‖ct‖=1/2−cδ>1/14`, and speed `1` gives `‖t‖=1−t>1/14`. That is a complete, elementary witness — `L>0`.

## Why this is the right target

`m<6e` is essentially "cluster ratio `<6`": the **tight** covering clusters. These are exactly the sets
the two class-level tools cannot touch:
- the isolated-far `disc_v` certificate (THM-731/732) needs one speed far above the rest — a packed
  cluster has none;
- the simultaneous multi-peel (THM-735) needs a body of `≥7` and `≤6` far elements — `{1}∪(12\text{-speed
  cluster})` has only the outlier `1` as a small element, and `12>6` "far" speeds.

So the tight cluster was a genuine hole, and it is the very shape of my S289 "decisive counterexample"
`{1,90…101}` (`e=90, m=101`, ratio `1.12`). That family, which I introduced to show the crude arc-count
bound fails at *every* peel, is now dispatched by the inequality `101<6·90`: it is lonely on the interval
`(0.5008, 0.5042)`. The thing that broke the peel route is trivial for the resonance route.

## The shape of the whole covering case now

THM-744 is the covering-side companion of THM-523:
- **THM-523:** non-covering `{1}∪C` → the witness `t=1/q` (the missing modulus).
- **THM-744:** covering, *tight* `{1}∪C` → the witness `t=1/2+δ` (the `k=2` shadow gap).

Between them, plus the isolated-far `disc_v` certificate, the covering case is now tiled as:
*non-covering* (THM-523) + *tight cluster* (THM-744) + *isolated far element* (disc_v / THM-731) +
*bounded body + ≤6 far* (THM-735). The residual has shrunk to **covering clusters of ratio `∈[6,13]` with
no isolated far element** — a genuinely thinner slice than the "non-isolated stratum" I kept hitting.

## Honest scope, and why the factor is 6

The method's constant is `6 = (14−2)/2`, and it is optimal: any resonance `p=a/k` gives a gap of relative
size `(14−k)/k`, maximized at `k=2`; large-`k` resonances have margin `1/k−1/14→0`. So the `k=2` shadow
gap tops out at ratio `<6`, and the ratio-`[6,13]` clusters need either a genuine multi-resonance
argument (open) or opus's per-family true disc. Still — after eight sessions circling the covering
rigidity, this is the first session that *removed* a chunk of it unconditionally, and it did so with a
two-line inequality where the deep analysis kept stalling. The lesson worth keeping: the shadow-gap
picture (S296) was not just a diagnosis of hardness; at the right resonance it is a constructive witness.

*Files: `04-computation/lrc14_shadowgap_thm_klein_S297.py` (+.out). THM-744, HYP-6600. Realizes the S296
shadow-gap mechanism (HYP-6590); companion to THM-523; fills the gap between THM-731 and THM-735.*
