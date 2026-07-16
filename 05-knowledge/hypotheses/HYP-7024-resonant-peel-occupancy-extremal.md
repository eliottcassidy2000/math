# HYP-7024 — the resonant-peel occupancy extremal

**Status:** CLAIMED / OPEN (codex-2026-07-16-S17).  The reduction is exact; the proposed
uniform inequality is proved only on bounded boxes so far.

For a six-offset core `E={0,e1,...,e5}`, let `N(x)` be the number of missed inner
sectors among `{1,...,6}`, and write `p_j=meas{N=j}`.  At the binding owner-resonance
class `a = 1 mod 7`, `THM-891` reduces the signed two-scale coefficient exactly to

`C_1(E) = -(6 p_1 + 2 p_2)/49 = -2(3p_1+p_2)/49`.

The sharp candidate is

`3p_1+p_2 <= 8/7`,

with equality at `E={0,1,2,3,4,6}` (`p_1=1/4`, `p_2=11/28`).  Equivalently, if `C`
counts colliding runner pairs, the target is `3 Pr(C=0)+Pr(C=1)<=8/7`.  A tempting
first-moment shortcut is false: pair-collision mass is not uniformly `1/7`
(`(a,b)=(1,8)` gives `1/4`, `(2,9)` gives `2/9`), so the proof must retain arithmetic
relations.  Exact scans through diameter 20 find no larger coefficient in any of the
six nonzero resonance classes.  A proof would replace a factor-17.6
constant-propagation loss by the actual cross-section cancellation and could collapse
much of the remaining bounded-ratio band.

Missing: a global moment/rearrangement certificate for the collision process and a
proof that residue `1` is the worst of the six classes for every core.  The computation
reserved for this claim is
`04-computation/lrc14_resonant_cross_section_cancellation_codex_S17.py` with output
`05-knowledge/results/lrc14_resonant_cross_section_cancellation_codex_S17.out`.
