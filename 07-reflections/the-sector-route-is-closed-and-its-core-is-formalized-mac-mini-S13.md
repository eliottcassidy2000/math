# The LRC(14) sector route is closed, audited, and its core is formalized

**Session:** mac-mini-2026-06-21-S13 (overnight, rigor-first, concurrent with kps and codex).

This session crossed a threshold: the LRC(14) S3 sector route — the last open case after the
covering-set reduction — went from "verified at every link" to **closed, adversarially audited
ALL-PASS, with its core bound machine-checked in Lean.** Not a finished gap-free proof of LRC(14),
but the nearest the project has stood to one.

## The two lanes met

The cover bound `measS7(E) ≤ cap_k` (S3 ⟺ this) splits exactly as `measS7 = P2(B) + R`:

- **P2 — the decorrelated plateau** is the Delsarte/moment-LP main term. The per-shape bound
  `measS7(E) ≤ L_y(E)` is a Delsarte linear program whose dual is Krawtchouk-nonnegative (HYP-2726),
  and it is now **formalized in Lean, sorry-free** (`LRCFactorialAtom.lean`): `delsarte_bound_k8/k9/
  k11` and the dual-feasibility `gK8/gK9/gK11_dominates`, covering every binding row k=8..13, each
  depending only on `[propext, Quot.sound]`.
- **R — the resonance correction** is bounded by the L1 cell-discrepancy `D_{p,q}` of the `(q,p)`
  torus geodesic, and kps closed it elementarily: `D_{p,q} ≤ 14/p` by Koksma on equally-spaced
  points (gcd(p,q)=1 forces the q sub-arc starts `{pm/q} = {0,1/q,…}` equally spaced; `Var(h_j)=2/7`;
  `|μ(i,j)−1/49| ≤ 2/(7p)`; sum the 49 cells). I re-verified it independently: 0 violations,
  `sup(D·p)=20/7` in the window (24/7 over all ratios) — both far under 14, so the bound is robust.

`measS7 = P2 + R` with `P2` Delsarte-bounded and `|R| ≤ 14/p`, plus a finite exact atlas for the
small denominators and `THM-546` for finite `f1`, closes L7 — the sole gap — uniformly in the number
of clusters (`r≥3` splits pairwise). All seven links L1–L7 now stand.

## The audit found no gap

A four-thread adversarial audit came back **ALL PASS**: the caps `cap_8..cap_13` are exactly
`2243/5880, 1979/4004, 55/91, 66/91, 6/7, 1` (a stale `cap_11=25/91`/`cap_12=cap_8` copy-paste lives
only in non-load-bearing scratch files); over all **11432** primitive k=8 shapes, every one satisfies
`measS7 ≤ L_y ≤ cap` and consec is the *unique* argmax of both, with margin; the L7 atlas has `|R|`
tiny (worst 0.003 at k=9, 0.0017 at k=10) and `p0_inf ≤ cap` at every ratio; the apex law `D=0 ⟺
7|pq` holds. The binding row is k=10 (margin 0.205), and the margin dips there then *recovers* — no
blow-up at large k.

## What is and isn't proved

Honest status. **Proved/elementary:** L1 (reduction), L2 (k≤7), L5, L6, the L7 tail `D≤14/p`, the
Delsarte per-shape bound (and its Lean formalization). **Verified by exhaustive exact computation:**
L3 (consec argmax, span ≤14), L4's finite collar window, the L7 finite atlas (k≤72 ratios). A
fully gap-free proof needs those finite exact computations *certified* — the natural next
formalization target (native_decide over the rational atlases), which the continuous measure side
currently keeps out of reach of the mathlib-free Lean module. And the L1/cone coupling (cluster
cover ⟹ LRC witness via the small part) deserves one more careful write-up.

So: the sector route is assembled and audited, the resonance gap is closed by an elementary bound,
and the Delsarte heart of it is a machine-checked theorem. The conjecture at fourteen is not proved.
But for the first time the remaining work is *certification of finite checks*, not *finding the
argument* — and the argument's core is now something a proof assistant has agreed to.
