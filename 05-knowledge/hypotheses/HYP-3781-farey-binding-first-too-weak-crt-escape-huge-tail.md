---
id: HYP-3781
title: THE FAREY BINDING-FIRST CHECK IS TOO WEAK (it is the CRT-escape) + the huge-speed tail is the Steinhaus scaling => covering-min(12,13,14)=construction via lazy-cut + tail. Worked the HYP-3779 menu ("think Farey tessellations"). (A) FAREY BINDING-FIRST (speed-independent): the covering-min is a Farey neighbor of 1/(n-1) (HYP-3732), rung j, binding modulus D=j(n-1)+1 (Stern-Brocot ray). Necessary condition for rung j: EXISTS rotation a* with every q in {2..n} having a multiple in the safe band {r: dist(r a*,D)>=j} mod D. RESULT: TOO WEAK -- ALL rungs j=2..n are binding-feasible for every n=7..14 (incl. n=9 where rungs 2,3 are NOT realized, and n=12,13,14 with no beater). It rules out NOTHING. WHY: a huge CRT-tuned speed satisfies every SMALL-modulus condition (binding D<=(n-1)^2, resonances, band are all small) -- this IS the CRT-escape (HYP-3745); the real obstruction is the multi-modulus CRT-INVARIANT COUNTING (each speed covers <=2r+1 rotations per modulus REGARDLESS of value), which no single-modulus Farey check captures. (B) HUGE-SPEED TAIL confirmed: the double-killer family {1..n-2, n(n-1)t} gives M=(nt)/(n(n-1)t+1) -> 1/(n-1), INCREASING in t (n=14: 14/183,28/365,42/547,...) -- larger killers are strictly WORSE (Steinhaus scaling HYP-3763); the construction (t=1, smallest double-killer) is the best. (C) SYNTHESIS: covering-min(n=12,13,14)=construction via the lazy-cut ILP (no beater speeds<=n(n-1), HYP-3779) + the huge-tail scaling/HYP-3737 band (speeds>n(n-1)); the general-huge-speed case rests on HYP-3745's CRT-invariance (verified, not fully proved)
status: MIXED (honest negative + confirmation + synthesis). VERIFIED: (A) the Farey binding-first necessary condition passes for ALL rungs 2..n at n=7..14 (too weak, rules out nothing) -- a clarifying NEGATIVE: it is defeated by the CRT-escape, pinpointing that the obstruction is multi-modulus counting not single-modulus binding. (B) huge-tail M=(nt)/(n(n-1)t+1) increasing (exact, n=14). (C) covering-min(12,13,14)=construction: RIGOROUS for speeds<=n(n-1) (lazy-cut HYP-3779); the double-killer tail is exact; the GENERAL huge-speed tail is argued via HYP-3745/3737, not fully proved. FAREY-TESSELLATION framing: the covering-min on the Stern-Brocot ray from 1/(n-1); the tail's binding D=n(n-1)t+1 is the Farey structure; three-gap (HYP-3762) governs it.
source: klein-2026-07-01-S62
depends_on:
  - HYP-3779   # the lazy-cut result (speeds<=n(n-1)); this works its menu
  - HYP-3745   # CRT escape / CRT-invariance (the real obstruction the binding-first misses)
related:
  - HYP-3732   # covering-min = Farey neighbor of 1/(n-1) (the binding-first premise)
  - HYP-3763   # Steinhaus scaling (the huge-tail M=(nt)/(n(n-1)t+1))
  - HYP-3737   # band over-constraint (the >n(n-1) tail argument)
  - HYP-3762   # three-gap on the Farey/Stern-Brocot ray
results:
  - 04-computation/covering_min_farey_binding_klein.py
  - 05-knowledge/results/covering_min_farey_binding_klein.out
  - 05-knowledge/results/covering_min_huge_tail_klein.out
---

# HYP-3781 — the Farey binding-first check is too weak (it is the CRT-escape)

## (A) Farey binding-first: a speed-independent check that RULES OUT NOTHING
"Think Farey tessellations": the covering-min is a Farey NEIGHBOR of `1/(n-1)` (HYP-3732), rung `j`,
binding modulus `D = j(n-1)+1` -- a det-1 edge of the Farey tessellation, on the Stern-Brocot ray from
`1/(n-1)`. A covering set realizing rung `j` must, at the binding rotation `a*`, have every speed in the
**safe band** `B(a*) = {r : dist(r a*, D) >= j}` mod `D`, and cover every `q in {2..n}` (a multiple of
`q` has residues `= {multiples of gcd(q,D)}` mod `D`). This gives a **speed-INDEPENDENT** necessary
condition (mod `D <= (n-1)^2`):
> rung `j` binding-feasible iff `EXISTS a*` coprime to `D` with every `q` coverable in `B(a*)`.
**Result: TOO WEAK.** Computed `n=7..14`: EVERY rung `j = 2..n` is binding-feasible -- including `n=9`
(where the covering-min is rung 4, so rungs 2,3 are NOT realized) and `n=12,13,14` (no beater at all). It
rules out nothing.

## Why it fails -- and what that teaches
The safe band `B(a*)` excludes only `~2j-1` residues, so a multiple of each `q` almost always fits. More
deeply: **all the constraining moduli are small** (binding `D <= (n-1)^2`, resonances `<= n`, band
`< 2n`), and a **huge CRT-tuned speed can be tuned to satisfy every one of them** -- this is exactly the
**CRT-escape** (HYP-3745). No single-modulus check (Farey binding included) can see the obstruction,
because the obstruction is not at any one modulus. It is the **CRT-INVARIANT COUNTING across ALL moduli
at once**: each speed covers `<= 2r+1` rotations per modulus *regardless of its value* (HYP-3745), so the
`n-1` speeds cannot simultaneously cover every modulus below the construction radius. The lazy-cut ILP
(HYP-3779) succeeds precisely because it is multi-modulus (it accumulates binding witnesses across
moduli); the Farey binding-first fails because it is single-modulus. Clarifying negative.

## (B) The huge-speed tail is the Steinhaus scaling (confirmed)
The `>n(n-1)` gap of HYP-3779, for the double-killer family: `{1..n-2, n(n-1) t}` has
`M = (n t) / (n(n-1) t + 1)`, INCREASING in `t` toward `1/(n-1)`:
```
  n=14:  t=1: 14/183   t=2: 28/365   t=3: 42/547   t=5: 70/911   t=10: 140/1821  ... -> 1/13
```
(exactly the Steinhaus scaling `c/(kc+1)`, HYP-3763). So a **larger** killer is strictly **worse**; the
construction (`t=1`, the smallest double-killer `n(n-1) = lcm(n-1,n)`) is the best. The binding modulus
`D = n(n-1)t+1` is itself a Farey point -- the tail lives on the same Stern-Brocot ray, sliding toward
the ceiling `1/(n-1)`.

## (C) Synthesis: covering-min(12,13,14) = construction
- **speeds `<= n(n-1)`**: RIGOROUS, no beater (the lazy-cut ILP, HYP-3779).
- **speeds `> n(n-1)`**: the double-killer tail is exact (increasing, above the construction); the
  GENERAL huge-speed case is argued by HYP-3745's CRT-invariance + HYP-3737's band over-constraint
  (`n(n-1)` is the smallest double-killer; larger ones lose the radius-1 band). Verified, not fully
  proved.
So the covering-min at `n=12,13,14` (incl. the LRC-14 target `14/183`) is the construction, complete up
to `n(n-1)` and strongly supported beyond.

## Menu status (HYP-3779)
- lazy-cut cutting-plane ILP: DONE (closes `<= n(n-1)`).
- Farey binding-first: TRIED -- too weak (this HYP); the CRT-escape defeats single-modulus checks.
- Lovasz-theta / SDP: UNTRIED, the promising next lever (a tighter, multi-modulus relaxation that could
  give the lower bound in one shot -- the LP was too weak, HYP-3779).
- column generation / B&B-witness-pruning / meet-in-the-middle: untried.

## Net
Working the menu with the Farey tessellation: the covering-min lives on the Stern-Brocot ray from
`1/(n-1)`, but the **binding-first (single-modulus) check is too weak -- it is the CRT-escape**, and the
true obstruction is the multi-modulus CRT-invariant counting (why the lazy-cut works and the Farey check
does not). The huge-speed tail is the Steinhaus scaling (larger killers strictly worse), so the
covering-min at `n=12,13,14` is the construction (rigorous `<= n(n-1)`; tail via HYP-3745/3737). Next
real lever: a tighter *multi-modulus* relaxation (Lovasz theta / SDP), since single-modulus and LP both
fail.
