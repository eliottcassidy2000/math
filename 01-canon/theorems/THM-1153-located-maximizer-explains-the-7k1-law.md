---
id: THM-1153
title: BOXEPH'S LOCATED-MAXIMIZER FORMULA EXPLAINS MY 7k+1 LAW — and corrects its scope. (I) CROSS-AGENT VERIFICATION: boxeph-S120's theorem M = |vᵢaⱼ − vⱼaᵢ|/(vᵢ+vⱼ) at a straddling active pair holds exactly on every clustered family in my thread, and the active pair is always **(core minimum, active killer)** with aᵢ = 0, so **M = a_k·v_min/(v_min + k)**. Verified: {2..12}∪{169,182} gives 8/57 with pair (2,169), sum 171, 24 = 12·2; {2..12}∪{182,196} gives 13/92 with (2,182), sum 184, 26 = 13·2; {1..11}∪{312,364} gives 26/313 with (1,312), sum 313, 26 = 26·1; {1..11}∪{168,338} gives 28/339 with (1,338), sum 339, 28 = 28·1. (II) IT EXPLAINS MY cont.55 "7k+1" LAW: for the {2,…,12} core the exact value is M(core) = 1/7, and adding killers perturbs it to the **Farey neighbour n/(7n+1)** — 8/57 has 7·8+1 = 57, 13/92 has 7·13+1 = 92. The 7 is not a coincidence of the covering construction; it is 1/M(core). (III) AND IT CORRECTS THE LAW'S SCOPE, which I had not noticed: the 7n+1 form holds **only for v_min = 2 cores**. For cores containing 1 the core value is M({1,…,11}) = 1/12, and the law becomes **n/(12n+1)** — 26/313 has 12·26+1 = **313** exactly. My cont.55 statement of the law as universal was wrong; 7·26+1 = 183 ≠ 313. (IV) THE GENERAL FORM: the killers push the core's own exact value 1/D to a Farey neighbour **M = n/(Dn + c)** with D = 1/M(core) and c small — c = 1 in the four cleanest rows, c = 3 for {1..11}∪{168,338} where 12·28+1 = 337 against the actual 339
status: (I) VERIFIED exactly on six families in rational arithmetic — a genuine cross-check of boxeph-S120's theorem against an independently computed set of exact optima. (II),(III) PROVED in the sense that the arithmetic identities are exact and checked; (III) is a **correction to my own cont.55 law**, whose universality I had asserted without testing a v_min = 1 core. (IV) is a measured pattern with one row already deviating (c = 3), so the general form is a shape, not a theorem. Uniform r=5 remains OPEN — none of this bears on that directly
source: kind-pasteur-2026-07-18-S128 (cont.81–82; owner: look for creative closed forms, pull from other agents, explore their ideas)
depends_on: []
related: [THM-1152, THM-1151, MISTAKE-173]
script: 04-computation/closed_form_hunt_kps_S128c81.py, closed_form_synthesis_kps_S128c82.py (+ .out)
---

# THM-1153 — the located maximizer explains 7k+1

## (I) boxeph's formula, checked against my exact values

boxeph-S120 proved that the optimum sits at a **straddling active pair** with

> M = |vᵢaⱼ − vⱼaᵢ| / (vᵢ + vⱼ).

My thread had independently computed exact optima for several clustered families, so this is
a sharp cross-check. It holds every time, and the active pair is always **(core minimum,
active killer)** with aᵢ = 0, giving

> **M = a_k · v_min / (v_min + k).**

| family | M | active pair | v_min+k | M·(v_min+k) |
|---|---|---|---|---|
| {2..12} ∪ {169,182} | 8/57 | (2,169) | 171 | 24 = 12·2 |
| {2..12} ∪ {169,196} | 8/57 | (2,169) | 171 | 24 = 12·2 |
| {2..12} ∪ {182,196} | 13/92 | (2,182) | 184 | 26 = 13·2 |
| {1..11} ∪ {312,364} | 26/313 | (1,312) | 313 | 26 = 26·1 |
| {1..11} ∪ {168,338} | 28/339 | (1,338) | 339 | 28 = 28·1 |

(My first search missed these because I capped |a| ≤ 3; the correct a_k are 12, 13, 26, 28.)

## (II) It explains the 7k+1 law

In cont.55 I found empirically that the covering families satisfy **M = n/(7n+1)** and read
the 7 as a feature of the covering construction. It is not. The core {2,…,12} has exact
value **M(core) = 1/7**, and adding killers perturbs it to the **Farey neighbour** just
below:

> 8/57 with 7·8 + 1 = 57,  13/92 with 7·13 + 1 = 92.

The 7 is simply 1/M(core).

## (III) And it corrects the law's scope

I stated the 7n+1 law without testing a core containing 1. It fails there:
{1..11} ∪ {312,364} has M = 26/313, and **7·26 + 1 = 183 ≠ 313**. But M({1,…,11}) = 1/12,
and

> **12 · 26 + 1 = 313** exactly.

So the law is **n/(Dn+1) with D = 1/M(core)**, not 7n+1 universally. My cont.55 statement
was scoped to v_min = 2 without my realising it.

## (IV) The general shape

> **M = n/(D·n + c)**, D = 1/M(core), c small.

c = 1 on four of the five rows; {1..11} ∪ {168,338} gives 28/339 against 12·28+1 = 337, so
c = 3 there. The shape is right, the constant c is not yet pinned.

## Why this is worth having

It converts a numerology I had been carrying since cont.55 into a statement about the core's
own optimum, and it locates the mechanism: **the killers do not set the scale — the core
does, and the killers only push M to the neighbouring Farey fraction.** That also explains
why the clustered families were never close to 1/14: the core value 1/7 or 1/12 is already
comfortably above, and a Farey step down from it cannot cross the threshold.

## Fleet notes

- **death-star-S58b independently refuted my centre-hitting criterion** (MISTAKE-173), with
  117 non-proportional hitters against my 96 — they counted permutations of the centre,
  which I had not. Their salvage is a "six-ray sojourn-max", which is the right shape for the
  invariant THM-1152 left unidentified.
- codex-S78 has proved universal five-comb dual noncoverage, which is further than the
  four-comb line I was pushing.

## Named next
- Pin c in (IV): it should be computable from which killer is active and the core's Farey
  neighbours, since the located pair is known.
- Apply the located-maximizer formula to the r=5 and r=6 clustered families directly — if the
  active pair is always (v_min, killer) there too, M has a closed form in those strata and
  the threshold comparison becomes arithmetic rather than computational.
