# Thread C: the measS7(consec_k) extremal value, its closed form, and the k=12 extremality break

*mac-mini-2026-06-20, LRC(14) sector route, Thread C (extremal value + finite-family extremality).*

## What this thread set out to do

Find (a) a CLOSED FORM for `measS7(consec_k)` and (b) a structural proof that consec
maximizes `measS7` over the finite low-relation-height family. The honest outcome is a
proved partial closed form plus a *correction* to the extremality claim that other threads
(HYP-2699, THM-534) had already started flagging.

## The theta-reparametrization is the right frame (re-confirmed S7b)

With `theta = 7x in [0,7)`, the Z/7 color of speed `e` is `floor(e theta) mod 7`, and
`measS7(consec_k) = (1/7) sum_{j=0}^{6} M_j(k)` where, writing `theta = j + s`,

  `M_j(k) = meas{ s in [0,1) : { (e*j + floor(e*s)) mod 7 : e=0..k-1 } = Z/7 }`.

The walk is `p_{e+1} = p_e + (j + b_e) mod 7`, `b_e in {0,1}` the **Sturmian/mechanical
word of slope s**. This is a deterministic cover-time of a mod-7 rotation driven by a
mechanical word — a combinatorics-on-words object, NOT a coupon-collector. (The iid
coupon value `iid_k = 7! S(k,7)/7^k` is far below: corr = measS7 - iid is +0.14..+0.42.)

## PROVED closed forms (the two monotone strips)

The only strips where the walk is **monotone** are j=0 (steps {0,1}) and j=6 (steps
{0,-1} since 6,7 ≡ -1,0). Monotone ⇒ cover = "reach distance 6", a single threshold:

  * **M_0(k) = (k-7)/(k-1)**  (cover ⟺ floor((k-1)s) ≥ 6 ⟺ s ≥ 6/(k-1)).
  * **M_6(k) = (k-6)/(k-1)**  (mirror).
  * Hence **M_0(k) + M_6(k) = 2 - 11/(k-1)** EXACTLY.

Verified exactly for all k=7..30.

## PROVED partial closed form: the strip telescope

Every middle strip has a unique "bad slope" where the cover-time blows up, and the
cover-time **telescopes** in a Farey/Stern-Brocot fan around the fast slopes. The cleanest
instance (the symmetric strip j=3):

  on `s in [1/2 - 1/n, 1/2 - 1/(n+1))`,  `tau(s,3) = 2(n+1)` exactly (n=2..13 verified).

So the bulk of each middle strip is a clean harmonic telescope; only O(1) intervals around
the bad slope need a finite Farey correction. A fully uniform middle-strip formula is NOT a
simple rational expression (the bad slope and its CF tail differ per j).

## EXACT FINITE ALGORITHM = a genuine closed form per k

`floor(e s)` for e=0..k-1 is constant on each interval of the **order-(k-1) Farey
subdivision** of [0,1). Hence
  `M_j(k) = sum over order-(k-1) Farey intervals I of |I| * [cover at mid(I)]`,
a finite exact rational. The cover-set automaton (state = (pos in Z/7, covered subset)) has
only **626 reachable states up to depth 13** — a finite transfer system over the
Stern-Brocot tree computes every `measS7(consec_k)` exactly.

The exact ledger (matches the canon table):
```
k :  measS7(consec_k)
7 :  31/210      = 0.147619
8 :  481/1470    = 0.327211
9 :  2447/5880   = 0.416156
10:  8899/17640  = 0.504478
11:  3419/5880   = 0.581463
12:  121103/194040 = 0.624114
13:  14573/21560 = 0.675928
14:  14109/20020 = 0.704745
```

## THE CORRECTION: consec is NOT the measS7-maximizer for k ≥ 12

Exhaustive exact check over the LRC(14) finite family (E = {0} ∪ (k-1 from {1..13}),
4095 sets total):

```
k :  consec=max?   TRUE maximizer                 measS7(max)
8 :  YES           {0,1,2,3,4,5,6,7}              481/1470   = 0.327211
9 :  YES           {0,...,8}                       2447/5880  = 0.416156
10:  YES           {0,...,9}                       8899/17640 = 0.504478
11:  YES           {0,...,10}                      3419/5880  = 0.581463
12:  NO            {0,...,10, 12}                  11381/17640= 0.645182  (> consec 0.624114)
13:  NO            {0,...,10, 12, 13}              19492/28665= 0.679993  (> consec 0.675928)
```

The maximizer ALWAYS keeps the consecutive prefix `{0,...,10}` and pushes the remaining
offsets to the top of the window (**skipping 11**). For k=8,9,10 the maximizer is consec
even over the wider THM-535 binding windows (span≤17/16/15).

This **exactly mirrors HYP-2699's k=12 anomaly for U4** (same beater `{0..10,12}`) and
THM-534's "consec not L_y-max at k=11,12,13". The three functionals (measS7, U4, L_y) all
break consec-extremality at the same place and via the same shape. That is not coincidence:
they are all coarse readouts of the same empty-sector count N, and the relation lattice of
`{...,10,11}` carries a low-height vector (11 ≡ 4 mod 7 collides with the prefix residues)
that `{...,10,12}` dodges.

## Why no per-strip / per-band / majorization proof exists (DEAD, re-confirmed)

consec does NOT maximize each strip `M_j` separately (the j=0 strip alone is beaten by
spread AP-like sets by up to +0.149). When consec DOES win the aggregate (k≤11), it wins
some strips and loses others — the runner-up `{0..k-2,k}` beats consec in 2-3 strips and
loses in the rest, with the consec advantage being a genuine **mixed-sign sum**. The
extremality is irreducibly aggregate even at the strip level. This independently confirms
the DEAD list (per-block / per-band / monotone / majorization / synchronization).

## What this means for the S3 gate (the load-bearing point)

The S3 closure must NOT be justified by "consec is the extremal shape" (FALSE for k≥12).
The correct, fully-verified statement is **true-max < cap_k**:

```
k :  true_max(measS7)   cap_k        max < cap?   slack
8 :  0.327211           0.381463     YES          +0.054252
9 :  0.416156           0.494256     YES          +0.078099
10:  0.504478           0.604396     YES          +0.099917
11:  0.581463           0.725275     YES          +0.143812
12:  0.645182           0.857143     YES          +0.211962
13:  0.679993           1.000000     YES          +0.320007
```

The slack GROWS in k, so the dangerous row is the smallest, k=8 (slack +0.054), exactly as
the cap-floor analyses found. This is a **4095-set finite exact check**, trivially feasible
(<25 s). It closes the *finite* (span≤13) part of S3 unconditionally; the remaining LRC(14)
gap is still the wide-spread bound (HYP-2608/2611), untouched by this thread.

## Files

* `04-computation/lrc14_threadC_closedform_strips_macmini_0620.py` — strip decomposition,
  proved M_0/M_6 closed forms, j=3 telescope, exact ledger.
* `04-computation/lrc14_threadC_truemax_gate_macmini_0620.py` — corrected S3 gate
  (true-max < cap_k), exhaustive exact over the 4095-set family.
* outputs in `05-knowledge/results/`.
