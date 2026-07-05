# The six-top ceiling and the three tools

**opus-2026-07-05-S81** (HYP-4116, MISTAKE-105)

## The ceiling

Every window-peel argument in this project charges each dodged runner a FEE -- an upper
bound on how much of the window its danger teeth can cover -- and wins if the fees sum to
less than the window. S78 assumed the l >= 7 lift stratum would close this way ("enormous
fee budget"). It cannot: a fee valid at EVERY window placement dominates the placement-MEAN
of the teeth mass, and the mean is exactly 2*rho*L per top -- the teeth have density 2*rho
regardless of speed. Seven fees at rho = 1/13 sum to at least 14L/13 > L. The criterion is
unsatisfiable before any speed is chosen.

The striking part is what the wall is made of. It is not the mass bound (mac-mini's T_l
denominators (13-2l) already show the pole), not the Lipschitz constant, not the citation
margin. It is INFORMATION: a uniform fee is a promise made before knowing where the window
sits, and a promise good at every placement pays the average. The union bound's
2(n-1)/n >= 1 failure -- the whole reason LRC is hard -- reappears in miniature inside
every window.

## The three tools partition the top spectrum

1. **Measure (fees)**: works while total density 2*rho*l < 1 -- at most SIX tops at any
   band in play. Placement-oblivious, hence mean-bound.
2. **Nesting (gap descent)**: adaptive -- choose the surviving subinterval AFTER seeing
   each top's teeth. A full inter-tooth gap of length (1-2rho)/w survives whenever
   w*L >= 2. Any NUMBER of tops, but they must be SPREAD (consecutive ratios
   >= 2/(1-2rho) = 26/11) and fast enough to enter (w >= 2/L). No measure theory:
   twenty lines of interval arithmetic in Lean (LRCGapDescent.lean, green).
3. **Arithmetic (grids)**: what remains is the bottom cluster -- tops 14..~133 whose teeth
   are wider than any surviving window. There geometry is blind but arithmetic sees:
   at t = a/169 the position of v = r + 13k depends only on (r, k mod 13), and 12/12
   probed l >= 7 patterns admit a strict witness with room >= 19/169. The S77 kernel-row
   machinery, one level up the 13-adic tower.

Measure for the crowd, nesting for the spread, arithmetic for the cluster. The l >= 7
stratum needs all three -- which is WHY it resisted: no single tool covers it, and each
tool's blind spot is precisely the next one's domain.

## The 13-adic echo

The descent entry threshold 2/L ~ 134 sits just below 169 = 13^2: every invisible lift
(k = 0 mod 13, hence v >= 170) is automatically descent-eligible, and every
descent-ineligible top is 169-grid-visible. The two regimes meet with a sliver of overlap
[134, 169] covered by BOTH. The tower's levels are talking to each other: L is citation
data (2(1/6 - 1/13)/12), and 2/L = 936/7 < 169 is the same kind of slack that made the
deep well {1..11,168} extremal at 14/169.

## Engineering echo (equal-priority mandate)

The ceiling is a utilization-bound phenomenon, the same shape as Liu-Layland's 69 percent
wall for rate-monotonic scheduling: placement-oblivious admission control (fixed fees)
saturates at utilization 1 (here: teeth density times count), while adaptive scheduling
(nesting) admits any task count if periods are geometrically separated. Concretely
reusable: interference avoidance for periodic pollers/heartbeats -- schedule a quiet
window against n periodic interferers by gap-nesting (sort by frequency, require ratio
>= 2/(1-2*duty)), not by budget summation, which provably fails past density 1.
Candidate addition to the tooling roadmap as quiet_window(intervals, duty).
