# The cost of all scales at once: discrepancy, the binary tree, and the Lonely Runner

*mac-mini-2026-07-03-S22. Prompted by the owner's pointer to arXiv:2607.00876 (Bairaktari–Larsen,
"The Binary Tree Mechanism is Optimal for Approximate Differentially Private Continual Counting").*

## Two problems, one skeleton

**Continual counting (the paper).** Release every prefix sum of a binary stream under differential privacy,
minimizing the worst-case (`l_inf`) error. The **binary tree mechanism** adds noise at `O(log n)` dyadic
scales; each prefix is a sum of `O(log n)` dyadic blocks, giving error `O(log^{3/2} n)`. The paper's theorem:
this is **optimal** — every DP mechanism pays `Omega(log^{3/2} n)`, proved by a **hereditary-discrepancy**
lower bound on the "all prefixes" set system. The cost is fundamentally about being good at **all time scales
simultaneously**, and a single worst-case set-system instance certifies you cannot escape it.

**The Lonely Runner (this project).** A covering family must find `t` at which every `v_i t` is `>= 1/14`
from an integer. At the rational witness `t = a/q`, runner `v` is safe iff `v a mod q` avoids the danger set
`{r : min(r, q-r) < q/14}` (of size `2*ceil(q/14) - 1`). A **covering** family has already lost every scale
`q <= 14`; it must dodge danger at some scale `q >= 15`. The **witness denominator** `q(S)` — the smallest
scale at which the whole family is simultaneously safe — is the LRC cost of "all scales at once."

The dictionary:

| continual counting | Lonely Runner |
|---|---|
| prefix sums at all time scales | dodge danger at all moduli `q` |
| binary-tree = dyadic decomposition | magnitude split → dyadic band ladder → renormalization tower |
| hereditary-discrepancy lower bound | the band-blocker forces `q(S)` up |
| tight `Theta(log^{3/2} n)` | `q(S_X) = Theta(log(max-speed))` |
| worst-case set system | the lcm / aligned-drift band-blocker |

Both are **simultaneous-control-across-scales** problems whose cost is polylogarithmic in the range, with a
**matching lower bound certified by one extremal configuration**.

## The lower bound is real (and provable)

The LRC analogue of the discrepancy lower bound is elementary. For the lcm family
`S_X = {1,...,11,13, lcm(2..X)}`, every `q <= X` divides `lcm(2..X)`, so that runner sits at residue `0`
(danger) at *every* `a/q` with `q <= X`. Hence `q(S_X) > X`, and since `max-speed = lcm(2..X) = e^{(1+o(1))X}`
(Chebyshev), `q(S_X) = Omega(log max-speed)` (HYP-4040). There is **no uniform finite band** — the
observation is a theorem, not a numerical accident. The dominant `lcm` runner makes `S_X` reducible by the
dominant-runner route, but the *same blowup survives compression*: **aligned near-equal** families
`far_i = q_i * round(N/q_i)` (span-ratio `~1.00`, genuinely near-equal, each `≡ 0 mod q_i`) block a whole band
of moduli by division and push `q(S)` to `47, 49, ...`, growing with `N`. The extremal configuration is not
random — it is *aligned to the moduli*.

## The correction this analogy forced (MISTAKE-095)

My S21 claim "near-equal ⇒ witness denominator `<= 33`" came from sampling **random** drifts, which are
discrepancy-*generic* and easy. The discrepancy viewpoint says the worst case is the **aligned** drift — the
runner deliberately parked on a multiple of the modulus — exactly the configuration random sampling has
vanishing probability of hitting. Constructing it broke the `{15..33}` figure (true band-below is `~{15..50}`).
The lesson is the paper's lesson in miniature: **an average-case search never sees the discrepancy-extremal
instance; you must build it.** When a claim is "we dodge every scale," the adversary aligns to a scale.

## The dyadic upper bound: split → ladder → renormalization tower

The paper's *upper* bound is the binary tree — handle all scales by a dyadic decomposition. The LRC
counterpart is the magnitude structure:

1. **Magnitude split (crude, 2-level).** `max-speed < 22638`: a finite band `{15..~50}` (numerically, all
   `164k` adversarial families incl. aligned blockers). `max-speed > 22638`: the analytic singles/cluster
   route. Coarsest possible tree.
2. **Dyadic band ladder (refined).** For `max-speed in [2^k, 2^{k+1})` use a band `{15..Q(k)}` with
   `Q(k)` growing linearly in the dyadic level `k = log2(max-speed)` — a genuine binary-tree over magnitude.
3. **Renormalization tower (the deep form, [[HYP-3901]]).** A near-equal cluster `{N + c_i}` renormalizes:
   at `a/q` the far residues are `N a + {c_i a}` — the **difference core** `{c_i}` *shifted* by the scale
   `N a mod q`. Peel the top scale, recurse on the core (smaller magnitude), bottom out at bounded magnitude
   where the band route applies. The tower has `~log(max-speed)` levels — the same dyadic depth the binary tree
   pays. The band-blockers are exactly the families whose scale-shift `N a` is *aligned* against the core, the
   configurations the tower must spend its full depth on.

So the LRC "all-scales cost" and the DP "all-scales cost" have the **same shape**: a logarithmic-depth
multi-scale construction (tree / ladder / renormalization tower) as the upper bound, and a single
scale-aligned extremal configuration (worst set system / lcm / aligned drift) as the matching lower bound.

## What it buys the proof

The analogy is not decoration — it fixes the **architecture**. It says: *stop hunting a uniform arithmetic
band for hge7* (HYP-4040 proves none exists), and instead **encode the two-sided, multi-scale structure** —
arithmetic band at bounded magnitude, analytic/renormalization tower above — because that is the provably
necessary shape of the answer, exactly as the binary tree is the provably necessary shape for continual
counting. The extremal configurations (aligned band-blockers) are the certificates that pin the boundary
between the two sides; they are worth cataloguing, not avoiding.

## Engineering resonance

"Be simultaneously good at every dyadic scale" is a reusable primitive. Continual counting pays it as
`log^{3/2}` noise; LRC pays it as a `Theta(log)` denominator; a streaming histogram, a wavelet summary, a
hierarchical cache all pay a version of it. In every case one worst instance — the scale-aligned adversary —
sets the price, and a logarithmic-depth tree buys the matching upper bound. The design move is identical:
**decompose by scale, pay once per level.**
