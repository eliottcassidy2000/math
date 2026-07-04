# Seven dormant threads into the covering-min core — an idea-generation pass

*mac-mini-2026-07-04-S39. Owner: work the remaining core AND mine prior work (related or unrelated) for
ideas. The remaining core (after S30–S38): prove `primitive covering ⟹ M ≥ 14/183`, equivalently the
extremal LRC config is a `{kα}` three-gap progression. Below: seven repo threads aimed at it, each with a
concrete proposal and an honest verdict. Two verified facts anchor the pass.*

## Two anchors (verified this session)
- **Φ₆ iterated = Sylvester's sequence.** `2 → 3 → 7 → 43 → 1807 → …` (`a_{k+1}=a_k²−a_k+1=Φ₆(a_k)`). So
  `covering-min = n/Φ₆(n)` lives on the **greedy Egyptian-fraction** tower (Sylvester = greedy expansion of
  `1 = 1/2+1/3+1/7+1/43+…`), the multiplicative twin of the **Ostrowski/Zeckendorf greedy CF** ladder
  `M_k=[0;n−1,k]` (HYP-4078). The covering-min is a **double greedy fixed point**.
- **Ladder quantization FAILS (0%).** Generic covering families have off-ladder `M` (`5/24, 2/21, 3/17,…`);
  only the *extremizers* sit on `{k/((n−1)k+1)}`. So the Ostrowski ladder is the **extremal locus**, not a
  universal quantization — the honest form of "M-uniqueness" (HYP-3739).

## Thread 1 — Covering radius / Ramanujan spectral slack (opus, 2026-06-30)
`M(S)` is a *dynamic* covering radius; the even block realizes the equally-spaced cycle `C_n` (radius `1/n`);
the primitive covering-min = `C_n` perturbed by a forced **killer-defect**. The apex gap = `C_p`'s slack
from the Ramanujan bound, `2+λ_min(C_7)=4sin²(π/14)=0.198`. **Verdict:** clean framing, but the *dynamic*
(flow, not rigid rotation) nature blocks static covering-radius/sphere bounds (sphere gives only `M≥1/26`).

## Thread 2 — Ostrowski ⊗ Sylvester greedy duality (mac-mini S38/S39) — NOVEL
The covering-min is the unique **double-greedy** value: `[0;n−1,n]` (CF) and `n/Sylvester-iterate` (Egyptian).
**Proposal:** recast the covering constraint (a multiple of every `q≤n`) as an **Egyptian/covering-system**
condition, and derive `M≥n/Φ₆` from **Sylvester-greedy optimality** (the greedy Egyptian expansion is the
`smallest`/`most-efficient` — Kellogg/Curtiss-type bounds). **Verdict:** the most novel handle; explains
M-uniqueness; the missing piece is the precise "covering-system ↔ Egyptian-greedy" dictionary (cf. HYP-3724).

## Thread 3 — Delsarte LP / Toeplitz-PSD hyperbolicity cone (opus, 2026-07-01) — MOST PROMISING (rigorous)
Covering ⟺ a **nonnegative trigonometric polynomial** (Toeplitz moment matrix PSD, a hyperbolicity cone,
barrier `−log det`). The `M≥c` bound is a **Beurling–Selberg / Fejér certificate** (opus's magic function
`F_7`). **Proposal:** solve the SDP for the *sharp* covering-min certificate — a nonneg trig poly `p=Σĉ_m e(mt)`
with the danger-arc pairing forcing `M≥14/183` for every covering `S`. This is the rigorous form of the
measure route klein/opus already push (THM-515 singular series). **Verdict:** the likeliest actual proof route;
the honest obstacle is that LP/Delsarte natively gives *upper* covering-radius bounds — the LRC needs the
*lower* (deep-hole-exists) direction, which is the dual **existence** SDP (feasibility of a hiding measure).

## Thread 4 — Heptagon field `Q(cos 2π/7)` Galois rigidity (kps) — SPECULATIVE
`14=2·7`; the apex-7 cubic `x³+x²−2x−1` (min poly of `2cos 2π/7`) governs the extremal family. **Proposal:**
the extremal config is **Galois-stable** under `Gal(Q(ζ_7))`; rigidity = "fixed by the heptagonal Galois
action." **Verdict:** the apex-7 is fundamental (7-Fourier-zeros, `F_7`, Φ₆ genus-1 `X₀(14)`), but no concrete
field-theoretic lever yet — a lead, not a mechanism.

## Thread 5 — Erdős–Turán equidistribution band-barrier (opus, 2026-07-01) — GROUNDED, coarse
Disorder is penalized by Weil `√`-cancellation of `Σ_v e(mv/D)`; the covering-min = the **lowest-discrepancy**
covering family (the interval `{1..n−2}`+killer, band-gap 3). **Proposal:** the sharp all-moduli Erdős–Turán
discrepancy (`q=2..2n` simultaneously) collapses beaters onto the interval. **Verdict:** opus's honest
character-sum path; grounded mechanism, coarse metric — needs the full multi-modulus discrepancy sum.

## Thread 6 — Stern–Brocot / Farey geodesic recursion (opus + mac-mini) — RECURSION
The CF rung descent is a **geodesic in the Farey tessellation** = the `(2,3,7)` hyperbolic plane = the
self-concordant barrier cone. **Proposal:** an **inductive** covering-min proof — `covmin(n)` reduces to a
smaller rung via the Stern–Brocot mediant step. **Verdict:** the recursion is real (the kernel `K(p,q)` is a
three-gap/Stern–Brocot function); whether it bottoms out to a finite base is untested.

## Thread 7 — Tournament bridge: LRC-extremal = the `χ²` orbit (repo, "unrelated") — WILD
The repo's other half (Rédei/tournaments) has `chi-separates-regular-tournaments…LRC-extremal-is-the-χ²-orbit`,
`lrc-orbit-functor-rigidity`, `lrc-dual-burnside-orbit`. **Proposal:** the LRC tight locus = a distinguished
**orbit** under the metagraph/Burnside action, and *orbit-rigidity* (a fixed-point/functor argument) forces
uniqueness. **Verdict:** genuinely orthogonal; if the LRC extremizer really is a single group orbit, orbit
uniqueness is a categorical rigidity worth a dedicated look.

## Net — where I'd put the chips
The core is open (LRC(14)). Ranked by proof-likelihood: **(3) Delsarte/Beurling–Selberg certificate**
[rigorous, the measure route's endgame] > **(2) Ostrowski⊗Sylvester greedy** [novel, needs the Egyptian
dictionary] > **(5) all-moduli Erdős–Turán** [grounded, coarse] > (6) Stern–Brocot induction > (1) covering
radius > (7) orbit rigidity > (4) heptagon Galois. The unifying picture: the covering-min is the **double
greedy (CF+Egyptian) fixed point**, certified by a **Fejér/Beurling–Selberg positive polynomial**, with the
apex-7 as the arithmetic and Steinhaus as the rigidity. See HYP-4079.
