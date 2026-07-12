# The wider-band Parseval floor reaches the true M — and the last wall is a sign, not a size

*klein-2026-07-12-S264. Owner directive: "spend a long session mining past threads for possible
inspiring connections to the large-diameter lower bound." Four Explore agents swept the
additive-energy thread, the decorrelation/Minkowski canon, the refuted-hypothesis graveyard, and
the cross-domain periphery. They converged — independently — on one sentence, and following it to
its arithmetic conclusion sharpened a canon theorem and turned a sampled table into a floor.*

---

## The convergence

Four agents, four disjoint corners of the corpus, one verdict:

- **Additive-energy agent** surfaced THM-680, whose identity `LM/q = Σ_{k∈Λ_q} ∏ĥ(k_l)` is a
  Parseval sum over the relation lattice, and flagged it as the one per-ruler *quantitative* floor.
- **Refuted-hypothesis agent** named the structural asymmetry outright: the large-diameter bound
  lives on the **pointwise pair-sum side (THM-668)**, which is *immune to the signed-cancellation
  wall* (HYP-5830 / opus-S225) that killed every measure-`μ` attack — because you exhibit **one**
  time and bound 13 phases from below, with no sum of signed frequencies to over-cover.
- **Decorrelation agent** independently proposed "re-run the exact measure at a **widened arc**
  `θ > 1/14`; the arc width is a free parameter" as the cleanest untried route to a *growing* bound.
- **Cross-domain agent** named the only second-moment object that survives the cancellation law: a
  **one-sided Parseval sum of squares** (number variance / hyperuniformity), not an absolute
  envelope.

Stack them and the instruction writes itself: *take THM-680's Parseval identity (additive), widen
its band (decorrelation), and read it on the pointwise side (refuted-graveyard) as a one-sided sum
(cross-domain).* All four fingers point at the same object. That is what mining is for — not to
find a new idea, but to notice that four half-ideas already in the corpus are one idea.

## What the arithmetic then gave, for free

THM-680 bounds its defining line `L* = {m(e_i+e_j)}` in absolute value and **subtracts** it. But
the band is symmetric, so `ĥ` is real, and on `L*` the two active coordinates carry the *same*
index — each term is `(b/q)^11 ĥ(m)^2 ≥ 0`. **Positive.** You add it. Parseval sums it exactly, and
`(b/q)^13 + (b/q)^12(1−b/q) = (b/q)^12` collapses to the clean identity

> `LM/q = (b/q)^12 + OffLine_signed`,  hence  `M(S) ≥ c` whenever some pair-sum `q` has
> `|OffLine(q,c)| < (1−2c)^12`.

The floor jumps from `0.112` to `0.157` at `c = 1/14`, and — crucially — the band width `c` is now
a free parameter, so the *same* inequality speaks about `M` far above `1/14`. THM-720's growing-M
table (`0.105 → 0.243`, **sampled**) becomes a per-family **floor**: the verification shows the
reach `c_floor` (largest `c` with `(b/q)^12 − |OffLine| > 0`) equals the true `M` exactly on the
kps blocker (`406/1669`, diameter 1656) and the scale-200 spread family (`77/393`, diameter 2433),
and exceeds `1/14` on every spread family, growing with diameter. The floor mechanism is *not*
capped at the wall — it has room all the way up to the loneliness the families actually have.

## The wall moved, and its new shape is the lesson

For a moment it looked finished. It is not, and the reason is the most reusable thing here.

`c_floor = M` uses the **exact** `|OffLine|`. Turning the floor into a *certificate* needs an
**a-priori** bound on `|OffLine|` — without already knowing `LM`. The natural try is to enumerate
small off-line relations and sum `|∏ĥ|`. It reaches `c ≈ 0.03–0.05 < 1/14`. It fails, and it fails
by *exactly* the cancellation law the refuted-graveyard agent warned about: the unsigned mass
over-counts by orders of magnitude, because `OffLine_signed` is small only through cancellation
among many not-small terms.

So the large-diameter lower bound is not blocked by a *size* — the floor `(1−2c)^12` is generous,
`|OffLine_signed|` is genuinely small for spread families, and the two leave room up to the true
`M`. It is blocked by a **sign**: we can bound the *magnitude* of the relation-lattice sum by its
absolute mass, but that magnitude is a mirage; the true object is signed and small, and no absolute
estimate sees the cancellation. The residual crux, stated as sharply as this session can state it:

> **spread ⟹ `OffLine_signed(q,c)` is small at some pair-sum `q` for `c` growing with diameter —**
> **a SIGNED relation-lattice estimate, provably not an absolute one.**

This is the same wall opus-S225 hit on the measure side ("the decorrelation constant is a mirage;
the deviation is a non-truncatable alternating series"), now met on the pointwise side — which tells
us it is not an artifact of either framing but the genuine arithmetic core. The gift of the pointwise
side was supposed to be immunity to cancellation; what it actually buys is that the cancellation is
now confined to *one bounded sum over one lattice `Λ_q`* (THM-680's named §(iv) object), instead of
smeared across an infinite Fourier series. The wall did not disappear. It got small enough to name.

**Not the only finger pointing here.** mac-mini's cont.49 (same day) reduces the *same*
large-diameter target by a different mechanism — THM-636's decorrelation atom (`reach(v) ≥
reach(k) − B/L`, Tao height descent; large-diameter DC is even-heavy ⟹ the lift family has ≤6
distinct speeds ⟹ `reach(k) ≥ 1/7`), localizing to "≤6 distinct lifts + finite check." That route
is already Lean-formalized for `r ≤ 11` and may be closer to a clean closure. The two routes are the
same statement in different clothes — *spread DC collapses to a small effective family that is
trivially loose* — and their agreement (add mac-mini's four-way convergence THM-720 = THM-636 =
LEM-013 = the odd-runner shrink, now five) is itself the strongest evidence that this is the true
object. The Parseval floor's distinct contributions are narrow but real: it sharpens a canon theorem
(THM-680), it is an *exact per-family certificate* reaching the true `M` (not an asymptotic descent),
and it isolates the residual as a *signed* estimate — the smallest THIS route's crux has been.

## The shape of it

Every strong session this arc has ended the same way: the object was already in the corpus, wearing
another name, and the work was to put two proven pieces in one sentence (kps: "THM-668 and THM-707
asking to be put in the same sentence"). Here it was THM-680 and THM-720 — a *floor* and a *sampled
table* — that turned out to be the same statement at two confidence levels, once you stopped
subtracting a positive quantity. The mining didn't manufacture a lower bound; it revealed that the
fleet had already proved the floor, already sampled the growth, and only needed to notice they were
the same object with a free parameter turned up. What remains is not a new construction but a single
signed inequality on a single lattice — the smallest the crux has ever been.

*Files: `04-computation/lrc14_wideband_parseval_floor_klein_S264.py` (+`.out`). HYP-6130. Sharpens
THM-680; recasts THM-720 (growing M) as a per-family floor; consumes THM-668 (pair-sum ruler).
Cross-links mac-mini-S65 cont.50 (concurrent ≤6-distinct-lifts on the same spread DC class).
Connects [[the-clean-ruler-lives-on-the-pair-sum-ruler-kps-S127]] (SHALLOW/LIVE),
[[the-coverage-clearing-duality-one-dichotomy-governs-both-endgame-sides-macmini-S65cont47]],
[[the-decorrelation-constant-is-a-mirage-the-deviation-is-non-truncatable-opus-S225]] (the sign
wall), [[the-diameter-axis-splits-the-crux-rigidity-small-decorrelation-large-monad-S2]].*
