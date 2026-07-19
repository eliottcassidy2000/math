# The centering witness closes the spread case: `M` of *every* arithmetic progression, and the full AP rigidity for LRC

*boxeph-2026-07-18-S118. Owner: work a new creative angle on the LRC(14) open math.
Result: a single **inverse-scaling witness** `t = d⁻¹/(2a+11d)` that reduces **any** odd-`d` AP to a
perfectly-centered 12-consecutive-block mod `q`. This **closes the spread case (`d≥17`) that S117 left
open** — the AP-rigidity face of HYP-4382 is now fully proved — and, as a bonus, pins down the **exact
loneliness of every arithmetic progression** (verified exhaustively). Verified S118.*

## The new idea (in one line)

For a primitive AP the whole loneliness problem is *the same* "center a run of consecutive residues mod
`q`" problem — you just have to pick the right modulus `q` and rescale by `d⁻¹`. S117 saw this only for the
consecutive case `d=1` (where `d⁻¹=1` is invisible). The move that unlocks the spread case is to **rescale
by the inverse of the common difference.**

## The witness

Let `C = {a, a+d, …, a+(n-1)d}` be primitive (`gcd(a,d)=1`), `n=12`, and suppose `d` is **odd**. Put

> `q = 2a + 11d`,  `p ≡ d⁻¹ (mod q)`,  and take `t = p/q`.

Because `d·p ≡ 1 (mod q)`, each speed reduces to
`(a+dk)·p ≡ ap + k (mod q)` — a **run of 12 consecutive residues** starting at `s := ap mod q`. And with
`q = 2a+11d` one computes `2ap ≡ -(11) (mod q)`, so (as `q` is odd, since `11d` is odd)

> `s = (q-11)/2`,  giving residues `{(q-11)/2, (q-9)/2, …, (q+11)/2}` — **12 consecutive integers,
> symmetric about `q/2`.**

The residue nearest an integer multiple of `q` is an endpoint, at distance `(q-11)/2`. Hence

> **`min_k ‖(a+dk)·t‖ = (q-11)/(2q) = 1/2 − 11/(2q)`,  so  `M(C) ≥ 1/2 − 11/(2(2a+11d))`.**

*(Verified exactly for all `a ≤ 59`, odd `d ≤ 59`: `s=(q-11)/2` and the witness `=(q-11)/(2q)` every time.)*

**Example (the case S117 could not decide).** `a=2, d=41`: `q = 2·2+11·41 = 455`, `p = 41⁻¹ mod 455 = 111`,
`s = 222`, residues `{222,…,233}` centered at `455/2 = 227.5`, witness `= 222/455 ≈ 0.488 > 1/13`. So
`M({2,43,…,453}) ≥ 222/455`: **not tight.** The maximizer sits at `q = 455 = 5·7·13` — exactly the
"non-obvious modulus" S117 flagged — because `455 = 2a+11d`.

## The corollary: AP rigidity, now fully proved (closes S117's spread case)

The even-`d` case is trivial: `gcd(a,d)=1` with `d` even forces `a` odd, so **every** term `a+dk` is odd,
and `t=1/2` gives `min_k ‖·‖ = 1/2`. So:

> **Theorem (AP rigidity, PROVED).** Every primitive 12-term AP has `M ≥ 1/2 − 11/(2(2a+11d))` (`d` odd) or
> `M = 1/2` (`d` even). This is `> 1/13` for **every** primitive AP except `{1,…,12}` (where `2a+11d=13`
> gives exactly `1/13`). Hence the only 12-term APs with `M = 1/13` are the dilates `c·{1,…,12}`.

*Proof.* Reduce to primitive by dilation-invariance of `M`. `d` even ⇒ `M=1/2`. `d` odd ⇒
`M ≥ 1/2 − 11/(2q)`, and `1/2 − 11/(2q) > 1/13 ⟺ q > 13`. Since `q = 2a+11d ≥ 13` with equality **iff**
`a=d=1`, every primitive AP other than `{1,…,12}` has `M > 1/13`. ∎

This subsumes S117's necessary conditions (ii)–(iv) for odd `d` into a **single** witness, and — crucially —
supplies the spread-case (`d≥17`) argument S117 was missing. The `a=d` rigidity for APs is complete.

## Bonus: the exact loneliness of *every* AP (verified)

The exhaustive exact search (maximize over all `t=j/Q`, `Q ≤ 2·max(speed)`, which provably contains the
maximizer) returned `M = witness` **exactly** in every case — not merely `≥`:

> **`M({a, a+d, …, a+(n-1)d}) = 1/2 − (n-1)/(2q)`  (`d` odd, `q=2a+(n-1)d`);  `= 1/2`  (`d` even).**

Confirmed for `n ∈ {4,6,8,10,12}` and many `(a,d)`. At `a=d=1` this is `1/2 − (n-1)/(2(n+1)) = 1/(n+1)`, the
LRC tight value — so **`{1,…,n}` is the unique tightest primitive `n`-AP for every `n`**, and the S117
consecutive formula `M({a,…,a+n-1}) = a/(2a+n-1)` is just the `d=1` slice.

## Lean (kernel-pure)

`LRCAPCentering.lean` (`namespace LonelyRunner`), built, both `[propext, Classical.choice, Quot.sound]`,
no `sorry`:
- **`centered_block_far`** — the integer core: `q ≥ 13` and `(q-11)/2 ≤ r ≤ (q+11)/2` ⟹
  `q - 11 ≤ 2·|r − q·k|` for every `k` (two-case split on `sign k`, `nlinarith`+`omega`).
- **`centered_block_witness`** — the real witness: if every speed's reduced numerator `N i` lands in the
  centered band mod `q`, then `∀ i m, (q-11)/(2q) ≤ ‖N i/q − m‖`. (In the AP application `N i = (a+d·k)·d⁻¹`.)
Companion to `LRCMod13Blocking` — same shape (omega-checked integer inequality lifted by `abs_div`).

## Status, honestly split

- **PROVED (elementary, exact; Lean kernel-pure):** the *lower* bound `M ≥ 1/2 − 11/(2q)`, hence the AP
  rigidity and the closure of the spread case. This is the owner-facing win.
- **VERIFIED, not yet proved:** the *upper* bound (that the witness is the global maximizer), i.e. the exact
  closed form. Exhaustively confirmed; a proof handle exists (below) but is not finished.

### The upper-bound handle (reflection pairing)

Pair runner `k` with runner `11-k`: their speeds sum to `(a+dk)+(a+d(11-k)) = 2a+11d = q`. So for **any**
`t`, `y_k + y_{11-k} ≡ qt (mod 1)` where `y_k = (a+dk)t`. If all `y_k` lay in the safe arc `(μ, 1-μ)` with
`μ = (q-11)/(2q)` (length `11/q`), then each pair-sum `qt ∈ (2μ, 2-2μ) mod 1`, forcing `‖qt‖ < 11/q`
(**proved** necessary condition — and the witness saturates it, `‖q·p/q‖=0`). Completing the upper bound is
a three-distance argument on the coupling between `at` and `dt` (`gcd(d,q)=1`); left open (HYP-7710).

## Why this is the right frame ("everything is the triangle")

`q = 2a + (n-1)d` is exactly **twice the mean speed of the AP** (`= n·mean` for the symmetric run). The
witness places the AP's *centroid* at `1/2` on the circle and reads off the endpoint distance — the runners
tile the safe arc symmetrically, the extreme runner is the lonely-runner bottleneck. The reflection
`k ↔ n-1-k` (sum `= q`) is the horizontal-leg/complement symmetry of the staircase applied to the speed
line: the AP is self-paired about its center, and `q` is the fixed axis.

Cross-links:
[[the-a-equals-d-rigidity-consecutive-case-proved-and-the-consecutive-loneliness-formula-boxeph-S117]],
[[mod13-blocking-formalized-and-the-exact-n12-tight-locus-is-homogeneous-boxeph-S116]],
HYP-4382 (n=12 tightness), HYP-7705 (S117 spread case, now resolved), HYP-7710 (exact-formula upper bound),
`lrc14_ap_centering_witness_boxeph_S118.py`, `LRCAPCentering.lean` (kernel-pure).
