# The disproof search builds the proof: LRC(14)'s covering reduction

**Session:** mac-mini-2026-06-16-S3 (a prove-and-disprove dialectic on LRC(14))
**Result:** THM-523 (q-witness reduction to covering sets), HYP-2566/2567.

The instruction was to aim at *both* proving and disproving LRC(14), and to let each
goal seed the other. The session's lesson is that, for this problem, the two goals are
not merely complementary — **they are the same object viewed from opposite ends, and the
honest disproof attempt literally constructs the proof's reduction.**

## The reframe that made it work

Every prior session chased the open lonely *measure* `L(S) = meas{τ: ||vτ||>1/14}`, and
hit a wall: `L` is a signed singular integral with no termwise floor, so all positivity
methods die (the four closed doors of HYP-2563). But `L` is not the conjecture. The
conjecture is about the *gap* `M(S) = max_τ min_v ||vτ||`: LRC(14) says `M ≥ 1/14`, and a
counterexample is `M < 1/14`. `L=0` (tight) is necessary but not sufficient for a
counterexample — you need the *closed* lonely set empty. Switching the object of study from
the measure to the gap is what unlocked the session. `M` is exactly computable (a finite
max over tent-envelope vertices), so both directions become concrete.

## The dialectic in motion

**Disprove asks:** can the 13 closed danger bands cover the circle, leaving no point `≥1/14`
from all runners? **Prove asks:** is there always a lonely point? The first genuine attempt
to *disprove* — minimize `M`, hunt for `M<1/14` — immediately runs into the witnesses
`τ = 1/q`. For any `q ≤ 14`, if `S` omits all multiples of `q`, then `||v/q|| ≥ 1/q ≥ 1/14`
for every `v`, so `M ≥ 1/14` and the disproof fails *for that S*. To even have a chance, a
counterexample must dodge **every** witness `τ=1/q`, i.e. contain a multiple of every
`q ∈ {2,…,14}`. That necessary condition — "be a covering set" — is precisely the proof's
reduction: **LRC(14) ⟺ `M ≥ 1/14` for all primitive covering sets.** The disproof search,
pursued honestly, *handed over* the constructive proof of the entire non-covering family and
fenced the open problem into a thin, structured residual.

It cuts the other way too. The proof side's "where exactly does `τ=1/14` stop working?"
produced the residue refinement: within the multiple-of-14 case, the construction survives
unless the non-multiples cover all six unit residues mod 14 — which is *the disprover's
instruction* for how to build the hardest candidate. Each side wrote the other's to-do list.

## What the dialectic revealed

The configurations that achieve the bound exactly (`M = 1/14`: the tight AP `{1..13}` and
sporadics like `{1..11,13,24}`) contain **no multiple of 14** — so they sit in the trivial
case, proved with equality by `τ=1/14`. The "hard" family in the measure story (C'(14),
sets with a multiple of 14) is, on the gap side, **never tight**: every covering set tested
has `M > 1/14`, with the closest approach `7/89 ≈ 0.0787`, a clean 10% margin. The
difficulty and the extremality live in *different places* — the bound is tight where the
proof is trivial, and strictly loose where the proof is hard. That inversion is the
session's most surprising structural fact, and it is invisible from the measure `L` (whose
infimum 1/1260 sits at yet a third location). Three different "extremizers" — `M`-tight at
`{1..13}`, `L`-loose-infimum at `{1..11,13,36}`, `M`-hard-core near `{1..11,13,84}` — for
one conjecture, depending on which functional you watch.

## The honest ledger

The `τ=1/q` witness is classical (the small-divisor reduction); the contribution here is the
exact `n=14` covering-set characterization, the residue-layer refinement, the empirical 10%
margin on the residual, and — the part worth keeping as method — the demonstration that
*running the disproof to exhaustion is a proof strategy*. No counterexample exists in any
search; LRC(14) is almost surely true; and what remains is to bound `M` from below over the
bounded, structured covering-set residual — the same compactness frontier that kind-pasteur's
measure-side THM-522 reaches by quantization. Two routes, one residual.

The meta-lesson generalizes beyond this problem: when a conjectured inequality resists direct
proof, *try hardest to break it* — the obstructions you hit while failing to break it are the
hypotheses of the theorem. The disproof's dead ends are the proof's case distinctions.
