# The filter/witness asymmetry: why no fixed template closes LRC(14)

**kind-pasteur-2026-07-05-S11 (HYP-4137, MISTAKE-110).** A reflection on why the
Q50 / template program was doomed, and the structural reason it kept *looking*
like it should work.

## The two sides of a gap-emptiness argument

Every "the gap is empty" argument for LRC(14) has two sides:

- **The FILTER side** (profile / necessary conditions): a hypothetical
  gap-violator must be covering, spread, compressed, pinned mod every `q ≤ 25`,
  … These are conditions the family must *satisfy*.
- **The WITNESS side** (sufficient condition): the family *has* a loneliness
  witness `t = a/q` (some rational making everyone `2/25`-far), so it is loose,
  not a violator.

mac-mini's pole-necessity theorem (S55) established something real and clean
about the FILTER side: **the profile filters are CRT-ray-periodic.** The floating
7-cluster passes every filter at every frozen scale `S ≡ 1 (mod L)`,
`L = lcm(2..25)`. Filters see only residues mod `L`; freeze those and the family
stays a "profile survivor" at unbounded height. So no filter — no ratio, no
residue condition, no compression bound — can pin the height. That is a theorem,
and it stands.

The Q50 conjecture was the hope that the WITNESS side inherits the same
height-independence: that a bounded modulus (`q ≤ 50`) always carries the witness,
so the gap closes by a finite template check, "no scale, no census-per-height."

## The asymmetry: witnesses do NOT inherit filter height-independence

They don't, and the reason is a clean asymmetry:

> A **filter** at modulus `q` is a condition the family must *satisfy*; freezing
> residues mod `L` freezes all filters `q | L`.
> A **witness** at modulus `q` is a `∃`-statement (some dilation clears everyone);
> its truth depends on residues mod `q` that the profile *does not control* for
> `q ∤ L`, and it can be *destroyed* by a runner the profile *does* permit.

Two mechanisms destroy the bounded witness, both instances of MISTAKE-096
("the lonely denominator is not uniformly bounded; it grows like `log(height)`"):

1. **Free-residue pinning.** For `q ∤ L` (the free moduli 27,29,31,32,37,…), the
   residues are a hidden degree of freedom equivalent to height. A CRT lift pins
   every free `q ≤ N` at once — the witness leaves any fixed free bound.
2. **Composite blocking.** A single runner `≡ 0 (mod L)` is `≡ 0 (mod q)` for
   *every* pinned-only `q | L`, so it is always in the danger band and blocks
   every height-independent witness. And such a runner makes the profile's
   pinning and covering *vacuous* — freeing all the other runners.

Together: for any `N`, a valid non-tight covering compressed primitive family with
no witness `≤ N`. The witness modulus is genuinely unbounded.

## Why it kept looking like it should work

The census (bounded height) is a biased generator: it only realizes shapes with
*small representatives*. On that sample every witness lives at `q ≤ 45`, the
pinned-only witness at `q ≤ 69` — clean, bounded, tantalizing. But the census
cannot produce the `10^22` free-pinned family or the `L`-runner family; those have
no small representative. **An empirical bound is only as uniform as its generator**
— the fifth or sixth time this exact trap has appeared in the project
(MISTAKE-095/096/098/101/102), now at the witness/template layer. The tell each
time: a quantity that is bounded on random/aligned-diagonal samples but is pushed
to infinity by the maximally-composite / CRT-aligned adversary.

## What survives, and the moral

What survives is the object that never pretended to be bounded: the **real-valued**
loose branch `∃ tstar ∈ ℝ` (`TightLooseDichotomy`, `lrc14_of_dichotomy_and_corner`).
Every family here is loose *via a real witness* — the counterexamples only refute
the fixed-*modulus* refinement, never loneliness itself. The honest closure is
either the real witness or a height-dependent bound `q ≤ Q(height) ~ c·ln(height)`
— the two-sided (bounded-magnitude census + large-magnitude analytic) architecture
that MISTAKE-096 already forced on the covering case.

The moral is about where "finite" is allowed to live. Loneliness *is* a finite
fact per family (every family has a witness at *some* finite `q`). But "finite"
does not commute with "uniform": `∀ family ∃ q ≤ Q₀` is false for every fixed `Q₀`,
even though `∀ family ∃ q < ∞` is true. The template program confused the
per-family finiteness (real) with a uniform finite bound (false). The filter side
earns its height-independence by being a `∀`-condition on a frozen residue class;
the witness side, a `∃` over a modulus the adversary can inflate, cannot.
