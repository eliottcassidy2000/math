---
source: claude-2026-06-06-S679
status: clean sieve reduction for prime n=19 (proved) + a cross-thread reframe (apex = involution fixed point); LRC@19 not proved (sieve-covered core open)
tags: [LRC, n19, prime-frontier, divisibility-sieve, negation-involution, fixed-point-residual, signed-lrc, redei, reframe, exploration]
---

# LRC@19, the prime frontier: a clean sieve reduction, and the reframe it points to

Prompt: spend a long session on LRC@19; if you hit a wall, explore creatively
until the connections give a reframe. That is exactly what happened.

## Why n=19, and the clean structure

`n=19` (prime) with pair-sum modulus `C=2n−1=37` (prime) is the "no side-channel"
frontier (S678): bigger than `n=14` but *cleaner* — there are **no sporadic
walls** (`3∤37`, so no `V*` doubling; verified the AP `{1..18}` is the unique
single-swap-tight). A proof mechanism that can't say something sharp here is
overfit to `Res₂₇`. The AP's tight witnesses are **all** of `t=j/19`
(`j=1..18`), each binding the antipodal pair `{a, 19−a}`.

## The clean sieve reduction (proved)

At `t=a/q` with `gcd(a,q)=1`, `q∈{2..19}`: `‖v·a/q‖ = dist(va mod q)/q`, and since
`q ≤ 18 < 19`, this is `< 1/19` iff `va ≡ 0 (mod q)` iff `q ∣ v` (`a` a unit). So

> **`t=a/q` is a lonely witness (all `‖vt‖ ≥ 1/19`) unless some speed is divisible
> by `q`.**

Hence any LRC@19 counterexample is **sieve-covered**: for every `q∈{2..19}` some
speed is divisible by `q`. Verified with zero failures over ~34k non-sieve-covered
random 18-sets (each has an explicit lonely `t=a/q`). This is Rosenfeld's
divisibility sieve, at its cleanest for prime `n` (every `j` invertible). The 8
primes `≤19` each force their own multiple among the 18 speeds.

## The wall

The sieve-covered family is infinite, and clearing it is exactly what the 2026
`n≤13` proofs needed the **polynomial method** for. A direct proof of LRC@19 is
out of reach in one session. So — explore.

## The reframe (the gem in the connections)

Pulling on the prime structure and three repo threads (signed-LRC S674, the
sin/cos even–odd wall S655, my own `κ`-even / negation finding S605) gave a clean
reframe:

- At prime `n`, the AP's binding pairs `{a, n−a} = {a, −a} (mod n)` are exactly
  the **orbits of the negation involution** on `𝔽_n`.
- A sieve-covered config must contain a **multiple of 19** (the apex), i.e. a
  speed `≡ 0 (mod 19)` — and `0` is the **unique fixed point of negation**.
- [verified] the apex **kills every division-point witness** (`‖19·j/19‖ = 0`
  for all `j`), and the antipodal pairs `{a,−a}` prevent local recovery (any
  `ε ≠ 0` drops one member of the binding pair below `1/19`). So loneliness must
  be recovered **off the `j/19` grid** (e.g. `V={2..19}`: `M = 2/21` at
  `t = 20/21`).

So the **hard core of LRC@prime-n is the negation-involution fixed-point
residual** — the apex sitting at `0 ∈ 𝔽_n`, breaking the clean `±`-pairing.

### Why this is a *gem*: it is the same residual everywhere

The fixed-point-of-an-involution residual recurs across the whole complex:

| problem | involution | hard residual = its fixed points |
|---|---|---|
| LRC@prime-n | negation `a ↦ −a` on `𝔽_n` | the apex (`0 mod n`, the sieve multiple) |
| Rédei / H-gaps | reversal `T ↦ T^op` | self-converse tournaments (`{7,21}` residual) |
| `V*` (signed-LRC S674) | sign gauge `v ↦ −v` | the nonunit seam (`3+24=27` zero pair-clock) |
| resonance sum (S605) | `c ↦ −c` on `L(V)` | `κ`-even ⇒ sign-preserving (no cancellation) |

In each case the "easy bulk" is the paired orbits (killed by the involution's
symmetry), and the **hard part is the fixed-point set**. And the repo's own
solved face (Rédei) shows the cure: a *secondary, sign-reversing* involution on
the fixed-point set (the twisted-involution program, S606). For LRC@prime-n this
says: **find a second involution acting on the apex/sieve-covered residual** —
the apex runner against the `±`-paired rest — that recovers an off-grid witness
uniformly. That is the well-posed next object, and it is the *same* object the
H=21 proof needs (a twisted involution on the self-converse residual).

## Honest status

- **Proved:** the prime-`n` sieve reduction (counterexamples are sieve-covered);
  the apex = negation-fixed-point and its killing of the division witnesses.
- **Not proved:** LRC@19 itself (the sieve-covered/off-grid recovery — the
  polynomial-method core).
- **Reframe delivered:** the hard core is a *fixed-point-of-involution residual*,
  unifying LRC@19, Rédei/H-gaps, signed-LRC `V*`, and the `κ`-even obstruction,
  with the twisted-involution (S606) as the shared cure.

## Next

1. **Build the twisted involution on the apex residual.** For a sieve-covered
   `V = V' ∪ {19w}`, pair times/speeds by a second involution (e.g. the unit
   `u` with `u² ≡ 1 mod 37` composed with negation) and try to produce an
   off-grid lonely witness uniformly. If it works, it works for n=14's
   self-converse residual too (the shared cure).
2. **Off-grid witness density.** Quantify the measure of the lonely set once the
   apex is removed (`V'` has 17 runners at gap `1/19 < 1/18`, lots of room); show
   the apex's tiny blocking arcs cannot cover it — a measure/positional bound,
   the prime-clean analogue of the even-fold (HYP-2065).
