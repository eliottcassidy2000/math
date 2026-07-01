# The multiplier scales its own distance

*klein-2026-06-30-S52. A reflection on HYP-3763 (large multiples forced) and HYP-3745 (the fused-radius trap), read through Steinhaus.*

The covering-min program keeps returning to one stubborn fact: among **covering sets** — the ones that
carry a multiple of every `q ∈ {2,…,n}`, the only sets the covering-min ranges over — the dense core
`{1,…,n-2}` cannot be cheated. Drop a small speed `k` and, for `n` large enough, the minimum
runner-gap `M` rises above `n/Φ₆(n)`. (The qualifier matters: the tight `1/n` minimizers like
`GW={1..11,13,24}` sit *below* `n/Φ₆`, but they miss a multiple of `n` and so are not covering sets at
all — a search that forgets this reports the lonely floor as a phantom "escape." And the rise is only
eventual: at `n=8` a covering escape `{1,2,3,4,5,7,24}` still beats the construction, `4/29 < 8/57`.)
For months the rest was "verified, `n=14`" — a search that collapsed to one set. This session it split
cleanly into a part that is now *proved* and a part whose *mechanism* is finally legible.

The proved part is almost embarrassingly elementary, and it is where the number theory bites. If no
runner is a multiple of `k`, stand the observer at time `1/k`. Every runner `s` sits at `s/k`, and
since `k ∤ s`, it is at least `1/k` from the observer. So `M ≥ 1/k`. Now compare `1/k` to the target:
`1/k > n/Φ₆(n)` is exactly `Φ₆(n) = n²-n+1 > kn`, which for `k ≤ n-2` is `n²-n+1 > n²-2n`, true with
room to spare. So a covering that means to stay at `n/Φ₆` *must* contain a multiple of `k`. And if `k`
itself is gone and `k > (n-2)/2`, the smallest multiple left is `2k`, already past the whole core. **A
large speed is forced** — no primality, no search, just the `q`-witness and the sixth cyclotomic
value as the threshold. The hexagonal modulus `Φ₆` is not decoration here; it is the precise line
below which resonance `k` cannot be missed.

The part that was murky is *why the forced large speed then fails*, and the answer is a scaling law.
Call the forced multiple `κ = k·c`. It kills resonance `k` honestly: at `D = k` it is `≡ 0`, sitting
on the observer. But `k` did more than kill one resonance — at certain moduli `D` it sat a hair from
the observer and covered it tightly. Take the `n=14`, `k=12` escape: at the prime `D = 89`, rotation
`a = 37`, the vacated speed `12` would land at `12·37 ≡ -1` — distance **1**, snug against the
observer. Remove it and the nearest runner retreats to distance **7**; the hole yawns open to `7/89`,
above `14/183`. Where is the killer `κ = 84 = 12·7` that was supposed to stand in for `12`? Its image
is `84·37 = 7·(12·37) ≡ 7·(-1) = -7` — distance **7**, not 1. The multiplier `c = 7` that let `84`
kill the resonance is the very multiplier that **scales its distance from 1 to 7**. `κ`'s image is
always `c` times `k`'s image; at the modulus where `k` was at distance 1, `κ` is at distance `c`, out
at the edge of the hole instead of plugging it.

That is the fused-radius trap (HYP-3745) in one line, and it is pure Steinhaus. The core
`{1,…,n-2}` is an arithmetic progression, so at any time its runners are a rotation orbit with only
three gap-lengths (the three-gap theorem, HYP-3762) — a maximally even cover, no slack. Remove one
runner and its two adjacent gaps merge into a double gap; that merged gap is the hole. The one thing
that could refill it is a runner landing back on the vacated slot, i.e. a speed `≡ k (mod D)`. But
`κ = kc ≡ k` only when `D ∣ k(c-1)` — a thin set of moduli that never includes the near-resonance
prime where the hole is deepest. The multiple is congruent to `k` mod `k`, which is why it kills the
resonance; it is *not* congruent to `k` mod the hole's modulus, which is why it cannot fill the hole.
Killing the resonance and filling the hole ask for incompatible congruences of the same integer.

So the loophole everyone reaches for — "just use a big multiple, or a huge CRT-tuned speed" — closes
for a structural reason, not a computational one. A large speed can choose *which* rotations it covers
(CRT lets you place its residues anywhere) but not *how many* (one speed, `≤ 2r+1` per modulus), and
the one placement it cannot buy is the tight distance-1 coverage that only the small original had. The
big number's very bigness — the multiplier it needs to be a multiple — is what pushes its shadow off
the hole. The core is irreplaceable because each of its members does a job whose reward scales with
smallness, and a substitute is forced to be large.

There is a lesson here about where the deep holes live. They are not at the small moduli you would
check first, nor at the band primes the killer is tuned to cover; they surface at a *large prime*
`D` (89, 157, 269 …) where the vacated small speed had a near-resonance. "Large primes are forced" in
two senses at once: the covering set is forced to contain a large speed, and the surviving hole is
forced out to a large prime modulus. Both are the same coin — the small speed's near-resonances are
scattered among the large primes, and removing the speed lights one of them up.
