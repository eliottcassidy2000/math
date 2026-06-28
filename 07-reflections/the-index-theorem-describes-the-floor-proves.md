# The index-theorem describes; the floor proves — testing the topological frame to its honest edge

*kind-pasteur-2026-06-28-S256. The owner: boldly predict how the abstract frames continue, test them at
length, and look for the bigger picture. I did — and the bigger picture is a clean division of labour
between a topological description that unifies everything and an analytic floor that does the actual
work. The boldest bridge between them failed the test, and that failure is itself informative.*

## What converged

My Chebyshev equioscillation frame and mac-mini's Borsuk–Ulam index frame turned out to be one object,
verified from both ends:

- the AP's safety function `f_S(t)=min_s‖st‖` touches its max `1/14` at exactly the `φ(14)=6` **unit
  points** `a/14`, in **3 antipodal pairs** (the complement involution);
- that "3" is the **saddle index `(p−1)/2`** (mac-mini), the **equioscillation count** (me), the
  **#QR mod p**, and the **quadratic Gauss sum** `i√7` (since `7 ≡ 3 mod 4`, it is *imaginary*);
- `p mod 4` is the parity of the index and the reality of the Gauss sum: `p≡1` (real) ↦ Brouwer/SOS,
  `p≡3` (imaginary, n=14) ↦ Borsuk–Ulam.

These are not three coincidences; they are one number `(p−1)/2` seen through equioscillation, topology,
and the cyclotomic field. As a *description* of why the AP is the extremal, the index-theorem frame is
the most unifying thing the project has found. **Predictions P1 (index value) and P3 (Gauss sum) hold.**

## Where it stops being a proof

The frame claims "LRC(2p) ⟺ an index `(p−1)/2` ≠ 0." I tested the forcing (P2), and it has a real gap:

> The index `(p−1)/2 = 3` is **ambient** — it depends only on `n=14`, not on the speed set `S`. But
> the lonely-set measure `meas{f_S ≥ 1/14}` is **S-dependent** (0 for the AP, 0.012 for the loose
> `12→26`, 0.005 for the covering `{1..11,13,84}`). A single ambient integer cannot distinguish a
> counterexample from the AP.

The naive Borsuk–Ulam map (`t ↦` nearest-runner signed position, odd under `t↦−t`) confirms it: its
forced zero is a point where a *runner sits on the observer* (`S=0`) — the **opposite** of lonely. So
the topological forcing, as stated, certifies the AP's own equioscillation saddle (which we already
knew is `1/14`); it does not obstruct an arbitrary `S`. "LRC ⟺ `χ(nerve_S) ≠ 0`" is a faithful
*restatement* — the Euler characteristic of the danger cover is S-dependent and nonzero exactly when a
hole (a lonely point) exists — but proving it nonzero for **all** `S` *is* the conjecture.

The index-theorem is a **"what," not a "how."** It names the target with great precision — the extremal
is the `(p−1)/2`-dimensional equioscillation saddle on `(ℤ/2p)*`, of Gauss-sum type — and it explains
*why n=14 resists the real/SOS methods*: its Gauss sum is imaginary, so the saddle is antisymmetric and
no positive (SOS) certificate can exist. That is genuine insight. It is not a proof.

## The bold bridge that failed (and why the failure matters)

The most exciting prediction was that the topology and the analysis are the *same* dichotomy: that
`p≡3 mod 4` (imaginary Gauss sum) would force the floor's resonance deviation to be **sign-oscillating**
— which would *explain* why my absolute Schur envelope failed (MISTAKE-078) and the signed spectrum sum
was required (Borsuk-Ulam ↔ signed), while `p≡1` (real) would give a sign-definite deviation (SOS ↔
absolute). It would have welded the two halves of the proof into one principle.

It does not hold. The cancellation ratio `|signed SPEC| / |abs SPEC|` came out `0.71, 0.62, 0.22, 0.056`
for `n=10,14,22,26` — *decreasing with p*, and worst at `n=26` (`p=13 ≡ 1 mod 4`), the opposite of the
prediction. The sign-cancellation in the floor is a **general** phenomenon (it grows with the number of
speeds / the resonance complexity — the HYP-2606 `F3` loss), not a signature of the imaginary Gauss
sum. The topological parity and the analytic sign-cancellation are *not* the same thing.

The failure is informative: it says the two frames are **complementary, not identical**. The topology
governs the *symmetry* of the extremal (real vs imaginary, Brouwer vs Borsuk–Ulam); the floor's
sign-cancellation governs the *quantitative* decorrelation, and it is generic. One should not expect a
single arithmetic invariant (`p mod 4`) to carry the analytic content.

## The bigger picture, honestly

> **LRC(14) = [the index-theorem describes the goal] ⊕ [the decorrelation floor does the work].**
> The index-theorem pins the extremal — the `(p−1)/2`-saddle on the unit group, the imaginary `i√7`
> Gauss sum, the antisymmetric Borsuk–Ulam structure — and tells you the proof must be topological/
> signed, not SOS. The actual lower bound `M(S) ≥ 1/14` for every `S` is the S-dependent floor: the
> `u=14t` lift, the Gaussian wide-V decoupling, the multi-far `R'≥c`, and the census. The index is the
> *target*; the floor is the *instrument*. They are not bridged by `p mod 4`.

So the route to close n=14 is unchanged by this session — it is still the uniform `R'≥c` floor — but
the *understanding* is sharper: we now know exactly what we are aiming at (a Gauss-sum-type saddle of
index 3) and why the easy tools can't reach it (imaginary, hence no SOS). The frame is the map; the
floor is the road.

## The pointer beyond

If the index is the imaginary Gauss sum `i√7` and the obstruction is in `ℚ(√−7)` — the Hurwitz prime
the project keeps meeting (the spectral bridge: `29 = 1² + 7·4` is *forbidden* in `ℚ(√−7)`; the
`{2,3,7}` Hurwitz triple) — then the right home for the floor's resonances may be the class group /
the ideal structure of `ℚ(√−7)`, not the rationals. The conjecture to chase: the multi-far floor's
resonance sum, reorganized over `ℚ(√−7)` (where `7` ramifies and `−1` is a non-residue), is *sign-
definite* in the right basis, turning the signed-SPEC obstacle into a positivity — the imaginary Gauss
sum made real by passing to its own field. That would finally weld the topological description to the
analytic floor, through the apex prime's number field rather than through `p mod 4`.

— Related: [[lrc14-thread]], HYP-3249 (predictions), HYP-3251/3252 (tests), HYP-3246/3247 (Chebyshev),
HYP-3132 (multi-far floor), HYP-3125 (Gaussian wide-V), mac-mini S78/S79 (the index frame), THM-250
(the spectral bridge / ℚ(√−7)); `lonely-runner-as-chebyshev-equioscillation.md`,
`the-tournament-spectrum-is-the-object.md`.
