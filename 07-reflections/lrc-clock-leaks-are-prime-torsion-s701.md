---
source: opus-2026-06-06-S701 (user seed: the n=14/n=15 leaks are prime-torsion)
status: CONFIRMS + SHARPENS the user's observation into THM-421. The LRC clock-witness leaks at
  composite n are LOCKED to the prime-torsion of ℤ/n: a runner that defeats one prime-cofactor's
  clocks but survives the complementary prime sits at a p-torsion element x∈T_p=⟨n/p⟩, CRT-decomposing
  as (nonzero mod p, 0 mod cofactor) — verified n=12..50. n=14: half-turn 7=(1,0)/(2,7)=2-torsion gen;
  n=15: 5,10=3-torsion, killed in the 5-base. Geometric: at t=1/n a p-torsion runner sits at the EXACT
  order-p rotation j/p, margin 1/p. DICHOTOMY: squarefree p ⟹ socle survives mod p, plugged by the
  prime clock t=1/p; prime-power p ⟹ socle ≡0 mod p too, needs the deeper p-adic clock t=1/p^a — the
  prime-power hardness, SAME face as THM-420's caveat (2n-1=27=3³). Hard core (all primes leak) =
  multiples of rad(n) = {0} for squarefree n = exactly THM-398's multiple-of-n residual.
tags: [lonely-runner, LRC, CRT, prime-torsion, clock-witness, half-turn, p-adic, prime-power,
  THM-420, THM-398, THM-403, cyclotomic, worry-set, squarefree, radical, torsion-tower]
---

# The LRC clock-leaks are the prime-torsion of ℤ/n

**Prompt (user):** "that n=14 half-turn leak at coordinate 6 with residue 7 is exactly 1 (mod 2) and
0 (mod 7) — it lives precisely in the 2-torsion subgroup that projects to zero in the 7-runner base.
And for n=15, those order-3 leaks at residues 5 and 10 are exactly the 3-torsion subgroup projecting
to zero in the 5-runner base. The leakage isn't random — it's completely locked into the algebraic
torsion of the composite divisors."

The user read the leak through CRT, and it is exactly right. This turns the *where* of the LRC
obstruction (left open by THM-420's k-clock witness) into a clean arithmetic statement — **THM-421**.

## 1. The CRT reading (confirmed)

For `n = ∏_p p^{a_p}`, CRT gives `ℤ/n ≅ ∏_p ℤ/p^{a_p}`. The **`p`-torsion** `T_p = ⟨n/p⟩` (order
`p`). Every nonzero `x ∈ T_p` is **0 in the cofactor base** `m_p = n/p^{a_p}` and **nonzero in the
`p`-base** — i.e. it defeats the clocks of every prime *but* `p`, while `p` still sees it. The
user's instances, verified:
- **`n = 14 = 2·7`:** the half-turn `7 ↔ (7 mod 2, 7 mod 7) = (1, 0)` — the generator of the
  2-torsion `{0,7}`, dead in the 7-base.
- **`n = 15 = 3·5`:** `5 ↔ (2,0)`, `10 ↔ (1,0)` mod `(3,5)` — the nonzero 3-torsion `{0,5,10}`, both
  dead in the 5-base.

"Locked into the algebraic torsion of the composite divisors" is literally true: the leak set is
`⋃_p (T_p\{0})`, the single-prime-surviving residues (n=14: `{2,4,6,7,8,10,12}`; n=15:
`{3,5,6,9,10,12}`).

## 2. Why the torsion is *safe*, not dangerous — the geometric margin

The crucial point: a `p`-torsion runner is exactly the one the **complementary prime plugs**. At the
full clock `t = 1/n`, `x ∈ T_p\{0}` sits at `x/n = j/p` — an **exact order-`p` rotation**, distance
`≥ 1/p` from the origin. So the leak through the cofactor's clocks is *caught* by the `p`-direction
with margin `1/p`. CRT turns "leak in one base" into "controlled in another." This is the structural
reason THM-420's k-clock witness works on composite `n`: each prime's blind spot is another prime's
clear shot.

## 3. The dichotomy that names the hardness: squarefree vs prime-power

Computing the *least controlling clock* of the socle generator `n/p = p^{a_p−1}m_p` splits cleanly
(verified n=12..50):
- **squarefree `p` (`a_p=1`):** `n/p = m_p` is coprime to `p` ⟹ survives mod `p` ⟹ plugged by the
  **prime clock `t=1/p`**. Clean. (n=14, 15, 21, 30, 33, 35.)
- **prime-power `p` (`a_p≥2`):** `n/p = p^{a_p−1}m_p ≡ 0 (mod p)` ⟹ the socle *also* dies mod `p`, so
  `t=1/p` sends it to the origin; control needs the **deeper `p`-adic clock `t=1/p^{a_p}`**. (n=12,
  18, 20, 45, 50 — socle `survive mod p = False`, verified.)

> This is the **same prime-power obstruction** as THM-420's honest caveat. There, the shell modulus
> `2n−1` for `n=14` is `27 = 3³`, a pure prime power with no coprime cofactor to plug the shell-leak.
> Here, on the *clock* side `ℤ/n`, prime-power factors force descent down the `p`-adic tower
> `C_p ⊂ C_{p²} ⊂ ⋯`. **Both faces of LRC — the multiplicative clock `ℤ/n` and the additive shell
> `ℤ/(2n−1)` — hide the difficulty in prime-power torsion towers.** Squarefree moduli are the easy
> regime (every leak coprime-plugged); prime powers are where the worry-set lives.

## 4. The hard core recovers THM-398

The runners that leak at *every* prime at once = `{x ≡ 0 (mod rad n)}` (multiples of the radical).
For **squarefree `n` this is `{0}` alone** — a runner is all-prime-dangerous iff `v ≡ 0 (mod n)`,
which is exactly **THM-398's multiple-of-`n` residual**. So the THM-398 residual *is* the CRT hard
core, and it collapses to the single origin precisely when `n` is squarefree. For prime-power `n` the
hard core grows (n=12: `{0,6}`), seeding extra worry beyond the origin.

## 5. Tournament reading (standing directive)

The order-`p` rotation clock `t = j/p` is the **`C_p` cyclotomic round tournament** (THM-403); `T_p`
is the lattice of its sub-clocks. The leak-localization says the worry-set's cyclotomic witnesses are
indexed by the **prime-torsion lattice of `ℤ/n`**, and the prime-power depth (§3) is the height of
the cyclotomic round-tournament tower `C_p ⊂ C_{p²} ⊂ ⋯`. The squarefree/prime-power split is the
tournament statement that a product of distinct-prime round tournaments is "transversal" (CRT-split,
each leak plugged), while a prime-power round tournament is an indecomposable tower.

## 6. Honest status

- **PROVED:** the CRT localization (1), the geometric margin `1/p` (2), the squarefree/prime-power
  dichotomy (3), the hard-core = multiples-of-`rad n` = THM-398-for-squarefree (4). All elementary
  CRT once seen; the content is the *identification*, which the user supplied.
- **VERIFIED:** n = 12, 14, 15, 18, 20, 21, 30, 33, 35, 45, 50 (`…s701.py`): geometric localization
  holds for all; squarefree clean-plug holds; prime-power socle dies mod `p`.
- **What it gives the program:** a precise map of *where* the k-clock witness (THM-420) can leak, and
  a clean statement that LRC's residual difficulty is concentrated in prime-power torsion towers on
  *both* the clock and shell faces — unifying THM-420's caveat, THM-398's residual, and THM-403's
  cyclotomic witnesses under one CRT/torsion picture.
- **Not a new bound on LRC(14):** this localizes and explains the obstruction; it does not resolve
  the prime-power core. **Next:** the `p`-adic clock descent for prime-power factors (does
  `t=1/p^{a}` + the shell witness cover the prime-power hard core?) — the concrete `n=14` question is
  whether the `3`-adic tower `C_3⊂C_9⊂C_27` (the `2n−1=27` shell) plus the clock `ℤ/14` torsion close
  the residual.

**Artifacts:** `04-computation/lrc_torsion_localization_s701.py` (+`.out`). Theorem: **THM-421**.
Builds on THM-420 (k-clock witness, S700), THM-398 (multiple-of-`n` residual), THM-403 (cyclotomic
witnesses); relates to HYP-2281 (shell-transversality) and the S708 `3`-tower work. New: **HYP-2285**.
