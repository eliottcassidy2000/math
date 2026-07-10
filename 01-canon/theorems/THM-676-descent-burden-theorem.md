# THM-676 — The descent-burden theorem: parity forces ≥ 11 distinct half-sum descent moduli on every 13-set, the burden is minimized EXACTLY by arithmetic progressions, and composite pair sums cannot escape

**Status:** PROVED ((i)–(v), elementary, proofs below) + VERIFIED (0 violations: Freiman
lower bound + equality characterization over 200k random 7-sets; composite proper-divisor fact
exhaustively to 20,000; real covering sets at cap 250 carry 19–62 forced moduli, far above the
floor). Plus the HYP-5765 stress-test result (Part B): the conjecture SURVIVED adversarial
hunts at caps 150/250/400 with the dodging deficit GROWING in scale (min open 3 → 5 → 7).
**Source:** mac-mini-2026-07-09-S65 (cont. 4).
**Depends on:** THM-668 (C2), THM-672 (window conditions), THM-674 (domination form).

## Theorem

Let `S` be any set of 13 distinct positive integers.

> **(i) (parity pigeonhole)** Some parity class `A ⊆ S` has `|A| ≥ 7`, giving ≥ C(7,2) = 21
> pair sums that are EVEN.
>
> **(ii) (Freiman burden)** The pairwise sums of any 7 reals take ≥ 2·7 − 3 = 11 distinct
> values; hence the half-sums `h = (v_i + v_j)/2` over the majority parity class take ≥ 11
> distinct integer values, each a PROPER divisor of its own pair sum `q = 2h`, with `h > 14`
> whenever `q > 28`.
>
> **(iii) (the burden)** Consequently every 13-set whose majority-parity pair sums exceed 28
> carries ≥ 11 distinct proper descent moduli `h > 14`; a full descent-dodger must BLOCK every
> one, paying THM-672's torsion conditions (`h ≤ 28`) or THM-674's ledger/domination
> (`h ≥ 29`: ≥ ⌈m/j⌉ occupied ±-classes in a T-dominating pattern; prime `h` exactly `T·D = G`).
>
> **(iv) (AP rigidity of the minimum)** The burden equals 11 iff the majority-parity class is
> an ARITHMETIC PROGRESSION; every non-AP majority class forces ≥ 12. (And for an AP with
> difference `d` — necessarily even — the 11 moduli themselves form an AP of difference `d/2`.)
>
> **(v) (no-primality escape)** Every COMPOSITE pair sum `q > 196` has a proper divisor
> `> 14` (write `q = pb` with `p` its least prime factor; then `b ≥ √q > 14`). Even pair sums
> escape for no scale (`h = q/2 > 14` once `q > 28`); only ODD PRIME pair sums carry no proper
> descent modulus.

*Proofs.* (i) 13 = 7 + 6. (ii) For `a_1 < ⋯ < a_7`: the chain
`a_1+a_2 < a_1+a_3 < ⋯ < a_1+a_7 < a_2+a_7 < ⋯ < a_6+a_7` exhibits 11 distinct sums; halves of
even sums are integers; `h = q/2` divides `q` properly and `h > 14 ⟺ q > 28`. (iii) THM-668's
C2 + THM-672/674. (iv) Classical Freiman equality: 11 distinct sums force
`a_1 + a_{i+1} = a_2 + a_i` for all `i` (every sum must lie in the chain), i.e. constant
differences; conversely an AP achieves 11. For the moduli: `h_{ij} = a_1 + (i+j−2)(d/2)` ranges
over an AP. (v) `q` composite ⟹ `q = p·b`, `p ≤ √q ⟹ b ≥ √q > 14` for `q > 196`; `b < q`
proper. ∎

## Part B — the HYP-5765 stress test (run BEFORE attempting its proof; honest method note)

Before attempting to prove HYP-5765 I hunted for counterexamples to it AS STATED, because its
clause (b) covers only PRIME moduli ≥ 29 and the cont-3 dodger had been captured at the
composite 49 — a suspected formulation gap. RESULT: **the hunt failed to refute it, with a
growing margin** — adversarial covering sets at caps 150/250/400 could not push the count of
open [(a) window ∪ (b) prime ≥ 29] descents below 3/5/7 respectively. The dodging deficit
GROWS with scale: larger sets have more prime divisors ≥ 29 among their 91 sums, and blocking
at such a prime demands THM-674 domination, whose incidence collapses (13.8% at 29 → 0.25% at
41 for random covering sets). The suspected composite gap is real in principle but empirically
never needed — composite `k ≥ 29` moduli are extra insurance on top of (b), not a repair.

**Where the conjecture's proof must come from (the quantified base camp):** by (ii)+(iv), any
dodger's burden is ≥ 11 blocked moduli, achieved only by parity-AP structure; real covering
sets carry 19–62. Each blocked `h ≥ 29` commits ≥ ⌈m_h/j_h⌉ ≈ 7 of only 13 residue classes to
a rigid dominating pattern (THM-674, tight at 31/37/41), and the patterns at distinct moduli
are CRT-independent demands on the same 13 integers. The three rigidity handles now available —
E3-global (LEM-015/LRCSchurRigidity), E_H-spectral (klein THM-673 C2-stability), parity-Freiman
local (this theorem's (iv)) — all say the only way to afford the burden is to BE a near-dilated
interval, and primitive covering sets cannot be exact dilated intervals (klein-S211:
primitivity thins the branch; the AP itself is non-covering). **The missing piece is a
STABILITY version: near-minimal burden ⟹ near-AP ⟹ explicit perturbation family ⟹ finite
check.** That is Freiman-3k−4 territory — the parked BSG lead (HYP-5682), now with a concrete
target: majority-parity classes with ≤ 12–13 distinct pair sums.

**Verification:** `04-computation/lrc14_hyp5765_test_and_burden_macmini_S65cont4.py` (+ .out).
**Related:** THM-668/672/674, HYP-5765/5730, LEM-015, klein THM-671/673, kps LRCSchurRigidity +
LRCPairSumDispatch, opus-S186 (Lemma A retired — hfloor now proved-route), HYP-5682 (3k−4).
