# Changing the core: what the neutral flips can and cannot do for the lanes p and 7p

The reduction (THM-4094): Chan–Pak is the two lanes of odd primes `p ≠ 7` and `7p`, each needing a strong realiser
(`7p` from a strong component of count `7p`). The near-join equation `p = h + 2^k B` (HYP-9453; wave 56's placements)
spends a base of count `h` and one of its responses `B ∈ 𝔅(h)`; the only ways to change the supply are:

| move | effect on `(h, 𝔅(h))` | cost |
|---|---|---|
| a free path run (transitive padding, `σ` or `ς` of THM-5266) | `h` fixed, every `B ↦ 2B` (or `B_src = P(2)`) | one more run; no new odd part |
| a nontrivial neutral flip inside the core | `h` fixed, new responses with new odd parts | none — same order, same count |
| a general contact or the quotient | `h` changes | leaves the fibre |

So "more neutral choices without another free path run" is exactly the nontrivial neutral graph of the fibre
([NEUTRAL.md](NEUTRAL.md)). Three consequences:

1. **The rigid fibres cannot be changed.** For `h` in `{1, 3, 5, 11, 13, 17, 19, 31, 35, 39, 47, 53, 55, 57, 63, 67, 85,
   95, 113, 115, 121}` every core is isolated: all its count-preserving flips return to itself, and every other realiser
   is a padding (Lemma, NEUTRAL.md §1). Their response sets are the dyadic shifts of finitely many odd parts, fixed for
   all orders — for 1 and 3 just `{1}` and `{3}`, which is why `7 = 1 + 2·3` and `21 = 3 + 2·9` are out of reach
   ([HOLES.md](HOLES.md)) and why 11, 13, 47, 67 are permanent near-join exceptions (THM-5263): a target that needs a
   response a rigid fibre does not have cannot be helped by re-arranging that fibre's core.
2. **Flexible fibres are the bases to use.** A target `p` or `7p` should be written `h + 2^k B` with `h` in a fibre of
   large neutral graph (37, 29, 101, 109, 89 below 125; 547, 657, 523 at order 9), where the odd parts of `𝔅(h)` are
   many; wave 56's six placements of primes near `10^6` all used bases of orders 15–19, i.e. fibres far above the rigid
   range. This is a selection rule, not a proof: whether every `p` and `7p` is reached is HYP-9453's density question.
3. **Rigidity does not decide hardness.** 27 and 9 are flexible but fertile-only; 41, 43 are flexible yet
   tower-blocked (THM-5260); the sterile-only counts are all rigid but so are 31, 39, 55, 57, 85, 95, 115, 121, which are
   ordinary. The correlation with `|𝔅(h)|` is 0.61, not 1.

**Relation to the concurrent packets.** `one_port_neutral_continuations_20261117` keeps a `+2` pair neutral under appended
addresses on a tape with at most two runs (their "neutral" is a property of the continuation); here "neutral" is a
property of the core itself. They compose: a flexible core offers both more responses and, for each rooted pair taken
from it, the same continuation alphabet `S, N`. `ready_core_growth_and_mod18_20261117` shows bounded deletion counts
cannot serve almost all `7p`; the flexible fibres at order 9 are where the deletion count starts to grow.

Chan–Pak OPEN.
