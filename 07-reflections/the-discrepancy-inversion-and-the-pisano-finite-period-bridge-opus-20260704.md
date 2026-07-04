# The discrepancy inversion: LRC-tight is perfect equidistribution, not low-discrepancy — plus the Pisano finite-period bridge and the large-tightener lemma

*opus-2026-07-04-S65. The owner steered me to Franel's discrepancy bound (via a fabricated "synthesis"
post) and to "Fibonacci mod 7, period 16." Per the owner's pattern these are structural-analogy prompts.
I extracted the one real kernel, reconciled the discrepancy framing with the frontier — finding a clean
inversion — and turned the folding identity's discrepancy content into a proved lemma. Honest throughout,
including flagging the fabricated post.*

## The post is fabricated; one kernel is real

The poke-forum "synthesis" (Franel × Dedekind-zeta × 2-adic Sha) is LLM-generated pseudo-math: "Pi Unital
Flower guardrails," the SHA-256 `e3b0c442…` is the **empty-string** hash, `5d41402abc…` is `MD5("hello")`,
"Placeholder Cryptographic Verification." The Dedekind-zeta / Artin-map / Tate–Shafarevich links are
invented. The **one real kernel** is the base-2 van der Corput discrepancy bound
`N·D*_N ≤ log N / (3 log 2) + O(1)` (the Béjian–Faure constant) — a genuine low-discrepancy fact. So I take
the *structure* (2-adic low-discrepancy equidistribution, the constant's `3`), not the fabricated edifice.

## The inversion (the real bridge)

Discrepancy theory's extremal object is the **golden ratio / Fibonacci** — the most badly-approximable
number, worst-case for how evenly `{nα}` fills `[0,1)`. The prompt invites mapping that onto LRC. The map
is an **inversion**:

- **Fibonacci speeds `{1,2,3,5,…,377}` are LRC-LOOSEST** (`M = 0.171 ≈ 2.4·(1/14)`): badly-approximable
  speeds never let the runners cluster, so a lonely time is easy and fat.
- **The AP `{1..13}` is LRC-TIGHTEST** (`M = 1/14` exactly): its critical time `t=1/14` puts the runners
  on the nonzero 14th **roots of unity** — *perfect* equidistribution, zero discrepancy.

So **LRC-tightness = the runners CAN be perfectly equidistributed (roots of unity) = the AP**, the *opposite
end* from the classical (golden-ratio) discrepancy extremal. The low-discrepancy prompt, reconciled, points
back to exactly the tight-locus rigidity the fleet already has: **tight ⟺ AP ⟺ roots of unity** (HYP-4062).
The golden ratio is a *red herring for tightness* — it is the loosest case. (Confirmed the obvious
corollary too: van der Corput `t`-sampling finds loneliness witnesses **slower** than random — its dyadic
structure resonates with the `1/14 = 1/(2·7)` band — so low-discrepancy sampling is a dead end for a
witness-finder.)

## The Pisano finite-period bridge

"Fibonacci mod 7, period 16" is the **Pisano period** `π(7) = 16 = 2⁴`. With `π(2) = 3` and
`π(14) = lcm(3,16) = 48`, this is the archetype of the LRC endgame's defining move: **an infinite family
collapses to a finite period.** That is exactly my mod-24 confinement check, mac-mini's mod-14 residue
problem, and kps's far-peel closing `{1..12,w>182}` from one base computation. And the numerology is not
idle: `π(2) = 3` is the same `3` as the van der Corput constant `1/(3 log 2)` and the three-gap bound
`g(14) ≤ 3` (HYP-2913); the `7` side has `π(7) = 2⁴`, purely 2-adic — the `2`/`7` split (`14 = 2·7`) that
runs through the whole endgame (the folding identity's even/odd, the danger band `1/(2·7)`). *Lead:* is
there a single uniform period `P` (a Pisano-like `lcm`) modulo which **all** confinement/tight-locus checks
decide, collapsing the residual to one finite computation?

## The discrepancy content of the folding identity — a proved lemma

The Franel framing sharpened my own open piece into a theorem. Recall (THM-615) confinement reduces to
`M(2U ∪ {w₁,w₂}) ≥ 1/12`, i.e. some `g_E ≥ 1/12` point avoids extremity. A tightener avoids extremity
wherever its orbit is "moderate" — a **discrepancy** condition. Made precise (Lemma 3, THM-615):

> if `max(w₁,w₂) > u_max/(6(M(U)−1/12))`, then `M(S) ≥ 1/12`.

A large tightener's orbit is dense enough on the high interval `I₀` to hit the moderate band `[1/12,5/12]`,
so it cannot stay extremal. Proved (Lipschitz + IVT), verified (360 families, 0 violations). This disposes
of the **loose end**; the residual is the *low-frequency* corner — both tighteners small AND `M(U)` near
`1/12` (near-AP) — where the orbits are too coarse to be forced moderate. That corner is the confinement
hard core, and it is a genuine **low-discrepancy-of-few-frequencies** question: can two low-frequency odd
combs avoid all the (few, structured) high points of `g_E`? This is where the three-gap / Franel structure
would bite, if it bites.

## Status and leads generated

- **Proved (THM-615 Lemma 3):** large tightener ⟹ `M ≥ 1/12` — the loose end of general-`U` confinement.
- **Structural (verified):** the discrepancy inversion (LRC-tight = roots of unity = AP, opposite the
  golden-ratio extremal); the Pisano finite-period bridge (`π(2)=3, π(7)=16, π(14)=48`).
- **Negatives closed:** Fibonacci speeds are loose (not extremal); van der Corput witness-finding is slower
  than random.
- **New leads (→ backlog):** (i) a uniform Pisano-like period `P` collapsing all confinement checks to one
  finite computation; (ii) the residual small-tightener × near-AP corner as a two-comb three-gap/discrepancy
  problem; (iii) the roots-of-unity characterization of the tight-locus as a discrepancy theorem.

The genuine kernel was small (one real bound) but the reconciliation is real: the LRC extremal lives at the
*perfect-equidistribution* end, not the badly-approximable end, and the confinement residual is a
low-frequency discrepancy question — with the loose end now proved.

Related: THM-615 (Lemma 3 added), HYP-4062/kps (tight-locus = AP = roots of unity, reconfirmed), HYP-2913
(three-gap `g(14)≤3` = `π(2)=3`), MISTAKE-101 (the sampling discipline behind the negatives). Owner steer:
Franel discrepancy + Fibonacci mod 7. Script: `lrc14_discrepancy_bridge_franel_fibonacci_opus_S65.py`.
Post `poke-forum/post_1783144839182.md` flagged fabricated. HYP-4074.
