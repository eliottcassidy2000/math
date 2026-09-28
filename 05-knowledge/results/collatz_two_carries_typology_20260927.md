# The two carries and the pointwise form of memorylessness (D24), the dual obstructions (D25), in-degrees (D26), and a typology of arithmetic iterations by invariant, drift and memory, with STICKY proposed as a control

**Session:** opus, `collatz-posets-zeta5-20260927` (tenth note), 2026-09-27.
**Owner's directive:** "pursue D24 and D25 toward a Collatz proof and consider
D26 and other deep principles you can discover as they arise."
**Inherits (cited):** the ninth note (the Lehmer five; drivers lock in for
free, shadows are priced; the measured transition matrices), the fifth note
(price sheet), the first note (Apéry forms `A_j = 3^j n + S_(j-1)`,
`v_2(A_j) = d_j`), the fourth note (family dichotomy), the barrier atlas
[`collatz_procgen_20260922_barrier_atlas.md`](collatz_procgen_20260922_barrier_atlas.md)
(controls SHEET, DRIFT, DEFECT, INTEGRAL, UNIFORM, DIMENSION; "every
mechanism that overcomes DRIFT is blind to SHEET"), the Kuratowski note's
strategy square, the parallel spine note (descent-tree in-degree `c_D`).
Classical facts from memory: Terras 1976, Erdős 1976 (aliquot growth
persists over any fixed horizon for almost all abundant `n`), the
lifting-the-exponent lemma, Ducci/Kaprekar/look-and-say, the Juggler and
Lychrel problems.

**Status: PROVED elementary (Propositions 1–4, with the Mersenne free-
regeneration lemma the one exact new item) + FINITE-EXACT (all identities
re-checked on ranges; persistence statistics on both sides; in-degrees to
632; Juggler and reverse-and-add statistics) + CITED + DIRECTION. Collatz
OPEN. No proof step: D24 and D25 resolve into exact statements that locate
the pointwise difficulty precisely and a typology that explains, without
proving, why the open problems are open. Independent audit OWED.** Script
`04-computation/experiments/collatz_two_carries_typology_20260927.py`,
output beside it.

## 0. The answer in one paragraph

The two maps' opposite memories have one algebraic source each, and both
are one-line identities. Collatz: `d_j = v_2(3^j n + S_(j-1))` with
`v_2(3^j n) = 0`, so every bit of 2-adic content in the orbit's Apéry form is
*created* by the carry `S_(j-1) = sum 3^(j-1-t) 2^(d_t)`; without the `+1` the
map is `3x`, valuation-free and trivially divergent. Aliquot: `v_2(s(2^a m))
= min(a, v_2 sigma(m))` when the two differ (and `>= a + 1` when equal), so
the carry `-n` *caps* the 2-adic content of `sigma` at the driver; without
the `-n` the map `sigma` has exploding 2-adic content. Created content is
fresh each step (memoryless); capped content persists (sticky)
(Propositions 1–2). The pointwise form of memorylessness is then exact
(Proposition 3): after a run of `K` ones the orbit is at `2·3^K t - 1` where
`t = (m+1)/2^(K+1)` is the free cofactor, and the next valuation is `1 +
v_2(3^(K+1) t - 1)`; regeneration to depth `J` needs `t = 3^(-(K+1)) mod
2^(J-1)`, whose least positive member is `1` exactly when `J <= 3 + v_2(K+1)`
(`K` odd) or `J <= 2` (`K` even), so the Mersenne number `2^(K+1) - 1`
regenerates for free only to depth `3 + v_2(K+1)`, logarithmic in the run,
and every deeper regeneration is paid in cofactor bits. That is what
"memorylessness for every integer" means: the orbit's memory is the
cofactor, the cofactor is the integer's own higher bits, and the integer
has only `log_2 m` of them. D25 is the same statement read from both
sides: Collatz divergence for `n` says the valuations of the Apéry forms
stay below the critical slope for ever, aliquot divergence says `v_2
sigma(m_k) > a_k` for ever; the growth event persists with probability
`0.52` on the Collatz side (Terras's `1/2`) and `0.80, 0.72, ..., 0.59` over
two to six steps on the aliquot side (Erdős), and the average-case
theorems on the two sides reach opposite conclusions by the same
averaging method, while neither side has a pointwise technology.
In-degrees (D26): the aliquot graph below 632 has 52 leaves, mean in-degree
`7.0` and a long tail (up to 41 preimages), the Syracuse graph has leaf
density exactly `1/3` and infinite in-degree elsewhere, compressed by the
descent tree to `c_D ≈ 1.67`; bookkeeping, no principle. The typology
(section 5): among the classical arithmetic iterations, provability tracks
the presence of an exact invariant (finite state space, linear algebra),
and among the invariant-free maps the conjectured fate tracks the sign of
the drift and the memory of the driving quantity; Collatz and Juggler sit
in the corner "negative drift, memoryless, no invariant", the aliquot map in
"positive drift in sticky classes", `5n+1` in "positive drift, memoryless",
reverse-and-add in "positive drift, carries". The deep principle for the
barrier atlas is a new control, STICKY: an argument for Collatz termination
that would apply verbatim to a map whose valuations persist is refuted by
the aliquot map, so every valid argument must use memorylessness, i.e.
Terras's bijection, for the specific integer, and the only place that
bijection is available pointwise is the cofactor identity of Proposition 3.

## 1. The two carries (PROVED)

**Proposition 1 (Collatz creates its 2-adic content).** For odd `n` with
Syracuse orbit `m_j` and `A_j = 3^j n + S_(j-1)`, `S_(j-1) = sum_(t<j)
3^(j-1-t) 2^(d_t)`: `v_2(A_j) = d_j` while `v_2(3^j n) = 0`. Hence the
valuation sequence is `d_j = v_2(3^j n + S_(j-1))`, the 2-adic valuation of
a fixed odd number perturbed by a carry whose own 2-adic structure is the
past; and `S_(j-1) = 3 S_(j-2) + 2^(d_(j-1))` adds a term of valuation exactly
`d_(j-1)` at every step. Without the `+1` the map is `3x` on odd numbers,
valuation-free and divergent.

*Proof.* First note, Proposition 3; the valuation of `S_(j-1)` is `d_0 = 0`
(its `t = 0` term is `3^(j-1)`), and `v_2(A_j) = d_j` is the identity `A_j =
2^(d_j) m_j`. Checked on all odd `n < 4000` along their whole orbits. ∎

**Proposition 2 (the aliquot map caps its 2-adic content).** For `n = 2^a
m` with `m` odd and `s(n) = sigma(n) - n = (2^(a+1) - 1) sigma(m) - 2^a m`:
`v_2(s(n)) = min(a, v_2 sigma(m))` when the two differ, and `>= a + 1` when
they are equal. Without the `-n` the map `sigma` has `v_2(sigma(n)) =
v_2(sigma(m))`, unbounded above along iterates.

*Proof.* `2^(a+1) - 1` is odd; compare valuations. Checked on even `n <=
60000` (4548 equal cases). ∎

Reading. A created valuation depends on the current cofactor's relation
to the carry and is fresh each step; a capped valuation depends on the
driver and on `v_2 sigma(m)`, which is at least the number of odd-exponent
prime powers of `m` and grows with size. The ninth note measured the two
behaviours; this is their algebra.

## 2. D24: the pointwise form of memorylessness (PROVED)

**Proposition 3 (the free cofactor).** Let `m = 2^(K+1) t - 1` with `t >= 1`
(`m` odd, `K >= 1`). Then:
(a) the orbit follows `v = 1` for exactly `K` steps when `t` is odd, and
`U^K(m) = 2·3^K t - 1`;
(b) the next valuation is `1 + v_2(3^(K+1) t - 1)`;
(c) the next valuation is at least `J` iff `t = 3^(-(K+1)) mod 2^(J-1)`; the
least `m` realizing depth `J` after a run of `K` ones is `2^(K+1) rho_J - 1`
with `rho_J` the least positive residue of `3^(-(K+1))` modulo `2^(J-1)`;
(d) (free regeneration) `rho_J = 1` iff `2^(J-1) | 3^(K+1) - 1` iff `J <= 2`
(`K` even) or `J <= 3 + v_2(K+1)` (`K` odd); so the Mersenne number `2^(K+1)
- 1` continues after its run of `K` ones with valuation exactly `2` (`K`
even) or `3 + v_2(K+1)` (`K` odd), and no deeper regeneration is free.

*Proof.* (a),(b): `U(2^(K+1) t - 1) = 3·2^K t - 1`, iterate; then `3(2·3^K t -
1) + 1 = 2(3^(K+1) t - 1)`. (c): `v_2(3^(K+1) t - 1) >= J - 1` iff `3^(K+1) t = 1
mod 2^(J-1)`. (d): lifting the exponent, `v_2(3^(K+1) - 1) = 1` for `K + 1`
odd and `2 + v_2(K+1)` for `K + 1` even. Checked for `K <= 11, t < 60` (a, b),
`K, J < 8` (c) and `K < 200` (d). ∎

What this says for a proof. The "memory" of a Collatz orbit at the end of
a run is the cofactor `t`, which is the integer's own higher bits; Terras's
bijection is the statement that over `t` uniform modulo `2^(J-1)` the next
valuation is geometric; pointwise, the next valuation is decided by the
2-adic distance of `t` from `3^(-(K+1))`, and (d) shows the distance is
never small for free beyond `log_2 K` bits. A divergence proof must use
this for the specific integer: it must show that the cofactors met along
the orbit are never, infinitely often, 2-adically close to the inverses of
the powers of three that the run lengths dictate. That is a statement about
the 2-adic expansion of one integer against the 2-adic expansions of `3^(-k)`
for the `k` the orbit chooses, i.e. against the 2-adic digits of `1/3`,
`1/9`, ... — the same object that the first note's `F_n(1/3) = -3n` and the
sixth note's fixed point `Q(1/3) = 1/3` keep pointing at.

## 3. D25: the dual obstructions (PROVED reading + FINITE-EXACT)

Collatz divergence for `n`: `3^L > 2^(d_L)` for all `L` up to the carry
correction (no coefficient descent), i.e. the created valuations stay below
the critical slope `log_2 3` for ever. Aliquot divergence for `n` (via
drivers): `v_2 sigma(m_k) > a_k` for ever along the odd cofactors, i.e. the
capped valuations stay at the driver for ever. Both are "a 2-adic valuation
of an arithmetic expression along a deterministic sequence stays on one
side of a threshold for ever"; the first is a universal statement whose
failure is the desired event, the second an existential-forever statement
whose success is the desired event.

FINITE-EXACT: the growth event persists with probability `0.520` on the
Collatz side (`P(v_(j+1) = 1 | v_j = 1)` along orbits of odd `n < 2·10^5`;
Terras: `1/2`) and, on the aliquot side, `P(s_2 > s_1 > n) = 0.80`, `P(s_3 >
s_2 > s_1 > n) = 0.72`, ..., `P(six increases) = 0.59` over 3000 random
abundant `n <= 10^6` (Erdős's persistence). The average-case theorems on
the two sides (Terras/Tao for descent almost everywhere; Erdős for growth
almost everywhere among abundant numbers) are proved by the same averaging
over residues or prime factorizations, and reach opposite conclusions
because of Propositions 1–2. Neither side has a pointwise technology, and
the pointwise statements are dual in the sense above.

## 4. D26: in-degrees (FINITE-EXACT)

Aliquot, exact for `n in [2, 632]` (preimages `m <= n^2` scanned): 52 leaves
(the untouchables, `8.2%`), in-degrees `1: 137, 2: 102, 3: 42, 4: 13, 5: 21, 6:
17, ...` with a tail to 41 preimages, mean in-degree `7.03`. Syracuse: leaf
density exactly `1/3` (multiples of 3), infinite in-degree elsewhere, and the
parallel note's descent tree compresses the graph to mean in-degree `c_D in
[1.67, 1.70]`. The aliquot in-degree of `1` is infinite (the primes), like
the Syracuse in-degree of every non-multiple of 3. No principle beyond the
bookkeeping; the leaf densities (`0.08` here, `~0.17` asymptotically by
Pollack–Pomerance's heuristic; `1/3` exactly for Syracuse) are the one
comparable number.

## 5. The typology, and STICKY as a control (FINITE-EXACT + CITED + DIRECTION)

| map | drift | memory of the driving quantity | exact invariant | conjectured fate | status |
|---|---|---|---|---|---|
| Collatz `3n+1` | negative (`-0.415` bits per odd step) | memoryless (`0.52`) | none | all terminate | OPEN; Terras density 1; Tao |
| `3n-1` | negative | memoryless | none | three cycles, all bounded | OPEN |
| `5n+1` | positive (`+0.16`) | memoryless | none | most diverge | OPEN; no divergence proved |
| aliquot `s(n)` | `-0.35` in class 1, `+0.13..+0.47` in classes `>= 2` | sticky (`0.90, 0.75, 0.57`) | none | many diverge (Guy–Selfridge); the Lehmer five | OPEN |
| Juggler `floor(n^(3/2))`, `floor(sqrt n)` | negative on `log log n` (`-0.32` measured, `-0.21` fair coin) | nearly memoryless (parity persists `0.56`) | none | all terminate | OPEN |
| reverse-and-add | positive (`+0.41` digits per step) | carries | none | Lychrel numbers never palindromic (`196`; 6091 candidates below `10^5`) | OPEN |
| Ducci, length `2^k` | — | — | linear over `GF(2)`, nilpotent | terminates | PROVED |
| Kaprekar, fixed digit count | — | — | finite state space | terminates | PROVED |
| look-and-say | positive (Conway's `1.3036`) | — | linear over 92 atoms | diverges | PROVED |

The principle: provability tracks the presence of an exact invariant;
among the invariant-free maps the conjectured fate tracks the sign of the
drift and the memory of the driving quantity. Collatz and Juggler share the
corner "negative drift, memoryless, no invariant", the corner in which
termination is expected and nothing certifies it; the aliquot map is the
sticky corner in which divergence is expected. In the barrier atlas's terms,
DRIFT is the drift coordinate and SHEET the sign of the intercept; memory is
a third coordinate the atlas does not have.

**STICKY (proposed control).** The aliquot map has the same 2-adic engine as
Collatz with sticky valuations and is conjecturally divergent from its own
smallest open cases (the Lehmer five). Hence any argument for Collatz
termination that never uses the memorylessness of the valuations, i.e.
never uses Terras's bijection or an equivalent, would apply to the aliquot
map and prove the Catalan–Dickson conjecture, which the community expects
to be false. So: every valid divergence-half argument must use
memorylessness for the specific integer, and Proposition 3 shows the only
pointwise form of it: the cofactor's 2-adic distance to `3^(-(K+1))`. The
published mechanisms that overcome DRIFT (Terras, Korec, Tao; the
predecessor counts of Krasikov–Lagarias) all use residue averaging and are
therefore not refuted by STICKY, as they should not be; the cycle-exclusion
mechanisms do not address divergence. STICKY does refute any "Lyapunov
function of the valuation word" (S13's Proposition 4 already did), any
"stopping-time rank from local data" (the parallel note's section 6.1),
and any argument of the shape "drift is negative, so orbits descend" that
does not certify the memorylessness the drift computation assumes.

## 6. Directions (DIRECTION; none pursued)

* **D27.** Turn Proposition 3 into a lemma about the full orbit: at every
  time, the next valuation is `v_2` of `3^(j+1) n + S_j` written as `3^(j+1)
  (n - x_w) + (stuff)` for the cycle point `x_w` being shadowed; the
  cofactor is `(m - x_w)/2^K`. A pointwise divergence proof would have to
  show that the sequence of cofactors is never, infinitely often, 2-adically
  close to the required inverses; whether this can be phrased as a
  statement about the 2-adic expansion of `n` alone (it can: the cofactor at
  time `j` is a function of `n` and the word) is the D24 form of the
  transversality statement.
* **D28.** Enter STICKY into the barrier atlas as a control with the
  aliquot map as witness, and re-type the published mechanisms against it;
  a Collatz-like map on `N` with provably sticky valuations and provable
  divergence would be a cleaner witness than the aliquot map (whose
  divergence is only conjectured); none was found here.
* **D29.** The typology suggests a fourth column, "does the map's growth
  regime strengthen with size?" (aliquot: yes, `v_2 sigma(m)` grows with
  the number of prime factors; Collatz: no, scale-free; reverse-and-add:
  carries scale with digits): a quantitative version would be the
  size-dependence of the transition matrices measured in the ninth note.

## 7. Reproduction

    cd 04-computation/experiments
    python3 collatz_two_carries_typology_20260927.py > collatz_two_carries_typology_20260927.out

Needs sympy for the aliquot sample; about eight minutes.
