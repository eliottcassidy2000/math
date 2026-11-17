# Circulant Hadamard {1, 4} against the holes {7, 21}: three finite exceptional sets from one kind of constraint

Session klein-2026-11-16-circ (2026-11-16). Owner: *"Consider that a circulant Hadamard matrix only exists if the order
is {1, 4} and how this relates abstractly analogously to {7, 21}. Also consider related to the packing-polynomial
lattices."* The concurrent `packing_hadamard_flatness_20261117` (other machine, same words) already records the
identity `P_T(q) = G_T(2q, q)`, the exact Hadamard defect of a regular tournament as two coefficients of `G`, the
orbit-phase reduction of `G` modulo a prime-order automorphism, and the verdict that 1, 4 are allowed *orders* while
7, 21 are excluded *counts*; it cites the 2026 proof of the circulant Hadamard conjecture. Nothing below repeats
those; this note states the abstract analogy as the owner asked, with the three finite sets side by side, and the
specialisation of the orbit-phase reduction to circulants with an exhaustive check.

## 1. Three finite exceptional sets

| statement | infinite family | constraint | survivors |
|---|---|---|---|
| circulant Hadamard (proved 2026; cited by the concurrent packet) | orders `n` | `Z_n`-invariance of a `±1` matrix with flat autocorrelation (`|ĥ(j)|² = n` at every frequency) | `{1, 4}` |
| vertex-transitive tournaments with a prime count (THM-4861, THM-4747) | orders `n` | `Z_n` acts freely on Hamiltonian paths, so `n | H`; strong floors | `{3}` (the 3-cycle, `H = 3`) |
| Chan–Pak holes (THM-1370; completeness open) | odd counts | strong floors `f(n)` and the finite kernel of orders ≤ 6 | the excluded `{7, 21}` |

The first two are the same mechanism — a cyclic symmetry turns an existence question into congruences, and the
congruences leave finitely many survivors; the third is the opposite sign: a growth constraint (Busch's floors) and a
finite check leave finitely many *non*-survivors. So `{1, 4}` and `{7, 21}` are not two coincident pairs; they are the
survivor set of a symmetry constraint and the exception set of a growth constraint. The honest tournament analogue of
"circulant Hadamard only at 1 and 4" is "a cyclic tournament has a prime count only at order 3" — the count of a
circulant tournament of order `n` is `n` times the number of rotation orbits of its Hamiltonian paths, and only `C_3`
has a single orbit (its floor `f(n) > n` for `n ≥ 4`). In that reading the holes are the single-orbit counts that
orders 7 and 21 would have if they could have one orbit: `7 = 7·1` and `21 = 21·1`, both forbidden by the floors
exactly as `C_3`'s `3 = 3·1` is forced.

## 2. The circulant congruence (the orbit-phase reduction specialised; checked on every connection set)

For a circulant tournament of prime order `p` the rotation acts on the spanning odd-cycle packings; packings with a
cycle of length `< p` lie in orbits of size `p`, and the only fixed packings are the all-singleton one and the
`(p − 1)/2` constant-step Hamiltonian cycles (one per element of the connection set). Hence, coefficientwise,

    G_T(x, y) ≡ y^p + ((p − 1)/2)·x   (mod p),      so   P_q(T) = G_T(2q, q) ≡ q^p − q ≡ 0   (mod p)

for every `q` (Fermat), and in the coefficients `a_m` of `P_q = q^p + Σ a_m q^m`: `a_1 ≡ −1`, `a_m ≡ 0` for `1 < m < p`.
This is the circulant case (`F = ∅`, one moved orbit) of the concurrent ORBIT_PHASES.md. Checked on every connection
set of orders 3, 5, 7, 11, 13 ([receipts/circulant.out](receipts/circulant.out)): 908 congruence checks, 0 failures (H-1).
The circulant counts are `3` (order 3), `15` (order 5), `175, 189` (order 7: the `{1, 2, 3}` circulant and Paley), `92411, 93027,
93467, 95095` (order 11, Paley last) and six values from `3,669,497` to `3,711,175` at order 13; the orbit numbers `H/p` are
`1, 3, 25, 27, 8401–8645, 282269–285475` — one orbit only for `C_3` (H-4). The coefficients also show the concurrent packet's
flatness: `a_{p−2} = 2·(number of triangles)` and `a_{p−4}` are the same for every circulant of a given order (110 and 3168 at
order 11), and the Hadamard defect first appears at `a_{p−6}`.

**Why this is the Hadamard analogue in substance.** A circulant Hadamard row is a `±1` sequence whose *every* Fourier
coefficient has modulus `√n`; a circulant tournament is a `±1` skew sequence whose packing polynomial is congruent,
modulo `p`, to the polynomial of the empty tournament plus `(p − 1)/2` copies of one Hamiltonian cycle: the symmetry
flattens the packing data to two terms. Both flattenings are what kill primes: a prime count would need `H ≡ 0 (mod p)`
with a single orbit.

## 3. Typing

- PROVED: the circulant congruence (Burnside + Fermat), `n | H` for vertex-transitive tournaments (free action), the
  single-orbit statement (floors).
- KNOWN: the circulant Hadamard classification (cited through the concurrent packet, not re-derived here).
- ANALOGY: the table of §1; no identification between the Hadamard orders and the holes is claimed.
