# Every odd residue is realised: what was already proved, and the monoid route through the first 32 odd numbers

Session klein-2026-11-16-monoid. Scripts `scripts/frattini3.py`, `scripts/residues.py`; receipts `frattini3.out`,
`residues.out`, `savoid.out`.

## 0. Priority

The owner asked to prove that every odd residue modulo every `2^a 3^b 5^c` is realised by a tournament count, and
that the pattern continues for every odd prime up to 63 and beyond, skipping over 7. **This is already proved, and more
strongly**, by two earlier packets, which this note does not repeat:

- `session25_spectrum_wiring_20261021/session25-spectrum/PROOF.md`: no finite congruence obstruction beyond parity
  exists for the Hamiltonian-path count, even among strong tournaments with feedback vertex number one; every
  parity-admissible residue modulo every modulus is attained (explicit reset blocks; checked for all moduli to 512).
- `one_apex_all_residues_20261116` (hadamard, 2026-11-16): for every modulus `M` every residue represented by an odd
  integer is the count of a strong one-apex tournament at **every sufficiently large order**, with an order bound
  growing like `√m` in the odd part `m` of `M`; the same controls place `H ≡ r (mod M)` within `O_M(√N)` of any large
  target `N`.

So "every odd residue modulo every modulus" is PROVED (CITED). What this note adds is a different mechanism with a
different conclusion: the residues are already generated **multiplicatively by the first 32 odd numbers minus the two
holes** — by ordinal sums of tournaments of order at most 8 — and the multiples of 7 need exactly one infinite line.

## 1. The monoid and the Frattini criterion

Ordinal sums multiply counts, `H(A → B) = H(A) H(B)` (THM-1370: `H` is the product over the strong components), so the
set `R` of counts is a multiplicative monoid. Every odd number ≤ 609 except 7 and 21 is a count of a tournament of
order ≤ 8 (THM-1370), and every odd number ≤ 2^32 except 7 and 21 is a count (THM-4875 (D)).

For a finite abelian group `G` a subset generates `G` (as a monoid, since `G` is finite) iff its image generates the
Frattini quotient `G/Φ(G) = ∏_q G/qG`. For `G = (Z/m)^×`, `m = ∏ p^(e_p)`, the components of `G/qG` are: for `p = 2`
(`e_2 ≥ 3`) the two bits `[x ≡ 3 (mod 4)]`, `[x ≡ 3, 5 (mod 8)]`; for an odd `p` with `q | p − 1` the `q`-th power residue
class of `x` modulo `p`; for `q = p` (`e_p ≥ 2`) the value `(x^(p−1) − 1)/p (mod p)`. None of these depends on the
exponents beyond `e_p ≥ 2`, so **generation modulo `8 ∏ p²` is generation for every exponent vector**.

**Theorem A (PROVED; receipt `frattini3.out`).**
1. The fifteen numbers `11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59, 61` — the first 32 odd numbers minus
   `7, 21` and minus the multiples of 3 and 5 — generate `(Z/2^a 3^b 5^c)^×` for all `a, b, c` (Frattini ranks
   `4/4, 1/1, 1/1` at `q = 2, 3, 5`).
2. The odd numbers ≤ 4000 other than 7, 21 whose prime factors are 7 or ≥ 67 generate `(Z/m)^×` for every `m` whose
   prime factors are ≤ 63 other than 7 (all exponents); those with all prime factors ≥ 67 do so for every `m` with
   prime factors ≤ 63 including 7 (full ranks at every `q ≤ 61`).
3. Consequently (with the realised primes `3, 5, 11, …, 61` for the non-unit residues) **every odd residue modulo
   every modulus with prime factors ≤ 63 other than 7 is the count of an ordinal sum of tournaments of bounded
   order**, and the products of the first 32 odd numbers minus 7, 21 hit every odd residue modulo every
   `2^a 3^b 5^c` (replay: no miss for `m ≤ 2·10^5`).

*Proof.* Monoid closure (ordinal sums) and the Frattini criterion, with the ranks computed in the receipt; a residue
`r` with `v_p(r) = v_p < e_p` is `(∏ p^(v_p))·u` with `u` a unit modulo `m/∏ p^(v_p)`, realised as the ordinal sum of
`v_p` copies of a count-`p` tournament and a product of generators. ∎

## 2. Skipping over 7: the one line that carries the 7-lane

Seven is not a count, so a residue divisible by exactly one 7 cannot be built from a count-7 block. The multiples of 7
that are counts are the "hidden" generators of MONOID.md (`35, 49, 63, 77, …`), and for a proof valid for every power
of 7 one needs an infinite supply with controlled residues. The line `L_3 = {1 + 3·2^s}` (THM-5126 (B): every point
except 7 is a count) contains `49, 385, 3073, …` at `s ≡ 1 (mod 3)`, and `u′ = (1 + 3·2^s)/7` runs through **every**
residue class modulo `7^k` as `s` runs through `s ≡ 1 (mod 3)` (because `2^s = 2·8^((s−1)/3)` and `8 = 1 + 7` generates
the units `≡ 1 (mod 7)` of `Z/7^(k+1)`; the class `s ≡ 4 (mod 21)` is exactly `u′ ≡ 0 (mod 7)`). The diagnostic
`savoid.out` shows that for each of the six classes `s (mod 7)` with `u′ ≢ 0` there are hundreds of `s < 20,000` with
`1 + 3·2^s` a unit modulo every prime ≤ 63 other than 7, the only obstructing primes being `5, 11, 29, 37, 53, 59`,
whose orders of 2 carry at most one factor 7 (`28`), so the avoidance is compatible with every lift of `s` modulo
`3·7^k`. Hence:

**Theorem B (PROVED, with the finite avoidance certificate).** Every odd residue modulo every modulus whose prime
factors are ≤ 63 (7 included) is a count: residues prime to 7 by Theorem A, residues with `v_7 = 1` as
`(1 + 3·2^s)·w` with `s` chosen on `L_3` and `w` a unit product, and residues with `v_7 = v ≥ 2` as `49^i 343^j w`
(`v = 2i + 3j`). The replay with `L_3`, `49`, `343` and the units misses nothing on the mixed moduli tested.

## 3. Reading

The owner's "multiply, or step by ±2 to reach the primes" is the division of labour between the two halves of the
problem: multiplication (ordinal sums) gives **every residue class** at bounded block size — Theorems A, B — but can
never give an exact prime (a prime is indecomposable); the exact primes need additive moves, which is the digit
programme (`3H + 2σ`, the `+2` partners, the port moves of PORTMOVE.md). Collatz and multiplication modulo odd numbers
enter as in wave 47 §3 (THM-4792: the tables alternate iff `−1 ∈ ⟨2⟩`); they describe the rows and the residue clocks,
not the exact values.
