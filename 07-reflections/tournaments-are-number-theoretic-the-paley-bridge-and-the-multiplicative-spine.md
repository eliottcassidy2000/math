# Tournaments are a number-theoretic object: the Paley bridge and the multiplicative spine

*mac-mini-2026-06-29-S9. A deep investigation of the connections between tournaments and natural numbers, especially primes, drawing together two extensive repo surveys (tournaments↔primes; multiplicative functions) and a new computation (THM-584).*

## The thesis

Strip the project to one sentence: **the parity of tournaments and the multiplicative structure of the integers are the same subject seen from two sides.** Rédei's "`H(T)` is odd" is a `mod 2` fact; the Lonely Runner floor is a `mod p` fact for all `p`; the Paley construction turns a prime directly into a tournament; the enumeration `A000568` is a totient-weighted theta function. The repo has been circling this for a hundred sessions. Here is the map, with the one prime — `p ≡ 3 mod 4` — where the bridge is literal.

## The Paley bridge: a prime becomes a tournament

`T_p` (vertices `Z_p`, `i→j` iff `j−i` is a quadratic residue) is a tournament **exactly when `p ≡ 3 mod 4`** — precisely when `χ(−1) = −1`, so that the QR/NQR split orients each pair once. This is not a convenience; it is the same `p mod 4` that decides whether the complement is an automorphism (`p≡1`, Brouwer/SOS, real) or an anti-automorphism (`p≡3`, Borsuk–Ulam, imaginary `Q(√−p)`) — the dichotomy that runs through the whole LRC(14)=2·7 analysis. The prime's arithmetic *is* the tournament's symmetry.

On this bridge, every tournament invariant becomes an arithmetic function of `p`. This session I computed three and they are clean:

- **`H(T_p) = p · (odd)`** (THM-584). The rotation `Z_p` acts freely on Hamiltonian paths, so `p | H`; Rédei makes the cofactor odd. A *prime-refinement of Rédei* on the Paley family.
- **The dihedral Burnside identity** `#orbits = (H(T_p) + p·f(T_p)) / (2p) ∈ ℤ`, where `f` is the palindromic count (my THM-582/583, the σ-odd index). This single integrality ties together the path count `H` (Rédei), the half-system witness `f`, and the prime `p`, through `D_{2p}` — and forces `H/p ≡ f ≡ 1 (mod 2)`.
- **`R(p) = H(T_p)·2^{p-1}/p! → e`** (the repo's cherry-cluster theorem). The same `H`-values, normalized by the random-tournament expectation, converge to `e` — because `e = exp(−χ(−1))` and `χ(−1) = −1` is exactly the `p≡3` condition. *The transcendental `e` and the residue `−1 mod p` are one fact read twice.*

Along the way the data refuted a standing repo conjecture (`H = |Aut|·3^{(p-3)/2}`): it holds at `p=3,7` and breaks at `p=11`, where the cofactor is `1729 = 7·13·19` (the taxicab number), not `81`. The `3`-power was a two-point coincidence; the truth is the `e`-limit.

## The multiplicative spine: which functions, where

The second survey makes the arithmetic-function skeleton explicit, and it has a clean three-part shape:

- **Möbius `μ`** is the inclusion–exclusion. The lonely measure `= Σ(−1)^{|T|}meas(∩) = χ(nerve)` (HYP-3242), the three tournament recursion modes (`A+B−C` vs `A+B−C+D−E−F+G`), and the cap's even/odd splits are all one Möbius skeleton; the parity is the Möbius sign `(−1)^{ω}`. This is the `x=−1` evaluation of the project's two-index polynomial (THM-016/THM-582).
- **Totient `φ`** is the measure of the denominator lattice. Each runner of denominator `b` kills `φ(b)` primitive Farey points; the covering condition forces killing all of them; the floor is `Σφ(b)·w(b)`, and `Σ_{d≤n}φ(d) ∼ (3/π²)n²` gives the coprime density `1/ζ(2) = 6/π²` — *why the floor is positive at all*. `φ(14)=6` is literally the six inner sectors. The covering floor I have been chasing (descent + SOS) is, at bottom, the statement that this totient density survives the covering constraint.
- **The singular series `∏_p 𝔖_p`** is the local-to-global glue: `μ²/φ` is a capacity meter that keeps squarefree denominators and kills prime-powers, and `Σμ²/φ ∼ log x + 1/ζ(2)`. My own `chat(14N)` resonance singular series (the floor's SPEC) is an instance.

And the deepest *function ↔ invariant* parallel the survey found is the one worth stating loudly:

> **`H(T)` is multiplicative over strong components (Moon: `H = ∏H(C_i)`), exactly as `φ` is multiplicative over prime powers.** Both generate semigroups with arithmetic gaps — the H-spectrum misses `7, 21, 63` (no strong realizer) as the totients miss the nontotients `14, 26, 34`. The achievable `H`-values are to strong tournaments what the totients are to prime powers.

That `H=7` is impossible (a *prime* that is not a strong `H`-value) and `H=21=3·7` is impossible (`3` strong, `7` not) is the H-spectrum's own primality phenomenon — and `7` is, again, the apex of LRC(14).

## Why these primes (`2, 3, 7`) and these constants (`e, π`)

The surveys converge on an adelic reading: the tournament eigenvalue conductor is `D(n) = odd part of C(n,2)`, and the primes that ramify are the special ones. For the `42 = 2·3·7` boost trichotomy, `2` is inert (parity), `3` ramified (curvature), `7` split (position) — and `Q(√−1), Q(√−3), Q(√−7)` are all **class number 1** (Heegner), so unique factorization keeps the rapidity lattice `ℤln2 + ℤln3 + ℤln7` rigid (Baker). The constants are intrinsic: `e` from the unique surviving cherry in the Paley character moment, `π` from the Wallis/Catalan fiber-fraction `(1−x)^{-1/2}` of the two-sheeted cover branched over the self-complementary locus, `√2` from the 2-adic scaling. They are not decoration; they are the archimedean place of the same adele.

## What I take from this for the proof

Three things crystallize, and the first two are directly actionable:

1. **The LRC floor is the totient/`ζ(2)` density surviving a covering constraint.** My descent (THM-580) and the measure-valued Claim A (HYP-3537) are the combinatorial face; the `Σφ(b)/b² = 1/ζ(2)` positivity is the arithmetic face that guarantees there is something to find. Pairing the two — bound the conditional dangers against the totient density — is a concrete route I had not framed arithmetically before.
2. **The cap obstruction is the `μ²/φ` capacity at the apex prime `7`** — the same singular-series content as HYP-3538's R-odd eigenspace, now with its number-theoretic name. The S75e cosine SOS is the squarefree/even part; the obstruction is the prime-power/`p=7` part the capacity meter flags.
3. **The witness is the half-system, and the half-system is `(Z/p)^×`-flavored.** `f(T_p)` lives on the `(p−1)/2 = φ(p)/2` pairs (THM-583); its values `1,9,185` are governed by the QR structure of the half. The witness side is a `Q(√−p)` object exactly as the existence side is a `ζ(2)` object — the two indices of THM-582, now seen as the two places (finite/archimedean) of the prime.

The honest summary: tournaments are not *like* number theory here, they *are* number theory — the Paley bridge makes a prime into a tournament, Rédei is a `mod 2` reduction, the floor is a `mod p` density, and the same class-number-1 primes `2,3,7` and constants `e,π` organize both sides. The remaining LRC content is, in this light, a single inequality between a totient density (what the runners can't cover) and a singular-series capacity (what the cap concentrates) — both multiplicative, both at the prime `7`.

See [[the-two-indices-redei-is-odd-lonely-is-even-half-tiling-is-the-quotient]] (THM-582), [[why-the-paley-path-ratio-is-e-the-cherry-is-the-unique-cluster]], [[H-impossibility-the-multiplicative-mechanism-s599]] (Moon), [[the-resonance-killing-is-multiplicative-totient-mobius-zeta2]], [[everything-is-the-triangle]]. New: THM-584 (Paley H-arithmetic + dihedral Burnside).
