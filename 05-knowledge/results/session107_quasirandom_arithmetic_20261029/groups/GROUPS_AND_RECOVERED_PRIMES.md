# Monster/Baby Monster, prime labels, and two incompatible meanings of quasirandomness

Session 107 — 2026-10-29. Status: recovered prior results; cited external structure; elementary proofs below. Independent audit is recorded in `../audit_structure/AUDIT.md`.

## Useful conclusions first

1. The exact prior source of `{4,2,5,3,7,13}` is a **valuation chain**, not a set of six primes. From the Monster/Baby Monster order ratio it is `13→2→5→3→7→4`; identifying 13 with 4 closes a five-state graph cycle. That identification is not arithmetic equality.
2. There is a genuine Monster-local-to-Paley action in the prior work: `N_M(<7A>)/He ≅ AGL(1,7)`, with orientation-preserving subgroup `7:3`. **Neither affine group is quasirandom in the minimal-complex-representation sense:** each has a nontrivial cyclic quotient, so its minimum degree is 1. The local bridge forgets, rather than transfers, Monster quasirandomness.
3. Even more strongly, neither the whole Monster nor the whole Baby Monster can act nontrivially by ordinary automorphisms of a finite tournament. This does not forbid local subgroup actions, actions on other objects, or representation-valued statistics.
4. Their actual group-quasirandomness parameters are `d(M)=196883=47·59·71` and `d(B)=4371=3·31·47`. These are ordinary complex character degrees, not the smallest modular dimensions 196882 and 4370 in characteristic 2.
5. Ogg's fifteen primes have a precise modular-curve meaning, and the missing Baby primes have a precise **level-structure** meaning. Neither theorem identifies prime lists with quasirandomness parameters.

## A. Exact recovery and uncertainty ledger

### A1. Monster/Baby order quotient (recovered and arithmetic verified)

The ATLAS orders are

```
|M| = 2^46 3^20 5^9 7^6 11^2 13^3 17 19 23 29 31 41 47 59 71
|B| = 2^41 3^13 5^6 7^2 11 13 17 19 23 31 47.
```

[ATLAS M](https://brauer.maths.qmul.ac.uk/Atlas/spor/M/), [ATLAS B](https://brauer.maths.qmul.ac.uk/Atlas/spor/B/).

Put `N=2^5 3^7 5^3 7^4 11 13^2=39046339332000`. Then
`|M|/|B|=N·29·41·59·71=194478922284018372000`.
The common unchanged primes are `{17,19,23,31,47}`; new primes are `{29,41,59,71}`; increased common valuations occur at `{2,3,5,7,11,13}`.

Recovered files:
- `session32-monster/MONSTER_BABY_PRIME_PARTITION.md`;
- `session32-root/VALUATION_CYCLE_AND_QUOTIENT.md`;
- `math/05-knowledge/results/monster_levels_square_sum_7_21_20261022.md`.

The arrows `p→v_p(N)` give `13→2→5→3→7→4` and `11→1`. Only after the additional set identification `13~4` is there a cycle. Modulo 9 does not make the whole map descend: `11≡2` but `1≠5`. Modulo 3, `2≡5` but their images `5,3` differ. Every modulus identifying 13 and 4 divides 9, so no nontrivial quotient ring realizes the full arrow map. Moreover the actual 2A centralizer is `2.B`, and replacing the order ratio by its conjugacy-class index replaces N by N/2; the displayed quotient cycle then has period two, not five.

### A2. The exact requested cycle `{2,3,7,5,2}`

**RECOVERED by the parallel dynamics lane:** this is the orbit of 2 under Φ6(x)=x²−x+1 modulo 19:
`2→3→7→5→2 (mod19)`. Over the integers the same orbit starts `2→3→7→43→1807`; reducing these terms modulo19 closes the four-cycle because 1807≡2. Thus 19, 13·139, and the four-cycle do share a precisely specified polynomial orbit. They are distinct from the Monster valuation chain, whose map is p↦v_p(N). The companion dynamics lane investigates the lift modulo powers of19; no lift theorem is independently asserted here.

The primary Euclidean-prime-generator search also found an unrelated sequence beginning `2,3,7,5`; its matching prefix is not needed and supplies no evidence for this modular orbit.

### A3. The requested pair 13 and 139

The exact recovered connection for the requested pair is `1807=13·139`: the Sylvester/Φ6 orbit `2,3,7,43,1807,...`. See `math/07-reflections/bernoulli-fixed-point.md`, `tower-crossings-synthesis.md`, and `math/05-knowledge/results/prime_sorting_s116m.out`, section VII. Here Φ6(x)=x²−x+1, so Φ6(43)=1807. This gives an actual cyclotomic mechanism: for either prime q=13 or139, the residue 43 has multiplicative order 6 modulo q. It does not relate 139 to the Monster order; 139 is absent from Ogg's list.

The historical Bernoulli note correctly computes the fixed denominator 1806: the primes p with p−1 dividing 1806 are exactly 2,3,7,43. But its later claim that `{2,3}` is self-selecting or that 6 is fixed is false, since 7−1 divides 6. The verified chain is `6→42→1806→1806`. Do not inherit the erroneous sentence or its broad unproved slogan about arbitrary Sylvester trajectories.

The number 19 is both the exact modulus closing the Φ6 four-cycle and an Ogg prime/unchanged common divisor of |M| and |B|. The two predicates are individually exact. No mechanism equating them has been established.

## B. What the fifteen primes actually mean

The list is
`2,3,5,7,11,13,17,19,23,29,31,41,47,59,71`.
All supersingular j-invariants in characteristic p are rational over F_p exactly at these primes; equivalently X_0^+(p) has genus zero. This is not the assertion that these are the only characteristics admitting supersingular curves. The genus-zero prime levels of **X_0(p), without the Fricke quotient**, are only `2,3,5,7,13`.

Aricheta generalizes the rationality property to level structures; his Table 1 identifies the odd Baby Monster prime divisors with the rationality primes for X_0^+(2). Thus `29,41,59,71` retain rational supersingular j but fail rationality after adding the coarse unordered 2-isogeny structure. [Aricheta, Supersingular Elliptic Curves and Moonshine](https://arxiv.org/html/1809.07421), Theorems 1.1 and 1.3 and Table 1.

The prior Session32 packet went further with exact polynomial and independent finite-field calculations: at each of these four primes exactly one pair of level-two points is exchanged by Frobenius. This is inherited computational evidence, not rerun in this worker; its scripts and independent audit are available in `session32-monster/` and `session32-monster-audit/`. It is a real group/modular bridge but not a bound on representation degrees or a theorem about tournament spectra.

For Fermat primes F_k=2^(2^k)+1, the intersection with the Ogg list is exactly `{3,5,17}`: these are F_0,F_1,F_2, whereas every F_k for k≥3 is at least257>71. This needs no assumption about whether infinitely many Fermat primes exist. The same exact intersection holds for Baby Monster prime divisors. The prime2 is not a Fermat number in this indexing. Set intersection is a proved arithmetic fact, not an explanation of either classification.

## C. The rigorous quasirandomness channel

Let d(G) be the minimum dimension of a nontrivial irreducible complex representation. High d(G) forces mixing of products of dense subsets. If A,B,C have densities α,β,γ in G, then

```
| E_{x,y} 1_A(x)1_B(y)1_C(xy) − αβγ | ≤ sqrt(αβγ/d(G)).
```

Consequently αβγ>1/d(G) guarantees some product ab=c; a product-free subset has density at most d(G)^(−1/3). This is [Gowers, Quasirandom groups](https://arxiv.org/pdf/0710.3877), Theorem 3.3, with normalized counting measure. The conclusion concerns the group multiplication table and dense Cayley relations, not any arbitrary graph bearing its prime divisors.

For M and B the standard character tables give d=196883 and 4371. Primary verification: [Breuer–Magaard–Wilson, Baby Monster table](https://arxiv.org/pdf/1902.07758), and [Monster table](https://arxiv.org/html/2412.12182v1). The [computational supplement](https://www.math.rwth-aachen.de/~Thomas.Breuer/ctbllib/doc/ctblm.pdf), section 2, explicitly lists the only irreducible degrees of 2.B up to196883 as `1,4371,96255,96256`; the degree4371 character is trivial on the central involution, so descends to B. Hence this also certifies the Baby minimum. The Monster minimum is standard ordinary-character-table data; this worker did not recompute its 194 characters.

A particularly concrete representation bridge is restriction of the Monster's degree196883 character to its 2A centralizer:

```
196883 = 1 + 4371 + 96255 + 96256.
```

This is an actual decomposition of a specified representation, not just a numerical partition. The last constituent detects the central involution; the first three do not. It is a much stronger relationship than a common prime label. The restriction creates a fixed vector and much smaller constituents, showing directly why subgroup restriction does not preserve the ambient minimum degree.

Useful exact arithmetic: d(M)/d(B)=4189/93, and the sufficient product-free density bounds are about 0.01718952703 and 0.06116078691, respectively (decimal approximations only). The prime factorization of these degrees is exact; no claim is made that prime factorization alone determines irreducibility or mixing.

## D. Two firewall proofs

### D1. Tournament automorphisms prohibit whole-group actions

An involution σ acting on a finite vertex set, if nontrivial, swaps some pair x,y. In a tournament exactly one of x→y and y→x holds; σ reverses that arc. Thus a tournament automorphism group contains no involution. By Cauchy's theorem its order is odd. A nontrivial action of an even-order simple group has trivial kernel and would embed that group into this odd-order group, impossible. Since M and B are even-order simple groups, every ordinary action of either on a finite tournament is trivial.

This elementary obstruction applies regardless of tournament size. It also excludes a nontrivial whole-group action by automorphisms and global anti-automorphisms together: the orientation sign would give a homomorphism from the nonabelian simple group to C2, hence would be trivial, reducing to the ordinary case. For arbitrary finite sets (without oriented edges), a weaker general obstruction is that any nontrivial permutation action of a group with minimum degree d has at least d+1 points: remove the constant vector from the permutation representation and a nontrivial constituent remains. Neither statement rules out subgroup actions.

### D2. The genuine Paley bridge has minimum degree one

The recovered bridge is
`N_M(<7A>)/He ≅ 7:6 = AGL(1,7)`;
its index-two orientation subgroup is `7:3=Aut(Paley7)`.
The external local structure is checked in the Monster computational supplement, §4.13 (subdirect product of 7:6 and He.2), also recorded with detailed proof and sources in `session32-paley/PALEY_MONSTER_BRIDGE.md`.

For H=F7 ⋊ K, with K=F7* or its square subgroup, projection `(b,a)↦a` surjects onto C6 or C3. Composing this projection with a nontrivial complex character gives a nontrivial one-dimensional representation. Hence d(H)=1 exactly. Taking a subgroup of M and then a quotient is not the same as taking a quotient of M itself. Quotients cannot reduce d; subgroups can.

Paley graphs/tournaments over additive finite fields can nevertheless be spectrally quasirandom: the ambient additive group is abelian, so d=1, but the **chosen connection set** has small nontrivial Fourier coefficients. Group-quasirandomness says all large connection sets mix because low-dimensional representations are absent; Paley pseudorandomness says one specially constructed connection set has cancellation. These are different, compatible mechanisms. Ogg's finite prime list does not supply an asymptotic Paley family; arbitrary large admissible prime powers do.

## E. Completion and bounded uncertainty

Completed: local recovery of the valuation chain, modular-level bridge, Paley local quotient, 13/139 factor pairing, primary-source checks, elementary no-go proofs, and exact arithmetic script. No external state changed. No denied action. No group matrices or full character tables computed. The exact `{2,3,7,5,2}` transition rule was recovered by the parallel lane before completion. The factor-pair connection is established; no literal ratio interpretation is needed.
