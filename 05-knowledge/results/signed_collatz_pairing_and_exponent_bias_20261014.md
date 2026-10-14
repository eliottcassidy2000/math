# Signed Collatz pairing, exponent parity, and prime-power bias

Date: 2026-10-14. Status: PROVED elementary identities and density statements; FINITE-EXACT checks; OPEN fixed-sign convergence and complete cycle classification. No novelty claim.

## Main findings

The proposed pairing is exact at the one-halving level: for every odd integer n, the numbers (3n−1)/2 and (3n+1)/2 are consecutive and lie in residues 1 and 2 modulo 3, respectively. What changes under full acceleration is the number of powers of two removed. Retaining that exponent gives a precise signed coupling and an exact bias.

For a positive odd source n, write

    a = v₂(3n+1),  b = v₂(3n−1),
    u = U₊(n) = (3n+1)/2ᵃ,  v = U₋(n) = (3n−1)/2ᵇ.

Then:

1. min(a,b)=1 and max(a,b)≥2.
2. 2ᵃu−2ᵇv=2.
3. u≡(−1)ᵃ and v≡−(−1)ᵇ modulo 3.
4. The residue pair (u mod 3,v mod 3) is one of (1,1), (2,1), (2,2), never (1,2). Each allowed pair has relative natural density 1/3 among positive odd starting integers.
5. For every j≥1, each residue u≡2 mod 3 modulo 3ʲ has U₊-output density 2/3ʲ, and each residue u≡1 mod 3 has density 1/3ʲ. U₋ reverses these weights. Multiples of 3 have output density zero.
6. For n>1, exactly one of U₊(n),U₋(n) is smaller than n and the other is larger. Choosing the descending sign at each step proves convergence to 1 for that *controlled-sign system*. It does not prove convergence for either fixed-sign system.

These are starting-input densities, not frequencies along an individual orbit, a stationary distribution of the actual orbit, or densities of cycle basins.

## 1. Map types and the modulo-12 table

It is essential to distinguish three integer maps and one partial map. For ε∈{+1,−1}, define:

- Cε(n)=3n+ε for odd n, and n/2 for even n: the full, unshortened map.
- Tε(n)=(3n+ε)/2 for odd n, and n/2 for even n: the shortcut map.
- Sε(n)=(3n+ε)/2, defined only on odd integers: its odd branch.
- Uε(n)=(3n+ε)/2^v₂(3n+ε), defined on odd integers: the fully accelerated odd-to-odd map.

Cε and Tε are self-maps of the integers and preserve the positive integers. Sε is not an odd-to-odd self-map. Its output can be even, so Sε cannot simply be iterated on its declared domain. Uε is an odd-to-odd map. The valuation is taken on the absolute value for negative arguments; no integer odd source makes 3n+ε zero.

Since Sε(n+12)=Sε(n)+18, the exact modular table is:

| n mod 12 | S₊ mod 18 | S₋ mod 18 | v₂(3n+1) | v₂(3n−1) |
|---:|---:|---:|---:|---:|
| 1 | 2 | 1 | ≥2 | 1 |
| 3 | 5 | 4 | 1 | ≥2 |
| 5 | 8 | 7 | ≥2 | 1 |
| 7 | 11 | 10 | 1 | ≥2 |
| 9 | 14 | 13 | ≥2 | 1 |
| 11 | 17 | 16 | 1 | ≥2 |

Thus all six odd source classes go to 2 mod 3 under S₊ and 1 mod 3 under S₋. The source's original residue modulo 3 is immaterial for this one step.

Writing n=2k+1 gives the particularly transparent bijection

    (S₋(n),S₊(n))=(3k+1,3k+2).

Over all integer k, these pairs partition the integers not divisible by 3. Over k≥0, they partition the positive integers not divisible by 3. They do not include multiples of 3.

## 2. Which modulo-3 arrows are actually valid?

For the shortcut self-map Tε, the parity of the source must also be known. Its exact one-step transition relation modulo 3 is:

| Source residue | T₊ possible residues | T₋ possible residues |
|---:|---|---|
| 0 | 0 or 2 | 0 or 1 |
| 1 | 2 | 1 or 2 |
| 2 | 1 or 2 | 1 |

Proof: even sources map to 2r mod 3, while odd sources map to −ε mod 3. A residue modulo 3 by itself does not specify parity.

Consequently, “0,1→2→2” is not a deterministic description of T₊. For example, 6→3 stays at residue 0, and 2→1 goes from residue 2 to 1. Likewise, “0,2→1→1” does not describe T₋: 6→3 stays at 0 and 4→2 goes from 1 to 2.

For the unshortened map Cε, the odd branch instead goes to ε mod 3. This is another reason to state the convention before drawing arrows.

For Uε, the correct formula is

    Uε(n) ≡ ε·2^(−a) ≡ ε(−1)ᵃ (mod 3),
    a=v₂(3n+ε).

No output of Uε is divisible by 3. Either nonzero residue is possible from any fixed source residue modulo 3, depending on a. Even modulo 12 is insufficient: U₊(1)=1 but U₊(13)=5; U₋(3)=1 but U₋(15)=11.

### No fixed source modulus determines the accelerated output mod 3

This failure is not repaired by choosing one larger fixed modulus. Let M≥1 and L=lcm(M,2)=2ᴷm, where m is odd and K≥1. Choose an odd r satisfying 3r+ε≡0 mod 2ᴷ. Along n=r+Lt,

    (3n+ε)/2ᴷ = (3r+ε)/2ᴷ + 3mt.

Since 3m is odd, t can make the right side 1 mod 4 or 2 mod 4. The corresponding positive sources are congruent modulo M but have valuations K and K+1, hence opposite accelerated residues modulo 3. Taking sufficiently large t preserves positivity.

Valuation parity is therefore an indispensable side coordinate for this output observable. Knowing it decides the current output residue, but does not determine the valuation at the next iterate.

## 3. Exact signed coupling and its densities

The numerators 3n+1 and 3n−1 differ by 2. Both are even, and exactly one is divisible by 4. This proves min(a,b)=1 and max(a,b)≥2. Subtracting their factorizations proves the lifted identity

    2ᵃu−2ᵇv=2.

It has the two explicit forms

    a=1:  u=2^(b−1)v+1;
    b=1:  v=2^(a−1)u−1.

These relate the two maps at the *same original source*. They do not identify their subsequent orbits.

### Complete one-step coupling distribution

For either sign and every k≥1, the condition v₂(3n+ε)=k is exactly one odd residue class modulo 2^(k+1), since 3 is invertible modulo powers of two. Its relative density among odd n is 2^(−k). For k≥2 the opposite sign necessarily has valuation 1. Therefore

    density(a=1,b=k)=2^(−k),
    density(a=k,b=1)=2^(−k),    k≥2.

These pairs exhaust all odd sources. The tails beyond k>K have total relative density 2^(1−K). Truncation followed by this tail bound justifies summing the infinite collection in natural density, rather than merely in a formal probability model.

If a is even then a≥2, so b=1 and (u,v)≡(1,1) mod 3. If b is even then b≥2, so a=1 and (u,v)≡(2,2). If the high valuation is odd, the pair is (2,1). Hence

    density(1,1)=Σ_{k≥2 even}2^(−k)=1/3,
    density(2,2)=Σ_{k≥2 even}2^(−k)=1/3,
    density(2,1)=2Σ_{k≥3 odd}2^(−k)=1/3.

This immediately gives

    density(U₊≡2 mod 3)=2/3,
    density(U₋≡1 mod 3)=2/3.

The two signed outputs are not independent. In the quadratic character χ₃, with χ₃(1)=1 and χ₃(2)=−1, the paired signs are (1,1),(−1,1),(−1,−1), each of mass 1/3. In particular χ₃(U₊)≤χ₃(U₋) pointwise. Their means are −1/3 and 1/3, their product mean is 1/3, and their correlation is 1/2 under the starting-value distribution.

### Lifting to all powers of 3

Fix j≥1, k≥1 and a sign ε. On the exact valuation-k cylinder, the source residue modulo 3^(j−1) is uniformly distributed by the Chinese remainder theorem. Multiplication by 3 and then by the inverse of 2ᵏ sends these 3^(j−1) source classes bijectively onto the output residues modulo 3ʲ satisfying

    u≡ε(−1)ᵏ mod 3.

Thus, conditioned on k, every such residue has relative mass 1/3^(j−1). Summing 2^(−k) separately over even and odd k gives

    density(U₊≡u mod 3ʲ) = 1/3ʲ if u≡1 mod 3,
                             2/3ʲ if u≡2 mod 3,
                             0     if u≡0 mod 3.

For U₋ the first two weights are reversed. The same shrinking dyadic tail proves existence of each natural density. The higher ternary digits are conditionally uniform; the bias sits entirely in the first nonzero ternary digit.

This is a pushforward of uniform odd *sources*, including multiplicities of the map. Counting each distinct output value once is a different measure.

## 4. A genuine theorem for the combined-sign graph

For n>1 odd, the valuation-1 branch exceeds n:

    (3n+1)/2 > n,
    (3n−1)/2 > n.

The valuation-at-least-2 branch is at most (3n+1)/4<n. Therefore precisely one signed accelerated branch descends. The selector is

    ε(n)=+1 when n≡1 mod 4,
    ε(n)=−1 when n≡3 mod 4.

Let G(n)=U_{ε(n)}(n), and G(1)=1. Then every positive odd n reaches 1 by strict descent. Moreover,

    G(n)−1 ≤ (3/4)(n−1).

For odd n>1, it suffices to take any integer t with (3/4)^t(n−1)<2, because a positive odd state other than 1 is at least 3. A simple valid bound is t≤1+floor(log((n−1)/2)/log(4/3)). For general positive integer sources, first divide out their initial power of two.

This proves that the two-sign directed graph contains a path to 1 from every positive integer. It does not show that the path can use the plus sign throughout, nor that the minus sign throughout reaches a specified cycle. Choosing between signs changes the dynamical system. At n=1 both U₊ and U₋ equal 1; this is the sole equality of their positive odd outputs.

## 5. Signed conjugacy and the cycle boundary

Negation J(n)=−n gives exact conjugacies

    C₊∘J = J∘C₋,
    T₊∘J = J∘T₋,
    U₊∘J = J∘U₋.

For example, 3(−n)+1=−(3n−1), and v₂ is unchanged by sign. This identifies the positive 3n−1 system with the negative part of 3n+1. It is not a conjugacy between the two systems on positive integers.

The familiar accelerated positive cycles, directly verified here, are:

- U₊: 1→1.
- U₋: 1→1.
- U₋: 5→7→5.
- U₋: 17→25→37→55→41→61→91→17.

The last cycle has valuations (1,1,1,2,1,1,4). Under negation the three U₋ cycles become the three familiar negative U₊ cycles. Under T their periods are 1,3,11; under C they are 2,5,18. The positive plus cycle has periods 1,2,3 under U,T,C, respectively. Comparing cycle lengths without a common map convention mixes types.

These are known cycles, not a proof that the list is exhaustive. Neither “exactly one positive plus cycle,” “exactly three positive minus cycles,” nor “every start enters one of these cycles” follows from the residue diagrams. On all integers C₊ and T₊ also have the fixed point 0. Zero is not in the odd-only domain.

### What the affine cycle equation does explain

For a word of exact positive valuations a₀,…,a_{m−1}, let A=Σaᵢ and

    B=Σ_{i=0}^{m−1} 3^(m−1−i)2^(a₀+…+a_{i−1}) > 0.

The empty prefix sum is zero. Composing the branch formulas gives

    Fε(n)=(3ᵐn+εB)/2ᴬ.

A positive periodic point on that native word must satisfy

    n=εB/(2ᴬ−3ᵐ).

Thus the plus sign requires 2ᴬ>3ᵐ, whereas the minus sign requires 2ᴬ<3ᵐ. This gives an exact sign-sensitive cycle constraint. For the four displayed cycles (m,A,B) are (1,2,1), (1,1,1), (2,3,5), (7,11,2363). In the last case 2¹¹−3⁷=−139 and 2363=17·139.

A candidate from the affine equation still requires integrality, positive odd states, and verification of the declared exact valuations at every step. An affine fixed point alone must not be relabelled a native Collatz cycle without that check. The equation does not count all cycles or show eventual entry into them.

### Three meanings of “cover all integers”

1. Signed domains: positive integers, negative integers and {0} partition the integers. This is exact and says nothing about eventual behavior.
2. Paired one-halving images: S₋ and S₊ partition the *nonmultiples of 3*. Multiples of 3 are excluded from this image. Once an odd step has been taken, neither sign can create a new multiple of 3.
3. Basins of the displayed fixed-sign cycles: covering every integer by these basins is an unproved global claim. The controlled-sign path theorem in section 4 is a different, proved statement.

## 6. Squaring and the precise prime-power connection

The correct squaring diagrams are

    mod 3: 0→0, 1→1, 2→1;
    mod 4: 0→0, 1→1, 2→0, 3→1.

In particular, 0 mod 3 does not square to 1. On the unit classes in either modulus, squaring erases a quadratic-character sign:

    χ(x²)=χ(x)²=1.

The signed Collatz and prime-power computations share the exact character law χ(xᵏ)=χ(x)ᵏ, but have different sources, maps, and weights:

- Collatz: χ₃(Uε(n))=εχ₃(2)^(−a)=ε(−1)ᵃ. The exponent is an exact dyadic valuation, and its input-density weight is 2^(−a).
- The prime race modulo 4: χ₄(pᵏ)=χ₄(p)ᵏ. The exponent is the prime-power multiplicity; the logarithmic prime-power counting function uses weight 1/k and a cutoff pᵏ≤x.

For clarity, write Dπ(x)=π(x;4,3)−π(x;4,1), let πodd count odd primes, and set

    DJ(x)=−Σ_{p odd, k≥1, pᵏ≤x} χ₄(p)ᵏ/k.

Splitting by exponent parity proves the exact finite identity

    DJ(x)=Dπ(x)−½πodd(√x)
           +Σ_{odd k≥3}Dπ(x^(1/k))/k
           −Σ_{even k≥4}πodd(x^(1/k))/k.

The companion [prime-power bias and primitive-orbit note](prime_power_bias_and_primitive_orbit_bridge_20261014.md) keeps the prime-power weights and analytic scope explicit.

The prime-square term has a fixed sign because every odd prime square is 1 mod 4. This is the exact algebraic bridge to the prime-power bias discussion. An assertion about how often the prime race leads still requires analytic information about DJ and its oscillations. The identity alone does not supply that distribution.

Connection contract:

- Source objects: signed accelerated Collatz edges with their valuation labels; prime powers with their prime and exponent labels.
- Common target: the two-element group of quadratic-character signs.
- Maps: χ₃ of the Collatz output and χ₄ of the prime power.
- Preserved information: exponent parity and its effect on the sign.
- Destroyed information: source magnitudes, higher residues, integer trajectories, which prime supplied a power, and all counting weights unless kept explicitly.
- Required sidecars: exact valuation and source for Collatz dynamics; prime, exponent, cutoff and weight for the prime-power count.
- Decisive tests: the forbidden Collatz residue pair (1,2), the exact 1/3 coupling masses, and the negative square contribution −½πodd(√x) to DJ.

This is a common character calculation with explicitly different measures. No weight-preserving map from prime powers to Collatz orbits, dynamical isomorphism, or transfer of a cycle/convergence theorem has been established.

## 7. Inheritance, corrections, and source limits

The inherited local audit `collatz-recovery/collatz_wilson_audit.md` was read. It already stresses that a proved local dependency is not a global ROOT certificate, that affine word guards must remain native and source-preserving, and that arithmetic residue coincidences cannot be upgraded to convergence.

Bounded repository reconciliation recovered the relevant operating contract, navigation, active guardrails, research protocol, and corrections:

- MISTAKE-554 in `01-canon/MISTAKES.md` explicitly separates an explained modular coincidence from a genuine conjugacy and warns against mixing shortcut/non-shortcut codes.
- MISTAKE-605 in the same file warns that an affine fixed-point equation needs the native-word argument, and that matching a known numerical mechanism is not a novelty result.
- `05-knowledge/results/collatz_mod6_20260921_g_negatives_joint_carry.md` already inherits signed negation, the affine carry and cycle sign gate from `arithmetic_braids_20260917_collatz.md` and `arithmetic_braids2_20260917_signed_cycles.md`. These are predecessors, not new claims here.
- The repository also contains `04-computation/experiments/collatz_parity_mod3_characters_20261013.py` and `collatz_parity_mod3_roots_20261013.py`. Their paths were identified, but their contents could not be read in the bounded pass because the lazy-object network approval was cancelled. No claim of priority or repository novelty is made for the density or character results.

External primary literature checked:

- Jeffrey C. Lagarias, *The 3x+1 Problem: An Overview*, published 2010, arXiv version 2021: https://arxiv.org/html/2111.02635v1 . It treats periodicity and cycle finiteness as conjectural problems, distinct from local probabilistic models. Its survey's mathematical cutoff is 2010, not 2026.
- Lagarias's official bibliography identifies *The set of rational cycles for the 3x+1 problem*, Acta Arithmetica 56 (1990), 33–53: https://websites.umich.edu/~lagarias/3x+1.html . The publisher record confirms DOI 10.4064/aa-56-1-33-53. The available PDF download returned 403, so no unread theorem from it is used as a dependency here.

The elementary statements in sections 1–6 are proved directly above. The known-cycle list is verified by exact arithmetic rather than inferred from an inaccessible source.

## 8. Reproduction and checks

Repository artifacts:

- `04-computation/experiments/signed_collatz_pairing_exact_check_20261014.py`: standard-library-only exact checker, with explicit exceptions rather than removable assertions.
- `05-knowledge/results/signed_collatz_pairing_exact_check_20261014.json`: deterministic output.
- `05-knowledge/results/signed_collatz_pairing_exact_check_20261014.optimized.json`: independent replay with Python optimization enabled, byte-identical to the normal output.

Commands from the repository root:

    python3 04-computation/experiments/signed_collatz_pairing_exact_check_20261014.py > 05-knowledge/results/signed_collatz_pairing_exact_check_20261014.json
    python3 -O 04-computation/experiments/signed_collatz_pairing_exact_check_20261014.py > 05-knowledge/results/signed_collatz_pairing_exact_check_20261014.optimized.json
    cmp 05-knowledge/results/signed_collatz_pairing_exact_check_20261014.json 05-knowledge/results/signed_collatz_pairing_exact_check_20261014.optimized.json

The checker passes 1,244,038 checks. Its main signed-pair universe is every positive odd n<2¹⁸, totaling 131,072 sources. The observed residue-pair counts are (1,1):43,691; (2,1):43,691; (2,2):43,690. Exact density 1/3 is established by the proof and rational sums, not by treating these finite counts as exact thirds.

It also checks signed conjugacy on −2048≤n≤2048, the four displayed accelerated cycles, complete joint-valuation cylinders through k=14, all sign/valuation fibres for j=1,…,5 and k=1,…,8, hostile source-modulus witnesses for every M≤64, and strict controlled-sign descent for every positive odd source below 10,000. The largest observed controlled-sign step count in that last finite universe is 20, first attained at 9175.

Canonical output-payload SHA-256, excluding the hash field itself:

    9593a267d19216526bc14368c171959c85a787ebcc26e59b348e9a5f1758734b

Stopping boundary: local coupling, its coherent 3-power lift, and the controlled-sign theorem are proved; the fixed-sign global problem is unchanged. A useful next test of a proposed stronger connection is whether it preserves the ordered native valuation word and the original source, rather than only the residue sign.
