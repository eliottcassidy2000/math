# Wieferich primes as failures of tournament-lift rigidity

Date: 2026-10-14. Status: elementary theorems with exact finite controls; one explicitly open hypothesis.

## Result in brief

There is a precise combinatorial distinction at a base-2 Wieferich prime: a natural family of cyclic tournaments on `Z/p²Z` acquires new independent orientation choices which are absent at an ordinary prime. At an ordinary prime every tournament in the family is a lexicographic product across the residue classes modulo `p`; at a Wieferich prime explicit members break those modules. For 3511 the number of labelled connection sets jumps from the ordinary formula `4` to `2^3512`. For 1093 the correct odd-order multiplier is 16, and the analogous count is `2^6564` versus the ordinary formula `2^12`.

Separately, the repository's odd-cycle collection formula contains the Fermat quotient as an exact first-order correction when one compresses a prime-sized orbit of cycles. This explains where base 2 can enter a Hamiltonian-path calculation at prime-square precision.

Neither statement derives the numerical pair 1093,3511 from the forbidden Hamiltonian-path values 7,21. Those gaps have a known structural proof. The new bridges identify a shared carrier and a testable next question, without pretending that the two pairs classify one another.

## 1. What 7 and 21 actually mean

Throughout, `H(T)` counts directed Hamiltonian paths as ordered vertex lists, with every consecutive arc in its forward direction. Reverse lists are not identified. Directed cycles below are identified under cyclic rotation, not reversal, and have length at least three.

The relevant canon is [THM-1370](../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md), not a congruence excluding multiples of seven. Its proof uses:

1. `H(T)=∏ H(C)` over the strong components, in condensation order.
2. The minimum for a strong tournament of order `n` is

   `f(n)=3·5^((n−3)/3)` for `n≡0 (mod3)`,

   `f(n)=5^((n−1)/3)` for `n≡1 (mod3)`,

   `f(n)=9·5^((n−5)/3)` for `n≡2 (mod3)`.

3. The strong spectra at orders 3,4,5 are `{3}`, `{5}`, `{9,11,13,15}`. At order 6 they omit 21. Already `f(7)=25`.

Thus 7 fits neither a strong atom nor a product of smaller permitted strong atoms. The same argument excludes 21: a factorization `3·7` inherits the first obstruction, while a strong atom of value21 cannot fit the small spectra or the floor beyond order6. This is a positive-enumeration and factorization obstruction, not a residue-class obstruction. For example, 63 is attainable in the canon.

The floor is the Moon–Busch theorem; [Busch's primary paper](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/) gives the component product and the floor proof. The small spectra through order6 are independently reproduced by the attached script. The fact that 7 and21 are missing, together with the candidate spectrum of all other positive odd values, was already discussed in [the 2016 MathOverflow question and Gordon Royle's computations](https://mathoverflow.net/questions/232751/the-number-of-hamiltonian-paths-in-a-tournament). Global completeness remains a conjecture; the note makes no novelty claim for the two omissions.

## 2. Exact arithmetic, with the correct carrier

For an odd prime `p`, write

`q_p(2) = (2^(p−1)−1)/p mod p`,

`d = ord_p(2) = 2^s δ`, with `δ` odd,

`e=2^s`, `a=2^e`, and `m=(p−1)/(2δ)`.

The normalization removes precisely the even part of the multiplicative order. It does not replace the base-2 Wieferich condition: `q_p(a)=e q_p(2)`, and `e` is invertible modulo `p`. Put

`w=v_p(2^(p−1)−1)=v_p(a^δ−1)`.

The standard order-lifting identity is

`ord_(p^r)(a)=δ p^max(0,r−w)`.

Indeed, `a^δ=1+p^w c`, with `p∤c`. The binomial theorem, or the elementary odd-prime lifting lemma, gives `v_p((a^δ)^(p^j)−1)=w+j`. An exponent returning to1 modulo `p^r` must first be divisible byδ, and this valuation gives the remaining factor.

Exact checks:

| p | ord_p(2) | ord_(p²)(2) | ord_(p³)(2) | a | δ | m | w |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 7 | 3 | 21 | 147 | 2 | 3 | 1 | 1 |
| 13 | 12 | 156 | 2028 | 16 | 3 | 2 | 1 |
| 31 | 5 | 155 | 4805 | 2 | 5 | 3 | 1 |
| 1093 | 364 | 364 | 397852 | 16 | 91 | 6 | 2 |
| 3511 | 1755 | 1755 | 6161805 | 2 | 1755 | 1 | 2 |

The next nonzero Fermat digits are

`(2^1092−1)/1093² ≡487 (mod1093)`,

`(2^3510−1)/3511² ≡51 (mod3511)`.

Consequently both known examples have depth exactly2, rather than an untested higher lift. The corresponding digits at the shorter order exponents are891 and1781.

A crucial hostile control: multiplication by2 cannot be a tournament automorphism on `Z/1093Z`. Its order364 is even, so one of its powers is `x↦−x`, which exchanges the endpoints of an arc. More generally, a nonidentity involution cannot be a tournament automorphism: a transposed pair would reverse its unique arc. Multiplication by16 has odd order91 and avoids this obstruction. At3511 multiplication by2 already has odd order1755.

[Conrad's notes](https://kconrad.math.uconn.edu/blurbs/ugradnumthy/wieferich-primes.pdf) give the standard definition and the two known examples. [Dorais–Klyve's published search](https://cs.uwaterloo.ca/journals/JIS/VOL14/Klyve/klyve3.pdf) certifies the historical bound `6.7·10^15`; that is not asserted to be the latest search boundary. The live [PrimeGrid results](https://www.primegrid.com/stats_ww.php), checked for this note, report no zero-quotient discoveries in that project's results. None of these observations proves that the known list is complete.

## 3. A complete tournament-lift characterization

### Theorem A: orbit choices encode the full Wieferich depth

A Cayley tournament on `Z/p^kZ` is specified by a set `S` of nonzero residues satisfying

`S ∩ (−S)=∅`, `S ∪ (−S)=(Z/p^kZ)\{0}`,

with arc `x→y` iff `y−x∈S`. Require additionally `aS=S`.

The number of such **labelled connection sets on this fixed cyclic group** is `2^B_k`, where

`B_k = m Σ_(r=1)^k p^min(r−1,w−1)`.

These are not counts up to isomorphism.

**Proof.** The residues of additive order `p^r` number `φ(p^r)`. Their multiplication-by-a orbits all have length `ord_(p^r)(a)`. That length is odd, so no orbit meets its negative. Choosing one orbit from each negative pair determines the tournament and exhausts the possibilities. The contribution of this layer is

`φ(p^r)/(2 ord_(p^r)(a)) = m p^min(r−1,w−1)`.

Summing over the layers proves the formula. ∎

For an ordinary prime, `B_k=mk`. At depth2, `B_k=m[1+p(k−1)]` for `k≥2`. Thus the entire order-lifting depth has a precise combinatorial readout as the growth of the number of independent orientation bits.

### Theorem B: Wieferich is exactly failure of residue-module rigidity

With the normalized multiplier above, the following are equivalent:

- `p` is base-2 Wieferich;
- some a-invariant Cayley tournament on `Z/p²Z` has a residue class modulo `p` which is not a module;
- some such tournament is not a lexicographic product with respect to the fixed residue-class partition.

A module means that every vertex outside the set has the same arc direction to every vertex inside it. The qualification about the **fixed partition** matters: no assertion is made that every such example lacks all other decompositions.

**Proof in the ordinary case.** Here `ord_(p²)(a)=pδ`. A unit orbit reduces to a δ-element orbit modulo p and has sizepδ, so it is the entire inverse image of that orbit. Its orientation is therefore constant on every nonzero residue class of differences modulo p. Cross-fibre arcs are uniform. Differences divisible byp specify an independent a-invariant tournament on each fibre, the same in every fibre. Hence the tournament is exactly `B[C]`, for two a-invariant Cayley tournaments B,C on p vertices. There are `2^m·2^m=2^(2m)` choices.

**Proof in the Wieferich case.** Unit orbits now have sizeδ. Each unit orbit downstairs has p separate lifts. Begin with any lexicographic choice whose connection set downstairs contains `<a> modp`. Let `K=<a> modp²`. Reverse the orientations associated with K and−K, leaving all other differences unchanged. This preserves both the tournament condition and a-invariance.

Reduction `K→<a> modp` is a bijection. Hence1 belongs toK but neither `1−p` nor `1+p` belongs toK. After the flip, `1→0`, whereas `p→1` still holds. The vertices0 andp lie in the same residue class; vertex1 sees them in opposite directions. That class is not a module. ∎

Exactly `2^(2m)` of the `2^(m(p+1))` connection sets in the Wieferich case remain lexicographic over that partition. Consequently there are exactly

`2^(m(p+1)) − 2^(2m)`

non-module examples in this specified labelled family.

At3511 the counts are `2^3512` total and4 lexicographic. At1093 with multiplier16 the counts are `2^6564` total and4096 lexicographic. The script verifies the explicit flipped-orbit witness using only91 or1755 modular powers, rather than materializing tournaments with more than a million vertices.

This is a structural characterization of an arithmetic predicate, not an algorithm which makes its classification easier. It also concerns the ring `Z/p^kZ`, not the finite field `F_(p^k)`.

## 4. The actual Hamiltonian-path consequence

Every tournament in Theorem A admits the affine automorphism group

`G_k=(Z/p^kZ) ⋊ <a>`.

This group acts freely on ordered Hamiltonian paths: an automorphism fixing an ordered path fixes every vertex. Therefore

`δ p^(k+max(0,k−w)) | H(T)`.

At prime-square level an ordinary prime forces `p³δ | H(T)`; a Wieferich prime forces only `p²δ | H(T)` by this particular group. “Only” here describes what the displayed group proves, not an exact valuation of the integer H. Other symmetries and arithmetic cancellation may supply more factors.

The difference has a direct structural meaning. Ordinary lifting builds a new order-p multiplier symmetry and forces uniformity between residue fibres. Wieferich lifting fails to build that symmetry, freeing p separate unit-orbit choices. Theorems A–B say exactly where that freedom resides; the displayed divisibility says what is lost in the compulsory path-orbit factor.

### Open hypothesis H-W: minimal forced p-valuation is attained

For every odd prime p and integer a with `gcd(a,p)=1`, `ord_p(a)` odd and `ord_(p²)(a)=ord_p(a)`, there exists an a-invariant Cayley tournament T on `Z/p²Z` with

`v_p(H(T))=2`.

The two base-2 targets are consequently `p=1093,a=16` and `p=3511,a=2`. A counterexample must show that **all** admissible connection sets have an additional p factor. No computation of H at either target order is claimed.

There is a useful genuinely small positive/negative control. At `p=3,a=10`, the order is1 modulo3 and9, so the 16 cyclic tournaments on9 vertices are all admissible. Their H-values and multiplicities are

- `3159=3^5·13`, four connection sets;
- `3267=3^3·121`, six connection sets;
- `3357=3²·373`, six connection sets.

Thus the proposed minimal valuation occurs, but certainly not for every orientation. For the ordinary normalized base `a=4`, only the four lexicographic choices are admissible, and all haveH3159, exceeding even their forced affine p-valuation3. This prevents turning the group argument alone into the conjecture's proof. The base10 example is a toy for the general lift mechanism, not an additional base-2 Wieferich prime.

## 5. Fermat quotients occur inside the odd-cycle collection formula

The repository's [THM-002](../../01-canon/theorems/THM-002-ocf.md) is

`H(T) = Σ_C 2^|C|`,

where C runs over collections of pairwise vertex-disjoint directed odd cycles, including the empty collection. Its external foundation is the [Grinberg–Stanley Rédei–Berge formula](https://arxiv.org/abs/2307.05569), also derived in matrix-algebra form by [Irving–Omar](https://arxiv.org/abs/2412.10572). The weights here are genuinely the repo's binary weights; no artificial encoding of a Fermat congruence is needed.

Let τ be an automorphism of T of odd prime orderp. It acts on cycles, and therefore on their disjoint collections. A τ-fixed collection C decomposes uniquely into:

- `u(C)` cycles fixed individually byτ;
- `t(C)` orbits of cycles, each containingp distinct cycles.

Define the integer sums

`F = Σ_(C fixed) 2^(u(C)+t(C))`,

`D = Σ_(C fixed) t(C) 2^(u(C)+t(C))`,

`O = Σ_(free collection-orbits [C]) 2^|C|`.

The last sum uses one representative per free orbit; the weight is invariant, so O is independent of the choice of representatives. All sums are finite and canonical for the pair(T,τ).

### Theorem C: the first-order OCF compression identity

`H(T) ≡ F + p[q_p(2)D + O] (mod p²)`.

**Proof.** A fixed collection has actual weight

`2^(u+pt) = 2^(u+t) (2^(p−1))^t`

`≡ 2^(u+t)(1+p·t·q_p(2)) (modp²)`.

Each free collection-orbit contributes exactly `p·2^|C|`. Summation gives the identity. ∎

The factorq is thus the first-order failure of replacing p equal binary cycle choices by one binary choice. At a Wieferich prime that correction vanishes. The free-orbit termO and the compressed termF remain, so the identity is not an unconditional detector `p²|H ⇔ p Wieferich`.

A polynomial version is useful if cancellations hide the correction. Put

`P_(T,τ)(x)=Σ_(C fixed) 2^u(C) x^t(C)`.

Then the fixed-collection contribution is `P_(T,τ)(2^p)`, while

`P_(T,τ)(2^p) ≡ P_(T,τ)(2) + 2p q_p(2) P'_(T,τ)(2) (modp²)`.

In particular `D=2P'(2)`. This is the exact derivative channel to retain in any proposed modular OCF compression.

### Hostile controls where every correction is visible

The script freezes two nine-vertex tournaments, with τ equal to translation by3 on labels0,…,8. It independently checks H by subset dynamic programming and by full odd-cycle collection enumeration. A mask uses pairs `(0,1),(0,2),…,(7,8)` in that order, with bit1 meaning `i→j`.

| mask | H | F | D | O | F mod9 | 3qD mod9 | 3O mod9 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 22024714438 | 2715 | 39 | 4 | 888 | 3 | 3 | 0 |
| 56041042850 | 1557 | 39 | 4 | 502 | 3 | 3 | 3 |

Here `q_3(2)=1`. The first example fails if the Fermat correction is dropped. The second fails if the free-orbit contribution is dropped, and it has `9|H` even though3 is not Wieferich. These are literal finite obstructions to the most tempting overclaims.

Both fixed-collection polynomials are `P(x)=35+2x`; thus `F=39` and `D=4` are not accidental unreduced summaries. Their free collection-orbits differ, precisely the information a compressed quotient would otherwise lose.

## 6. Where the original gap question now stands

Three distinct layers should remain visible:

1. The permanent omissions7,21 come from permitted positive strong-component counts and an extremal floor.
2. Wieferich primes mark loss of a multiplicative orbit-lifting constraint, creating new tournament orientation freedom.
3. The OCF formula puts the same base-2 Fermat quotient into the derivative of a natural cycle-orbit weight compression.

The third statement is the closest direct link to the mechanism behind Hamiltonian-path counting. It gives a concrete research instruction: retain both the derivative channelD and the free-orbit channelO, and ask whether the new non-module orientations from Theorem B can make the path quotient `H/(p²δ)` nonzero modulo p. That is precisely H-W, not a classification of Wieferich primes or a claim that7,21 cause1093,3511.

For comparison, the nearly transitive tournament R_n satisfies `H(R_n)=1+2^(n−2)`: its nontrivial paths are indexed by the subset of interior vertices lying before the reversed extreme arc. Consequently the first n−2 with `H(R_n)≡2 (modp^r)` is `ord_(p^r)(2)`. This is an exact useful calibration, but Theorems A–C supply more structural information than that direct exponential readout alone.

## Reproduction and scope

Run `python 04-computation/experiments/wieferich_tournament_lifts_20261014.py`. The paired JSON records:

- exact orders modulo p,p²,p³ and nonzero second digits;
- explicit one-orbit-flip witnesses at1093 and3511;
- all labelled tournament spectra through6 vertices;
- near-transitive path counts through10 vertices;
- all16 cyclic connection sets at order9;
- the two OCF controls, including full adjacency matrices and independent H checks.

The proofs are elementary deductions from established order lifting, Cayley orientations, and OCF. No exhaustive prime search, giant tournament path count, or novelty survey for the general Cayley enumeration formula is claimed. Repository foundations used here were read from checkpoint `0f7b1b200b1b96a6ce61bbe304bb92c603515bbf`; only the new staged files are intended for publication.
