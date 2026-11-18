# Collatz branch clocks, multiplication modulo odd primes, and physical count renewal

**Status:** exact algebraic bridge and scope separation, with elementary derived marking obstruction; no Collatz convergence or Chan–Pak coverage proof. Most ingredients are recovered classical/repository material, not claimed new.

## Sources retained

* `modular_789_bridge/PALEY_COLLATZ.md`: the valuation-marked mod-seven branch action, native finite-word cylinders, and the distinction from unmarked integer dynamics.
* `math/05-knowledge/results/arithmetic_clock_tournaments_20261021/session24-arithmetic/README.md`: the affine trunk conjugacy, conductor-stratified skew sections, and Hamilton-path divisibility by a free group action. Its341 example demonstrates why unit-only reasoning fails for composite moduli.
* `arithmetic_reformulation/ARITHMETIC_RENEWAL_SELECTOR.md`: ordered strong-component renewal uses Dirichlet convolution in Hamilton counts; it does not act by affine changes of integer labels.
* Classical quadratic reciprocity and its supplementary laws: [Mathlib's formal theorem documentation](https://leanprover-community.github.io/mathlib4_docs/Mathlib/NumberTheory/LegendreSymbol/QuadraticReciprocity.html). The mod24 table below is a direct consequence, not an empirical rule.

## 1. Which Collatz information the multiplication table retains

Fix b=+1 or−1. On positive odd integers the accelerated branch with native exponent k≥1 is

    F_k(n)=(3n+b)/2^k,    k=v₂(3n+b).

For an odd prime p>3, its marked reduction is the affine permutation

    f_k(x)=2^(−k)(3x+b) modp.                         (1)

The same formula is meaningful for any odd modulus coprime to3, but the prime hypothesis is retained for the Paley construction below. A residue x does NOT determine k, so(1) is a family of marked charts, not a deterministic quotient of Collatz.

Every finite exponent word(k₁,...,k_L) occurs on exactly one odd residue class modulo2^(K+1), K=sum k_i. Proof: start with an odd terminal residue1 mod2. Working backwards through one letter, if the next residue is r mod2^s, solve

    3n+b ≡ 2^k r (mod2^(s+k)).

Since3 is a unit, this gives one residue n, necessarily odd, and forces exact valuation k. Iterate. Large positive representatives avoid the finite low boundary for b=−1. CRT then supplies infinitely many such finite-word starts in EVERY residue modulo p.

This does not supply a positive integer for an arbitrary infinite word. An infinite compatible collection determines a2-adic integer, which need not be an ordinary positive integer. Nor does a modular periodic point certify an integer cycle.

For a word, the affine slope is exactly

    a=3^L 2^(−K) modp.

Its translation retains the branch order; retaining only L,K loses that information.

## 2. An exact tournament orientation character

Let p≡3 mod4 and p>3. The Paley tournament on F_p has x→y precisely when y−x is a nonzero quadratic residue. Write χ for the quadratic character. Sinceχ(−1)=−1, each pair has exactly one direction.

For an affine map f(x)=ax+c,

    χ(f(y)−f(x))=χ(a)χ(y−x).

Therefore each branch(1) preserves or reverses EVERY arc according to

    χ(3)χ(2)^k,

and a whole word has orientation signχ(3)^Lχ(2)^K. The translation and b do not enter this orientation sign, although they still affect the actual permutation and fixed points.

The four Paley-prime classes give:

| p mod24 | χ(2) | χ(3) | one-letter orientation sign |
|---|---|---|---|
|7|+1|−1|−1|
|11|−1|+1|(−1)^k|
|19|−1|−1|(−1)^(k+1)|
|23|+1|+1|+1|

Hereχ(2)=(−1)^((p²−1)/8), and quadratic reciprocity givesχ(3)=+1 for p≡±1 mod12,−1 otherwise. The earlier mod-seven result is the first row of this general table. This is a precise interaction of mod8, mod3, and tournament converse, not a relation between these residue labels and desired Hamilton counts.

## 3. The generated group and the distinction between two symmetries

The branch permutations generate

    F_p ⋊ <2,3>  <= AGL(1,p).                         (2)

Indeed f₁ f₂^(−1) is doubling D(x)=2x. Then D f₁ is U(x)=3x+b. The maps DUD^(−1) and U differ in their translation, and (DUD^(−1))U^(−1) is translation by b. This generates all translations; the available slopes generate<2,3>. Conversely every branch lies in that affine group.

The subgroup<2,3> need not equal the entire group of square multipliers. We do not assume equality. Independently, the Paley tournament has the entire affine-square automorphism subgroup

    G={x↦ax+c: χ(a)=+1},   |G|=p(p−1)/2.

An automorphism fixing an ordered Hamilton path fixes every vertex, hence is the identity. Consequently

    p(p−1)/2 divides H(Paley_p).                      (3)

At p≡23 mod24 all Collatz branches belong to this arrow-preserving group. At other rows some reverse arcs. An anti-automorphism is not automatically a free action on ordered paths: combining it with path reversal can have fixed paths. Thus one must not double the divisor in(3) just because the signed branch group contains arc reversals.

Literal controls give H(Paley₇)=189 divisible by21, and H(Paley₁₁)=95095 divisible by55. At p=23 the branch affine group has order253, but this packet does not enumerate its Hamilton paths.

This explains a limitation of the analogy: importing arithmetic multiplication as an actual unbroken symmetry forces count divisibility. For p≥7 the Paley count cannot itself be prime, because the divisor p(p−1)/2 is composite. This says nothing against realizing prime counts by other, asymmetric tournament constructions.

## 4. How many labels must a physical perturbation remember?

Suppose a construction attaches new vertices to a Paley core using contacts invariant under the stabilizer of a marked core vertex. Its residual affine-square stabilizer has size(p−1)/2: every square multiplier has a unique translation fixing that vertex. If those symmetries extend to the output, they still act freely on its ordered Hamilton paths, so

    (p−1)/2 divides H(output).

For p7,11,23 these residual divisors are3,5,11. A one-point-marked construction cannot simply forget this divisor while retaining that symmetry.

Two DISTINCT individually fixed core labels eliminate the affine stabilizer: subtracting their fixed-point equations gives a=1, then c=0. Thus a perturbation using genuinely two-point-dependent contacts can remove this particular group obstruction. It does not prove that the full automorphism group is trivial, that the chosen contacts actually remove every symmetry, or that any desired count results. It identifies the information that a potentially useful physical repair must retain.

This is the constructive reading of the modular clock: use its predictable arithmetic to organize the core, but distinguish enough contacts to prevent its divisibility from surviving the read. The earlier equal-gain return repairs do exactly this kind of nonuniform reading, although no claim is made that those small cores are Paley cores.

## 5. Why ordered tournament renewal is not integer Collatz iteration

Strong components compose by ordinal sum and multiply their Hamilton counts. Their isomorphism-class inventory obeys the divisor recurrence

    a_n(N)=sum_k sum_{d|N} s_k(d)a_{n−k}(N/d).

The factors are actual core counts. Collatz branch clocks instead act affinely on LABELS modulo p, with valuation data needed to lift them back to integers. A clock orbit is not a list of strong components, and its modular period is not a Hamiltonian path count. Neither equation transforms the other without an additional physical construction.

The analogy still points to a precise target: find a target-decreasing physical core move whose continuation coordinates play the role of the native valuation k. A reduction of numerical H alone is no more adequate than reducing Collatz modulo an odd number while discarding its valuation. The missing theorem is simultaneous preservation of legality and exact target height, not another modular clock computation.

## Reproduction

`python collatz_modular_replay.py` checks the orientation character for every Paley prime7≤p<200, both signs and exponents1..6; verifies all exponent words of length≤4 with letters1..3, exact cylinders and CRT lifts; computes literal Hamilton counts for Paley7/11; and records the branch multiplier groups at7,11,19,23. Normal and optimized outputs agree. All-order claims above have separate algebraic proofs.
