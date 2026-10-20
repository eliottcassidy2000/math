# Independent audit: generator-free circulant naturality

2026-10-20. **PASS for the stated four theorems and their stated scope.**

The audit checks the mathematical arguments, a genuinely separate exact implementation, the inherited prime-order identification, the q=9 hostile example, and the cited primary circuit lower bound. It makes no historical novelty determination and performs no broad census. No producer, repository-checkout, or Git files were edited by this reviewer.

The final reviewed `GENERATOR_FREE_NATURALITY.md` has SHA256 `36a1b2fd5ad0910828cdf3bd0a3fea07dcc6f5b90661fab6fbe61852c91449fd`. The producer applied both precision corrections before this final hash was verified.

## 1. Verdict and corrections

All four main claims pass:

1. For odd m, the parity projection is a generator-free, C-equivariant surjection with precisely complementary two-element fibers and an Aut(C)-equivariant type bijection.
2. For even m at least 4, the automorphism with exponent m+1 fixes every target but moves a source type, prohibiting an equivariant injection. The m=2 exception is genuine.
3. In even-degree finite-field extensions, Frobenius has no fixed source multiplier type, while the empty target type is fixed. This obstructs a Galois-equivariant bijection in the multiplier/even-necklace model.
4. The particular parity projection is not uniformly first-order definable over all primes congruent to 3 modulo 4. This does not prohibit other first-order bijections.

Two wording corrections were requested from and made by the producer. The Reynolds averaging operator is J, while P=I+J is its complementary augmentation projector. The group of translations has order p^r; an individual nonidentity translation has order p. Neither correction changes a theorem or any test result.

## 2. Proof audit

### A. Odd-order projection

For odd m, H and the involution z are characteristic, and C=H times the order-two subgroup. The map P(s)=s+(sum s)1 has image ker(sum), because m=1 in F_2, and kernel precisely the constant functions. The target primitive parity equals total parity because every repetition factor divides odd m. Consequently every target has exactly two lifts, exchanged by z. This also proves injectivity on C-orbits: after translating two outputs to agree, their lifts differ only by z. Reindexing H proves both multiplier and automorphism naturality.

For p=2m+1 prime, the inherited Sylow argument is sound. The translation subgroup C_p is Sylow in S_p, and conjugating Sylow subgroups inside the relevant automorphism group reduces graph isomorphisms to affine maps. For sign detection, the sign kernel has index two and contains the odd-order translations; conjugating inside that kernel moves an odd-sign automorphism into the affine normalizer without changing its sign. Thus the distinction between multiplier types and actual graph types disappears for these prime-order objects. For odd m, |A| divisible by four is equivalent to even primitive-root weight.

### B. Even-order abstract-cyclic obstruction

For even m, gcd(m+1,2m)=1. The automorphism is identity modulo m and hence fixes every symmetric connection word. On f=0^m1^m it toggles precisely the odd positions. The cyclic transition count becomes 2m-2 instead of 2, so its source type moves for m at least 4. Equivariance and injection are incompatible. The stronger no-fixed-source-type assertion when 4 divides m also passes its separate odd/even-shift argument.

The obstruction forgets field addition. Its nontrivial cyclic automorphism is generally not a prime-field automorphism. Therefore it cannot be used to exclude arbitrary ring-defined constructions over prime fields.

### C. Even-degree Frobenius obstruction

A fixed multiplier orbit would supply a semilinear permutation L(x)=a x^p preserving the Cayley tournament. Its image in Gal(F_(p^r)/F_p) has order r, so its order is divisible by even r. An involutive power of L would transpose a vertex pair and reverse that pair's directed edge. This is impossible in a tournament. The empty graph multiplier type is fixed, proving the fixed-point mismatch. No counting equality is needed for this contradiction.

A parameter-free first-order output map commutes with Frobenius on the structures with the input predicate supplied. If it induces the stated multiplier-type bijection, the induced bijection also commutes with Frobenius and contradicts the mismatch. Prime-field parameters remain fixed and do not change the argument. Arbitrary extension-field parameters are explicitly outside the claim.

### D. Nondefinability of this map

Choose arbitrary Boolean data on the m nonzero squares H, extend to a skew selector by complementation on -H, and set S(0)=false. At x=1, the displayed output is exactly parity of the m-1 data bits in H other than 1. There are arbitrarily large primes congruent to 3 modulo 4, so these input sizes are unbounded.

Unrolling a fixed first-order formula over a universe of p elements gives bounded-depth, polynomial-size AND/OR/NOT circuits. Every ring operation and every extra input-independent relation, constant, or field-dependent parameter value can be hardwired. An atomic S-term becomes a single Boolean input, its negation, or a constant after the skew-input substitution. This contradicts the nonuniform parity lower bound even though the relevant input lengths form a subsequence. The proof needs no uniformity condition on the hardwired circuit family.

The conclusion concerns the parity projection itself. It does not establish nonexistence of some other first-order prime-field bijection. A supplied primitive generator makes the different local derivative map first-order definable, which is fully consistent with this result.

## 3. Independent exact controls

`independent_controls.py` imports no producer code. It uses integer subset masks and image permutations for exhaustive orbit calculations, GF(2) row elimination for source affine fixed words, and direct Möbius inversion for target exact-period parity. These differ from the producer's tuple/precomposition representation, signed-component method, and period-recursion implementation.

The independent normal and optimized runs are byte-identical and pass **176,751 explicit checks**. The producer checker was also independently replayed normally and under optimization; both outputs agree with its submitted JSON and retain **183,074 explicit checks**.

Controls cover every m from 1 through 13 and every unit multiplier modulo 2m. They verify all fixed-type counts, total counts, every odd-m projection input, exact complementary fibers, multiplier/automorphism equivariance, type bijectivity, the output bit used in the parity reduction, the even-m moved-source witness, and the m=2 exception. For m at most 8, every individual affine fixed-word count is also compared with direct exhaustive enumeration, before Burnside averaging.

The separately computed prime-power counts are:

| q | Common total types | Frobenius-fixed source | Frobenius-fixed target |
|---:|---:|---:|---:|
| 9 | 2 | 0 | 2 |
| 25 | 172 | 0 | 46 |
| 27 | 316 | 16 | 16 |
| 49 | 349536 | 0 | 3372 |
| 81 | 13743895360 | 0 | 2522 |

The q=27 odd-degree case is a material control, not evidence of an all-odd-degree theorem.

### Independent formulas

For any finite C-set Z and automorphism alpha normalizing its C-action,

    |Fix(alpha on Z/C)| = (1/|C|) sum_(c in C) |Fix(c alpha on Z)|.

Indeed a nonstable orbit contributes zero. In a stable orbit, each point contributes the size of its stabilizer to the sum over c, so the orbit contributes exactly |C|.

For the source, use m independent bits x_i and enforce f_(i+m)=1+f_i. Invariance under i -> a i+k becomes the linear equations

    x_i + x_(j mod m) = floor(j/m),  j=(a i+k) mod 2m.

An inconsistent system contributes zero; otherwise it contributes 2^(m-rank). This avoids relying on the producer's component-count procedure.

For the target, let c_e be the number of cycles of i -> a i+k modulo e. There are 2^c_e fixed words. For a length-e word repeated d/e times, the resulting length-d weight is automatically even if d/e is even; otherwise its weight-even condition is one nonzero linear constraint exactly when an odd cycle occurs. Thus define

    B(e,d) = 2^c_e, if d/e is even or all cycles have even length;
             2^(c_e-1), otherwise.

Direct inclusion-exclusion gives the number of fixed words of exact primitive period d and even primitive weight:

    E_d = sum_(e divides d) mu(d/e) B(e,d).

The target labelled fixed count is sum_(d divides m) E_d. Burnside averaging then gives its fixed-type count. This direct Möbius formula agrees with exhaustive enumeration and with the producer's different recursion.

The common total type count was also checked against

    T(m) = (1/(2m)) sum_(d divides m, d odd) phi(d) 2^(m/d).

One independent derivation differentiates an antiperiodic word once to obtain an odd-total-weight length-m word; complementary integrations are the same source orbit. Burnside on odd-total-weight words allows a rotation only when its repetition factor d is odd, leaving 2^(m/d-1) valid cycle-label assignments. This yields the displayed closed formula.

## 4. Prime-power hostile and scope boundary

The reviewer independently implemented F_9 as a+b t with t^2=2, using integer encodings. It verified that g=1+t generates the nonzero group and that exponent positions 2,3,6,7 give precisely {plus or minus t, plus or minus (1+2t)}. In the indicated basis the graph is the 3 by 3 rook graph, with 9 vertices and 18 edges.

The coordinate involution (u,v) -> (-u,v) preserves all edges and swaps the two endpoints of exactly three edges. Under the numerical reference orientation it reverses nine edges, also odd. Therefore the graph is not fully sign-even. This is a concrete counterexample to interpreting the even-root target as full sign-even graphs in extension fields.

There are two source multiplier types at q=9, and the explicit Frobenius action swaps them. They therefore already become one type under the larger semilinear isomorphism group. This directly confirms that enlarging the equivalence relation changes the source type count. Theorem C must not be promoted to a result about full-isomorphism types or cyclic groups of order p^r. The note keeps those boundaries explicit.

## 5. Primary-source verification

Verified directly on 2026-10-20:

- Johan Håstad, *Almost Optimal Lower Bounds for Small Depth Circuits*, STOC 1986, pages 6–20, DOI https://doi.org/10.1145/12130.12132. The ACM publisher record confirms title, author, venue, year, and DOI.
- Author-hosted full version: https://www.csc.kth.se/~johanh/largesmalldepth.pdf. Section 2.1 states the unbounded-fanin AND/OR/NOT model; Section 4, Theorem 1 and the following corollary on polynomial-size parity circuits provide the required lower bound. Its bibliography identifies the STOC 1986 paper as the earlier version.

Only the consequence excluding fixed-depth polynomial-size parity circuits is imported. No novelty or priority assertion is required by the audit.

## 6. Reproduction

    python3 independent_controls.py --output independent_results.json > independent_results.out
    python3 -O independent_controls.py --output independent_results.optimized.json > independent_results.optimized.out
    cmp independent_results.json independent_results.optimized.json

The producer replays are saved separately as `replayed_producer_results.json`, `replayed_producer_results.out`, and `replayed_producer_results.optimized.out`. Hashes of the exact reviewed producer files and the review artifacts accompany this report.

Finite computations validate the implementations and examples. The general assertions rely on the checked arguments, with Theorem D additionally using the cited circuit theorem.
