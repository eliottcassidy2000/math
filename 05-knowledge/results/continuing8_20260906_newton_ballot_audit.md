# Independent referee: separated Newton clusters and global ballot repair

**PASS — full analytic scope accepted; no mathematical repair requested.** The two producer candidates may be promoted in the precise scopes below. This audit distinguishes the all-profile proofs from their finite controls and the constructed factorial consumer from an unchanged coefficient row.

Audited frozen candidates:

| Artifact | Source SHA-256 | Output SHA-256 | Certificate SHA-256 |
|---|---|---|---|
| `continuing8_20260906_newton_clusters` | `29c4fdc5614380b11e1e69770fc141e831e27f35618b681aae0ced3bb439775e` | `3d980a902fd23299521e6a3d51b3878d787760853ebcefbe1cf09d7668282281` | `1d35470321fbc24289b62fb49955cafe9e10e7518f26b9ae17f59ee12145e17e` |
| `continuing8_20260906_ballot_repair` | `3d537f9a32a4b6f18c152d1986231a2fe4f5def78de0e072f558ad78b434902c` | `42f65a43f6e9a1ff3c074a2cd28cffbb759603eed16cc34a978a97dd7a043adc` | `151ddbd2240c210ccaf87f5f40f101439af2d58bf4fcebd765a6d6c06cc379d2` |

The proof-report hashes read were respectively `8f9f7962e637e135bf1ebde8330d96648f06bd1423a9b237dfdc89cc1c25a657` and `8fc96200857bd6c371735f71cf3b1dabd098ead3d187d484930b87d236e92fc6`. A later status/link-only filing promotion does not change the audited mathematics. The referee imports no producer implementation.

## 1. Uniform cluster theorem

The theorem fixes a finite profile `m_1,...,m_K`, with `K>=2` and every `m_j>=2`. It then constructs finite constants `T_0` and `epsilon_0` such that **every** choice of bases with adjacent ratio at least `T_0`, and **every** choice of root parameters in the stated relative-width boxes, has exactly `2K-3` zero-filtered circuit changes. This order of quantifiers is proved, rather than inferred from the fifteen-profile bank.

I independently derived the limiting ratios. Inside a block, with `A` roots before it, `B` after it, size `m`, and `k=A+l`, the limiting normalized Newton ratio is

`kappa_(A+l) = [(l+1)(m-l+1)/(l(m-l))] * [k(d-k)/((k+1)(d-k+1))`

`= f_A(l) f_B(m-l)`, where `f_A(l)=(l+1)(A+l)/(l(A+l+1))`.

For `A>0`, `log f_A` is strictly decreasing and strictly convex. A direct proof of the second fact uses the strictly decreasing function `h(x)=x^(-2)-(x+1)^(-2)`: the second derivative is `h(l)-h(l+A)>0`. This proof does not use the producer's symbolic positive-coefficient expansion. The `A=0` factor is exactly one. Therefore the first block's bounded ratios increase, the last block's decrease, and the interior circuit ratios in each middle block are strictly increasing.

The coefficient envelope retains the source normalization. For a size-`k` subset, every prefix of its cluster-count vector is bounded above by the top-filled count vector. Unless its vector is top-filled, at least one integer prefix deficit is negative. Expressing the monomial ratio as a product of adjacent base ratios therefore bounds it by `1/T_0`. There are no more than `2^d` subsets, and each product of within-cluster weights lies in `[1,(1+epsilon_0)^d]`. This proves the envelope for all coefficients, including `e_0` and `e_d`.

The explicit error budget is valid. With `epsilon_0=1/(6Md)` and `T_0>=3M2^d`, one may use the slightly sharper referee majorant

`(1+epsilon_0)^d (1+2^d/T_0) < (1+1/(3M))/(1-1/(6M)) < 1+1/M`.

The last inequality holds for every integer `M>=1`; it is not asymptotic. The producer's alternative product of two `1+1/(3M)` bounds is also sufficient. Thus each coefficient correction is in `[1,tau)`. A circuit correction has four positive and four negative total exponents, giving `(tau^-4,tau^4)`. Comparing adjacent circuits costs at most `tau^8`, exactly as retained in the margin set.

The finite margin set contains every needed comparison: first-band signs, last-band signs, and adjacent bounded-interior circuit ordering. Each member exceeds one. Consequently the dyadic search for `M` terminates for **every** allowed profile. The two spike requirements at each boundary use the same explicit `T_0>tau^4 S`, giving a strict positive circuit followed by a strict negative one. Unequal separation ratios do not invalidate these estimates.

## 2. Exact count, zero cases, and exclusions

The count is complete, not merely a lower bound from alternating spikes. There are `K-1` positive-to-negative boundary changes. Between successive boundaries there is a negative circuit, a strictly increasing bounded-interior portion, and a positive circuit. After deleting zeros this contributes exactly one further change. This argument includes empty and one-element interiors, hence cluster sizes two and three. It permits an exact zero at the change without counting that zero twice. The first and last portions contain no additional changes. The total is `(K-1)+(K-2)=2K-3`.

All indices are available because each multiplicity is at least two. The case `K=2` has one boundary and no middle block. When all multiplicities are two, the spikes exhaust the `d-2` circuits, giving `d-3` changes, the maximal possible count. Their strict inequalities persist on open sets of distinct parameters. No conclusion is drawn for arbitrary wide clusters or for an arbitrary number of distinct roots at arbitrary separation.

The source's named exclusions are genuine. The `(1,2,1)` profile loses the multiplicity hypothesis and has only one change. A single narrow cluster can have a change, so the `K>=2` scope must remain. The old `(1,1,3,3,8)` example still refutes a two-end classifier. These statements do not promote the wider finite-evidence claims in **THM-3004, `01-canon/theorems/THM-3004-circuit-sign-change-cluster-law-and-classifier-refutation.md`**.

## 3. Factorial consumer is correctly typed

I read the relevant supplier in **THM-3079, `01-canon/theorems/THM-3079-laguerre-pf-row-transform-and-strict-integer-mesh-terminal-minus-one.md`**, especially its falling-factorial/Laguerre coefficient polynomial. The consumer uses only its strictly negative real-rooted factors. It does not import the reciprocal-Gamma Hankel theorem outside its positive strip.

For a monic negative-root factor `F(n)=product(n+r_i)`, the coefficient Cauchy bound `r_i<B=1+max|nonleading coefficients|` is valid and strict. If `L>=B/epsilon_0` and `c=rho/L`, then

`c^m F(n/c+L)`

is monic and has root parameters `rho(1+r_i/L)` in the required narrow block. Translation and independent dilation genuinely change the original polynomial. The final product preserves negative real roots, but its Newton ratios are those of the new product. The producer explicitly retains this distinction.

The independent companion reconstructs the `m=2,3,4` factors by finite differences of `(x+1)_m`, not by the producer's factorial formula. Disjoint rational sign intervals certify every negative simple root, avoiding its Sturm path. Derivative evaluation at each translation point reconstructs every coefficient of the moved factors; a rational convolution recovers all ten coefficients of the degree-nine target. The circuit word is exactly `+,-,+,+,-,-,-`, with three changes. Thus the positive consumer is literal and nonvacuous, while any claim about untouched factorial/Laurent rows remains excluded.

## 4. Global ballot and bronze proofs

For integers `a,b`, the condition `k>=max(1,1-b,1-a+b)` makes all three adjacent binomial entries positive and the uncancelled denominator strictly positive. Factorial cancellation gives the producer's ratio. The independent source verifies the full identity in the polynomial ring `Z[k,a,b]` with a separate sparse arithmetic implementation:

`1-R_k=2Q_(a,b)(k)/D_k`,

`Q=k^2+(a+1-(a-2b)^2)k+[a(a+2)-(2a+1)(a-2b)^2]/4`.

This proves the all-shift quadratic identity before cancellation. If the reduced numerator still has degree two, monic normalization returns the same `Q`; a canceled linear factor would instead lower its degree.

For constant coefficient `-1`, put `delta=a-2b` and `D=2a+1`, a nonzero odd integer. The exact identity is

`D(4delta^2-D-2)=13`.

Its complete signed divisor list is `D=1,-1,13,-13`. The negative choices give `delta^2=-3`. The positive choices give `delta=+2` or `-2`, and hence precisely `(a,b)=(0,+1),(0,-1),(6,2),(6,4)`. Their monic quadratics are `k^2-3k-1` and `k^2+3k-1`, respectively. Only the first pair has positive metallic parameter. Their discriminant 13 is nonsquare, so there is no cancellation against the rational linear denominator factors. This proves the **global** classification, including the reduced-numerator formulation; it is not a finite rectangle result.

The integer-order Fuss formulas are valid for every integer order `p>=2`, by elementary factorial ratios. The old claim that rationality singled out `p=2` is false. The independent source also checks the displayed `p=4` quartic by an exact polynomial identity and evaluates its numerator at all denominator roots to exclude hidden cancellation. Rationality at higher order does not preserve a quadratic numerator.

Finally the exact `(1,1,3,3,9,9)` normalized row has a palindromic ratio vector and an antipalindromic, maximally alternating circuit. This refutes a general norm-plus-one/nonalternation assertion. It does not identify positive coefficient-polynomial root parameters with the signed characteristic roots of a metallic recurrence. The surviving two-element disjointness statement is correctly separated in the producer.

## 5. Independent reproduction and promotion basis

The adjacent `continuing8_20260906_newton_ballot_audit.py` pins both complete producer JSON certificates. It uses Newton power sums instead of root convolution, a max-plus coefficient recurrence and a separate factorial formula instead of the producer's top-coefficient construction, finite differences plus derivative transport instead of symbolic Laguerre substitution, and sparse polynomial arithmetic instead of SymPy for the universal ballot identities.

It reconstructs every one of the 30 declared cluster controls, every threshold and coefficient envelope, all first/middle/final signs, the exact reciprocal tie, all named scope hostiles, and the entire factorial transport. Its different ballot bank has 1,516 exact rows and 99 Fuss ratios. These computations corroborate the proofs; the proofs supply their universal quantifiers.

```text
python continuing8_20260906_newton_ballot_audit.py
python -O continuing8_20260906_newton_ballot_audit.py
```

Both modes pass **2,703 always-active exact gates**. Their actual output bytes are identical and contain no carriage returns. The audit source and output hashes are:

```text
source 6ef9b86c577a9fba03484880092a5b1f92065a903a72fbf25db7171c8e683088
output c1340a86b563eb0021219e3e166866a19811fc3889819b9c2836e54eb1029811
```

The promotion basis is the complete analytic proof review plus the independent exact reconstruction. No broad census, repository mutation, theorem-ID reservation, or proof of arbitrary-separation no-return was performed by this referee.
