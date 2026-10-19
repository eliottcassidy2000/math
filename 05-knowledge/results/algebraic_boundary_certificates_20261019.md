# Algebraic boundary certificates and exact threshold rigidity

**2026-10-19. PROVED scoped transfer statements; no universal 7-thin conclusion or priority claim.**

The Fibonacci irrational-ray obstruction does not require abandoning exact cones. Algebraic rays are natural for periodic contexts. This note shows why quadratic fields are not universally sufficient, proves a finite exact certificate for one-sided periodic contexts with a strict witness margin, and isolates the equality case where that certificate can fail.

## 1. A natural 4-thin transfer already requires a quartic Perron field

Let R_n be the poset on labels 1,…,n with

`i<j` exactly when `j−i≥3`.

Its incomparable pairs have label distance one or two, so its maximum incomparability degree is four for n≥5. At every sufficiently interior ideal rank k there are precisely four ideals:

- `S0={1,…,k}`;
- `S1={1,…,k−2} ∪ {k,k+1}`;
- `S2={1,…,k−1} ∪ {k+1}`;
- `S3={1,…,k−1} ∪ {k+2}`.

To see completeness, if an ideal contains label m, it contains every label through m−3. Among m−2,m−1 it may contain neither, one, or both. Fixing its cardinality gives exactly the four listed possibilities.

In that order the interior one-step matrix is

`T=[[1,0,1,1],[1,0,0,0],[1,1,0,0],[0,1,1,0]]`.

Each entry follows by testing whether the source ideal is contained in the target ideal of the next rank. Its characteristic polynomial is

`q(x)=x⁴−x³−x²−3x−1`.

Modulo two this is `x⁴+x³+x²+x+1`. It has no root in F₂ and is not divisible by the only irreducible quadratic `x²+x+1`. Thus q is irreducible over Q. The matrix T⁴ is strictly positive (also a consequence of the D-step positivity theorem at D=4).

Let λ be its Perron eigenvalue. Then λ lies between 2 and 3, is an algebraic unit of degree four, and has positive eigenvector

`v=(λ³,λ²,λ²+λ,2λ+1)^T`.

The first/second coordinate ratio is λ, so the exact Perron ray cannot have coordinates in any quadratic number field. This does not exclude rational or quadratic-coefficient *enclosing* cones; it only rules out the proposed universal restriction to quadratic fields for exact extremal-ray coordinates.

For periodic contexts with S states, the Perron eigenvalue of an integer period matrix is always algebraic of degree at most S. The useful algebraic field is therefore attached to the actual period matrix, rather than imposed as the golden-ratio field in advance.

## 2. Fixed rational marked rows and one-sided periodic contexts

Let M be a positive rational S×S matrix. Let b be a nonzero nonnegative rational terminal column. (The error bounds also hold for arbitrary real nonnegative b; rational data are needed for the promised exact checks of finite exceptions.) Let f,c be fixed rational rows satisfying

`f≥0`, `f≠0`, `0≤c≤f` entrywise.

Here f counts all prefixes ending at each frontier state and c counts prefixes with a specified event. Appending t repetitions of the same continuation gives

`p_t=c M^t b/(f M^t b)`.

Let v>0 be a Perron right eigenvector. Then

`p_∞=c v/(f v)`.

The rows c,f are fixed rational data. They are not stationary algebraic left vectors, and this is not yet a two-sided stationary-ratio theorem.

### Theorem 1. Strict algebraic margin gives a finite exact tail certificate

If `β<p_∞<1−β`, then p_t is β-balanced for all sufficiently large t, uniformly over every nonzero nonnegative boundary b. An explicit bound is obtained from

`|p_t−p_∞| ≤ tanh(Δ τ^(t−1)/4)`, for t≥1,

where Δ is M's Hilbert projective diameter and `τ=tanh(Δ/4)<1`.

**Proof.** Since `Mv=λv`, the projective distance between `M^t b` and v is at most `Δτ^(t−1)`: the first positive multiplication supplies the image-diameter bound and subsequent multiplications contract. Multiplying coordinates by f and normalizing does not increase the log-likelihood oscillation. The ratio p_t is the expectation of the [0,1]-valued observable c_i/f_i under those weights, with zero-f coordinates omitted. The likelihood-ratio bound gives the stated error. ∎

For an exact executable threshold, choose certified rational bounds `B≥Δ`, `τ≤ρ<1`, and

`0<ε≤min(p_∞−β,1−β−p_∞)`.

Every t satisfying `Bρ^(t−1)≤4ε` works, because tanh(z)≤z. The comparison can be done with integers. If τ=0, one positive step already aligns every boundary with the Perron ray.

For a D-step ideal kernel the preceding note supplies `Δ≤2 log C_D` and `τ≤(C_D−1)/(C_D+1)`, so no numerical eigenvector conditioning is needed for the error estimate. To prove the fixed witness for *all* lengths, check the finitely many smaller t exactly. A primitive matrix can be treated by taking a positive power and checking its finitely many residue classes.

## 3. Full-degree Perron eigenvalues eliminate nontrivial rational endpoint limits

### Lemma 2. Coordinate independence

Let T be an S×S rational matrix with eigenvalue λ of algebraic degree S, and let v≠0 be a corresponding eigenvector. Then the S coordinates of v are linearly independent over Q.

**Proof.** Let V be their Q-linear span. The eigenvector equation makes V invariant under multiplication by λ. For any nonzero a∈V, the S numbers `a,aλ,…,aλ^(S−1)` belong to V and are Q-linearly independent. Thus dim_Q V≥S. There are only S coordinates, so their span has dimension S and they are independent. ∎

### Corollary 3. Rational threshold rigidity

Under Lemma 2, for fixed rational rows c,f and rational β,

`c v/(f v)=β` implies `c=βf` entrywise, provided f v≠0.

Consequently `c h/(f h)=β` for every boundary vector h with nonzero denominator, not merely for a Perron limit.

This gives a complete finite decision procedure for the following specific question: is the fixed event β-balanced for all repetition lengths with positive denominator in a fixed one-sided periodic context whose primitive period matrix has full-degree Perron eigenvalue?

1. Compare the algebraic p_∞ with β and 1−β exactly.
2. If the limit is outside the closed interval, sufficiently long repetitions fail.
3. If it equals an endpoint, Corollary 3 makes every length exactly that endpoint.
4. If it is strictly inside, Theorem 1 gives a finite tail cutoff; inspect all smaller lengths with positive denominator by exact rational counting.

Both the Fibonacci transfer and the quartic transfer of Section 1 satisfy the full-degree hypothesis. The procedure addresses a fixed marked event and periodic context, not selection of a witness in an arbitrary nonperiodic 7-thin poset.

## 4. Why a two-sided stationary ratio needs a separate argument

An expression `l^T E v/(l^T N v)` contains products of algebraic coordinates. Those products can be rationally dependent even when the coordinates of l and v separately are independent. Corollary 3 does not apply to it.

A concrete integer marked-kernel control uses

`T=[[2,1],[1,1]]`, `K=[[0,1],[-1,0]]`,
`N=3T`, `E=T+K=[[2,2],[0,1]]`.

Here `0≤E≤N` entrywise. With the common Perron vector v of the symmetric T, skew symmetry gives `v^T K v=0`, so

`v^T E v/(v^T N v)=1/3`,

although E is not N/3 entrywise.

For finite boundaries e=(1,0)^T and nonnegative integers m,n, the discrepancy from 1/3 has numerator

`e^T T^m K T^n e`.

For m=0,n=1 it is +1, and for m=1,n=0 it is −1. More generally, because det(T)=1, the identity `T^r K T^r=K` shows that the sign depends only on n−m, and is positive when n>m and negative when m>n. Thus there are failures of the lower 1/3 threshold with both boundary distances arbitrarily large, despite exact stationary equality and strict positive mixing.

This is an abstract nonnegative marked-kernel example. A realization of E as a comparison event in a poset is not asserted. It demonstrates why a general transfer argument cannot infer finite-context equality from a two-sided stationary equality without using more of the comparison-event structure.

## 5. Complex subdominant modes can alternate across a rational threshold forever

Let J be the 3×3 all-ones matrix and let P be a cyclic permutation matrix. Put `M=J+P`, which is strictly positive. Its eigenvalues are 4,ω,ω², where ω³=1 and ω≠1. Its Perron vector is `(1,1,1)^T`.

Take `f=(3,3,3)`, `c=(2,1,0)`, and `b=e_1`. Then

`p_t=c M^t b/(f M^t b)`

equals 1/3 plus a repeating permutation of `1,0,−1`, divided by `3·4^t`. Indeed

`M^t b=4^t(1,1,1)^T/3 + P^t(b−(1,1,1)^T/3)`.

Thus the limiting ratio is exactly 1/3, but the finite ratios lie below and above it infinitely often. All entries, marked rows, and boundaries are nonnegative, with `0≤c≤f`. Again, this is a marked-kernel obstruction, not a claimed poset-comparison counterexample.

The algebraic issue is now exact: when the Perron degree is smaller than the state dimension, a nonzero rational deficit row may annihilate the Perron vector. Its sign is then controlled by the surviving subdominant modes. Exponential contraction controls magnitude but cannot determine that sign.

## 6. What this contributes to the 7-thin route

- Exact algebraic boundary rays repair the irrational-cone problem for a fixed periodic context.
- The field need not be quadratic, already at degree four incomparability.
- Strict margins, and full-degree rational-threshold rigidity, turn a one-sided periodic witness into a finite exact certificate.
- Two-sided equality and degree-deficient Perron limits require additional marked-event structure or an exact subdominant-mode analysis.
- Arbitrary legal nonperiodic products still need a global invariant-cone/witness-migration certificate. No reduction of all 7-thin posets to finitely many periodic types is proved.

The standard projective contraction theorem is cited and proved in `ideal_transfer_positivity_20261019.md`. All new statements here have explicit proofs. A comprehensive literature-priority check has not been completed.

Reproduce the small controls with `python3 05-knowledge/experiments/verify_algebraic_boundaries_20261019.py`; its output is `algebraic_boundaries_20261019.json`. It independently constructs the quartic ideal matrix at an interior rank, checks its characteristic polynomial and mod-two irreducibility test, verifies positivity of T⁴ and its exact maximum cross ratio 2, and checks both abstract threshold examples with integer arithmetic. The all-length claims rest on the displayed proofs, not the finite replay.
