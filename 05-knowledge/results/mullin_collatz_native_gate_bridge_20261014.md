# Mullin–Collatz native gates: a word compiler, ternary memory, and two rooted hypotheses

Date: 2026-10-14.

Status: **PROVED elementary statements** below, with **FINITE-EXACT** checks; the two final hypotheses are **OPEN / SPECULATIVE**. The signed-kernel rank application uses the independently checked companion [joint-kernel theorem](signed_collatz_joint_tower_kernel_20261014.md). No claim of literature novelty is made for the elementary identities or CRT techniques. No claim resolves Mullin's question or either fixed-sign Collatz problem.

## 1. The useful connection and its limits

There is an exact common operation: transport a native arithmetic cylinder through an affine word, then measure how much ternary source information survives in the next signed Collatz output pair. The word gates are different and remain explicit.

- Mullin: multiply the state by a prime, but only when that prime is the **least** divisor of the state plus one.
- Signed accelerated Collatz: apply `(3n+epsilon)/2^a`, but only when `a=v_2(3n+epsilon)` and `epsilon` is the declared sign.
- Both finite word problems compile to local residue conditions. The Mullin gate has a nontrivial backward-priority obstruction; the Collatz valuation cylinder is dyadic.
- These cylinders can be fused exactly by CRT. They are independent in the dyadic coordinate, but the Mullin cylinder leaves a measurable ternary signature in the signed output-pair kernel.

The output is a concrete construction and a pruning criterion, not a claim that two famous open problems become equivalent. In particular, arbitrary realizable seed cylinders do not imply that the canonical seed `P=1` visits them. Sampling many starting integers does not describe the time distribution along one orbit.

### Literature boundary, checked 2026-10-14

The sequence in the question is the first Euclid–Mullin sequence:

    P_0=1, p_(m+1)=lpf(P_m+1), P_(m+1)=P_m p_(m+1).

[Booker–Irvine, The Euclid–Mullin graph (2016)](https://arxiv.org/abs/1508.03039), §1 and Lemma 2.2, define the all-prime-factor graph and prove its finite-word CRT progression statement. Their graph allows every prime divisor; it is not the least-factor selector. Their graph is not a tree.

[Booker–Simon, A generalisation of the Euclid–Mullin sequences (2026)](https://arxiv.org/html/2601.21901v1), introduction, explicitly state that coverage and recursive membership for the first sequence remain open. They also recall Shanks' missing-prime residue-avoidance reasoning. The largest-factor sequence omits infinitely many primes. Their cyclotomic generalizations are prior work, not results claimed here.

[Booker, A variant of the Euclid–Mullin sequence containing every prime (2016)](https://arxiv.org/abs/1605.08929) proves coverage for a changed generator; it does not settle the original selector. The [OEIS A000945 data file](https://oeis.org/A000945/b000945.txt) supplies the 51-term prefix used in the bounded probes below. Those large factorizations are not repeated here.

Repository inheritance: the affine carry interface is [THM-4604, Collatz words are the reducible locus](../../01-canon/theorems/THM-4604-collatz-words-are-the-reducible-locus-carries-are-the-extension-cocycle-and-traces-only-hear-the-tuning.md). The closest new mechanism is the [signed joint ternary kernel](signed_collatz_joint_tower_kernel_20261014.md). The hostile boundary is that a controlled-sign descending path does not establish fixed-sign convergence. Here the additional hostile is the impossible least-factor word `(7,3)`, despite its perfectly valid all-factor CRT cylinder.

## 2. Exact deterministic Mullin word language

Work with an even seed `P_0=2n`, where `n` is positive and odd. This includes the canonical post-first-step seed when `n=1`, but the following theorem varies the seed.

Let `w=(p_1,...,p_k)` be distinct odd primes. Put

    A_0=1, A_i=p_1...p_i, M=max_i p_i,
    R=product of odd primes r<=M.

We ask whether the **native least-factor** map can follow `w`, that is,

    p_i = lpf(2 A_(i-1) n + 1),  1<=i<=k.                  (2.1)

### Theorem 1: priority obstruction, complete iff criterion

A positive odd seed realizing (2.1) exists if and only if, for every `i<j` with `p_i>p_j`,

    p_i p_(i+1)...p_(j-1) != 1 (mod p_j).                  (2.2)

When this holds, infinitely many positive odd seeds realize the word. In fact the exact relative natural density among odd `n` is an explicit positive rational number.

The obstruction means exactly this: a later smaller prime was already enabled when the word proposed choosing a larger prime. That word cannot be repaired by changing the seed while preserving all its divisibility gates.

**Proof.** For a selected prime `r=p_j`, its selection gate requires

    n = c_r := -(2 A_(j-1))^(-1) (mod r).                 (2.3)

For an earlier step `i<j` with `p_i>r`, least-factor priority additionally requires `r` not to divide `2 A_(i-1)n+1`. Substitution of (2.3) makes the prohibited equality

    A_(i-1) A_(j-1)^(-1) = 1 (mod r),

which is equivalent to the negation of (2.2). Earlier steps with `p_i<r` impose no priority exclusion for `r`. At every step after selection of `r`, the running product is divisible by `r`, so the next product plus one is automatically `1 mod r`.

For a prime `r<=M` not in the word, define

    F_r={-(2 A_(i-1))^(-1) mod r : p_i>r},
    S_r=F_r complement in Z/rZ.                          (2.4)

These are exactly its priority exclusions. Every forbidden residue is nonzero; consequently `0` always belongs to `S_r`. There is no further local emptiness obstruction. At a selected prime put `S_r={c_r}`, provided (2.2) holds. The factor `2` never divides `2A_(i-1)n+1`.

The complete solution set is therefore

    n mod r belongs to S_r for each odd prime r<=M.       (2.5)

The moduli are distinct primes, and parity is independent. CRT proves existence, infinitude, and the exact density

    delta(w)= product_(r<=M, r odd) |S_r|/r.             (2.6)

If any priority obstruction occurs, set `delta(w)=0`. These are sufficient gates because they exhibit the divisor `p_i` and exclude every smaller prime divisor at every step. ∎

Thus **realizability**, with primality of the supplied labels assumed or certified, needs only the quadratic list of modular suffix-product tests (2.2). Computing the full density or an explicit small witness may require a much larger prime wheel. This is a genuine distinction from factoring a new canonical Euclid–Mullin term.

Examples:

- `(7,3)` is impossible: `7=1 mod3`, so every seed enabling the later `3` enabled it before the proposed `7`.
- `(7,5,3)` is possible. Its density among odd seeds is `1/105`.
- Every strictly increasing word of distinct odd primes is possible; there are no reversed-priority pairs to check.

### Squarefree seeds are available too

The theorem is not an artifact of allowing repeated factors in the initial seed. Each nonempty residue component (2.5), combined with any odd dyadic cylinder, contains a positive density of squarefree integers `n`; then `P_0=2n` is squarefree.

To see this, its modulus has squarefree odd part and power-of-two part, while its residue is odd. For an odd prime dividing that modulus, a nonzero residue already excludes its square; a zero residue modulo that prime leaves a proportion `1-1/r` with valuation exactly one. At primes outside the modulus the squarefree sieve contributes `1-1/r^2`. The product is positive. Finite inclusion-exclusion followed by the convergent square-divisor tail proves the density assertion. Thus one may retain the usual interpretation of a seed as a product of distinct already-used primes.

## 3. Exact fusion with any signed Collatz valuation word

Fix signs `epsilon_1,...,epsilon_l` in `{+1,-1}` and positive valuations `a_1,...,a_l`. Let `S_i=a_1+...+a_i` and form the carries

    B_0=0,
    B_i=3 B_(i-1)+epsilon_i 2^S_(i-1).

Native composition gives

    2^S_i n_i = 3^i n+B_i.

### Theorem 2: native cylinder product

The exact signed valuation word occurs on one odd residue class

    n = (2^S_l-B_l) 3^(-l) (mod 2^(S_l+1)).             (3.1)

If the Mullin word `w` satisfies Theorem 1, then the two native words occur simultaneously for a set of odd starting `n` of relative density

    delta(w) 2^(-S_l) > 0.                              (3.2)

Both words also coexist at infinitely many squarefree odd sources.

**Proof.** Equation (3.1) makes the final iterate an odd 2-adic integer. Recovering earlier iterates backwards by `(2^a n_i-epsilon)/3` shows they are odd 2-adic integers and that each declared valuation is exact. Equivalently one may impose the valuations successively. Positive sufficiently large representatives make all integer intermediate states positive; in fact either signed map preserves positive odd integers. The source class has relative odd density `2^(-S_l)`. All Mullin constraints use the odd squarefree modulus `R`, coprime to (3.1). CRT gives (3.2), and the squarefree argument above applies. ∎

A concrete shared witness is `n=2201`:

    Mullin seed P=4402:
      P+1=4403=7*629, first prime 7;
      7P+1=30815, first prime 5;
      35P+1=154071, first prime 3.

    Signed Collatz:
      2201 --(+1, valuation 2)--> 1651
           --(-1, valuation 3)--> 619.

The joint density for these two words is `1/(105*32)=1/3360`.

This proves a useful limitation on analogies: no particular finite signed Collatz valuation word can be excluded merely because a compatible finite Mullin word occurred at seed `2n`. The native words are exactly independent in the dyadic direction. A shared lower-level multiply/add/divide presentation alone cannot supply a chronological conclusion.

## 4. The nontrivial bridge: ternary memory and exact rank loss

For odd `n` write

    u=U_+(n), v=U_-(n), a=v_2(3n+1), b=v_2(3n-1).

The inherited signed identity is

    min(a,b)=1,   2^a u-2^b v=2.

The full joint distribution of `(u,v) mod3^j` under uniform odd sources is a matrix of dimension `d=2*3^(j-1)` and rank `d`. This is the companion joint-kernel theorem, not a temporal Markov assertion.

### Theorem 3: the Mullin word's entire ternary influence at its source

Let `w` be admissible and let `S=S_3` from Theorem 1 (if the wheel has no `3`, set `S=Z/3Z`). Condition positive odd sources on realizing `w` from seed `2n`.

1. The dyadic coordinate remains Haar-uniform.
2. The ternary coordinate is Haar-uniform on `n mod3 in S`.
3. All other Mullin prime gates are independent of these two coordinates for every finite observation.
4. The conditional signed output-pair matrix has rank

       2,                         j=1;
       2 |S| 3^(j-2),             j>=2.                 (4.1)

Thus a native least-factor word can be invisible modulo `3` yet remove exactly two-thirds of the matrix rank modulo `9` and at every higher level. If the word contains `3`, `S` is a singleton and this rank loss is automatic.

**Proof of the conditional law.** The gate at `3` constrains only one ternary digit; all other gates are modulo primes different from `2,3`. CRT establishes the claim for every finite dyadic and ternary observation. The shrinking dyadic valuation tails justify passing to the accelerated maps.

**Proof of rank.** More generally fix `n=s mod3^h`, and let `j>=h+1`. On the branch `a>=2,b=1`,

    v=(3n-1)/2 = (3s-1)/2 (mod3^(h+1)).

On `a=1,b>=2`,

    u=(3n+1)/2 = (3s+1)/2 (mod3^(h+1)).

So the matrix is supported on `3^(j-h-1)` columns together with that many rows. Its rank is at most `2*3^(j-h-1)`. There are `3^h` possible source classes `s`, and the unconditioned full-rank matrix is their equally weighted sum. Rank subadditivity forces every singleton bound to be exact. Applying the same argument to a subset and its complement gives exact rank `2 |S| 3^(j-h-1)` for a uniform mixture over a subset of classes. For (4.1), take `h=1`.

At `j=1`, all three source classes have the original rank-two matrix, since the source's ternary residue does not affect the first output digit. ∎

**Overlap warning.** In the middle allowed rectangle a cell can receive contributions from both high-valuation branches. Conditioning must be applied to each branch's source separately. Multiplying the unconditioned matrix by the union-of-strips indicator is generally incorrect.

### Transport to the Mullin endpoint

After `w`, the odd state is `n'=A_k n`. Put `h=v_3(A_k)`, which is `0` or `1` because selected primes never repeat.

- If `3` is absent from `w`, the endpoint is uniform on `n' mod3 in A_k S`; its rank law is again (4.1).
- If `3` occurs, the source class `S={s}` is nonzero, and `n'` occupies one nonzero-multiple-of-three class modulo `9`, namely `A_k s mod9`. The endpoint's paired kernel has rank

       2*3^max(j-3,0).                                  (4.2)

The one-time selection of `3` therefore leaves one additional ternary digit of endpoint memory.

For comparison, a length-`l` native Collatz valuation prefix has affine numerator coefficient `3^l`. Its endpoint occupies one class modulo `3^l`, with Haar-uniform dyadic coordinate and higher ternary digits. The companion native-prefix theorem consequently gives rank

    2*3^max(j-l-1,0).                                   (4.3)

This is a shared cylinder-transport mechanism with a real distinction: the Mullin multiplier can acquire the prime `3` once; the Collatz affine multiplier acquires another factor `3` on every step. The divisibility depth, native carry, and exact gate are indispensable coordinates.

### Theorem 4: unconditional rank ceiling along the actual canonical orbit

Let `n_m=P_m/2` for the canonical sequence and `m>=2`. These odd integers have `v_3(n_m)=1`, because `3` is selected once and never again. Form any empirical matrix by placing arbitrary real weights at the pairs `(U_+(n_m),U_-(n_m)) mod3^j` for any finite selection of indices. No averaging or mixing hypothesis is imposed. Then

    rank <= 2,                         j=1,2;
    rank <= 4*3^(j-3),                 j>=3.             (4.4)

In particular the actual canonical product sequence is visible through a rank-at-most-four paired-output matrix modulo `27`, versus rank eighteen for unrestricted odd sources.

**Proof.** For `j>=3`, the possible source residues modulo `3^(j-1)` with valuation one number `2*3^(j-3)`. If the plus valuation is one, the row belongs to the image of these source residues under `(3n+1)/2`. If the minus valuation is one, the column belongs to their image under `(3n-1)/2`. At least one of those valuations is always one. Thus every entry is supported on `2*3^(j-3)` rows or the same number of columns, regardless of weights. The rank of a matrix with this support is at most the sum of those two counts. The two low levels give at most one such row and one column. ∎

Modulo `27`, the two rows are `u=5,23` and the two columns are `v=4,22`. The imported 51-term prefix, using indices `2..51`, already gives exact empirical ranks

    modulus:    3  9 27 81
    rank:       2  2  4 12.

These saturate (4.4) in this finite data. They do not establish a limiting distribution. The shell-mixing hypothesis below concerns the exact cell weights and asymptotic behavior; the rank ceiling needs no such hypothesis.

## 5. Frontier scheduler and an obstruction to seed-uniform coverage bounds

### Finite-veto lemma

Fix an unselected prime `q` at an arbitrary current product `P`. Every occurrence of `q | P+1` at which the least-factor selector chooses a different prime consumes a previously unselected prime smaller than `q`. Consequently the number of such vetoes is at most

    #{primes r<q : r does not divide the initial P}.      (5.1)

No blocker can be consumed twice, since a selected prime divides every subsequent product and cannot divide any later product plus one.

When `q` is the least missing prime, the budget is zero:

    next selected prime=q  iff  P=-1 modq.              (5.2)

Thus prime coverage is a target-hitting problem for the changing least-missing prime. This reduction is elementary and is already implicit in the Shanks reasoning recalled by Booker–Simon. It is not presented as a new solution.

### Theorem 5: no waiting bound depending only on q over all seeds

For every `L`, there exist squarefree even seeds for which `3` is the least missing prime and the next `L` least-factor choices all avoid `3`. They can simultaneously realize any prescribed finite signed Collatz valuation word at the odd source `n=P/2`.

**Proof.** Choose `L` increasing primes all congruent to `1 mod3`. Such arbitrarily long lists exist: infinitely many primes are `1 mod3`, for example by the standard elementary argument using a prime divisor other than `3` of `x^2+x+1`, with `x` divisible by `3` and any proposed finite list. The word is admissible by Theorem 1. All its prefix products are `1 mod3`, so its only exclusion at prime `3` is `n!=1 mod3`; choose `n=2 mod3` instead. Then `P=2n=1 mod3` throughout these `L` steps. The seed is even, not divisible by `3`, so `3` is its least missing prime. Theorems 1 and 2 and the squarefree sieve supply the claimed seeds and optional Collatz word. ∎

This is why a quantitative coverage hypothesis below must be rooted at the canonical initial product. Arbitrary-seed finite-cylinder abundance cannot justify a seed-uniform clock. There is no contradiction with possible eventual coverage from each individual seed.

## 6. Two sharp, deliberately ambitious hypotheses

These are research targets, not conclusions suggested by 51 terms. They ask for different missing information.

### H-M1: canonical quadratic frontier deadline

For the canonical sequence only, suppose prime `q` becomes the least missing prime at index `m`. Then `q` is selected in one of the next `q^2` steps:

    p_(m+t)=q for some 1<=t<=q^2.                        (H-M1)

This is much stronger than eventual coverage and has a single finite counterexample criterion. Theorem 5 rules out extending it to all seeds. A random-residue heuristic has typical waiting scale `q-1`, so `q^2` is a conservative scale; that is motivation, not arithmetic evidence.

If true, H-M1 gives effective prime coverage: every prime at most `y` has appeared by index at most `sum_(q<=y, q prime) q^2`, in particular at most `y^3`. This follows by charging only the successive frontier primes; primes selected ahead of their turn cost no additional frontier interval.

Cheapest test: record the activation and consumption of each actual frontier, without conflating a prime's first absence with its activation. The sourced 51-term prefix gives completed waits

    q:       2  3  5  11 17 19 31
    wait:    1  1  5   5  1 23 14.

All satisfy the proposed bound. The active frontier is `41` from index `50`; the prefix ends at `51`. The hypothesis predicts selection of `41` by index `1731`. The present data neither verifies this prediction nor provides the factorizations needed to approach that index. There is no promise to extend the sequence computationally here.

### H-M2: canonical dyadic–ternary squarefree-shell mixing

Put `n_m=P_m/2` for `m>=2`. These are odd and have `v_3(n_m)=1`, exactly, because `2` and `3` each occur once. For every fixed `K>=1,J>=2`, conjecture that the residues

    n_m mod (2^K 3^J)

are equidistributed, as `m` tends to infinity, among the `2^K 3^(J-2)` residues that are odd and have 3-adic valuation exactly one.

This retains the true arithmetic shell instead of guessing uniformity on all residues. It is a canonical-orbit assertion; the finite-seed CRT theorems above do not prove it.

Its measurable Collatz consequence is precise. Sampling `(U_+(n_m),U_-(n_m))` along canonical Mullin products would have the one-step joint law obtained from Haar odd sources conditioned on `n=3 or 6 mod9`. Equidistribution at each finite dyadic depth, with the exact geometric tail bound, justifies this implication. The resulting limiting matrix ranks would be

    rank=2 at j=1,2;
    rank=4*3^(j-3) at j>=3.                             (6.1)

Thus the predicted limiting rank at modulus `27` is `4`, saturating the unconditional ceiling in Theorem 4, whereas the unrestricted odd-source kernel has rank `18`. Both the forbidden cells and all rational cell weights can be compiled exactly before further Euclid–Mullin terms are computed. This would establish a nontrivial deterministic sampling law for Collatz observables, while saying nothing by itself about convergence of a Collatz trajectory.

Cheap descriptive test: at modulus `72`, the 50 sourced products indexed `2..51` give counts on the allowable residues

    residue:  3 15 21 33 39 51 57 69
    count:    7  5  9  5  5  6  9  4.

The uniform comparison is `6.25` per residue. This is a tiny descriptive sample, not a statistical confirmation. Systematic long-run concentration on a proper subset or a nonzero limiting character average would refute H-M2. Checking exact shell membership and finite-level joint weights is cheap; obtaining many additional canonical products is not.

H-M1 asks for chronological target recurrence at a moving prime. H-M2 asks for equidistribution after retaining fixed ramified coordinates. Neither follows from the other as stated. Together they suggest a more focused long-term program: combine the root-specific least-factor scheduler with the affine word's surviving residue coordinates, rather than replacing either gate by an iid model.

## 7. Reproducible bounded evidence and loss ledger

Files:

- [Exact test script](../../04-computation/experiments/mullin_collatz_native_gate_bridge_20261014.py)
- [JSON output](mullin_collatz_native_gate_bridge_20261014.json)

Commands from the repository root:

    python3 04-computation/experiments/mullin_collatz_native_gate_bridge_20261014.py
    python3 -O 04-computation/experiments/mullin_collatz_native_gate_bridge_20261014.py

Both modes produce byte-identical output; checks use explicit exceptions, not removable assertions.

The run contains **79,679 checks**. It exhausts all **64** distinct-prime words of lengths `1..4` over `{3,5,7,11}`; exactly **37** are admissible. Every word is checked over its entire odd-source prime wheel against direct integer least-factor gates. Every admissible word is fused with all **36** signed length-two Collatz words having valuations in `{1,2,3}`. Conditional kernel ranks are checked by exact rational elimination modulo `3,9,27`, including the singleton and squarefree-shell conditions. Eight endpoint-law checks independently impose actual Mullin word gates, multiply to the endpoint, compute native signed outputs, and compare every normalized matrix entry with the predicted ternary-cylinder law. The canonical prefix empirical matrices are also checked through modulus `81`, where their ranks saturate the unconditional shell ceilings. The imported 51-term prefix is checked for divisibility and distinctness, but its huge prime certifications and least-factor factorizations are explicitly **not** rerun.

Preservation and loss:

- The suffix test preserves native least-factor word existence over varying seeds; it loses the canonical root and seed height.
- The CRT product preserves both finite native words and exact density; it loses chronological sampling along the canonical Mullin orbit.
- The conditional kernel preserves the full same-source output pair at the specified finite ternary level; it loses the source's unobserved digits and both maps' future paths.
- Rank records the number of surviving residue cylinders in this particular kernel; rank alone is not target hitting, independence, mixing, or descent.
- The affine transport restores the exact ordered carry and ramification depth; passing to traces or just multiplier sizes would discard this information.

The reusable deliverable is the **native word compiler plus ramified cylinder transport**. The remaining global problem is to explain why the canonical root must supply the right cylinders in time.
