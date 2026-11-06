# Independent audit: bounded source-error obstruction

Outcome: **PASS, in the stated all-order scope.** No mathematical correction is required. This audits the frozen candidate, not an unrestricted target-occupancy claim. The twelve-template occupancy problem remains OPEN. No source files, repository files, branches, or GitHub state were changed.

## Inputs and isolation

`INPUTS.json` records original paths, frozen copies, and SHA-256 hashes. The candidate proof hash is `6b88ac774aa9b4bc320eb4d8ccf93c857b87b5931d010e5446bf8e6c2f00b4f9`. The author's manifest verifies. The referenced compact release `dependencies/CORE_PROOF.md` is byte-identical to `session201-central-arithmetic-release/session201/central-arithmetic/PROOF.md`, hash `49ad1123e80fbf751eac3e72ed0f5ad2ec39e7705f9ee89c77d4826ef5e1066c`.

The related positive approximation proof is frozen only as a dependency/context input. This report independently verifies its index bijection and run identity, the only portions needed here. It does not audit its discrepancy statements. The compact twelve-template PROOF is also frozen to identify the menu explicitly.

The shared math checkout was dirty at entry. Its status is recorded before and after; no startup/checkpoint/fetch/sync or other mutating repository command was run. Scratch-only auditing takes precedence over normal research-session publishing. AGENTS.md, navigation, relevant guardrails, and research protocol were consulted.

## 1. Actual source and index map

At the seed `(u,g,p,t)=(2,1,0,1/2)`, the error is zero; `p` here denotes parity to distinguish it from the fixed error threshold E. Direct substitution gives invariance of `g-ut` under O and increment `(t+1)/2-p` under Z. Before a maximal Z-run, p=0, either because it starts at the seed or follows O. For subsequent letters in that run p=1. This validates the first-versus-later increment split; no hidden resetting of error occurs.

For chronological bits epsilon_i, induction gives `u=2^(k+1)-sum epsilon_i 2^(k-1-i)` and `t=(sum epsilon_i 2^i+1/2)/2^k`. Reversal permutes all k-bit labels. Thus each length has exactly the interval `2^k<u<=2^(k+1)`, and those disjoint intervals partition integers u>=2. Equivalently, removing the initial 1 from the binary expansion of u-1 and changing 0 to Z and 1 to O yields its unique chronological word. The empty word has k=0 and u=2. The excluded singleton u=1 is never used. For u exactly a power of two, k=log2(u)-1, not log2(u); the candidate's summation through floor(log2 N) is a harmless upper bound and does not miscount at equality.

## 2. Run contribution and quantitative lower bound

If a run starts at s and has r>=1 Z letters, summing the first increment `(s+1)/2` and later increments `-(1-s)2^(-j)`, j=2,...,r, gives `eta=s+(1-s)2^(-r)`. This remains true for a final run without terminal O, including a one-letter run. Every eta is strictly positive.

Let o be the number of immediately preceding O letters since the preceding Z-run (or since the seed). The preceding base t is 1/2 initially and strictly above 1/2 after any Z-run. Halving o times proves `s>=2^(-o-1)`. The seed and the O separating successive runs prove `s<=1/2`. Therefore `eta>=2^(-o-1)+2^(-r-1)>=2^(-(o+r)/2)`, the last step being two-term AM-GM. This includes o=0 for the initial run. Trailing O letters are omitted, never double counted, and the sum of assigned o+r is at most k.

For R>=1 the function f(x)=2^(-x/2) is convex and decreasing. Jensen gives `sum f(o_i+r_i)>=R f(sum(o_i+r_i)/R)>=R f(k/R)`. Combining with the exact run decomposition proves the candidate's lower bound for EVERY word, not merely sampled words. Raising to the positive even power 2R gives the rational-power check because both sides are nonnegative. The R=0 case is separate, with e=0; no expression dividing by R is used there.

## 3. Uniformity and subexponential counting

Fix finite E>0 once and for all. For each delta>0, if R>=delta k then R>=delta k and `k/(2R)<=1/(2delta)`, so `e>=delta k 2^(-1/(2delta))`. Thus for every integer k strictly greater than `(E/delta)2^(1/(2delta))`, no permitted word has R>=delta k. This proves the supremum over ALL permitted words has R/k tending to zero. It is not a pointwise limit that permits an exceptional moving sequence. If E=0, positivity of each run forces R=0. For R>E, taking logarithms gives exactly `2R log2(R/E)<=k`; the restriction R>E is useful but not needed for the algebra when E>0.

There are exactly binomial(k+1,2R) words with R Z-runs: choose starts and ends from the k+1 integer boundaries, in strictly increasing order, and pair them in order. Strict separation of successive boundaries supplies at least one O letter between Z-runs. Conversely, every word has that unique boundary set. R=0 gives the single all-O word; when 2R>k+1 the binomial is zero.

Let r_k be the actual maximum allowable run count, an integer, which exists since all-O is always permitted. With n=k+1 and j=2r_k, eventually j<=n/2. Monotonicity of binomial coefficients in this range bounds the sum over even indices by `(j+1) binomial(n,j)`. For j>0 its base-two logarithm is at most `log2(j+1)+j log2(exp(1)n/j)`. Since j/n tends to zero and x log(1/x) tends to zero, this is o(k). If j=0 the count is 1. These observations prove subexponential word counts without treating a finite census as an asymptotic argument.

## 4. Summing lengths and handling constants

Write a_k for the number of bounded-error words of length k. For a requested eta>0, set alpha=eta/2. There exists k0(E,alpha) such that a_k<=2^(alpha k) for k>=k0. Let M=sum of a_k over k<k0. With K=floor(log2 N),

`C_E(N)<=M+sum_(k=0)^K 2^(alpha k)<=M+D N^alpha`,

where `D=2^alpha/(2^alpha-1)` is fixed. For N>=1 this is at most `(M+D)N^alpha`, hence at most N^eta once `N>=max(2,(M+D)^(1/alpha))`. This supplies the precise quantifiers and absorbs both the geometric-sum constant and the finitely many exceptional lengths. If E=0, the only indices are 2,4,8,... and `C_0(N)=floor(log2 N)` for N>=2. Nonempty counts make the standard notation N^{o(1)} legitimate. No uniformity in growing E is claimed or proved.

## 5. Product-bound inheritance and exact obstruction scope

Section 4 of the frozen reviewed core proves `Psi_W(u,v)>=uv/2` on actual nontrivial sources for every central word, hence for all twelve listed words. The independent checker reconstructs the physical matrices and cone coordinates, verifies `P^t Q P >= (1/2) z z^t`, z=(1,1,2,2), and verifies that O,Z,J have nonnegative cone-coordinate matrices and do not decrease cover z^t y. These finite 4-by-4 coefficient identities are all-order certificates: iteration preserves the inequalities for arbitrary words and all nonnegative cone vectors. The actual seed has cone coordinates (0,0,0,1), so induction places every actual binary source in the cone. Direct coefficient checks for the explicit twelve words provide an additional control, not the premise of a target census.

For u,v>=2 and value at most N, `uv<=2N` implies separately u<=N and v<=N. Each ordered source pair and one template determines one value. Thus the number of DISTINCT represented targets with BOTH errors <=E is bounded above by `12 C_E(N)^2`. Collisions only reduce this number. Allowing both source orders is important and already included. Given any beta>0, apply the source bound with exponent beta/4 and absorb the fixed factor 12 for sufficiently large N to obtain a target count <=N^beta. A tail of odd integers has linear counting growth, so it cannot be contained in this subpolynomial set.

This is not a target-by-target lower bound, not an exclusion of a bounded first error alone, and not a refutation of unrestricted occupancy. It does not cover thresholds growing with N. For two distinct fixed thresholds use their maximum. For any other fixed finite menu with bound `Psi>=c uv`, c>0, source cutoffs become at most N/(2c), a fixed scaling of N; subpolynomiality and the conclusion persist. Singleton sources without a comparable cutoff would require separate treatment and are expressly excluded.

## 6. Explicit hostile and scope controls

For (ZO)^m, the value before each pair obeys `t_(j+1)=(t_j+1)/4`, starting at 1/2. Hence `t_j=1/3+(1/6)4^(-j)`. An isolated Z contributes `(t_j+1)/2=2/3+(1/12)4^(-j)`; summing gives `2m/3+(1-4^(-m))/9`. This is valid also for m=0 and grows linearly with length 2m. It demonstrates unbounded errors, but the counting obstruction follows from the independent universal lower bound and entropy argument, not this example alone.

## 7. Reproduction and result scope

Author verifier: frozen source, normal and `-O`, all 32767 words through k=14, thresholds 0,1/2,1,2, hostile m=1,...,80. Both reruns match each other and the frozen supplied output byte-for-byte.

Independent verifier: `independent_check.py`, no author imports, derives g from the scalar integer recurrence and enumerates actual indices rather than trusting an error accumulator. It checks the actual-to-word inverse, reversal labels, run error against g(u)-ut, per-run and Jensen inequalities, exact binomial run counts, threshold counts, all actual cone memberships, explicit twelve-template product coefficients, zero-error cumulative counts at every N=2,...,4096, and the alternating family for m=0,...,100. Normal and optimized outputs must agree. Both independently regenerated threshold tables equal the author's.

The asymptotic results are **PROVED by the written reasoning above**. Exhaustive bounded ranges are **FINITE-EXACT controls**. Matrix coefficient identities are **exact algebraic certificates**. No target coverage census was performed or relied on. No numerical root finder, floating-point asymptotic inference, or external runtime dependency is used.

Reproduce from the audit folder using `python3 independent_check.py` and `python3 -O independent_check.py`; author runs are `python3 inputs/session201-bounded-carry-obstruction/verify.py` and its `-O` variant. The independent checker compares its table to the included AUTHOR-normal.json and verifies frozen-input hashes. The audit manifest covers all deliverables except itself.
