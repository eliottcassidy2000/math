# An alphabet-three shield, a dyadic extinction formula, and a refutation of the size-six logarithmic law

**Session:** 2026-10-19, Gilbreath / cellular-automaton lane.  
**Status:** PROVED elementary structural theorems below; PROVED finite-state all-length certificates at explicitly listed ages; FINITE-EXACT numerical tables; OPEN sharp decision horizon and complete asymptotic. No claim to settle Gilbreath's prime conjecture or to establish literature priority.

## 1. Inheritance and the changed question

The starting points are THM-4511, `01-canon/theorems/THM-4511-gilbreath-size-four-wall-theorem.md`, and HYP-9163, `05-knowledge/hypotheses/HYP-9163-gilbreath-size6-excess-dyadic-law.md`. The former supplies the `{0,2,4}` wall; the latter conjectures a logarithmic re-emission excess for a lone defect of size six. The historical warning that an infinite Pascal seed can hide an arbitrary future trace is retained from `05-knowledge/results/forest_20260926_gilbreath.md`. All objects here have finite state descriptions and finite exact certificates.

Divide the even tail by two. The automaton is

    a(t+1,i) = |a(t,i) - a(t,i+1)|.

The leading original 1 survives precisely while the first normalized tail cell stays in `{0,1}`. A normalized defect 3 is the original size-six defect. Write `p6(F)` for its eventual edge-destruction probability when it sits at column F, with F−1 independent fair binary cells to its left and an independent fair infinite binary sea to its right.

The live concepts were: absorbing alphabets; right-boundary information; the Pascal front transform; exact first-failure languages; the gap/depth staircase; and the distinction between random-sea and prime-specific hypotheses. The key change is to replace an absorbing Markov chain's entire state graph by a stopping threshold in a zero extension.

### Primary literature checked

Chase–Hunter–Tao, [arXiv:2607.08712v1](https://arxiv.org/html/2607.08712v1), proves a Cramér-model analogue and a deterministic obstruction criterion; its Lemma 3.8 immediately implies the maximum-rigidity fact used below. Its Proposition 4.1 is relevant to random depth bounds. These are cited dependencies, not new results here. Odlyzko's [1993 paper](https://www-users.cse.umn.edu/~odlyzko/doc/arch/gilbreath.conj.pdf) already studies prime-difference depth. Plouffe's [2025 preprint](https://arxiv.org/abs/2510.06688) reports a larger numerical verification. These references do not prove the actual prime conjecture.

## 2. Definitions and maximum rigidity

For a finite word w of nonnegative integers, let D(w) be its adjacent absolute differences. Call w **good** if the first entry of every nonempty word D^t(w), including t=0, is at most 1. A good word over `{0,1,2,3}` may contain a 3 away from its left edge; for example `(0,1,3,2)` is good.

**Maximum rigidity.** If all initial values of a finite difference triangle lie in `[0,M]` and its apex equals M, every initial value is in `{0,M}`. This holds for real values as well. One proof starts with the singleton apex and repeatedly pulls its `{0,M}` block up one row: the parent lies in a single coset modulo M, and the bound `[0,M]` forces the endpoint coset. Equivalently, apply the parentage argument of the cited Lemma 3.8; no probabilistic input is involved.

## 3. The alphabet-three extension theorem

**Theorem 1 (extension rules).** Let w be a nonempty good word over `{0,1,2,3}`.

1. The word w1 is good.
2. The word w2 is good if and only if w is not entirely zero.
3. If w is zero, or its rightmost nonzero letter is 1 or 2, then w0 is good.

Consequently, if a good word ends in 1 or 2, every finite or infinite binary continuation preserves its good left edge.

### Proof by simultaneous induction

For a one-letter good word, the only possibilities are 0 and 1; the three assertions are checked directly. Suppose they hold for shorter words. The word D(w) is good, and

    D(wb) = D(w) followed by |last(w)−b|.

For assertion 1, when last(w) is 0 or 2 the appended difference is 1, so use assertion 1 for D(w). When last(w)=3 it is 2; the exception in assertion 2 would require D(w) to be zero, making w constant 3 and hence not good. When last(w)=1 it is 0. The rightmost nonzero entry of D(w), if any, is the jump into w's final constant run of 1s, and that jump is at most 2. Assertion 3 applies.

For assertion 2, last(w)=1 or 3 gives a final difference 1 and invokes assertion 1. Last(w)=0 invokes assertion 2: D(w) is zero exactly when w is constant, here entirely zero. Last(w)=2 invokes assertion 3, since the last nonzero jump into a constant run of 2s is at most 2. For an all-zero w, appending 2 gives final apex 2, proving the exceptional direction.

For assertion 3, if the last letter is zero and w is nonzero, the rightmost nonzero entry of D(w) equals the rightmost nonzero entry of w. Apply assertion 3 again. If the last letter is 1, use assertion 1. If it is 2, use assertion 2; its constant-2 exception is incompatible with goodness. The all-zero case is immediate.

The last conclusion follows by appending bits one at a time: a new 1 becomes a rightmost nonzero marker, and appended zeros preserve the marker. ∎

### Independent geometric proof and boundary

An independent audit in `04-computation/experiments/gilbreath_alphabet_three_20261019/review/PREFIX_SHIELD_REVIEW.md` gives a backward-diagonal proof using maximum rigidity. In a supposed first bad apex beyond the marker, all cells on the relevant backward diagonal are at most 2: their cones either contain the interior marker or lie entirely in the binary tail. A first bad value must therefore be 2. Pulling backward, a row beginning with zeros followed by 2 forces the preceding row to begin with one more zero followed by 2; the alternative needs a 3 on that diagonal. Eventually the forced 2 enters a binary cone, a contradiction.

The alphabet restriction is sharp for this statement. The word `(1,2,3,4)` is good, but appending 1 produces edge `(1,1,0,0,2)`. A marker is not a universal shield at maximum 4.

## 4. Exact reduction of size-six risk to a zero-run threshold

For a binary word w of length F−1 put v=w3, and define

    tau(w) = min{r≥0 : v followed by r zeros is not good},

with tau=∞ if no such r exists. For tau=0, the original finite word is already bad.

**Theorem 2 (geometric zero-run formula).** Conditional on w, the probability that the original leading 1 is ever destroyed is

    2^(−tau(w)),     with 2^(−∞)=0.

Consequently

    p6(F) = 2^(1−F) sum over binary w of length F−1 of 2^(−tau(w)).       (1)

In particular p6(F) is a dyadic rational for every F.

**Proof.** Let R be the number of zeros in the right sea before its first 1. If R≥tau, the bad finite prefix has already appeared and destruction follows. If R<tau, the prefix v0^R is good; Theorem 1 says appending its first 1 preserves goodness and then shields it from every remaining binary continuation. Thus destruction is exactly `{R≥tau}`. Fair independent bits give probability `2^(−tau)`.

There are only finitely many left words, so (1) is a finite sum of dyadic rationals. No asymptotic or unproved decision-horizon bound is needed. For an algorithmic bound, evolve the F-cell strip with constant zero boundary and remember its states. A repeated state without a previous bad left entry proves tau=∞. There are at most `4^F` states, so the computation terminates. This coarse finite bound is enough for exact decidability; the much sharper observed dyadic horizon remains a separate question. ∎

This also explains why the Markov chain's mixed-probability core is acyclic: before the first right-boundary 1 there is only one potentially unresolved branch, the all-zero branch. The first 1 fixes the eventual verdict. An unresolved zero cycle would have no future bad zero prefix and hence would actually be safe.

## 5. Rotate the triangle: a small automaton for each exact delay

Let L=F−1 and reverse the left sea into v_j=w_(L−1−j). Define the front coordinates

    alpha_s = sum_{j=0}^s binom(s,j) v_j modulo 2,    0≤s<L.

This is the lower triangular Pascal transform over F2, with diagonal 1; it is a bijection on binary words. It is exactly the pure-binary cell encountered by the left-moving first defect front.

For zero right boundary, use s for distance leftward from the defect and k for age behind the first front. Put x_0(k)=3 for all k≥0. Then

    x_s(−1)=alpha_(s−1),
    x_s(k)=|x_s(k−1)−x_(s−1)(k)|,    k≥0.                              (2)

This is the original update in rotated coordinates: the physical cell is at column F−1−s and time s+k. In particular x_L(k) is the left-edge value at time L+k. All earlier left values depend only on the binary left sea and are safe.

For fixed r, scan alpha through the state `(x_s(0),...,x_s(r))`, starting at `(3,...,3)`. The transition on bit b is

    u_0=|b−v_0|,  u_1=|u_0−v_1|, ..., u_r=|u_(r−1)−v_r|.              (3)

The final state accepts exactly the left words with tau=r when its first r coordinates are at most 1 and its final coordinate is at least 2. This is a finite automaton, so exact threshold counts have rational generating functions in L.

The reproduction program compares (2)–(3) against direct shrinking triangles for every binary left sea through F=11, at ages 0 through 6: 2,047 seeds. It also verifies the Pascal-front bijection independently at every width.

## 6. A complete delay-two language and a rigorous refutation

**Theorem 3.** The binary front words with tau=2 are exactly

    (0000)*11       or       (0000)*1110*.                              (4)

Thus their number at distance F is

    A_2(F) = floor((F+1)/4).                                          (5)

**Finite-state proof, valid at every length.** For r=2, (3) has 23 reachable states. Its only states from which an accepting state can still be reached are

    A0=(3,3,3), A1=(3,0,3), A2=(3,3,0), A3=(3,0,0),
    B=(2,1,2), C=(1,0,2), E=(0,0,2).

A zero cycles A0→A1→A2→A3→A0. The surviving one-transitions are A0→B→C→E. The state E has a zero loop. States C and E accept. Every unlisted transition leaves this seven-state set permanently. These facts follow by direct substitution into (3); the full 23-state transition table is in `gilbreath_certificates.json`. They prove exactly (4). For length L, each j with `4j+2≤L` contributes one word: the terminal `11` when equality holds, otherwise `111` followed by zeros. This gives (5). ∎

There are exactly F words with tau=0: along the first front the value starts at 3 and loses one for each encountered alpha=1 until it enters the binary sea permanently. It remains dangerous at the left edge exactly when alpha has at most one 1. These F words are disjoint from (4).

Define the inherited scaled excess by

    c6(F) = 2^F [p6(F) − F 2^(1−F)].

Every delay-two word contributes 1/4 to the sum in (1). Therefore

    c6(F) ≥ (1/2) floor((F+1)/4),
    liminf p6(F)/(F 2^(1−F)) ≥ 17/16.                                (6)

**REFUTED:** HYP-9163's claims `c6(F)=1.4 log_2 F+O(1)` and `p6(F)/(F2^(1−F))→1` are false. The error was inferring a logarithmic asymptotic from one apparent doubling in a short, dyadically structured table. The missing coordinate is the multiplicity of distinct left words with a fixed short re-emission delay. That multiplicity grows linearly.

The front-only probability is still a valid lower bound. The shielding theorem and dyadic-denominator phenomenon survive and are now proved. No inference about actual prime failure probabilities is justified by this random-sea model.

## 7. Stronger all-length finite-state certificates

For each fixed age r listed below, the same exact DFA proves

    A_r(F) = C_r floor((F−1+h_r)/(2h_r)).                            (7)

The pairs `(r; h_r,C_r)` are

    (2;2,1), (4;4,3), (6;4,1), (8;8,7),
    (10;8,2), (12;8,3), (14;8,1).

These are all-F statements certified by finite integer algebra, rather than fits to a finite F table. To make the validation finite, discard DFA states that cannot reach acceptance, leaving m states. The actual counting generating function has denominator degree at most m, from `(I−zA)^−1`. The proposed generating function is

    C_r z^h_r / [(1−z)(1−z^(2h_r))].                                (8)

After subtraction, the numerator has degree at most m+2h_r. Checking its first m+2h_r+1 coefficients therefore proves equality identically. The program performs exactly these checks. The coaccessible state counts are respectively 7,18,15,52,34,50,31. It also proves no accepting states exist at odd ages 1,3,5,7,9,11,13. **Update, 2026-10-19:** absence at every odd age is now PROVED by the [paired-language invariant](gilbreath_even_delay_20261019/GILBREATH_EVEN_DELAY_THEOREM.md); the fixed-age certificates above remain unchanged.

Adding the nonnegative contributions in (7) yields

    liminf c6(F)/F ≥ 23533/131072,
    liminf p6(F)/(F2^(1−F)) ≥ 285677/262144
                               = 1.089771270751953125.              (9)

The elementary hand proof (6) already refutes the inherited asymptotic. Formula (9) is a stronger finite-state consequence with its own explicit certificate.

## 8. Exact computations and corrected old approximations

Two independent exact implementations agree through F=15:

- `probe_size6.py`: enumerate the finite absorbing chain, classify probability-zero and probability-one states by directed reachability, then solve the acyclic mixed core over exact rational numbers.
- `zero_extension_scan.py`: use (1), bitplane zero-boundary evolution, and an explicit repeated-state certificate for every word declared safe. This checks all binary left seas through F=23. It does not assume the conjectural sharp horizon.

Selected exact values:

| F | p6(F) | c6(F) |
|---|---:|---:|
| 12 | 212717/33554432 | 16109/8192 |
| 13 | 232429/67108864 | 19437/8192 |
| 14 | 248813/134217728 | 19437/8192 |
| 15 | 269293/268435456 | 23533/8192 |
| 16 | 285677/536870912 | 23533/8192 |
| 17 | 19796138989/70368744177664 | 1542527981/536870912 |
| 23 | 26993564653/4503599627370496 | 2297502701/536870912 |

The historical floating value `c6(13)≈2.372697` is replaced by the exact `2.3726806640625`; the corresponding F=15 and F=17 rounded values also require updating. This does not invalidate the historical exact rationals through F=11. The earlier floating algorithm is absent from the present checkout, so the cause of its numerical error has not been established.

## 9. Precise remaining hypotheses

1. **PROVED parity; OPEN sharp horizon.** Every finite positive tau is even, by the [paired-language theorem](gilbreath_even_delay_20261019/GILBREATH_EVEN_DELAY_THEOREM.md). The remaining conjecture is that tau at distance F is at most `2^ceil(log_2 F)−2`; the maximum equals this bound for F≥3 in the exact tested range. Neither parity nor this bound is needed for Theorems 1–3 or (9).
2. **OPEN all-age counting recurrence.** The block form (7) appears to continue: for `2^k≤r≤2^(k+1)−2`, even r, set h_r=2^k. The coefficients have the recursive pattern visible in the exact histograms. Only the seven ages in Section 7 are certified here.
3. **OPEN repaired asymptotic.** The evidence suggests `c6(F)=C F+O(1)` and a non-unit limiting ratio near 1.08978. Equation (9) is a lower bound, not a determination of C or a proof of convergence. Establishing the all-age recurrence and a summable tail bound would be a concrete next step.
4. **OPEN paired halving.** The excess halves from odd F to F+1 in the tested data. The all-age version has not been proved here; fixed-age formulas (7) have this property individually.

## 10. Reproduction and scope ledger

From `04-computation/experiments/gilbreath_alphabet_three_20261019/`:

    python verify_gilbreath_results.py > gilbreath_certificates.json
    python -O verify_gilbreath_results.py > gilbreath_certificates_optimized.json
    cmp gilbreath_certificates.json gilbreath_certificates_optimized.json
    python probe_size6.py 15 > probe_size6.out
    python zero_extension_scan.py 3 23 > zero_extension_scan.out

The primary certificate uses only the Python standard library and explicit failures that remain active under optimization. The fast all-word scan uses NumPy integer bit operations, not floating arithmetic. Every safe declaration in that scan has an exact cycle certificate; repeated observations without a repeated state are not accepted as safety.

**Connection contract.** Source: a lone normalized defect 3 and two binary seas. Target: a zero-run stopping threshold and a fixed-age finite automaton. Map: Pascal front transform plus rotated coordinates (2). Preserved predicate: the exact first dangerous left value and its age. Lost data: right-sea bits after the first 1. Restoration sidecar: none for this predicate, by Theorem 1; they remain relevant to other observables. Cheapest hostile: maximum 4, where the shield fails. The prime arithmetic and correlations have not been mapped to independent fair bits.
