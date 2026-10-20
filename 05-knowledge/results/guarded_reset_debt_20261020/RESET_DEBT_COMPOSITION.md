# Guarded reset-debt composition: exact closure and actual Collatz banks

Date: 2026-10-20. Status: **PROVED elementary identities and sufficient guarded-cylinder certificates**. The identities are standard affine algebra, not a literature-novelty claim. Source coverage, universal mean-current cancellation, and Collatz convergence remain **OPEN**. This is a staging artifact, not a repository promotion.

## 1. Inheritance, types, and the question actually answered

The three requested predecessor files do **not** define one common pair `(lambda, delta)`:

- `05-knowledge/results/collatz_reset_debt_20261015f.md`, section 3, proves a native signed-row norm recurrence `Q_j <= rho Q_(j-1) + (1-rho) R_j`. Its one-reset `R_j <= 7X` premise fails at scale `k=23`, though its finite discounted and paired repairs pass the declared tests.
- `05-knowledge/results/collatz_refinement_debt_20261016a.md`, sections 3–4, transports a signed Gram refinement remainder under a common native map. Its entire refinement charge is bounded once, not charged afresh at every reset. The unpaid coarse mean-current estimate is separate.
- `05-knowledge/results/reset_20260926_inequality.md`, sections 3–6, has a **multiplicative** split-product imbalance factor. Canonical all-ones first-collision episodes can grow without bound. That factor is not an affine additive delta without an additional change of carrier.

The closest prior affine mechanism is `05-knowledge/results/collatz_bank_reset_20261019.md`, section 1: contraction fixed-point envelopes compose over an arbitrary guarded bank. The present note identifies the exact reset-debt normalization, its inverse convention, exact valuation domains, and concrete grouped-reset banks. It does not relabel existing affine closure as a new Collatz breakthrough.

Anchor: source-priced affine debt. Niche: inverse branches and exact image guards. Wildcard: arbitrary finite prefixes versus a single infinite integer orbit. Live coordinates: chronological order; multiplier; debt units; residue address; literal first-exit clock; reset count.

## 2. Choose the forward convention before composing

For a nonnegative scalar resource use

    q_out <= (q_in + delta)/lambda,     lambda > 0.       (2.1)

A nonnegative **debt** additionally requires `delta >= 0`. This is not the same forward map as `q_out = lambda q_in - delta`.

For blocks 1 then 2 in chronological order,

    Lambda = lambda_1 lambda_2,
    Delta  = delta_1 + lambda_1 delta_2,                 (2.2)
    q_2 <= (q_0 + Delta)/Lambda.

Induction gives, with `Lambda_0=1` and `Lambda_j=product_(i<=j) lambda_i`,

    Delta_m = sum_(j=1)^m Lambda_(j-1) delta_j.          (2.3)

Positivity of lambda is what permits substitution in inequalities. Equality holds when each block is an exact affine map. No commutation or ordering hypothesis is needed beyond using the actual chronological order.

Set

    eta_j = lambda_j - 1 - delta_j.

Then

    Lambda_m - 1 - Delta_m
      = sum_(j=1)^m Lambda_(j-1) eta_j.                  (2.4)

Consequently `delta_j <= lambda_j-1` for each block implies `Delta_m <= Lambda_m-1` for **every finite m**, with no reset-count cutoff. The result also allows blocks which fail locally: a whole prefix or macroblock passes if and only if its weighted total on the right of (2.4) is nonnegative. Later slack can pay earlier negative slack, but must be proved on the actual concatenation.

For `lambda_j>1`, define normalized debt `d_j=delta_j/(lambda_j-1)`. The composite normalized debt is the convex combination

    Delta_m/(Lambda_m-1)
      = sum_j [Lambda_(j-1)(lambda_j-1)/(Lambda_m-1)] d_j. (2.5)

The coefficients are positive and sum to one. This statement needs `lambda_j>=1` (with zero-weight identity factors treated separately). It is false to call it a convex average when factors cross below one.

A tempting multiplicative defect is not multiplicative. Put `r_j=(1+delta_j)/lambda_j`. For two blocks,

    r_12 - r_1 r_2
      = delta_2 eta_1/(lambda_1 lambda_2).               (2.6)

If `delta_2>=0` and `eta_1>=0`, in particular when both blocks lie in the nonnegative unit-budget cone, the difference is nonnegative, not nonpositive. Nonnegative debts alone do not determine its sign. The correct quantity is the weighted additive slack (2.4).

### Signs and boundary cases

- `delta>=0` together with `delta<=lambda-1` forces `lambda>=1`.
- At `lambda=1`, the cone forces `delta=0`; this is a neutral resource bound, not strict progress.
- A forward physical multiplier `rho<1` corresponds to `lambda=1/rho>1` in (2.1).
- If signed delta is allowed, algebraic closure still holds for positive lambda, but signed delta is no longer a nonnegative discarded-debt quantity.
- Closure itself is not a progress theorem. The map `q -> q` satisfies it. Even infinitely many finite-prefix certificates do not prove source coverage.

## 3. Invariant region, descent, inverse map, and the native debt normalization

For exact maps satisfying the unit budget,

    q_out - 1 = (q_in - 1 - eta)/lambda.                 (3.1)

Thus the interval `q<=1` is forward invariant (intersect with the actual nonnegative domain), and a source `q>1` strictly decreases whenever `lambda>1`. The stronger condition `eta>0` is unnecessary for strict decrease above the barrier, but supplies an additional additive gain.

For a barrier `B>0` in original units `Q=Bq`, the condition is

    Q_out <= rho Q_in + e,  rho=1/lambda,
    e <= B(1-rho),                                      (3.2)

or, for a physical affine map `F(x)=(x+d)/lambda`, `d<=B(lambda-1)`.

The predecessor native recurrence translates exactly as

    q_j = Q_j/B,
    lambda_j = 1/rho_j,
    delta_j = (lambda_j-1) R_j/B.                        (3.3)

Therefore its unit debt budget is exactly `R_j<=B`, not an independent estimate. Equation (2.3), divided by Lambda, is precisely the usual discounted residual convolution. The `k=23` one-reset failure is not repaired by affine closure; it must be paid by verified macroblock slack or a proved discounted estimate. Refinement's nonaccumulating Gram remainder remains a separate carried sidecar.

For an exact forward block, the inverse real affine map is

    G(y)=lambda y-d.

`d<=B(lambda-1)` makes the ray `y>=B` invariant under G. Its chronological composition law is `d_12=lambda_2 d_1+d_2`, because these are now inverse maps applied in their own chronological order. Reversing a forward word reverses the inverse order. An inverse Collatz branch is **not** legal on the whole real invariant ray: its exact integer image congruence must also hold. Ray preservation alone proves no physical descent and creates no new orbit choices.

## 4. Exact odd-valuation words retain all addresses

Let `U(n)=(3n+1)/2^v2(3n+1)` on positive odd integers. A finite valuation word is `a=(a_1,...,a_k)`, each `a_i>=1`. Define

    A_0=0, C_0=0;
    A_j=A_(j-1)+a_j,
    C_j=3 C_(j-1)+2^A_(j-1).

Its affine branch is

    n_j = (3^j n + C_j)/2^A_j.                          (4.1)

For `A=A_k`, `S=3^k`, `M=2^A`, `C=C_k`, the word's lawful positive-odd cylinder is exactly

    n = r + 2M t,  t>=0,
    r = S^(-1)(M-C) mod 2M,  0<r<2M.                   (4.2)

The residue r is odd. Unlike plain final integrality modulo M, this modulus also enforces the **odd terminal** and therefore exact valuations.

Proof of sufficiency: the full congruence gives an odd integer terminal. For every `j<k`, write `C_k=3^(k-j) C_j+2^A_j R`. Modulo `2^(A_j+1)`, R is odd, since its first term is odd and all later terms contain at least one further factor two. Since `A_k>=A_j+1`, the final congruence implies `3^j n+C_j = 2^A_j mod 2^(A_j+1)`. Thus every intermediate value is an odd integer; (4.1) then proves each exact valuation. Necessity follows by the same recurrence.

The full shortcut-clock address of a reset step `a_j` consists of its odd shortcut followed by `a_j-1` literal halvings. Total shortcut time is A, not k and not the number of resets.

For ROOT-stopped dynamics require `n_j>1` at every proper odd-step prefix. These are explicit affine inequalities on (4.2), so they trim only a finite head of its lifts. Never continue through ROOT merely because the unstopped affine word allows repetitions of valuation 2.

The actual inverse-image domain also has an exact description. Let `z_0=(S r+C)/M`, an odd positive integer. Forward lifts give

    z=z_0+2S t,  t>=0.                                 (4.3)

Only these image values invert to the cylinder in (4.2), with additional ROOT-stop restrictions if applicable.

### Sharp cylinder-specific descent test

For any complete stopped cylinder with least allowed positive source `n_min`, all its sources strictly descend at the endpoint if and only if

    M>S and (M-S)n_min>C.                              (4.4)

For `M<=S`, C is positive for every nonempty odd word, so no positive source descends. For `M>S`, the endpoint displacement is strictly decreasing in source n, proving (4.4). Unit-barrier closure is sufficient but not necessary for (4.4).

All finite valuation words have positive integer sources, since (4.2) has arbitrarily large lifts; this does not provide any prescribed source with a lawful next word.

## 5. An infinite bank of actual first-reset blocks

In sections 5–7 the resource is the actual positive integer height, and a reset means an odd step of exact valuation at least two. This is not the precision-collision convention of the 2026-09-26 split-product note, nor one simultaneous terminal-graph closure of the native signed-row controller. These examples do not establish its `7X` norm bound.

A first-reset block has odd valuation word

    (1,...,1,a) = (1^h,a),  h>=0, a>=2.

There are exactly h valuation-one steps before the first valuation at least two. Direct induction gives

    F_(h,a)(x)
      = [3^(h+1) x + 3^(h+1)-2^(h+1)]/2^(h+a).         (5.1)

In original units with B=1,

    lambda=2^(h+a)/3^(h+1),
    delta=1-(2/3)^(h+1).                               (5.2)

The unit budget holds exactly when

    2^(h+a) >= 2*3^(h+1)-2^(h+1).                     (5.3)

For every h some integer a>=2 satisfies (5.3). Hence the bank includes arbitrary long valuation-one runs, not just isolated high valuations. Every such block has slope

    rho <= 1/[2-(2/3)^(h+1)] <= 3/4.                  (5.4)

For every finite sequence of these blocks, concatenate their complete valuation words and use (4.2) with the ROOT-prefix inequalities. The resulting nonempty cylinder satisfies the composed budget and every block lowers its input above 1. This is a genuine all-length certificate with explicit reset, valuation, residue and clock information.

The first shortcut descent within one such first-reset block is also exact. Before the reset, all h odd endpoints strictly exceed the original x. Put

    N=3^(h+1)(x+1)-2^(h+1).

At reset substep j, `1<=j<=a`, the actual value is `N/2^(h+j)`. The first descent is the least j satisfying

    N < 2^(h+j) x.                                    (5.5)

Condition (5.3) ensures it exists for x>1; it occurs at literal time h+j. This retains the distinction between first odd reset boundary and first literal descent.

**Coverage hostile.** The first-reset word `(1,2)` has map `(9x+5)/8` on `x=11 mod16`. Every such source grows across the reset, including `11 ->17 ->13`. It has `lambda=8/9`, `delta=5/9`, so no nonnegative unit-budget certificate can cover that block without grouping or changing the resource. The bank (5.3) is not a partition of all actual first-reset episodes.

## 6. A concrete grouped bank that pays an actual growing reset

Use these macroblocks:

    A: odd valuations (1,2,3), F_A(x)=(27x+23)/64,
       x=11 mod128, lambda_A=64/27, delta_A=23/27;
    B: odd valuations (1,3), F_B(x)=(9x+5)/16,
       x=19 mod32, lambda_B=16/9, delta_B=5/9.

They contain two resets and one reset respectively. Both satisfy the unit budget; their slacks are `14/27` and `2/9`.

For A's whole cylinder, write `x=11+128t`. Its odd endpoints are

    11+128t -> 17+192t -> 13+144t -> 5+54t.

The first reset `(1,2)` grows by `2+16t` at **every** lift, yet the following reset 3 pays the whole macroblock. This is not merely a bank in which each individual reset was assumed to pay.

The literal shortcut values are

    11+128t, 17+192t, 26+288t, 13+144t,
    20+216t, 10+108t, 5+54t.

Thus the first literal descent is uniformly at shortcut step 5, not at the macroblock's final step 6. For B, the values are

    19+32t, 29+48t, 44+72t, 22+36t, 11+18t,

so its first descent is step 4.

Every finite word over `{A,B}` has the unique complete cylinder (4.2); both the valuation guards and reset boundary counts are retained. Its composed unit debt certificate holds for any order and length. Arbitrarily large reset counts are certified. This proves a guarded bank theorem, not that an arbitrary positive integer can keep choosing A or B.

### Unit budget is not necessary even for an entire actual cylinder

The block E with valuations `(1,1,3)` has

    F_E(x)=(27x+19)/32, x=55 mod64.

Its unit debt is `19/27>5/27=lambda-1`. Nevertheless its fixed point is `19/5`, below its least source 55, so its entire lawful cylinder strictly descends; `55 ->83 ->125 ->47` is the first lift. Its literal first descent is step 5 for every lift. A barrier B=4 certifies E. The bank `{A,B,E}` therefore has a common barrier 4, uniform slope at most `27/32`, and all its legal block inputs are above that barrier. This is an explicit guard-aware enlargement of the unit bank.

## 7. What unbounded reset count does and does not quantify

The proved quantifier is:

    For every finite bank word w, every positive integer in its
    complete stopped cylinder satisfies the certificate for w.

It entails `for every depth there exists a source`, because each finite cylinder has arbitrarily large lifts. It does **not** entail `there exists one positive source that supports every depth`.

Two exact controls make the distinction sharp.

1. Repeating valuation 2 r times has sources `n=1+2*4^r t`; for t>=1 all r odd endpoints stay above ROOT. Its jth endpoint is `1+2*3^j*4^(r-j)t`. Every finite depth exists. The nested full cylinder intersection is just the integer 1; ROOT-stopped 1 supports no reset at all.
2. Repeating A gives

       F_A^r(x)=(27^r x+C_r)/64^r,
       C_r=23(64^r-27^r)/37.

   Every finite word has positive integer sources. Modulo `64^r`, its source address is `23/37`. The rational `23/37` is a fixed point of F_A and is `11 mod128` in the 2-adic integers, so it belongs to every full repeated cylinder, including the odd-terminal guard. Its nested cylinder intersection in the 2-adic integers is exactly `{23/37}`; the intersection among positive integer sources is empty. An infinite formal/2-adic branch can thus coexist with arbitrarily deep finite integer cylinders while no positive integer follows it forever.

If a **fixed** source had an infinite lawful sequence satisfying `lambda_j>=4/3`, the budget would give

    n_m-1 <= (3/4)^m (n_0-1).

Positive odd values above ROOT have `n_m-1>=2`, so such a nonterminal sequence is impossible. This only proves that the source must reach ROOT **or leave the prescribed bank**. Showing that the actual orbit always has a further certified macroblock is the missing source-coverage theorem. The analogous integer descent argument works without a uniform multiplier gap whenever every available nonempty block strictly decreases the same positive integer potential.

## 8. Exact validation and remaining target

`verify_reset_debt.py` checks rational composition and slacks, independently reconstructs cylinders by brute-force exact valuations, replays the actual literal shortcut trajectory, verifies the grouped-bank first exits, includes the growing `(1,2)` and nonnecessary-unit-budget E controls, and checks arbitrary short words plus much longer selected schedules. Every assertion is an explicit exception check, so `-O` does not remove validation. No census is promoted to an all-length theorem: the proofs above establish that quantifier.

The strongest remaining targets are not more affine algebra. They are:

- a source-preserving rule proving that an arbitrary unresolved positive integer reaches a next certified macroblock;
- or, for the native signed-row controller, a uniform weighted residual/discounted estimate that pays the actual original rows and common terminal labels;
- plus any application-specific first-exit, endpoint-balance and refinement sidecars.

Neither a positive-density cylinder, a closed scalar cone, arbitrary finite depth, nor inverse-ray invariance supplies those obligations.
