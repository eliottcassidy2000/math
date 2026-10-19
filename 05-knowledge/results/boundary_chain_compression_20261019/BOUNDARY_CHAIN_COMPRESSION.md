# Bounded chain exteriors, gap ratios, and shared degree budgets

19 October 2026. Status: PROVED and root-reviewed; exact finite application independently reconstructed by the supplied verifier. This is a necessary-boundary reduction for a fixed full root core. It does not prove seven-thin balance or weak Gold, and no historical novelty is claimed.

## Inheritance and the live obstacle

The preceding [degree-cap theorem](../boundary_increment_caps_20261019/boundary_increment_caps.md) (packet147) proves child increment <= exterior incomparability budget times parent increment. Its hostile example has a 408-extension core and a factorized formal weighting that cannot come from actual exteriors. The earlier one-slack split (packet145) shows why convexifying distinct exterior allocations can lose decisive information. Here the missing coordinate is the complete gap profile of one exterior extension. The main reduction retains that profile before taking any averages.

The live concepts are: full root cores; monotone gap words; bounded two-chain models; common lower/upper degree allocations; and fixed comparison/Gold inequalities. The map is from an exterior extension to terminal gap lengths. It preserves every core-extension lift count exactly, forgets the original exterior order beyond the terminal gaps, and restores averaging through exact extension-count weights.

## 1. Setting

Let P be a finite poset of maximum incomparability degree at most k. Suppose its vertices partition as Q, A, B, with two incomparable distinguished roots x,y in Q, such that:

- A is below both roots and B is above both roots;
- A<B;
- for a in A and q in Q, either a<q or a is incomparable with q;
- for q in Q and b in B, either q<b or q is incomparable with b.

These are the full critical-root-core hypotheses used in packet147. Let d_Q(q) be q's incomparability degree inside Q and delta(q)=k-d_Q(q). In particular

    s_A(q)+s_B(q) <= delta(q),

where s_A,s_B count actual exterior incomparable vertices. No probability or majority-order hypothesis is needed for the lemmas below.

## 2. Exact bounded-chain compression theorem

Fix a linear extension sigma=(a_1,...,a_m) of A. Put

    alpha_sigma(q)=m-max({i : a_i<q} union {0}).

Every A vertex after the last predecessor of q is incomparable with q. Therefore

    0 <= alpha_sigma(q) <= s_A(q).

Similarly, for a linear extension tau=(b_1,...,b_n) of B, put

    beta_tau(q)=min({j : q<b_j} union {n+1})-1.

Every B vertex before the first successor of q is incomparable with q, so beta_tau(q)<=s_B(q). Both profiles vanish at x,y. If q<r in Q, then alpha(q)>=alpha(r) and beta(q)<=beta(r). Their supports have size at most k: when A is nonempty, every q with alpha(q)>0 is incomparable with its final vertex a_m; use the first B vertex for beta.

Construct C_(alpha,beta) by retaining Q, adding a lower chain c_1<...<c_h with h=max alpha(q), and an upper chain d_1<...<d_j with j=max beta(q), with all lower-chain vertices below all upper-chain vertices. Declare

    c_i<q iff i<=h-alpha(q),
    q<d_t iff t>beta(q).

All other chain/core pairs are incomparable. The profile monotonicity makes these relations transitive and creates no additional relation inside Q. The roots have profile zero, so the lower chain lies below both roots and the upper chain above both roots.

Each core vertex has exactly alpha(q)+beta(q) new incomparable neighbors, hence degree at most k. Each lower-chain vertex is incomparable with a subset of {q:alpha(q)>0}, a set of size at most k; the upper argument is dual. Both chains are internally total and mutually comparable. Thus C_(alpha,beta) is k-thin and has the same full root core Q. Each chain has length at most max delta(q).

THEOREM 1. For every core linear extension w, let C_P(w) count full extensions of P whose restriction to Q is w. Then the following identity holds simultaneously for all w:

    C_P(w) = sum over sigma in LE(A), tau in LE(B)
                 C_(C_(alpha_sigma,beta_tau))(w).                 (1)

Proof. Condition on the exact orders sigma and tau. Because every A vertex precedes both roots, A can interleave only with the prefix of w before its first root. In this merged word, q must follow the last A predecessor of q, and that is its only A constraint. Consequently its permitted A gaps are exactly the final alpha(q)+1 gaps. The initial A segment before all those gaps contributes only one forced placement once sigma is fixed, and deleting it loses no interleaving multiplicity. Those initial vertices need not be universal predecessors in the original P: their forced placement here follows from the fixed sigma order and the last-predecessor constraints. The replacement chain has exactly those same permitted gaps. The upper argument uses the initial beta(q)+1 B gaps after the last root. A<B makes the two choices independent. Thus the number of lifts for each conditioned pair sigma,tau agrees with the replacement-chain number, for each w. Sum over all distinct conditioned exterior orders. This proves (1).

More explicitly, let G(t_1,...,t_r) count weakly decreasing integers z_1>=...>=z_r>=0 with z_i<=t_i, and put G(empty)=1. If pre(w) and suf(w) are the words before the first and after the last root, respectively, then the common conditioned lift count is

    G(alpha(pre(w))) G(beta(reverse(suf(w)))).                   (2)

Every gap assignment determines exactly one merged word, including when several core letters share a gap, because their order is fixed by w.

Summing (1) over w gives e(P)=sum_(sigma,tau) e(C_(alpha,beta)). Hence the normalized restriction law on Q is a convex mixture of the normalized laws of these bounded chain models, with weights

    e(C_(alpha_sigma,beta_tau))/e(P),

aggregated over repeated profiles. These weights sum to one. They are generally NOT uniform over pairs of exterior orders. Equation (1) preserves all joint comparison-event numerators and denominators simultaneously, not just individual pair probabilities.

In a full root core each q has at least one incomparable core neighbor. Thus max delta<=k-1 and |Q|<=2k, giving at most 4k-2 vertices in each replacement model. For k=7 the universal bound is 26 vertices. The sharper actual bound is |Q|+max alpha+max beta.

## 3. A finite profile universe and the exact scope of transfer

For a fixed Q, it suffices to enumerate integer profiles alpha,beta satisfying:

- 0<=alpha(q),beta(q), alpha(q)+beta(q)<=delta(q);
- alpha(x)=alpha(y)=beta(x)=beta(y)=0;
- q<r implies alpha(q)>=alpha(r), beta(q)<=beta(r);
- each positive profile support has size at most k.

Every such pair gives an actual k-thin two-chain model. Conversely every exterior-conditioned pair from P belongs to this universe. Thus this is an explicit finite model family containing every atom needed in (1).

Any homogeneous linear inequality in core-extension counts which holds in every model holds in every P with that core. In particular, a single fixed two-comparison Gold strategy transfers if all its answer-path inequalities hold in every model. A fixed nonnegative weighted comparison inequality likewise transfers.

Important boundaries:

- One replacement model need not preserve all comparison probabilities.
- The statement that each model has some balanced pair does NOT by itself imply that their mixture has a balanced pair. The pair or separating weighted inequality must be common, or another argument must retain the mixture structure.
- Nor does “each model has some Gold strategy” automatically supply one strategy for a mixture.
- Arbitrary convex mixtures of these models are not asserted to be realizable by an exterior poset. Actual profile multiplicities arise from the two exterior posets, with independent choices sigma,tau and fixed common degree allocations. The convex hull is a necessary relaxation.

## 4. Elementary stochastic monotonicity of bounded gap words

For positive cap vector t=(t_1,...,t_r), let Omega(t) consist of weakly decreasing positive integers z_1>=...>=z_r>=1 with z_i<=t_i. Let F(t)=|Omega(t)| and F(empty)=1. Under the uniform distribution on Omega(t), the last coordinate z_r is stochastically nondecreasing in each cap t_i.

Proof by induction on r. For r=1 this is uniform distribution on {1,...,t_1}. Increase just cap t_i by one. The new tuples, if any, form the layer z_i=t_i+1. Conditional on any fixed z_i=a, the prefix and suffix choices factor; the suffix is uniform on the shorter gap family whose caps are min(t_j,a) for j>i. By induction its last coordinate is stochastically nondecreasing in a. Thus the new layer's last coordinate dominates the conditional law at every old value z_i<=t_i; it dominates their mixture. Mixing the old distribution with the new layer cannot decrease it. If i=r the new last coordinate is simply larger than all previous ones. If there is no new layer nothing changes. Repeating unit increases proves the claim. This proves ratio domination, rather than relying on the insufficient fact that F itself is increasing.

## 5. Sharpened ancestor ratios

Use packet147's normalized lower interleaving counts L(w), with U_empty=1 and U_(wq)=L(wq)-L(w). For a fixed exterior extension sigma,

    U_w(sigma)=F(alpha_sigma(w))                               (3)

for nonempty w: subtracting the interleavings where the final core letter lies after all A leaves exactly those whose every core letter has at least one A letter after it. These are the positive gap sequences. Averaging (3) over sigma gives the normalized increment U_w.

Suppose w has positive total-slack cap vector s=(delta(q_1),...,delta(q_r)). Define

    rho_s(a)=F(s,a)/F(s),   rho_empty(a)=a,   rho_s(0)=0.

For every actual lower exterior and any appendable letter q,

    U_(wq) <= rho_s(s_A(q)) U_w
            <= rho_s(delta(q)) U_w.                            (4)

Indeed, condition on sigma. Given a positive ancestor gap tuple, the number of ways to append q is min(z_r,alpha_sigma(q)). Its mean is at most the mean of min(z_r,s_A(q)) in Omega(s), by Section4. This yields a bound with the SAME constant for every sigma, so it remains valid after averaging the numerator and denominator. No ratio of averages is replaced by an average of ratios. If some cap in w is zero, U_w and every descendant are zero, handled separately.

For equal ancestor and child cap s,

    F(s,...,s)=binomial(s+r-1,r),
    rho_(s,...,s)(s)=(s+r)/(r+1).                               (5)

Thus at depth r the old factor s improves to (s+r)/(r+1). For every nonempty ancestor and a>=2, rho_s(a)<a, since the positive gap family has a positive-probability tuple with last coordinate 1. At a=1 the factor is exactly 1. The dual suffix inequality uses the reversed suffix cap word.

The cap bound is sharp as a statement about prescribed gap caps: a chain exterior with terminal gaps exactly s realizes the full gap family whenever the core constraints permit those profiles. It is not claimed that every numerical slack vector is compatible with every fixed core.

## 6. Opposite siblings and the allocation envelope

Let lambda_(u,v)=U_u V_v be actual lower/upper increment coefficients, with u,v disjoint core words. Let q not occur in either. Assume the parent and whichever children are used are supported boundary generators. The children are (uq,v) and (u,qv), representing the same q on opposite sides, never q in both words of one generator.

The independent degree bounds imply the stronger SHARED FORK inequality

    lambda_(uq,v)+lambda_(u,qv) <= delta(q) lambda_(u,v).        (6)

This follows by adding U_(uq)V_v<=s_A(q)U_uV_v and U_uV_(qv)<=s_B(q)U_uV_v. Unsupported children can be omitted because their formal coefficients are nonnegative; no unsupported variable is silently equated to zero.

More generally let s be the slack cap word of u and t the reversed slack cap word of v. For any nonnegative rational c,d,

    c lambda_(uq,v)+d lambda_(u,qv)
      <= max_(0<=a<=delta(q)) [c rho_s(a)+d rho_t(delta(q)-a)]
         lambda_(u,v),                                       (7)

where a ranges over integers. This uses the fixed actual allocation a=s_A(q), with s_B(q)<=delta(q)-a, pointwise before any averaging. Thus it combines the quantitative gap ratio with the common lower/upper budget without allocation branching.

There is an explicit finite facet list. Put D=delta(q), x_a=rho_s(a), y_a=rho_t(D-a). For a=0,...,D-1 let l_a=x_(a+1)-x_a and r_a=y_a-y_(a+1). Both are nonnegative, l_a decreases with a and r_a increases, because rho_s(a)=E[min(Z,a)]. The nontrivial edge inequalities are

    r_a lambda_left+l_a lambda_right
       <= (r_a x_a+l_a y_a) lambda_parent.                    (8)

The consecutive profile points maximize this linear expression: its change at allocation b is r_a l_b-l_a r_b, nonnegative for b<a and nonpositive for b>a. Zero coefficients are allowed; identically zero rows are discarded. Together with coordinate bounds these describe the downward convex hull of the allocation rectangles, a necessary relaxation only.

Example: both ancestor words have one letter of cap3, and q has total slack3. Then rho(a)=(0,1,5/3,2) for a=0,1,2,3. Besides individual bounds the three facets are

    lambda_left+3 lambda_right <=6 lambda_parent,
    lambda_left+lambda_right <=(8/3) lambda_parent,
    3 lambda_left+lambda_right <=6 lambda_parent.

The middle row improves the unrefined shared budget3 to8/3. Repeated labels are excluded by hypothesis, and reversal of the suffix word is essential.

## 7. Fixed-core quantitative application

For the already-studied core

    predecessor masks (0,28,0,0,12,12,13,125,125)
    projected order (2,3,0,4,5,6,1,7,8),

the finite profile universe has 270 lower profiles, 270 upper profiles and 56,379 admissible pairs. Every ambient seven-thin extension of this full root core satisfies the following common four-comparison inequality:

    2 N(0<4)+2 N(6<1)+4 N(4<5)+5 N(5<6)
        <= (13/20)*13 E.                                     (9)

Here E=e(P), and N counts full extensions, so at least one of these four probabilities is at most13/20. The exact verifier reconstructs all atoms and checks

    20[2N(0<4)+2N(6<1)+4N(4<5)+5N(5<6)]-169E <= -472

on every atom. Summing (1) transfers strict negativity to P. This improves the previous fixed-order 2/3 exclusion quantitatively; it does not settle the 3/5 target. The core was already excluded at2/3, so this is supporting evidence for the new instrument, not a newly closed morphology.

For context only, numerical LP discovery gives relaxed maximum-minimum forward probabilities about0.657130 with old degree caps,0.653485 with gap caps, and0.643550 with the complete atom convex hull. These decimals and their numerical optima are not needed for (9). The deliverable retains only the exact four-term certificate and an independently reconstructing standard-library verifier, not large solver arrays.

## 8. Open next obligation

The finite chain family converts arbitrary-size exterior uncertainty into a finite collection of actual small posets for every fixed core. A global seven-thin proof still requires common separating comparisons or fixed Gold strategies across the relevant finite hulls, followed by structural coverage of all possible cores/orders. The new reduction neither supplies that coverage nor identifies arbitrary convex-hull mixtures with actual exteriors. Testing an unresolved morphology where the old capped-generator cone survives is the next targeted use, rather than replaying already closed scans.

## 9. Reproduction and frozen payload

From this directory run:

    python verify_boundary_chain_compression.py
    python -O verify_boundary_chain_compression.py

Both supplied receipts agree. The standard-library verifier SHA-256 is aa67a3a159a7b69accc118c02a47091b65b6393a492d88dd12f2220de8409a99. It reconstructs the finite profile universe and every coefficient of the four-term inequality from the displayed core; no binary matrix, numerical optimizer, or other packet payload is required. The finite checks support the elementary all-size proofs rather than replace them.
