# Session 116: finite phase-completion optimization and the two-reset obstruction

## Status and scope

**Concurrent precedence and provenance.** During this independent Session 116 investigation, commit `bdbe68566e8` reached `main` with [`entry_cosets_resonant_echoes_20261030/COLLATZ.md`](../../entry_cosets_resonant_echoes_20261030/COLLATZ.md). Its sections 2–3 already publish the same universal suffix-domination theorem, ultrametric phase carrier, and cyclic-rotation corollary, together with a more general single-outlier endpoint theorem for `a^L b`. We credit that incoming packet as the preceding repository publication; no priority or novelty claim is made for those overlapping results. The present independent proof and executable checks are retained as corroboration. Theorem B below explicitly supplies the **arbitrary primitive-controller maximizing-entry grammar**, followed by exact primitive-root reduction for nonprimitive controllers; the incoming packet states that grammar for its single-outlier family. Our two concrete two-reset counterexamples and conjugated-window description are additional scoped illustrations, not a new discovery of the failure of empty-entry optimality (already exhibited in incoming section 4). Its larger controller-clock certificate and finite-time discrepancy theorem are not repeated or superseded here.

**Independently proved and checked here, with the overlapping results credited above, conditional only on the exact periodic-mean theorem already proved in Session 114:** arbitrary positive valuation controllers admit a finite exact all-entry optimization. The one-reset empty-entry theorem does **not** extend to two resets. The replacement is a finite symmetric phase-depth matrix, not a universal empty-entry rule or a forced analogy with algebra dimension.

No assertion concerns convergence of an integer Collatz orbit. Positive finite words are separately realizable; their realizing source integers need not be the same. The entry is retained throughout, and no interchange of infinite layer sums with ordinary limits is introduced here.

## 1. Definitions and inherited ordinary mean

For a positive word x of length d and clock A, write

    F_x(X)=(P_x X+B_x)/Q_x, P_x=3^d, Q_x=2^A,
    B_x=sum_i 3^(d-i) 2^(a_1+...+a_(i-1)).

Concatenation is chronological: F_xy=F_y F_x. Content is
c(x)=v_5(gcd(P_x-Q_x,B_x)). For nonempty w of length d, prefix v_i=w[:i], entry u, define

    R(u,w)=lim_(N->infinity) (1/N) sum_(n=0)^(N-1) sum_(i=1)^d c(u w^n v_i).

Session 114, `collatz/PERIODIC_CONTENT_RATES.md`, proves this ordinary limit and an O_(u,w)(log(N+2)) error. Put M=F_w, C_e=<M mod 5^e>, T_e=|C_e|. A phase contributes at layer e precisely when F_(v_i u) belongs to C_e, with density 1/T_e. Define its membership depth D_i(u) in {0,1,...,infinity}. The exact result is

    R(u,w)=sum_i f(D_i(u)),  f(D)=sum_(1<=e<=D) 1/T_e.

Writing t=ord_5(P_w/Q_w), the inherited formula is

    r=v_5(P_w^t-Q_w^t)                    if t>1,
    r=min(v_5(P_w-Q_w),v_5(B_w))          if t=1,
    T_e=t*5^max(0,e-r).

Thus f(D)=D/t for finite D<=r; f(D)=r/t+(1-5^(-(D-r)))/(4t) for finite D>r; and f(infinity)=(r+1/4)/t.

The earlier proof of the ordinary limit handles infinite-depth phases by positive-word centralizers and logarithmic height control. We use that theorem, rather than assuming Haar densities alone justify an ordinary mean.

## 2. Universal completion domination

Write w=v_j s_j, j=1,...,d, so s_d is empty. These d suffixes are the **completion entries**.

**Theorem A.** Every entry u is dominated, at every layer simultaneously, by some completion s_j. Consequently

    max_u R(u,w) = max_(1<=j<=d) R(s_j,w).

If all D_i(u) are finite, the rate domination is strict.

**Proof.** Choose j with D_j(u)=D=max_i D_i(u). Fix e<=D. In the finite affine group, put U=F_u and V_i=F_v_i. Because U V_j belongs to C_e,

    U V_i belongs to C_e iff V_j^(-1) V_i belongs to C_e.

But F_s_j=M V_j^(-1), so the right side is also equivalent to F_s_j V_i belonging to C_e. Hence u and s_j have exactly the same successful phases through layer D. Above finite D, u has no successful phases, while s_j always has phase j since F_s_j V_j=M. If D=infinity, the equality applies at every layer. All weights 1/T_e are positive, with summable tails, giving the assertions. QED.

This argument is a general nested-coset fact, not peculiar to one or two resets. It also rules out a putative supremum approached only by ever-longer finite entries: a completion actually attains the maximum.

### Rotation corollary

Let rot_j(w)=s_j v_j. Then

    R(s_j,w)=R(empty,rot_j(w)),
    max_u R(u,w)=max_j R(empty,rot_j(w)).

Indeed s_j w^N=(s_j v_j)^N s_j. The two marked sums traverse the same infinite positive word, with their summation intervals shifted by |s_j| letters. The finitely many initial boundary charges are constant; the finitely many terminal boundary charges are O(log(N+2)) by the inherited fixed-phase height estimate. Dividing by N proves the equality. Thus entry optimization is exactly **cyclic-rotation optimization**; the clock is rotation-invariant, so the same statement holds for original-clock-normalized rates. This is an ordinary-limit statement, not just a finite-layer-density assertion.

## 3. Complete entry classification and finite exact formula

Define the depth matrix

    Delta_(j,i)=max{e: V_j^(-1) V_i belongs to C_e}
               =D_i(s_j).

It is symmetric, has infinite diagonal, and satisfies the ultrametric inequality
Delta_(i,k)>=min(Delta_(i,j),Delta_(j,k)). At each depth e its equivalence classes are exactly the phase cosets modulo C_e. The exact optimization is the finite row-score formula

    R_max(w)=max_j sum_i f(Delta_(j,i)).

**Theorem B (primitive controller).** If w is primitive, every off-diagonal Delta is finite. Let J be the rows attaining the largest score. Then all maximizing entries, and only those, are

    u=s_j w^k,  j in J, k>=0.

**Proof.** Infinite depth means F_(v_i s_j) is in every C_e, hence commutes with M over the rationals. Positive-word affine maps are faithful, so the words v_i s_j and w commute. The commuting-word theorem makes v_i s_j a power of primitive w. Its length is d+i-j, between 1 and 2d-1; only i=j can give a multiple of d. This proves finite off-diagonal depths.

By Theorem A an optimizer must have an infinite phase, or a completion would strictly improve it. An infinite phase similarly forces v_j u=w^b. Since j>=1, b>=1, so u=s_j w^(b-1). Multiplication by a controller power does not change any cyclic-subgroup membership, so its score is exactly the completion score. QED.

For completeness, faithfulness follows from equality of affine slopes: unique factorization forces equal word lengths and clocks. With fixed length d, B=3^(d-1)+2^a_1 B_tail and B_tail odd, so a_1=v_2(B-3^(d-1)) when d>1; recurse, with the last exponent determined by the total clock. Thus equality of maps gives equality of words. The standard word-commutation statement follows by cancellation/Euclidean descent on word lengths.

**Nonprimitive reduction.** If w=z^k with primitive root z, the marked sums satisfy S(N;u,w)=S(kN;u,z) exactly. Hence R(u,w)=k R(u,z). Apply Theorem B to z. This gives a complete finite exact optimizer classification for every positive controller, not just the two-reset family.

### Explicit arithmetic evaluation

For W=(P,Q,B), H=(R,S,C)=coeff(v_i s_j), put

    gamma=v_5((P-Q) C-B(R-S)).

If t>1, first require (R/S mod5)^t=1, otherwise Delta=0. Set beta=v_5(R^t-S^t). Then

    Delta=min(gamma, beta if beta<r else infinity).

If t=1, set s=v_5(P-Q), h=v_5(B), and alpha=v_5(R-S) if s<=h, otherwise alpha=v_5(C). Then

    Delta=max(0,min(gamma-r, alpha if alpha<r else infinity)).

These are precisely the Session 114 membership-depth formulas. Infinite valuation means an exactly zero integer, not a numerical cutoff. For primitive words only the diagonal is infinite, so d(d-1)/2 finite exact computations determine the entire classification.

## 4. Sharp two-reset counterexamples

### Small example: w=(3,1,2)

This is 1^L a 1^M b with (L,a,M,b)=(0,3,1,2), primitive. Its affine coefficients are (27,64,49), and (t,r)=(4,2). In prefix order 1,2,3,

    Delta = [[infinity,1,0], [1,infinity,0], [0,0,infinity]].

Thus the empty-entry rate is 9/16, whereas the exact all-entry maximum is 13/16. All maximizers are

    u=(1,2) w^k or u=(2) w^k,  k>=0.

The clock is 6, so the corresponding maximal original-clock rate is 13/96. This is an entry-optimization obstruction, not a counterexample to a clock-budget bound.

### Both runs nonempty: w=(1,2,1,1,2)

Here (L,a,M,b)=(1,2,2,2); coefficients are (243,128,287), and (t,r)=(1,0). Every off-diagonal depth is zero except Delta_(2,4)=Delta_(4,2)=1. Hence

    R(empty,w)=1/4,     R_max(w)=9/20.

Exactly the phases j=2 and j=4 win. All maximizing entries are

    u=(1,1,2) w^k or u=(2) w^k,  k>=0.

The original-clock maximum is 9/140. The failure therefore persists with two genuine nonempty one-runs and identical reset letters.

## 5. What replaces the one-window modulus

For w=1^L a 1^M b, take A=F_1 and D=F_(1^L a). Use the phase representatives before the final reset: replacing V_d=M by identity preserves membership. The phase set consists of two windows

    {A^j: 0<=j<=L} union {A^k D: 0<=k<=M},

with multiplicities retained. At layer e a selected entry asks which representatives belong to one left coset of C_e. In the first window the successful exponents, if nonempty, form a residue class modulo the index associated to <A> intersect C_e. In the second window they form a residue class for <A> intersect D C_e D^(-1).

Proof: two first-window representatives differ by A^(j'-j); two second-window representatives have ratio D^(-1) A^(k'-k) D. Each membership-difference set is therefore a subgroup of the cyclic exponent group. Counting a window amounts to counting one arithmetic progression in a finite interval.

Unlike the one-reset situation, the second window carries a **conjugated intersection and a relative residue offset**. Its favorable phase cluster need not contain the empty-entry phase. Merely reusing the first-window modulus loses precisely this cross-window compatibility. The depth matrix supplies an exact repair: it records the nested coset partition, and completion-row scores select compatible offsets at every layer simultaneously. A separate digit formula for special two-reset parameters could simplify this computation, but is not required for completeness of the finite algorithm.

## 6. Checks and remaining boundary

`check_phase_completion.py` is standalone. Its arithmetic primitives reproduce Session 114's coefficient/order/depth formulas; independent finite-group orbit enumeration checks the phase-depth formula modulo 5,25,125. Explicit exceptions keep all gates active under `python -O`.

Receipt: 127 controllers of clock <=7, 128 entries of clock <=7 per controller; 121,797 exact gates plus 2,016 independent modular checks. Normal and optimized runs emit identical receipts. Tests cover the rotation identity, symmetry, the ultrametric inequality, exact primitive-root reduction, primitive separation, simultaneous layerwise domination, strict improvement for finite maximal depth, optimizer word completion, and both exact counterexamples.

This completes all-entry optimization **for each fixed controller**. It does not optimize over all controller words, prove a universal sharp clock-normalized bound for stalled controllers, or prove Collatz convergence. The useful dimensional lesson is specific: adding a reset changes a single residue-window problem into coupled phase windows; the invariant that survives is their nested coset partition, not empty-entry optimality.
