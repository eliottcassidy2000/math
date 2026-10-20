# Two actual adjacent swaps: labelled sparse forms and the norm barrier

Date: 20 October 2026. Status: **PROVED structural and algebraic statements; FINITE-EXACT verification.** This is a bounded extension, not a two-swap integer-cycle exclusion. The one-swap exclusion and the negative cycle classification are unchanged.

## 1. Scope and inherited coordinates

Let p>=2, q>p, gcd(p,q)=1, and q=pk+r, with 1<=r<p. Start from the zero-intercept lower mechanical word of length q and weight p. Exactly two operations are performed, each an actual cyclic adjacent unequal-bit swap in the word present at that step. The critical-clock restriction used later is only an exploratory slice: the one-swap critical-clock reduction has not been proved for two swaps and is not inherited here.

Inherit B, D=2^q-3^p, g=gcd(|D|,B), and the Eisenstein norm mechanism from [the sharp one-swap norm theorem](../../collatz_one_swap_sparse_norms_20261020/SHARP_ONE_SWAP_CONTENT_BOUND.md). The critical legality mechanism is from [the critical-clock refinement](../../collatz_one_swap_critical_gap_20261020/CRITICAL_CLOCK_SQRT7_GAP.md).

The closest proved mechanism is the one-swap sparse norm. The near miss is treating all formal edit signs as physical at k=1. The hostile is the known -17 cycle, whose [necklace swap distance is exactly three](../../collatz_one_swap_negative_exclusion_20261020/MINUS17_THREE_SWAP_BOUNDARY.md), not two. The least-used sidecar is the labelled intermediate word and its cut transport. The concept board is: labelled gaps; moving cuts; denominator-cleared sparse norms; critical legality; norm-versus-content loss; asymptotic torus averages.

Use the increasing lifted source positions

    s_j=jk+ceil(jr/p)-1,   j in Z,   s_(j+p)=s_j+q,
    a_j=p ceil(jr/p)-jr=(-jr mod p).

The representative j=0 is at -1, hence at q-1 in the ordinary word. Every actual swap moves one labelled one by +1 or -1. Labels never pass each other. Let d_j be each label's total lifted displacement after the two operations, with d_(j+p)=d_j.

## 2. Complete two-operation shape classification

There are exactly the following possibilities:

1. Undo: all d_j=0; the second operation reverses the first.
2. Same one twice: a single d_j is +2 or -2.
3. Distinct ones: precisely two distinct d_j are nonzero, each in {-1,+1}.

Put G_j=s_j-s_(j-1), the gap before label j. An individual move of label j is initially legal exactly when

    +1: G_(j+1)>=2;       -1: G_j>=2.

After any displacement vector e, the gap before j is G_j+e_j-e_(j-1). Applying this formula after the first move gives the exact second-step test, including crossing the original cut.

A same-one +2 path exists exactly when G_(j+1)>=3; a same-one -2 path exists exactly when G_j>=3. Both intermediate moves are then legal. Therefore neither kind can occur at k=1, where every source gap is 1 or 2. This assertion deliberately has p>=2; the one-one length-two circle has its own coincident-edge convention and is outside this statement.

For a proposed vector moving two distinct ones once, existence of a legal ordering is equivalent to

    G_j+d_j-d_(j-1)>=1   for every cyclic j.             (2.1)

Necessity is preservation of labelled cyclic order. For sufficiency, a move blocked initially can only be blocked by the other moving one. Condition (2.1) rules out moving the pair toward or through one another, so that dependency must be a same-direction adjacent pair. Move the leading one first. Both moves cannot depend on each other: that would require p=2 and both circular gaps equal to 1, contrary to q>p. If neither is blocked, either order works unless the two moves share a destination; the latter violates (2.1). Thus (2.1) is an existence criterion, while the successive gap tests retain the actual order of a specified path.

In particular, overlapping swaps must not be discarded. The paths 110 -> 101 -> 011 and 011 -> 101 -> 110 move two distinct ones, not one one twice. In the first, the leading right one must move first; in the second, the leading left one must move first.

## 3. A moving-anchor formula, including p=2 and cut crossings

An unchanged label need not exist when p=2. Instead anchor the final word at any label t, at its final position s_t+d_t. Define, for 0<=m<p,

    eta_m=ceil((t+m)r/p)-ceil(tr/p),
    h=max(0, max_m[d_t-d_(t+m)-eta_m]),
    P(X)=sum_m 3^(p-1-m) 2^(h+eta_m+d_(t+m)-d_t) X^m,
    F(X)=2^r X^p-3^p.

The coefficient exponents are nonnegative, so P is integral, nonzero, and has degree p-1. Cyclic label order implies that its physical one offsets are exactly mk+eta_m+d_(t+m)-d_t in [0,q). Hence

    P(2^k)=2^h B(final word rotated to its moved anchor),
    g=gcd(|F(2^k)|,P(2^k)).                            (3.1)

Set z^p=2, A=a_t, x=3z^(-r), and S(z)=sum_(a=0)^(p-1)z^a=1/(z-1). Direct substitution gives

    P(x)=2^(h-d_t) 3^(p-1) z^(-A)
          [S(z)+sum_j(2^(d_j)-1)z^(a_j)].              (3.2)

Thus the universal rational sparse expression is exactly

    Q(z)=1+(z-1) sum_j(2^(d_j)-1)z^(a_j).              (3.3)

The factor 2^(-d_t) is essential when the anchor moves. Omitting it is a cut-transport error, even though it disappears after removing 2-primary factors.

Put ell=max(0,-min_j d_j), and Qhat=2^ell Q in Z[Z]. With Rhat=Res(Z^p-2,Qhat), the exact norm identity is

    |Res(F,P)|=2^[r(p-1)+p(h-d_t-ell)-A]
               3^[p(p-1)] |Rhat|.                    (3.4)

The exponent displayed on 2 can be negative; (3.4) is an equality of nonzero rational numbers whose left side is an integer. It must not be presented as an integral factorization without checking that exponent. Its 2,3-free consequence is unaffected.

As in the predecessor, gcd(p,r)=1 makes F irreducible: the fields generated over the rationals by 3z^(-r) and by z coincide, using a Bezout relation pu+rv=1. Since deg P<p, its resultant with F is nonzero. Formula (3.2), and the nonvanishing monomial and z-1 factors, show Rhat is nonzero too. An integer resultant identity evaluated at 2^k gives

    g divides the positive 2,3-free part of |Rhat|.    (3.5)

Alternatively construct z modulo g using the same Bezout relation. Then z^p=2, z^r=3/2^k, z-1 is a unit, and Qhat(z)=0. This proof works modulo the whole g, retaining prime-power depth. The norm gate forgets the clock z^r=3/2^k and does not identify g with the norm.

## 4. Exact finite family of integral representatives

The following canonical H has the same positive 2,3-free resultant as Qhat. Every equality used is in Z[Z]/(Z^p-2), and multiplying by Z has norm of absolute value 2.

- Undo: H=1.
- One +2 at residue a: H=1+3(Z-1)Z^a.
- One -2 at residue a: H=2Z^(p-a)-3Z+3. Here 4Q=Z^a H.
- Two +1 at distinct residues a,b: H=1+(Z-1)(Z^a+Z^b).
- Two -1 at residues a<b: H=Z^(p-a)-(Z-1)(1+Z^(b-a)). Here 2Q=Z^a H.
- A +1 at a and a -1 at b, a!=b:

      H=Z^(p-b)[1+(Z-1)Z^a]-Z+1.

  Here 2Q=Z^b H. Reduce powers >=p using Z^p=2 if desired.

These are the complete six shape classes, not arbitrary independent signed one-swap factors. In particular, taking a norm of the un-cleared rational Q and treating it as an integer modulus is invalid. The canonical minus and mixed forms remove the known dyadic monomial before numerical ranking.

## 5. Complete critical-band residue legality

Now let q=ceil(p log_2 3) and gcd(p,q)=1. There is no p=2 instance. For p>=3, k=1 and r>p/2. The latter follows from log_2 3>3/2, since 9>8. The individual initial legality conditions are

    + at a is legal iff a<r,
    - at b is legal iff b>=p-r.

There are no same-one +/-2 paths. For two distinct labels the following is an exact existence classification:

- ++: either both residues are <r, or the unordered pair is {u,u+r} with 0<=u<p-r. In the exceptional pair, residue u moves first and residue u+r follows into its vacated position.
- --: either both residues are >=p-r, or the unordered pair is {u,u+r} with 0<=u<p-r. In the exceptional pair, residue u+r moves first and residue u follows.
- Mixed, + at a and - at b: a<r, b>=p-r, a!=b, and b!=a+p-r. Both moves are initially legal. The excluded equation is precisely the consecutive source-one gap of length 2 in which the two ones would move into the same zero.

Proof: an initially blocked same-direction move can only follow its adjacent source one. The residue relation between consecutive labels is a_(j+1)=a_j-r mod p. The inequalities above follow immediately; r>p/2 makes the other move in each exceptional pair legal. Opposite moves cannot unlock one another. Two opposite initially legal moves are obstructed only by their shared destination, giving the excluded equation. This is the critical specialization of (2.1), not a replacement for the intermediate test when an ordered path is specified.

## 6. Exact barriers to a naive sub-three norm bound

### Small mixed barrier

At p=4,q=7,r=3 the actual path, with zero-based edges 0 then 3, is

    0101011 -> 1001011 -> 1000111.

It moves residue b=1 left and residue a=2 right. Its canonical polynomial is

    H=Z^3+Z^6-Z^5-Z+1
     =Z^3+2Z^2-3Z+1 mod (Z^4-2).

Exactly,

    Res(Z^4-2,H)=343=7^3>3^4=81.

This norm is already 2,3-free. The final word has B=331,D=47,g=1. Thus it is a norm obstruction, not an integer-cycle counterexample.

### Overlapping-chain barrier

At every critical coprime shape with p>=3, the pair

    a=p-1,   b=p-1-r

is a legal ++ exceptional chain: the b label moves first. Set n=p-r. Its polynomial is

    Q_p=3-Z^(p-1)+Z^n-Z^(n-1).                         (6.1)

At p=28,q=45,r=17, the source ones of labels 5 and 6 occupy positions 8 and 9; position 10 is zero. Move 9->10, then 8->9. The residues are 27 and 10. Exactly,

    |Res(Z^28-2,Q_28)|=23221006835657,
    3^28=22876792454961,
    difference=344214380696.

Again the norm is already 2,3-free, while D=12307579633871 and actual g=1. The receipt also gives independent examples at p=29 and p=31.

Consequently no estimate M<c^p with one c<3 can hold for every critical actual two-swap type using these exact norm moduli. A finite counterexample does not rule out an eventual estimate, a bound with a prefactor, or a smaller content estimate retaining the clock and common modular root. In particular this does not refute two-swap cycle exclusion.

## 7. A rigorous raw-norm asymptotic, with the stripping issue left open

For the legal family (6.1), along any unbounded sequence of critical coprime shapes,

    |Res(Z^p-2,Q_p)|^(1/p) -> 3.                       (7.1)

Such an unbounded sequence exists by taking p prime. This theorem concerns the raw norm, not its 3-free part.

Proof. Modulo Z^p-2, multiply (6.1) by Z to get

    H_p(Z)=3Z-2+(Z-1)Z^n,   n=p-r.

Its norm differs from Q_p's by absolute factor 2. Put rho_p=2^(1/p), beta_p=n/p, and beta=2-log_2 3. Then beta_p->beta, which is irrational, and rho_p^n->4/3.

At a conjugate Z=rho_p exp(i theta), write phi=n theta. The logarithmic integrand converges uniformly to

    f(theta,phi)=log|3 exp(i theta)-2
                    +(4/3)(exp(i theta)-1)exp(i phi)|.

There is uniform nonvanishing: for rho>=1,

    max_theta |rho exp(i theta)-1|/|3rho exp(i theta)-2|
       =(rho+1)/(3rho+2)<=2/5,

and rho_p^n<=2. Thus the ratio of the perturbation to the linear term is <=4/5. At the limit the sharper bound is 8/15.

The p points (2pi j/p,2pi n j/p) equidistribute on the two-dimensional torus as p grows through this sequence. Indeed, the empirical mean of character exp(i(h theta+l phi)) is 1 exactly when p divides h+ln, and is 0 otherwise. For each fixed nonzero integer pair (h,l), this divisibility eventually fails: if l!=0 it would force h/p+l beta_p to approach an integer while its limit l beta is irrational; if l=0 the claim is immediate. The character criterion, or elementary trigonometric-polynomial approximation, proves equidistribution.

Averaging f first over phi gives log|3 exp(i theta)-2|, because its constant term dominates the rotating term. Averaging this over theta gives log 3. Uniform convergence and equidistribution therefore give p^(-1)log|Norm H_p|->log 3, and the factor 2 between norms disappears in that limit. This proves (7.1).

For p>=5, n>=2 and Q_p has odd constant term 3, so its raw norm is odd (reduce the resultant modulo 2). However 3-adic valuations can be nonzero: p=13,r=8 gives raw norm 1052433=3^5*61*71, hence v_3 equal to 5. No o(p) estimate for v_3 along all critical shapes, nor an infinite 3-free subsequence with a controlled norm limit, is proved here. Thus (7.1) obstructs an eventual raw-norm sub-three bound but must not be promoted to the corresponding eventual stripped-norm or content claim.

## 8. Exact checks and stopping point

`verify_two_edit_symbolic.py` independently enumerates every ordered actual cyclic two-swap path for 2<=p<=10, p<q<=70, gcd(p,q)=1. It checks the labels, gap order, same-one gap criterion, moving-anchor carry polynomial, original-cut cases, and critical residue classification. It compares ordinary carry, cleared sparse, and canonical sparse resultants, including the exact power-of-two identity (3.4), on 486 labelled parameter types. Its critical converse check examines all two-distinct-label signed vectors through p=60 and verifies that the displayed residue classifier is equivalent to actual reachability.

The receipt records 48,796 actual paths and 82,940 critical signed-vector candidates, of which 29,732 are legal. Exact norm barriers are evaluated independently with integer SymPy resultants. Reproduce from this directory:

    python3 verify_two_edit_symbolic.py > receipt.replay.json
    python3 -O verify_two_edit_symbolic.py > receipt.optimized.json
    cmp receipt.json receipt.replay.json
    cmp receipt.json receipt.optimized.json

The verifier uses explicit check functions, not Python assert statements. It imports no repository mathematical implementation.

Stopping point: the sparse reduction and all critical legal shapes are now explicit. The cheap one-swap norm strategy does not extend unchanged: even its correctly stripped exact modulus can exceed 3^p on legal two-swap words, while the actual content remains 1 in the displayed hostiles. An exclusion proof would need to keep the clock/common-root sidecar, exploit cancellations beyond this modulus, or use a different invariant. No further unbounded exclusion is asserted.
