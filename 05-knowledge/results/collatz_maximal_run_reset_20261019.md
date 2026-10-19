# First-descent-safe expanding Collatz runs and exact source-sensitive resets

**Status, 2026-10-19:** PROVED a maximal-run decoder and three infinite
parameterized families of actual descent receipts. Two families have exact
first-stopping-time formulas. Universal entry into these families remains
OPEN. No finite residue census or uncomputed Collatz tail is used.

## 1. Calibration: the first expanding return was not an unresolved source

The [switching-height packet](collatz_switching_height_certificate_20261019.md)
exhibits an expanding adjacent return

    J(x)=(4782969x+10623730)/4194304,
    x=1817598 (mod 4194304).

Its lower source is even and already descends at its first step. Its odd
upper source also has an earlier descending prefix: at step 20 its branch
has multiplier 3^12/2^20 and fixed point 106493/103427<2. Every actual upper
source is therefore already smaller by that prefix. J remains a valid
obstruction to coarse *return-block* height ranks, but its long-run reset
families would not provide new first-descent coverage.

The present packet corrects that weakness by shifting to a phase whose
entire expanding word survives ordinary first-descent screening. It is a
return of the relation (3u+2,u), not an adjacent-pair return. The adjacent
odd-lower-source structural question is kept separate in section 7.

## 2. An expanding odd-source phase with no descending prefix

Put

    P=4782969=3^14, Q=4194304=2^22, C=5311865,
    d=117733, B=1062373, c=B/d.

The word obtained by moving J's initial even step to its end is

    w_R = 1111111110011011010000.

Its actual branch is

    R(u)=(Pu+C)/Q,
    u=908799 (mod Q),   p_R=-c,   lambda=P/Q>1.     (2.1)

Indeed C=5B and P-Q=5d. The branch coefficients and residue follow by the
ordinary parity-word recurrence. At the anchor 908799 the endpoint is
1036349. The companion source 3u+2 satisfies

    T^22(3u+2)=3R(u)+2

on this same cylinder. Its parity word is 1111111110110000000111, with odd count 14
and affine carry 14758265. Its source guard is checked separately in the
companion script. The descent results below concern the
odd driver u; they make no first-descent-free claim for its companion.

**No-prefix-descent lemma.** Every nonempty prefix of w_R strictly increases
its positive starting source. The same is true of every nonempty prefix of
any lawful repetition w_R^L.

**Proof.** A prefix of length t with a odd bits has map
(3^a u+C_t)/2^t, with C_t>=0. After the initial run of ones, the minima of
3^a/2^t occur at the ends of even runs. Their (t,a) pairs are

    (11,9), (14,11), (17,13), (22,14),

and all four satisfy 3^a>2^t. Initial all-odd prefixes also have multiplier
above one. Thus every prefix has value strictly greater than its block's
positive starting source. Each complete block increases that source,
so induction handles every prefix of R^L. ∎

In particular this is a genuinely first-descent-free expanding family of
odd sources. Its later reset is not made redundant by an earlier dip.

## 3. Compress the maximal run, retaining the actual exit shell

For a positive integer x in R's cylinder, decompose

    dx+B = 2^(22L+s) m,
    m positive odd, L>=1, 0<=s<22.                 (3.1)

The exact one-word valuation theorem gives

    L=floor(v_2(dx+B)/22).

This is the maximum number of consecutive complete R blocks. Their exit is

    y=R^L(x)=[P^L 2^s m-B]/d.                     (3.2)

It is a positive integer and has residual R rank exactly s. The identity
dR(x)+B=(P/Q)(dx+B) proves the formula. Integrality is also visible from
P=Q modulo d. The actual parity-word budget, rather than formal affine
integrality alone, establishes lawfulness.

The pair (y,s) remains unchanged as one moves through the run; the remaining
L drops by one. No intermediate orbit search is needed. All integers used
have bit length O(log(x+1)), since 22L<=log_2(dx+B).

Put epsilon=log_2(lambda)/22=0.0086125004589.... The run endpoint satisfies

    y+c=lambda^L(x+c)
       <=(dx+B)^epsilon(x+c)
        =d^epsilon(x+c)^(1+epsilon).               (3.3)

This bounds the endpoint after the whole maximal expanding run. It does
not claim that every intermediate Terras value is below that endpoint.

## 4. Source-derived reset capacity: a concrete exit grammar

Use the three contracting branches

    E(z)=z/2,                 word 0,         center 0;
    F(z)=(3z+1)/4,            word 10,        center 1;
    H(z)=(81z+101)/512,       word 110101000, center 101/431.

Their exact parity cylinders are respectively 0 modulo 2, 1 modulo 4,
and 43 modulo 512. The center separations from p_R are 0,1,2 in 2-adic
valuation. Formula (3.2) gives this source-specific grammar:

- s=0: y is even. Exactly

      k=v_2(y)=v_2(P^L m-1062373)                  (4.1)

  consecutive E steps are available.
- s=1: y=1 modulo 4. Its exact F rank is

      v_2(y-1)=1+v_2(P^L m-590053).               (4.2)

  Exactly floor(v_2(y-1)/2) complete F blocks are available.
- s=2: y=3 modulo 8. Its exact H rank is

      v_2(y-101/431)=2+v_2(431P^L m-117443449).    (4.3)

  Exactly floor(v_2(y-101/431)/9) complete H blocks are available. Having
  residual shell 2 alone does not guarantee a full H block.
- s>=3: y=7 modulo 8, so none of E,F,H is available as its next word.

These formulas come from clearing the odd denominators. In particular
431B+101d=469773796=4·117443449. None of the displayed valuations is infinite
at these actual exits: y>x>=908799, and the reset centers are smaller.

For a reset W(z)=q+mu(z-q), a lawful j-block prefix descends below x exactly
when

    mu^j(y-q)<x-q.                                (4.4)

All sides are positive and this can be checked by integer multiplication.
The required j exists as a real-affine comparison because mu<1; the grammar
above determines whether that many blocks are actually lawful. A decoder
returns the actual word and smaller endpoint only when both tests pass.
Otherwise it reports an uncertified exit, not a counterexample.

There are useful automatic cases. If s=0 and L<=5, even one halving is
already sufficient. If only one halving is available and L>=6, it is
insufficient. This follows from

    9/8<lambda<8/7,   c/908799<1/100000,
    (100001/100000)(8/7)^5<2<(9/8)^6.

Likewise, when s=1 a single F block suffices for L<=2; for L>=3 a single F
block cannot descend, because its multiplier already exceeds one. The
checker verifies the exact endpoint inequality at x=908799,L=2 and the
opposite multiplier inequality at L=3.

Whenever the even or F test succeeds, the decoder also supplies the exact
first stopping time. All R prefixes stay above x. Even resets decrease
monotonically; within F, the odd step increases its current input and the
complete two-step block contracts toward 1. Therefore the first descent
occurs at the first sufficient complete reset block. For H the decoder
supplies a receipt and an upper bound on first stopping time, not an
unproved assertion about its intermediate prefixes.

## 5. All-L optimal-reset cylinders with a uniformly bounded exceptional head

Define the expanding-run carry

    C_L=C(P^L-Q^L)/(P-Q).

Fix W in {E,F,H}, writing W(z)=(az+b)/q_2 where q_2 is a power of two.
Use symbols p_W for its center and ell_W for its word length. Let

    j_0(L,W)=min {j>=1 : a^j P^L < q_2^j Q^L}.    (5.1)

This is the *smallest possible* number of W blocks for descent: with fewer
blocks the multiplier is at least one and the positive affine offset
precludes any decrease. For j=j_0, put

    B_j=b(q_2^j-a^j)/(q_2-a),
    A'=a^j P^L,    M=q_2^j Q^L,
    C'=a^j C_L+Q^L B_j,    Delta=M-A'>0,
    rho=-C'(A')^(-1) mod M,  0<rho<M.              (5.2)

The full lawful cylinder is x=rho modulo M. Every positive source in it
has *maximal* initial R-run length L. Its residual rank is respectively
s=0,1,2 for W=E,F,H. Indeed a lawful W block puts its source closer to p_W
than the finite center separation from p_R, so the exit R rank is exactly
that separation. This also follows directly from the exit grammar.

**Optimal-reset tail theorem.** On the cylinder (5.2), the word R^L W^j
is a smaller-source receipt exactly when

    Delta x>C'.                                  (5.3)

Writing x=rho+tM, t>=0, the first successful lift is explicitly

    t_min=max(0, floor((C'-Delta rho)/(Delta M))+1). (5.4)

Every lift t>=t_min is successful. Uniformly in L,

    t_min<=10 for E,
    t_min<=11 for F,
    t_min<=10 for H.                              (5.5)

For E and F the first stopping time is exactly

    22L+ell_W j_0(L,W).                           (5.6)

For H this clock is a valid descent receipt, possibly later than the first
descending prefix.

**Proof.** The branch-word numerator gives the complete cylinder and the
endpoint (A'x+C')/M. This proves (5.3)-(5.4), including strictness at an
integer threshold. To bound the exceptional head, express the carry using
the centers:

    C'=cA'+p_W M-(c+p_W)a^j Q^L
       <(c+p_W)M,

because A'<M and the last subtracted term is positive. As Delta is a
positive integer,

    C'/Delta < (c+p_W)M.

Consequently every t>=ceil(c+p_W) satisfies (5.3). These ceilings are
10,11,10 for E,F,H. This argument remains valid even when the powers of
two and three in Delta are unusually close; no unproved gap estimate is
used. The prefix lemma handles R^L. Before the final E or F block, every
complete reset prefix has multiplier at least one by minimality of j_0,
and positive offset. E has no intermediate step; F's intermediate odd
step increases its input. Hence no earlier prefix descends, proving (5.6). ∎

For concrete controls, 5103103 follows R to 5819318 and then one even
step to 2909659; its first descent is exactly at step 23. The source 908799
follows R to 1036349 and then F to 777262; its first descent is exactly at
step 24. In contrast, the H receipt from 2005786111 ends at 361857301 at
step 31, while its actual first descent is already at step 27.

These are infinite families with arbitrarily long genuinely rising odd-source
excursions. Their source domains are explicit progressions with at most a
constant number of initial lifts removed. Distinct (L,W) families lie at
distinct initial valuation levels 22L+s, and are therefore disjoint.
The theorem is sufficient coverage of these exact domains, not a statement
that arbitrary sources enter them.

## 6. Whole cylinders are available with at most one extra reset block

If removing the finite initial head is undesirable, an all-positive-lift
version is also explicit. Let

    r_L=-C_L(P^L)^(-1) mod Q^L,
    y_L=(P^L r_L+C_L)/Q^L,

and choose the least j>=1 with

    mu^j(y_L-p_W)<r_L-p_W.                        (6.1)

Use that j in the cylinder formula (5.2). Every positive lift now descends,
because the composite affine endpoint divided by x decreases with x>0,
and every positive R^L source is at least r_L. Here the comparison at r_L
is merely an affine bound; no reset at r_L is assumed lawful.

This j is either j_0 or j_0+1. Indeed r_L>=908799>3c+4p_W, so

    (y_L-p_W)/(r_L-p_W) < (4/3)lambda^L.

Every reset multiplier is at most 3/4. One additional block beyond (5.1)
therefore forces (6.1). Thus even the full-cylinder version uses at most
one extra reset block. The asymptotic counts per expanding block are

    E: log_2(lambda)=0.1894750101...,
    F: log(lambda)/log(4/3)=0.4565...,
    H: log(lambda)/log(512/81)=0.0712271903....

All thresholds in the checker use exact integers or rationals, not these
decimal approximations.

## 7. Remaining source coverage and the adjacent-return question

The decoder can fail: s>=3 gives y=7 modulo 8 and no next E,F,H word;
s=2 may lack H's deeper digits; and either of the first two shells can
have insufficient reset capacity. These are explicit source conditions,
not cases declared solved by their having finite-state descriptions.
Repeated successful decoder applications strictly reduce the positive
integer, but future decoder entry has not been proved.

A separate, concrete structural question remains: must every *oriented
adjacent return with odd lower source* either contract, or have an earlier
descending prefix on that odd source? J did not test this question because
its lower source is even. The R family is first-descent-free but returns
the different relation (3u+2,u), so it does not refute the adjacent-return
statement. This is a next conjecture to investigate, not an inherited
fact or a theorem of the present packet.

## 8. Exact reproduction

The [checker](../experiments/collatz_maximal_run_reset_20261019.py) performs **70,732 explicit
checks**. It reconstructs R and its
actual phase companion, checks every prefix multiplier, proves the displayed
integer comparisons, and replays the source decoder and the complete
parameterized reset cylinders. It includes exact first-stopping-time
controls for E/F, ordinary receipts for H, and explicit uncertified exits.
The [normal receipt](collatz_maximal_run_reset_20261019.json) and
[optimized receipt](collatz_maximal_run_reset_20261019.optimized.json) must
agree byte-for-byte.

```sh
python3 05-knowledge/experiments/collatz_maximal_run_reset_20261019.py > /tmp/maximal-run-reset.json
python3 -O 05-knowledge/experiments/collatz_maximal_run_reset_20261019.py > /tmp/maximal-run-reset-O.json
cmp /tmp/maximal-run-reset.json /tmp/maximal-run-reset-O.json
cmp /tmp/maximal-run-reset.json 05-knowledge/results/collatz_maximal_run_reset_20261019.json
```
