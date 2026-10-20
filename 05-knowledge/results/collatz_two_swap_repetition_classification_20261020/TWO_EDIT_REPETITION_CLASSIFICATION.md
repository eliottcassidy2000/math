# Two adjacent edits of repeated words: exact cancellation and all nonprimitive mechanical shapes

20 October 2026. Status: **PROVED elementary repetition-factor classification; PROVED mechanical consequence using the audited one-swap theorems; FINITE-EXACT hostile controls.** This is an actual-content result, not a resultant estimate. The primitive two-edit frontier remains open.

## 1. Exact theorem

Let v be any binary word of length q0 and weight p0, with 0<p0<q0. No balance or primitivity of v is required. Put

    c=v^d,    d>=2,
    A=2^q0,    C=3^p0,
    D0=A-C,
    D=A^d-C^d=H D0,
    H=sum_(j=0)^(d-1) A^(d-1-j) C^j.

Use chronological carry

    B(w)=sum_(i:w_i=1) 2^i 3^(number of ones after i).

Obtain w by exactly two actual cyclic adjacent swaps of unequal bits, at zero-based edge starts i then j in the current word. Let q=d q0. The following is an exact classification:

**Theorem R2.** H divides B(w) if and only if either

1. the second swap undoes the first, equivalently j=i, so w=c; or
2. d=2 and the two edges are one base period apart, j-i=q0 modulo q.

In the second case the two swaps are disjoint corresponding edits in the two copies. Their directions agree, and

    w=(v')^2,

where v' is obtained by one actual cyclic adjacent swap of v. The statement uses the displayed repetition representation; v need not have least period q0.

Since integrality requires |D| dividing B(w), any non-undo two-edit path with d>=3 is nonintegral for either clock sign. When d=2, a non-undo integral output can occur only in the corresponding-copy case, and then integrality is exactly that of the one-edited base v'.

**Mechanical corollary.** Let c be the zero-intercept mechanical word of any noncoprime shape 0<p<q, and w any two-actual-swap output. Then w is an integer-cycle word if and only if

- the reduced shape (p0,q0)=(p/d,q/d), d=gcd(p,q), is (1,2) or (2,3); and
- the path is undo, or d=2 and the edges are q0 apart.

Every such integer output has necklace distance zero from c. In particular there are no integral outputs of exact cyclic-adjacent-swap necklace distance two at any nonprimitive mechanical shape, at any length, for either sign. This does not assert the same for primitive mechanical sources or arbitrary nonmechanical sources.

## 2. Inheritance, cuts, and the exact factor being retained

The closest mechanism is the [one-edit repetition obstruction](../collatz_one_swap_negative_exclusion_20261020/NEGATIVE_ONE_SWAP_EXCLUSION.md), Section 7: one smooth carry perturbation is invertible modulo H. Two edits can cancel, so a second unit cannot be discarded. The relevant hostile is a squared word obtained by performing the same edit in both copies. The new coordinate is the shorter cyclic separation between the two edited edges, together with its intervening one-count. The concept board is: repetition factor; smooth perturbations; cyclic cuts; shorter-arc binomials; exact equality cases; inherited base integrality.

The repetition identities are

    B(v^d)=H B(v),    D=H D0.

Here H>=A+C>=7 for d=2, and H>=7 for every d>=2. Also gcd(H,6)=1: reduction modulo 2 leaves C^(d-1), which is odd, and reduction modulo 3 leaves A^(d-1), which is nonzero.

For a one-place left rotation rho and first bit epsilon,

    2 B(rho(w))=3^epsilon B(w)+epsilon D.

Consequently gcd(H,B(w)) is invariant under rotation, because H divides D and 2,3 are units modulo H. This permits a cut chosen around the edited edges without changing the exact H-primary content, including all prime-power depths.

An internal swap at edge h contributes

    Delta B=epsilon_h 2^h 3^s,

where epsilon_h=+1 for 10->01 and -1 for 01->10, and s is the number of ones after the two swapped positions. The full word contains q>=4 positions, since the base is nonconstant and d>=2. Thus the exceptional coincident-edge convention on a two-position circle never occurs here.

## 3. Undo and overlapping affected positions

If i=j, the second swap is the inverse of the first, so w=c and H divides B(w).

Suppose instead the two edges are distinct but share an affected position. On a circle of length at least four, their union consists of exactly three consecutive positions. Choose a cut outside this block. Inspection of the actual intermediate word leaves only the following non-undo possibilities and their reversals:

    100 -> 010 -> 001,
    110 -> 101 -> 011.

Their three-letter carry differences are respectively

    B(001)-B(100)=4-1=3,
    B(011)-B(110)=10-5=5.

The unchanged outside context multiplies these by a power of 2 times a power of 3. Since the unedited carry is divisible by H, the exact conclusions are

    same-one displacement two: gcd(H,B(w))=1,
    two-one overlapping chain: gcd(H,B(w))=gcd(H,5).

Both are strictly smaller than H>=7. This covers either temporal edge order, both directions, and a three-position block crossing the original cut. The factor 5 must not be called a 2,3-unit; it can survive partially in H.

## 4. Disjoint edits: the two-unit sum becomes a shorter-arc binomial

Now the two affected pairs are disjoint. Their swaps commute, and each is already legal in the source word c: neither touches either bit used by the other.

Choose the shorter cyclic distance between their edge starts, and rotate so the first selected edge is 0 and the other is ell. Then

    2<=ell<=q/2,

and both pairs are internal to the displayed word. Write c again for this rotated repeated source. Define

    m=sum_(t=2)^(ell+1) c_t,
    s=sum_(t=ell+2)^(q-1) c_t,
    epsilon_0=c_0-c_1,
    epsilon_ell=c_ell-c_(ell+1).

The count m is over an interval of exactly ell positions. It is at least one, because it includes the second unequal pair. The two source perturbations are

    epsilon_0 3^(m+s),
    epsilon_ell 2^ell 3^s.

The other swap does not change either applicable suffix count: a disjoint pair is unchanged, and the earlier swap preserves the number of ones in its pair. Thus the exact total perturbation is

    B(w)-B(c)=3^s U,
    U=epsilon_0 3^m+epsilon_ell 2^ell.                 (4.1)

Since H divides B(c) and is coprime to 3, this gives the full arithmetic statement

    gcd(H,B(w))=gcd(H,U).                             (4.2)

In particular this is not merely a comparison of prime supports. Also U is nonzero: if the signs differ, equality of 2^ell and 3^m is impossible for ell,m>=1; if they agree, the sum is nonzero.

### At least three repetitions

For d>=3,

    ell<=d q0/2<(d-1)q0.

Every interval of length (d-1)q0 in a periodic repetition contains exactly (d-1)p0 ones. The interval counted by m can be extended to such an interval, so

    m<=(d-1)p0.

No mechanical balance is used. Therefore

    0<|U|<=2^ell+3^m
            <A^(d-1)+C^(d-1)
            <=H.                                     (4.3)

The last inequality follows because both endpoint terms occur in the positive sum defining H. Hence H cannot divide U or B(w).

### Exactly two repetitions

For d=2, the shorter arc has ell<=q0. The same interval-containment argument gives m<=p0. Hence

    |U|<=2^ell+3^m<=A+C=H.                            (4.4)

Since U is nonzero, divisibility by H requires |U|=H. Equality holds exactly when

    ell=q0,    m=p0,    epsilon_0=epsilon_ell.

Conversely, at ell=q0 periodicity gives m=p0 and the two source pairs are identical, so their directions agree and U=+H or -H. Thus H divides B(w) exactly when the edges are q0 apart. Editing corresponding positions preserves the period q0 and gives the square of the one-edited base. Together with Section 3, this proves Theorem R2.

## 5. Integrality and the mechanical specialization

For undo, B(w)/D=B(v)/D0. For the exceptional doubled edit, w=(v')^2 and

    B(w)/D=B(v')/D0.

Every other path loses some factor of H and is nonintegral. More quantitatively, if h=gcd(H,B(w)), then

    gcd(|D|,B(w))<=|D0| h,
    reduced denominator(w)>=H/h>1.

The exact h is supplied by Section 3 or (4.2).

For the mechanical specialization, take d=gcd(p,q). Then c=v^d and v is the reduced, coprime-shape zero-intercept mechanical word. A primitive mechanical word has unit carry content. One short proof is to compare its intercept-0 and intercept-1 mechanical representatives: they are rotations, and the floor-function difference changes exactly one adjacent 01 to 10. Rotation preserves the carry gcd with D0, while the difference of these two carries is a 2,3-unit. Thus gcd(|D0|,B(v))=1. This is also [THM-4716, balanced smooth ladders and freeness](../../../01-canon/theorems/THM-4716-balanced-collatz-cycles-are-three-smooth-ladders-and-are-free.md); only its elementary coprime-shape arithmetic assertion is needed here.

Therefore an undo output is integral exactly when |D0|=1. Under 0<p0<q0 the only possibilities are (p0,q0)=(1,2),(2,3):

- If 2^q0-3^p0=1, reduction modulo 3 makes q0 even. The coprime odd factors 2^(q0/2)-1 and 2^(q0/2)+1 are powers of 3 differing by 2, hence are 1 and 3. This forces q0=2,p0=1.
- If 3^p0-2^q0=1, then p0>=2,q0>=3, and reduction modulo 8 makes p0 even. The factors 3^(p0/2)-1 and 3^(p0/2)+1 are powers of 2 with gcd 2 and difference 2, hence are 2 and 4. This forces p0=2,q0=3.

For the non-undo doubled edit, apply the independently audited [positive one-swap exclusion](../collatz_one_swap_exclusion_20261020/ONE_SWAP_EXCLUSION.md) and [negative one-swap classification](../collatz_one_swap_negative_exclusion_20261020/NEGATIVE_ONE_SWAP_EXCLUSION.md) to v'. Their small unit-clock boundaries and p0=1 calculation give the same two reduced shapes. At those shapes every one-swap output is a rotation of v. Hence (v')^2 is a rotation of v^2 and has necklace distance zero.

This proves the full mechanical corollary without using any two-edit critical-clock reduction, norm-growth estimate, or finite weight cutoff. It does inherit the analytic dependency of the already proved primitive positive one-swap theorem for the exceptional d=2 base reduction.

## 6. Hostile examples that delimit the theorem

### Partial 5-adic retention in an overlapping chain

Take v=0011, d=2, so H=25 and D=175. The ordered edges 3 then 2 give

    00110011 -> 00101011 -> 00011011.

The final carry is 680, so

    gcd(H,B)=5,    gcd(|D|,B)=5.

Thus overlap need not destroy every prime of H. It destroys enough prime-power depth to prevent integrality. Replacing the exact gcd(H,5) formula by coprimality would be false.

### The doubled-edit exception can restore a genuine integer cycle

The theorem for arbitrary repeated words must retain its exceptional case. Take

    v=11101111000,
    c=v^2,

and swap corresponding edges 3 then 14. The output is

    w=(11110111000)^2.

This is two traversals of the genuine -17 cycle. Exactly,

    H=4235,
    D=-588665,
    B(w)=10007305=-17 D.

The output has exact necklace distance two from c: it is neither a rotation of c nor of a one-swap output. The base v is not mechanical. This explicit example rules out any inference that Theorem R2 excludes all nonprimitive two-edit integer outputs without the mechanical-source hypothesis.

The known distance-three of the -17 necklace from its mechanical source is unaffected; the source in this hostile example is a different, nonmechanical repeated word.

## 7. Verification and audit

`verify_repeated_two_edits.py` checks every nonconstant binary base of length 2 through 8, every repetition count 2 through 5, and every ordered actual cyclic two-swap path. It tests the exact gcd in every overlap and disjoint case, the cut transport, interval bounds, mirrored equality, and base-square identities. It also tests the mechanical corollary wherever the base is a reduced mechanical word. The completed receipt records 494 bases, 1,976 repeated sources, 400,784 ordered paths, and 2,282,595 explicit checks. Exactly 25,088 undo paths and 3,584 mirrored doubled-edit paths retain all of H. Normal and optimized receipts are identical. Checks are explicit runtime conditions and remain active in optimized mode.

An independent checker supplied by the parent reviewer uses positional carry sums only, without importing the perturbation classification. It covers all 240 nonconstant binary bases of length 2 through 7, repetitions 2,3,4, and 81592 paths. Exactly 8448 paths retain all of H: 6912 undo paths and 1536 mirrored doubled edits. Every mirrored output is directly checked to be a square. Normal and optimized independent receipts agree.

A [fresh independent audit](audit/TWO_EDIT_REPETITION_AUDIT.md) passed without corrections. Its separate mask-based verifier records 3,013,615 explicit checks across 437,288 ordered paths, including repetitions 6 and 7, all nine-bit squares, nonminimal bases, cut crossings, all q=4 paths, and primitive mechanical shapes through q=200. Its normal and optimized receipts match exactly.

The parent-checker files are copied unchanged into `audit/`; the original reviewer files remain untouched. All claims concern these explicit algebraic families, not arbitrary Collatz trajectories or historical priority.

Reproduce from this directory:

    python3 verify_repeated_two_edits.py > receipt.replay.json
    python3 -O verify_repeated_two_edits.py > receipt.optimized.replay.json
    cmp receipt.json receipt.replay.json
    cmp receipt.json receipt.optimized.replay.json
    python3 audit/check_periodic_two_edits.py > audit/receipt.replay.json
    python3 -O audit/check_periodic_two_edits.py > audit/receipt.optimized.replay.json
    cmp audit/receipt.json audit/receipt.replay.json
    cmp audit/receipt.json audit/receipt.optimized.replay.json

The remaining two-edit exclusion problem is now entirely primitive after discarding distance-zero integer returns. The exceptional arbitrary-base square construction above shows exactly why the primitive one-edit theorem, rather than repetition alone, is needed to close the mechanical d=2 branch.
