# All-length integer classification for two mechanical endpoint chains

20 October 2026. Status: **PROVED, using the cited Matveev estimate and independently audited exact certificates.** The arithmetic input is independently audited and frozen; both exact computational parts and the complete application have passed independent review. This packet does not classify arbitrary primitive two-edit words.

## 1. The two families and the conclusion

Let p>=3, q=p+r, p/2<r<p, gcd(p,r)=1, and let

    c_i=floor((i+1)p/q)-floor(ip/q),  0<=i<q.

Label the ones by source positions and residues

    s_j=j+ceil(jr/p)-1 modulo q,
    a_j=-jr modulo p,                       0<=j<p.

The endpoint-plus chain moves the one of residue p-1-r right by one position, then the one of residue p-1 right by one position. The endpoint-minus chain moves the one of residue r left, then the one of residue 0 left. These are actual overlapping adjacent swaps, including across the original cut. The order is part of each definition; reversing it is not a legal path.

For a final word w, write

    B(w)=sum_{i:w_i=1} 2^i 3^(number of ones after i),
    D=2^q-3^p,    g=gcd(|D|,B(w)).

The rational fixed point of this branch word under

    T(n)=n/2 if n is even, and T(n)=(3n+1)/2 if n is odd

is B(w)/D. It is an integer cycle exactly when g=|D|. Its sign is the sign of D, since B(w)>0.

**Theorem.** Neither endpoint family produces an integer cycle for any p>=3 in the displayed coprime band. This holds for both clock signs and for every rotation of the final word.

There is one legal coprime boundary outside the strict band: p=2,r=1,q=3. Its two endpoint outputs are genuine negative-cycle words:

    plus:  011 -> 110 -> 101, fixed point -7;
    minus: 011 -> 101 -> 110, fixed point -5.

Both describe rotations of the cycle -5,-7,-10. The classification retains this exception; it does not exclude it by a D>1 convention.

The whole theorem concerns these two specified endpoint families. No critical-clock restriction is assumed. Other primitive overlapping chains, disjoint edits, larger gaps, and the general two-edit problem are not classified here.

## 2. Inherited exact-content bound

The independently audited [endpoint-content theorem](../collatz_two_swap_endpoint_content_20261020/ENDPOINT_CONTENT_AND_SHORT_RELATIONS.md) proves, at full prime-power depth,

    g_plus =gcd(2*13^p-10^p, 3*13^r-2*10^r),
    g_minus=gcd(9^p-2*5^p, 2*9^r-3*5^r).

More importantly for closure, with N=ceil(sqrt(p)),

    g_plus<52^N,    g_minus<36^N.

This comes from retaining both common clock-root equations z^p=2 and z^r=3/2. The endpoint equations make z equal to 10/13 or 9/5 modulo the actual content. A pigeonhole relation ap+br=m of size O(sqrt(p)) then produces a nonzero small integer divisible by the full content. A norm alone does not retain this simultaneous root information.

For p>=25 the same theorem gives

    g/3^p<exp(-37p/300)<exp(-p/9).                     (2.1)

Suppose an endpoint output is integral. Set

    U=2^q/3^p-1=D/3^p,
    L=q log2-p log3=log(1+U).

Then |D|=g, so

    0<|U|<exp(-p/9).                                 (2.2)

The lower inequality holds because positive powers of 2 and 3 cannot be equal. For p>=25, the upper bound is less than 1/2. The elementary integral formula for log(1+U) therefore gives, for either sign of U,

    |L|<2|U|<2 exp(-p/9).                             (2.3)

It is essential to distinguish U from L: the cited Matveev estimate below applies directly to the multiplicative-minus-one expression U.

## 3. An effective global cutoff

Apply the real-field form of Bugeaud--Mignotte--Siksek, Theorem 9.4, to the nonzero expression U, with algebraic numbers 2 and 3, integer exponents q and -p, field Q, and height bounds A_1=7/10, A_2=11/10. Its parameter B is q, since p<q<2p. The theorem gives

    log|U|>-C(1+log q),
    C=1.4*30^5*2^4.5*(7/10)*(11/10)<602494200<10^9.

The theorem and its hypotheses were checked directly in the [primary paper, printed page 16](https://arxiv.org/pdf/math/0403046). The constants here are a deliberately loose specialization, not an improved logarithmic-form estimate.

Combining with (2.2) implies

    p/9<10^9(1+log(2p)).                              (3.1)

At p=10^12, the left side exceeds 31*10^9, while log(2*10^12)<30. The difference between the left side and right side has derivative 1/9-10^9/p>0 for p>=10^12. Consequently (3.1) forces

    p<10^12.                                         (3.2)

Every numerical comparison in this step is certified by rational series in `large-range/certify_endpoint_large_range.py`. In particular, log(2*10^12)=log2+12(3log2+log(5/4)) is bounded without floating-point arithmetic.

## 4. Legendre reduction and the complete large range

Put alpha=log3/log2, which is irrational: a rational equality alpha=q/p would imply 3^p=2^q. From (2.3),

    |q/p-alpha|<2 exp(-p/9)/(p log2).                  (4.1)

For every p>=100,

    4p exp(-p/9)<log2.

At p=100 this is certified by a positive partial sum of exp(100/9); the left side decreases for real p>9. Thus (4.1) is strictly below 1/(2p^2). Since gcd(p,q)=1, Legendre's continued-fraction criterion forces q/p to be an ordinary convergent of alpha.

The exact certificate isolates all convergents with

    100<=p<10^12.

There are exactly eighteen:

| q | p | sign of q log2-p log3 |
|---:|---:|:---|
| 485 | 306 | positive |
| 1054 | 665 | negative |
| 24727 | 15601 | positive |
| 50508 | 31867 | negative |
| 125743 | 79335 | positive |
| 176251 | 111202 | negative |
| 301994 | 190537 | positive |
| 16785921 | 10590737 | negative |
| 17087915 | 10781274 | positive |
| 85137581 | 53715833 | negative |
| 272500658 | 171928773 | positive |
| 357638239 | 225644606 | negative |
| 630138897 | 397573379 | positive |
| 9809721694 | 6189245291 | negative |
| 10439860591 | 6586818670 | positive |
| 103768467013 | 65470613321 | negative |
| 217976794617 | 137528045312 | positive |
| 1193652440098 | 753110839881 | negative |

For each row, the exact interval for L proves

    |L| > 2*2^(-min(floor(p/9),256))
         > 2 exp(-p/9).

The second comparison is deliberately coarse: e>2 and min(floor(p/9),256)<=p/9. Each row therefore contradicts (2.3), on either side of alpha.

The next convergent is

    q=8573543875303, p=5409303924479,

whose denominator exceeds 10^12. Convergent denominators increase, so the eighteen-row list exhausts the required large range.

### What the exact certificate checks

For t in (0,1), the code bounds the logarithm by

    S_N(t)=2 sum_{j=0}^{N-1} t^(2j+1)/(2j+1),
    S_N(t)<log((1+t)/(1-t))
           <S_N(t)+2t^(2N+1)/((2N+1)(1-t^2)).

It uses N=220, t=1/3 for log2 and t=1/2 for log3, and rational arithmetic throughout. At each continued-fraction step, both ends of the current interval have the same integer part. This proves the extracted coefficient belongs to the exact irrational alpha. Interval inversion then preserves the enclosure.

The producer performs 63 explicit checks. An independent implementation uses log3=(3log2+log(9/8))/2, only 112 series terms, and a different continued-fraction extraction, yet obtains the same eighteen rows and terminal convergent. It separately checks the absolute-gap direction and the coarse upper bound. Both implementations run identically with Python optimization enabled.

## 5. The complete small residual and the negative exception

`verify_small_endpoint_residual.py` enumerates every coprime shape with

    3<=p<100,    p/2<r<p,

for both endpoint chains, plus the boundary (p,r)=(2,1). There is no restriction to q=ceil(p alpha) or q=floor(p alpha), and no restriction on the sign of D.

This finite universe contains 1,502 shapes and 3,004 actual paths, with 1,247 positive-clock and 255 negative-clock shapes, including 98 paths crossing the original cut. Each carry is calculated both by affine composition and by the positional power sum. The exact gcd formulas are checked but are not used as the sole integrality calculation.

All 3,002 paths with p>=3 have g<|D|. The only integral outputs are the two p=2 boundary paths:

| family | source and intermediate word | output | B | D | fixed point |
|---|---|---|---:|---:|---:|
| plus | 011 -> 110 | 101 | 7 | -1 | -7 |
| minus | 011 -> 101 | 110 | 5 | -1 | -5 |

The script directly iterates their integer branch orbits and verifies every parity and closure, rather than merely accepting a divisibility test. Its 22,538 explicit checks have byte-identical normal and optimized receipts.

This residual includes p<25, where the analytic estimate (2.1) was not asserted. Combining it with Sections 3--4 excludes every p>=3 in the theorem and retains exactly the declared boundary exception.

## 6. Dependencies, reproducibility, and the remaining frontier

The proof uses:

1. the audited exact-content and subexponential theorem for these two endpoints;
2. the cited real Matveev estimate, with its hypotheses and constants verified;
3. the elementary Legendre criterion;
4. exact rational continued-fraction exclusion and the complete finite residual.

It uses neither a pre-existing global Collatz verification threshold nor an assumed critical-clock reduction. The earlier norm barrier remains correct: it says the norm shortcut cannot have a uniform exponential base below 3 on certain nearby shapes, not that actual content is comparably large.

Reproduce from this directory:

    python3 verify_small_endpoint_residual.py > small.replay.json
    python3 -O verify_small_endpoint_residual.py > small.optimized.replay.json
    cmp small.receipt.json small.replay.json
    cmp small.receipt.json small.optimized.replay.json
    python3 large-range/certify_endpoint_large_range.py > large-range/replay.json
    python3 -O large-range/certify_endpoint_large_range.py > large-range/optimized.replay.json
    cmp large-range/receipt.json large-range/replay.json
    cmp large-range/receipt.json large-range/optimized.replay.json

Independent audit commands and their exact comparison universes are in `audit/`.

The remaining primitive two-edit families generally leave a nonlinear common-root equation. This endpoint theorem does not supply a uniform content bound for those families and is not a resolution of the Collatz conjecture.
