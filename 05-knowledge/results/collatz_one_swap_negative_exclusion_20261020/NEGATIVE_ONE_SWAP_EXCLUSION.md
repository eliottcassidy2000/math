# Negative integer cycles in the one-swap mechanical family

Date: 20 October 2026. Status: **PROVED for the stated one-swap family, with independently audited exact residual certificates**. The argument concerns a prescribed parity-word family, not all negative cycles or the Collatz conjecture. The negative cycle through -17 is retained as a hostile control, not excluded.

## 1. Exact statement and sign convention

Use the shortcut map

    T(n)=n/2 for even n,    T(n)=(3n+1)/2 for odd n,

on integers of either sign. For a chronological binary word w of length q with p ones, let

    B(w)=sum_(i:w_i=1) 2^i 3^(number of ones after i),
    D=2^q-3^p,    g(w)=gcd(|D|,B(w)).

The rational fixed point is B(w)/D. Since B(w)>0, a negative clock D<0 gives a negative rational fixed cycle. It is an integer cycle exactly when g=|D|; all rotations share this denominator, and integrality forces the indicated parities.

Let c be the zero-intercept lower mechanical word

    c_i=floor((i+1)p/q)-floor(ip/q),    0<p<q,

and let w be obtained by one actual cyclic adjacent swap of unequal bits.

**Theorem N (primitive negative classification).** If gcd(p,q)=1 and D<0, then

    g(w)=|D|  iff  (p,q)=(2,3).

At this exceptional shape the zero-intercept source is 011 and the two distinct edited words are 101 and 110. Their fixed points are -7 and -5, respectively, and both belong to

    -5 -> -7 -> -10 -> -5.

Every primitive one-swap negative clock with |D|>1 has g<|D|.

**Theorem R (nonprimitive case, separate mechanism).** If gcd(p,q)>1, then no such one-swap word is an integer cycle of either sign. This follows from the repetition-factor argument in Section 7, not from applying the primitive norm theorem outside its hypotheses.

Thus the same negative classification holds without a coprimality assumption. The fixed negative cycle {-1} has a constant parity word and is outside the nonconstant, unequal-swap family; it is not being denied. The known cycle through -17 also remains a genuine integer cycle, outside the one-swap family.

## 2. The negative-clock norm gate has base sqrt(6)

Inherit the cut-safe polynomial/norm construction from [the sharp one-swap content theorem](../collatz_one_swap_sparse_norms_20261020/SHARP_ONE_SWAP_CONTENT_BOUND.md) and the legality calculation from [the critical-clock refinement](../collatz_one_swap_critical_gap_20261020/CRITICAL_CLOCK_SQRT7_GAP.md). These establish exact prime-power divisibility, not merely prime support.

Assume first gcd(p,q)=1 and D<0. Then

    p<q<p log_2 3<2p,
    q=p+r,    1<=r<p,
    r/p<log_2 3-1<3/5.                                 (2.1)

The last inequality follows from 3^5<2^8. Let a be the residue label of the moved one in the norm construction. At k=1, its complete legality conditions are

    plus move is legal iff a<r,
    minus move is legal iff a>=p-r.                     (2.2)

The two large-norm endpoint types a=p-1,+ and a=0,- are therefore forbidden. The other endpoint types are free, with content 1.

For interior types, put rho=2^(1/p) and z_j=rho exp(2 pi i j/p). The norm polynomials are

    Q_+(Z)=Z^(a+1)-Z^a+1,
    Q_-(Z)=Z^(p-a)-Z+1.

Their three exponents are distinct modulo p, and their absolute resultants bound g in divisibility. Orthogonality gives exact second moments.

For every legal plus interior, a+1<=r, so

    mean_j |Q_+(z_j)|^2
      =1+4^(a/p)+4^((a+1)/p)
      <1+2*4^(r/p)
      <1+2*4^(3/5)<6.                                  (2.3)

Here 4^(3/5)<5/2 follows by raising positive quantities to the fifth power: 4^3=64<(5/2)^5.

For every legal minus interior and p>=3, p-a<=r, so

    mean_j |Q_-(z_j)|^2
      =1+4^(1/p)+4^((p-a)/p)
      <1+4^(1/3)+4^(3/5)
      <1+2+5/2<6.                                      (2.4)

AM-GM and the norm divisibility give the uniform negative-clock gate

    g(w)<6^(p/2).                                       (2.5)

For p=2 the only negative shape with p<q is q=3; its two legal endpoint outputs have g=1, so (2.5) holds there too. There are no negative clocks with p=1 and q>p.

## 3. A negative cycle must lie immediately below the critical slope

Suppose g=|D|, and put

    R=|D|/3^p=1-2^q/3^p,
    delta=p log 3-q log 2=-log(1-R)>0.

The norm gate gives

    R<(2/3)^(p/2)<exp(-p/6),                            (3.1)

because log(3/2)>1/3. For p>=4,

    R<(2/3)^2=4/9<1/2.

Consequently

    0<delta<R/(1-R)<2R<2 exp(-p/6),                    (3.2)

and 2^q/3^p=1-R>1/2. Together with D<0, this implies

    p log_2 3-1<q<p log_2 3,
    q=floor(p log_2 3).                                (3.3)

The logarithmic bound and the lower critical length cannot be copied from the positive proof without the factor 2 and this small-R argument. For p=2 and p=3, the only possible negative lengths are q=3 and q=4, respectively, which already satisfy (3.3). These cases are retained in the direct finite certificate.

## 4. Matveev bounds p effectively

The sole imported transcendence input is the real-field Matveev bound in Bugeaud, Mignotte and Siksek, *Classical and modular approaches to exponential Diophantine equations I. Fibonacci and Lucas perfect powers*, Theorem 9.4, printed page 16 of [arXiv:math/0403046](https://arxiv.org/pdf/math/0403046). It explicitly bounds the absolute value of the multiplicative form, so a negative form is allowed.

Apply it to

    Lambda=2^q 3^(-p)-1=-R,

with alpha_1=2, alpha_2=3, b_1=q, b_2=-p, n=2, field degree 1, B=q, A_1=7/10, A_2=11/10. The rational interval certificate checks the height inequalities. The resulting prefactor satisfies

    C=1.4*30^5*2^(9/2)*(7/10)*(11/10)
      <602494200<10^9.

Therefore

    log R>-10^9(1+log q)>-10^9(1+log(2p)).              (4.1)

Combining this with (3.1) gives

    p<6*10^9(1+log(2p)).                                (4.2)

The function p/(1+log(2p)) is increasing for p>=1. At p=10^12, log(2p)<30 and 6*10^9*31<10^12, so

    p<10^12.                                           (4.3)

The proof uses |Lambda|=R, not delta, in Matveev. No factor from the logarithm conversion is lost. It needs no imported finite Collatz convergence computation or minimum-cycle bound.

## 5. Complete continued-fraction exclusion for p>=100

For p>=100, (3.2) yields

    0<log_2 3-q/p
      <2 exp(-p/6)/(p log 2)<1/(2p^2).                 (5.1)

The final inequality is equivalent to 4p exp(-p/6)<log 2. Its left side decreases for p>6; the certificate checks it at p=100 using a positive rational Taylor partial sum for exp(100/6) and a rational lower bound for log 2.

Because gcd(p,q)=1, Legendre's criterion forces q/p to be an ordinary convergent of log_2 3. The exact interval algorithm enumerates all 18 convergents with 100<=p<10^12, including both sides of the irrational number.

- Nine are above log_2 3 and have the wrong clock sign.
- The nine below it have (q,p) equal to

      (1054,665), (50508,31867), (176251,111202),
      (16785921,10590737), (85137581,53715833),
      (357638239,225644606), (9809721694,6189245291),
      (103768467013,65470613321), (1193652440098,753110839881).

For each lower convergent the certificate proves, by rational comparison,

    delta > 2*2^(-min(floor(p/6),256))
           >=2 exp(-p/6),                              (5.2)

contradicting (3.2). The terminal convergent has denominator 5409303924479>10^12, so no candidate denominator in the relevant interval is omitted. Intermediate fractions cannot be candidates because (5.1) invokes the stronger ordinary-convergent criterion.

## 6. Direct small-weight certificate, with exceptions preserved

For every 2<=p<100, use the necessary lower critical length q=floor(p log_2 3). The program calculates it exactly as bit_length(3^p)-1 and independently verifies that value with rational logarithm intervals. After the primitive filter, 59 shapes remain. All 3,348 cyclic unequal-adjacent-swap occurrences are generated, including the original-cut crossings.

Each carry is computed independently by affine recurrence and by its positional sum. The program records all integral outputs before checking their classification. Exactly these survive:

| p | q | edited word | D | fixed point |
|---|---|---|---|---|
| 2 | 3 | 101 | -1 | -7 |
| 2 | 3 | 110 | -1 | -5 |

Every other supplied word has g<|D|. As an independent hostile universe, the same program checks every primitive negative mechanical shape with q<=120, not just lower critical lengths: 47,892 actual swap occurrences, with the same two integral outputs.

There is no hidden primitive |D|=1 shape. If 3^p-2^q=1 with 0<p<q and D<0, then q>=3. Reduction modulo 8 makes p even, say p=2m. The two even factors (3^m-1)(3^m+1)=2^q differ by 2 and have gcd 2. Therefore they are 2 and 4, forcing (p,q)=(2,3).

The genuine -17 cycle is an explicit hostile control:

    -17,-25,-37,-55,-82,-41,-61,-91,-136,-68,-34,

with parity word 11110111000, p=7,q=11,D=-139 and integral carry. Its minimum Hamming distance from the mechanical source 01011011011 over all rotations is 4, not 2. Thus it cannot be a one-adjacent-swap output, and is correctly retained as an actual cycle outside the theorem's family.

The stronger [exact adjacent-swap distance certificate](MINUS17_THREE_SWAP_BOUNDARY.md) shows that this genuine cycle is exactly three cyclic adjacent swaps away, not two. This preserves a sharp three-edit hostile boundary without claiming a two-edit classification.

Sections 2-6 prove Theorem N without treating negative clocks as positive ones.

## 7. Nonprimitive words: a separate repetition-factor proof

For d=gcd(p,q)>1, the mechanical source is c=v^d with reduced shape (q_0,p_0)=(q/d,p/d). Write

    D_0=2^(q_0)-3^(p_0),
    D=H D_0,
    H=sum_(j=0)^(d-1) 2^[q_0(d-1-j)]3^(p_0j)>1.

Then H is coprime to 6. Choose a cut not separating the two edited positions. The unedited rotation is a repetition of a rotated base, with carry H B(v'). One adjacent unequal swap changes that carry by a single signed 2,3-smooth monomial. Hence gcd(H,B(w))=1, and

    gcd(|D|,B(w))=gcd(|D_0|,B(w))<=|D_0|<|D|.

Rotation preserves each gcd because every modulus used divides D and the rotation cocycle multiplies carries by units modulo those divisors. This proves Theorem R for either clock sign. No primitive sparse-norm assertion is being extended without proof.

## 8. Reproduction, status, and limits

`certify_negative_swap_exclusion.py` uses only Python's standard library. All proof decisions use exact integers or fractions; no floating-point value controls a branch. It records 58,247 explicit checks, the exact logarithm intervals, all 18 large convergents, all 98 small p-values, the direct small residual, broad hostile census, and the two actual negative-cycle controls.

The log interval uses the positive series

    log((1+t)/(1-t))=2 sum_(j>=0) t^(2j+1)/(2j+1)

at t=1/3 and 1/2, with 220 terms and the exact geometric tail bound. A continued-fraction digit is accepted only when both interval endpoints have the same floor; reciprocal updates reverse the endpoints. The full normal and optimized receipts agree.

Reproduce from this directory:

    python3 certify_negative_swap_exclusion.py > negative-receipt.replay.json
    python3 -O certify_negative_swap_exclusion.py > negative-receipt.optimized.replay.json
    cmp negative-receipt.json negative-receipt.replay.json
    cmp negative-receipt.json negative-receipt.optimized.replay.json

The [completed independent audit](audit/NEGATIVE_ONE_SWAP_AUDIT.md) passed with 1,044,780 exact checks and an expanded hostile census of 113,674 swaps through q=160. It verifies the source application, the signed norm/logarithm steps, the exception list, and the interval/continued-fraction completeness. The result classifies this one-swap mechanical family only. It does not classify all negative cycles, prove convergence of negative trajectories, or justify several-edit analogues. No historical-priority claim is made for Matveev, Legendre, the continued-fraction prefix, or the classical carry identities.
