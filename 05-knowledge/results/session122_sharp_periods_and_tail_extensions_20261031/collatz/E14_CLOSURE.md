# Closing the fourteenth finite content layer by restoring clock equality cases

**PROVED with finite certificates and independent replay.** The full E10,R2
potential and independently reconstructed near-tight graph pass normal and
NDEBUG replay. Arithmetic and multiplier controls pass normal and optimized
Python. The inherited E10,R1 rational certificate is a dependency, not rerun
here. See the companion `../audit_collatz/E14_AUDIT.md`.

## Result and exact scope

For every finite word w of positive valuation letters, with original clock
A=sum(a_i), define its affine prefix maps by

    F_u(X)=(3^len(u) X+B_u)/2^A(u),
    c(u)=v5(gcd(3^len(u)-2^A(u),B_u)),
    S14(w)=sum over nonempty prefixes u of min(c(u),14).

The sharp inequality is

    S14(w)<= (762939453/1220703125)*A.

For a word containing a non-one letter, the stronger bound is

    S14(w)<= (762939453/1220703125)*A
               -134277301/1220703125.

This closes the finite-depth E14 content obligation. The all-layer bound
and Collatz convergence remain open. The theorem does not turn a formal
infinite valuation word into a positive integer orbit.

## 1. A potential that is deliberately not a semantic proof by itself

The inherited height-adaptive graph retains z=t-m+1 modulo5^10 and m
modulo5^max(2,v5(z)), with heights truncated at10. Its compatible lifts
admit every actual word, but can admit artificial multiplier histories.
The new integer potential uses excess8*S10-5*A. Its verified output has

    245312500 states,
    minimum V=-128, maximum V=-10, terminal V(0,1)=-30.

The producer's40-letter proposal does not itself certify all letters: its
coarse omitted-entry test fails. The separate checker therefore tests all
letters1..64 and every compatible lift, all630 paid canonical entries,
and the65+ tails. The final range is more than sufficient for the65+ tails.
That full independent replay passes all18850000000 edge inequalities and
all630 entries; normal and NDEBUG receipts match.

Even a valid terminal value -30 gives only the abstract bound H10>=30/8,
where H10=5*A/8-S10. It does not by itself improve the previous actual
bound. We instead recover every path within two integer score units of it.

## 2. Complete near-tight paths, not one artificial extremizer

`NEAR_TIGHT_CLOCK_THEOREM.md` proves why this is a finite exact problem.
Use augmented states(state,b), where b is the remaining nonnegative
integer potential deficit,0<=b<=2. All inverse compatible lifts are
checked. Original letters add their own(a,1) to(clock,length).

An arbitrary initial all-one prefix is also retained. If its canonical
phase is s modulo T_ell=2*5^(ell-1), then changing it to s+b*T_ell loses
exactly

    2*s5(b mod5^(11-ell))+10*floor(b/5^(11-ell))

integer score units. At deficit at most two, the only extra phases are
s+T_e for ell<=e<=10. This is a proved digit-sum formula, not a phase cutoff
chosen from numerical experiments.

The complete independently reconstructed graph has3908743 augmented states and3908754
edges and is acyclic. Its exact output is:

| integer score 8H10 | abstract paths | possible(clock,length) pairs |
|---|---:|---|
|30|2|(342,326)|
|31|0|none|
|32|38|(A,A-16), A in {72,352,392,552,592,1592,6592,31592,156592,781592,3906592}|

Every one of these twelve pairs already fails modulo625 and hence fails

    3^length=2^clock mod5^10.

Ten pairs have exact multiplier congruence height2; the pairs(352,336)
and(552,536) have exact height3. Thus none can represent an actual non-one
word whose endpoint is identity
modulo5^10. Translation and native-return constraints need not be inferred
or tested, because the necessary multiplier condition already fails.
Consequently, for every actual word in that endpoint class,

    H10>=33/8.                                         (1)

The graph includes long all-one rails and prefixes; the largest clock is
3906592. Excluding only the short representative(342,326) would not prove(1).

## 3. Integer and rational prices combine on the same actual word

The inherited E10,R1 certificate, at source commit
0578f16726e4e535189fb6074cc2ddb024c95d82,
`precision_lifts_ramsey_deficits_20261030/COLLATZ.md` C5, gives

    lambda12*A-S10>=eta,
    eta=183105426/48828125,
    lambda_e=5/8-1/(8*5^(e-1)).

Combine this with(1), with weights1/25 and24/25. Since

    lambda14=(24/25)*(5/8)+(1/25)*lambda12,

we obtain, for the same actual endpoint class,

    lambda14*A-S10 >= eta14,
    eta14=5017089801/1220703125 >4.                       (2)

The two abstractions need not identify artificial paths. Both inequalities
are applied only to the actual word and its unchanged clock.

## 4. The native returns now pay every E14 family

Each actual non-one primitive11 excursion c has S11(c)=S10(c)+1.
Thus(2) gives its lambda14 fee at least eta14-1.

Cut a non-one primitive13 word into its actual12 excursions and then cut
the non-one constituents at their actual11 returns. Let B be the number
of non-one12 constituents and C the number of non-one11 cores. Then
C>=B>=1. Every all-one11 constituent pays31/125, and every all-one12
constituent pays6/25. There is one terminal12 charge per non-one12 block
and one terminal13 charge. Hence

    delta13 >= (eta14-1)*C-B-1 + nonnegative all-one fees
             >= eta14-3
             =1+alpha,
    alpha=134277301/1220703125 >0.                       (3)

This holds for every actual non-one primitive13 word, hence for every
specified tangent class. In particular both inherited E14 obligations
are paid: g13(0)>=1 and g13(-(4,4))>=4/5. The paid-return criterion in the
same source's C6 yields the displayed S14 inequality.

Here is the boundary argument explicitly. At actual13 returns, record the
next affine tangent(alpha0,beta0) in F5^2. Tangents add under composition
modulo5^14. Set Phi(alpha0,beta0)=((-alpha0) mod5)/5. Let I be the indicator
that the new tangent is zero, namely the extra14th-layer charge.
For every arbitrary tangent increment, I+DeltaPhi<=1: if I=1 the new
Phi is zero, and otherwise DeltaPhi<=4/5. For an all-one13 excursion the
increment is(4,4), and I+DeltaPhi<=1/5. The latter follows by the five-step
carry in((-alpha0) mod5). Its actual fee is exactly1/5. Thus(3) pays every
non-one increment with surplus at least alpha, and every all-one increment
without deficit. Summation telescopes Phi; the initial value is zero and
the final value is nonnegative.

An unfinished suffix has no13 return, so its S14 equals S12, in the native
clock because its entry map is identity modulo5^13. If non-one, the inherited
global S12 slack gamma12=61035156/48828125 at slope lambda12 pays it at the
larger slope lambda14; gamma12 exceeds alpha. If all-one, the exact floor
sum gives S12<=lambda12*A<=lambda14*A. Thus the unfinished suffix is not
dropped. Every non-one full word has either a completed non-one excursion
or a non-one suffix and retains at least alpha.

## 5. Sharpness and verification boundary

At T14=2*5^13=2441406250 consecutive ones,

    S14=(5^14-1)/4=1525878906,
    S14/T14=762939453/1220703125.

This proves sharpness of the coefficient by exact finite floor sums;
no billion-letter enumeration is claimed.

The producer is a finite certificate proposal. The all-letter checker,
complete near-tight graph replay, clock-pair exclusions and interpolation
are separate verification obligations. No Lean formalization or external
human peer review is claimed. All-layer control and convergence remain open.
