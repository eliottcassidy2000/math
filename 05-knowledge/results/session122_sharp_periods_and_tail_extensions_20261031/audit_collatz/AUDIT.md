# Independent audit: nested native-core grammar and certificate interface

**Grammar PASS, conditional on inherited certificates.**
**Update: the full E14 closure audit is complete; see [E14_AUDIT.md](E14_AUDIT.md).**
The grammar and its exact arithmetic do not prove E14 or Collatz convergence.

## Dependencies and native charges

The audited source is `precision_lifts_ramsey_deficits_20261030/COLLATZ.md`
C1-C3/C5/C7 at commit d345f993475dd7b5272ac6c3b4ab693899515a38,
together with session121's `collatz/QUANTIZED_NATIVE14.md`.
The inherited E9/R2 and E10/R1 large certificates were not rerun in this
small arithmetic audit. Their actual-word terminal consequences, not an
unrestricted assertion about every altered abstraction, are the inputs.

For H_e=5A/8-S_e, inherited E9/R2 gives non-one primitive10 cores

    H10 >= (2832031/781250-1)+A/(8*5^10).

The score H10=j10/8 is quantized. The only possible score below22/8 is
21/8, which forces A<=25. All300 pairs1<=L<A<=25 fail the necessary
multiplier equality3^L=2^A modulo5^10. Thus every such core has H10>=22/8.

A residual session121 primitive11 core w has score j/8 with j=23 or24.
Cut at its actual level10 identity returns. Lower-layer native charges
are additive because every entering map is identity in that quotient.
There is only one level11 return, at the end. An all-one primitive10 block
has clock T=2*5^9 and score1/4. Hence

    H11(w)=sum_nonone H10(c)+a/4-1.

Two non-one cores would force H11>=2*(22/8)-1=9/2>3. There is at least
one because w is non-one, so there is exactly one. Its score is
j10=j+8-2a. This yields0<=a<=4 for j23 and0<=a<=5 for j24.

The word has the necessary form1^(pT)c1^((a-p)T). Both p<5 and a-p<5.
The prefix restriction follows from the first-return property. For the
suffix restriction, a suffix of five all-one blocks is identity modulo
5^11. Since the complete word is also identity there, the prefix preceding
that suffix is identity too, a proper earlier return. This argument
requires the complete native endpoint and does not follow from a scalar
clock alone. The allowed placements count15+19=34.

## Endpoints, clocks, and exact branch census

An all-one block acts by H=(h,h-1) modulo5^14. Full identity of
H^(a-p) after c after H^p forces c=H^(-a), independently of p.
This quotient loses intermediate return and charge information, so the
placements are retained separately and are not asserted realizable.

The core clock is bounded by both its inherited fee and the total clock:

    A_c<=min((8*5^10)*(j10/8-gamma), C_j-aT),
    A_c=5j10 mod8.

The discrete logarithm of2 base3 modulo5^14 is2883827743 modulo
4*5^13. Thus the necessary length is its stated affine residue in A_c.
All branch caps are below this length modulus, so the floor-difference
count correctly tests precisely1<=L_c<A_c without multiple possible lengths.

`check_nested.py` independently computes the branch counts with a recursive
lattice-reflection floor-sum routine, rather than the producer's iterative
Euclidean routine. It passes17,100 signed small literal controls, all300
small-clock exclusions, exact affine compositions, clock bounds/residues,
and first/last original modular pairs for every branch. Its eleven branch
records agree with the producer's complete records:34 forms and423,868
branch-counted multiplier pairs. These are not actual words or a complete
realizability census. Normal and optimized receipts are byte-identical.

## Conditional convex-combination interface

Let an independently verified E10/R2 certificate at slope5/8 have terminal
fee at least33/8 for actual non-one words ending identity modulo5^10.
Combine this actual-word inequality with inherited E10/R1 fee eta10 at
lambda12, in proportions24/25 and1/25. The slope is exactly lambda14,
and the terminal fee would be

    (24/25)*(33/8)+eta10/25
      =5017089801/1220703125>4.

The same calculation with32/8 gives4870605426/1220703125<4 and does
NOT close the desired threshold. The integer probe is not itself a proof.
The nested level11/12/13 return argument can consume three terminal units;
if the displayed fee is certified, a uniform primitive13 fee greater than1
follows. The separate paid-return criterion and final unfinished-suffix
handling must still be cited exactly. No certificate value is inferred
from this conditional algebra.

## Independent checker review and completed small control

`audit_integer.cpp` in this audit folder adapts the pinned source's separate
z-major checker to integer weightsD8,N5 and64 explicit letters. It uses
repeated modular halving, no producer transition table, no maximization,
and the separately proved global-minimum edge bound described below. It explicitly rejects the producer's
unreached sentinel and checks both infinite-letter tails.

    c++ -O3 -std=c++17 audit_integer.cpp -o /tmp/session122-independent-audit64
    /tmp/session122-independent-audit64 7 /tmp/s122-bucket7.bin 2

The completed E7/R2 control checks1,812,500 states,131,750,000 semantic
edges (8,562,359 explicit and123,187,641 bounded) and630 paid entries. It reports terminal-42, minimum-104,
maximum-10, and passes both tails. Raising its terminal value to zero
is rejected on edge(z,m,a,m2)=(0,1,1,39064). These controls do not
substitute for an E10 certificate.

The producer-side separate replayer additionally checks all64 letters and
uses this valid shortcut. Let common=min(oldmod,newmod) and let b be the
unit target multiplier residue modulo common. If b!=1, v5(b-1) is below
the precision of common and equals v5(m2-1) for every compatible lift m2.
If b=1, the convention v5(0)=E is an upper bound. Thus
min(h_target,v5(b-1)) bounds every compatible reward. If the resulting
right-hand side is at most the minimum potential, all compatible target
inequalities hold. The number of such edges is exactly newmod/common.
This shortcut is an inequality certificate, not an uncounted prune.

The completed E10 audit and independent corridor reconstruction are recorded in
E14_AUDIT.md. Normal/NDEBUG receipts agree. The conditional interface above
is now discharged by the complete near-tight multiplier exclusions, rather
than by a raw abstract terminal value of33.
