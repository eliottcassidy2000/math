# A unique primitive10 core in every residual E14 candidate

**CONDITIONAL PROVED necessary structure.** This uses the inherited E9,R2
certificate in `precision_lifts_ramsey_deficits_20261030/COLLATZ.md` C7
and the session121 reduction. Neither inherited large certificate is rerun
by the small checker here. This precursor alone does not close E14; the later complete closure is in `E14_CLOSURE.md`.

Let H_e(w)=5A(w)/8-S_e(w). All clocks are original positive-letter clocks.
For a non-one primitive10 excursion c, the inherited terminal bound gives

    H10(c) >= gamma + A(c)/D,
    gamma=2832031/781250-1=2.62499968, D=8*5^10=78125000.

H10(c)=j10/8 with integral j10. Thus j10>=21. At j10=21 the
clock is at most25. The finite check excludes 3^L=2^A modulo5^10 for
all300 pairs 1<=L<A<=25. Therefore

    H10(c)>=22/8.

Take a session121 surviving primitive11 core w, with H11(w)=j/8,
j in {23,24}, which is also primitive13 and has full identity endpoint
modulo5^14. Cut at its actual10 returns. If b is the number of non-one
constituents and a the number of all-one constituents, exact additivity gives

    H11(w)=sum_nonone H10(c_i)+a/4-1.

Here an all-one primitive10 block has length and clock T=2*5^9=3906250
and H10=1/4. At least one constituent is non-one. Two would give
H11>=2*(22/8)-1=9/2>3. Hence b=1.

All other letters are ones, and the complete word has the necessary form

    w=1^(pT) c 1^((a-p)T),
    j10(c)=j+8-2a.

For j=23, 0<=a<=4. For j=24, 0<=a<=5. Native primitive11 status forces
p<5 and a-p<5: five consecutive all-one10 blocks already return at11
from an actual10 endpoint and, at either end of this displayed word, would
produce a proper initial11 return. Thus there are15+19=34 placement forms.
This does not assert their realizability; the entire word must still have
all its required first-return properties.

## Retained clocks and endpoint map

Let C23=48828467 and C24=97656592 be session121's total-clock caps.
For every form,

    A(c)<=min(C_j-aT, D*((j+8-2a)/8-gamma)),
    A(c)=5*(j+8-2a) mod8.

Put q=5^14 and h=(3/2)^T modulo q. An all-one block has affine map
(h,h-1). Full identity of w forces the *actual* core map to be

    (m_c,t_c)=(h^(-a),h^(-a)-1) mod q.

This endpoint formula is independent of p, but intermediate charges and
first returns are not. Placements therefore remain separate.

With M=4*5^13 and log_3(2)=2883827743 mod M, the necessary multiplier
condition becomes

    L_c=log_3(2)*A_c+aT*(log_3(2)-1) mod M,
    1<=L_c<A_c.

The exact floor-sum checker counts these pairs along the retained modulo8
clock progression and extracts each branch's first and last pair. It checks
all affine compositions, clock caps and score identities,300 excluded small
pairs, and37440 small literal floor-sum controls. Normal and optimized
outputs match.

The summed423868 pairs are across the eleven core-score/padding branches,
not distinct original words or distinct original clock pairs. They must not
be compared as a smaller replacement for the earlier scalar pair census.
Translation, positive-word realization, native10 status of c, and native11/
13 status of w remain required. The refinement is a structural grammar,
not a finite word exhaustion or a proof of full E14.
