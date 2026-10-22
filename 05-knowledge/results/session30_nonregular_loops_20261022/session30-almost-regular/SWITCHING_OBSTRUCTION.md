# A complete almost-regular switching family

Status: elementary PROVED family description; FINITE-EXACT H values through n14. No spectral asymptotic claim.

Let n=2m and let S be the transitive skew matrix S_ij=−1 for i<j. For signs d_i=±1 put B=DSD. Every such B is a tournament skew matrix and has exactly the same characteristic polynomial as S. Fix d_1=1 to identify the duplicate D and−D.

Write p_i=Σ_(j<i)d_j and a=Σ_jd_j. The i-th signed score is

(B1)_i=d_i(2p_i−a)+1.

Thus B is almost regular iff each d_i(2p_i−a) belongs to{0,−2}. At i1 this forces a∈{0,2}.

If a0, the prefix walk starts at0 and each return to0 may take either ±1 step, after which the next step is forced to return. Hence d=(1,−1, ε_2,−ε_2,...,ε_m,−ε_m), with2^(m−1) choices.

If a2, the first step is+1. Every internal excursion from1 to0or2 is followed by a forced return to1, until the final+1 step ends at2. Hence d=(1,ε_1,−ε_1,...,ε_(m−1),−ε_(m−1),1), also2^(m−1) choices.

Therefore precisely2^m different labelled tournaments in this switching orbit are almost regular. Distinct normalized sign choices give different matrices because DSD=D'SD' forces all ratios d_i/d'_i equal.

All have characteristic polynomial ((λ+1)^n+(λ−1)^n)/2. This last polynomial identity is independently checked by exact Newton identities at n≤14; a standard determinant recurrence for the transitive skew matrix supplies the general identity. The spectrum invariance itself needs no polynomial formula, since D²=I.

Already n6, paired signs(1,−1,1,−1,1,−1) giveH41, whereas(1,−1,−1,1,1,−1) giveH45. Both H values are verified by enumerating all720 vertex orderings. This is an exact obstruction to any spectrum-only formula for H on almost-regular tournaments.

Complete fixed-switching-orbit H histograms:

- n2:{1:2}.
- n4:{5:4}.
- n6:{41:6,45:2}.
- n8:{629:8,653:8}.
- n10:{14937:10,15313:10,15461:10,15565:2}.
- n12:{513669:12,521869:12,526425:24,530293:12,531141:4}.
- n14:{24104937:14,24369345:14,24542997:28,24603877:14,24698969:14,24753381:14,24781785:14,24815125:14,24855901:2}.

These frequencies count labelled members of one fixed switching orbit, not all labelled tournaments or isomorphism classes. They invite a binary-word orbit description but no such classification is claimed here. Exact bit strings and signs are in switching_family.jsonl; computation is switching_family.py.

## Stronger ordered-score obstruction

In the a0 paired family, each first vertex of a pair has prefix0 and score+1; the second has prefixε and sign−ε, hence score−1. Thus the entire ordered signed-score vector is always(+1,−1,...,+1,−1), independent of the pair choices. The n6 H41/H45 witnesses above therefore have the same exact skew spectrum AND the same ordered score vector. Neither spectral data plus score data nor scorevariance can determine H. The missing information is their alignment: quantities such as δᵀf(M)δ distinguish how the degree-defect vector interacts with the orientation operator.

## Polynomial formula proof

Write U for the strictly upper-triangular all-ones matrix, so S=J−I−2U. Then λI−S=B−11ᵀ with B=(λ+1)I+2U. For λ≠−1, solving Bx=1 backwards gives x_i=(λ+1)^−1 r^(n−i), r=(λ−1)/(λ+1). Hence1ᵀB^−11=(1−r^n)/2, and the determinant lemma yields det(λI−S)=(λ+1)^n(1+r^n)/2. Polynomial continuation includesλ=−1. This proves the displayed formula for every n.
