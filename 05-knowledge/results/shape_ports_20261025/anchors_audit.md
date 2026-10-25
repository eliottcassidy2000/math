# Independent audit: the nonzero-defect quartic/sextic exclusion

2026-10-25. **Verdict: PASS in the explicitly stated polynomial normal form.**
The proof excludes the beta!=0 branch; together with the cited, separately
audited beta=0 predecessor it removes terminal (4,6) in the inherited
primitive monic genus-six setting. It proves neither planar JC nor a new
universal conductor bound. No high-dimensional JC claim is used.

## Reviewed input and independent method

The reviewed author draft is `anchors_quartic_sextic.md`, SHA-256
`5cb63771b337a0be877eb485d5e624cb20734f2641084efdf032b3e334ee207e`.
The author then displayed the optional constant-a decomposition explicitly.
That clarification was checked and preserves the PASS verdict; the clarified
draft hash is `8d61d7258d9d519e240332bfec781cbac9a47c60f33eb1b970061d53e13ef2ad`.
The final author version changes the audit status/link and states the direct
degree-spectrum corollary checked below. Its hash is
`954c9e945be8acc3e41a02f08e1a35c1cf7f4e1a3a23006db2602961d45a1017`.
The accompanying independent checker is `anchors_independent_audit.py`,
SHA-256 `aee16a26db9a7ce6faa75a6108c6e0c6dd9af3debd488a29f2b57100d6fe3eff`.
The checker imports no author or repository code. It reconstructs the complete
Jacobian in a formal polynomial jet ring, uses a denominator-cleared
substitution rather than the author's resultant calculation, and separately
tests comparison charts with finite poles. Sympy version: 1.14.0.

The proof itself, rather than finite computation, supplies the all-degree
exclusion. I independently read the two inherited packets at repository
revision `2b96df316ec324bb8785c45d672f5562f5f2abfa`:

- `05-knowledge/results/genus_six_source_frontier_20261020/FRONTIER.md`,
  particularly sections 6--8;
- `05-knowledge/results/conductor_ten_source_exclusion_20261020/CONDUCTOR_TEN_SOURCE_EXCLUSION.md`,
  particularly sections 4--7.

Their scope matches the use made here: generic codimension at most six forces
the missing sextic y^5 coefficient, while the old first-integral argument
closes the beta=0 normal form. The new argument does not assume beta=0 or
silently globalize a low-degree semigroup generator.

## Adversarial branch audit

1. **Types and normalization.** The source is k[x,y], char(k)=0, with the
   displayed monic y-degrees and absent y^3/y^5 coefficients. The top equations
   force beta to be a constant kappa. Dividing by rational integers is valid;
   no algebraic-closure assumption enters the new proof. The only discarded
   constants are the legal target shear Q -> Q-C P-F. The scalar terms d,f
   are retained throughout.
2. **Whole determinant, not only selected rows.** Independent expansion gives
   exactly `J=I' y^2+K' y+h d'-b f'` after the five integrated upper rows.
   Thus the final scalar coefficient is precisely the stated consumer.
3. **Elimination direction.** In the chart r=2b+kappa!=0, independently solving
   I=L and substituting into K-M gives `N=384 r^2 (K-M)`. This agrees with
   `Res_d(I-L,K-M)=3N/512`. Only necessity is used; zeros of N are never called
   Keller pairs. The omitted r=0 case is handled separately.
4. **b constant, including b=0 and r=0.** If a moved, the coefficient
   `9 kappa^2 a^4` would be the unique highest term. Hence a is constant.
   When r!=0, I=L forces d constant; when r=0, K=M is a quadratic equation in
   d with nonzero leading coefficient. Either gives J=0.
5. **a constant, including a=0.** Rewriting `4b+kappa=2r-kappa` makes d a
   constant plus a constant divided by nonconstant r. Polynomiality makes the
   numerator vanish. K'=0 then forces a=H=0. The remaining scalar is
   `-(3/8)b(2b+kappa)b'`, of degree `3 deg(b)-1>=2`. Its leading coefficient
   cannot vanish in characteristic zero. There is no dropped zero branch.
6. **Positive-degree forcing.** Once A=deg(a)>0 and B=deg(b)>0, every monomial
   of N other than a^4 and ab^4 has strictly smaller degree than
   `max(4A,A+4B)`. Their coefficients are respectively `9 kappa^2` and
   `-1152`, both nonzero. Hence `3A=4B`, so A=4m,B=3m,m>=1.
7. **Fraction-field chart.** z=a/r has degree at infinity m; it is never
   asserted polynomial. The correction numerator after subtracting the two
   main terms has weight at most 12m and denominator weight 13m, giving
   `deg_inf(b-b0)<=-m`. Consequently the a and d errors are at most 0 and
   4m. Finite poles are allowed and do not affect these infinity inequalities.
8. **Final cancellation audit.** The h and f errors are at most 3m and 8m.
   Expanding all six terms of `J-J0` bounds each by 11m-1. Independent
   differentiation gives the stated factored J0, of degree 12m-1, with leading
   coefficient `(15/8192) kappa^7 m lc(z)^12`. This is nonzero. The strict
   one-weight separation prevents the remainder from canceling it.

In particular the largest separate products h0 d0' and b0 f0' have degree
15m-1 and do cancel at their leading orders. The argument correctly computes
their difference first and uses its actual degree 12m-1; it does not assume
that either product's leading term survives.

## Supporting controls and outcome

The independent checker passes **28 exact checks**. The universe consists of
formal polynomial-jet identities and four rational comparison charts:

    z=x, x^2+1, (x^2+1)/x, (x^4+1)/(x+1), with kappa=2.

The final two controls have finite poles. All four reproduce the complete
source Jacobian and the `12 deg_inf(z)-1` degree. Normal and optimized Python
runs produce byte-identical JSON, and explicit exceptions enforce all gates.

No mathematical blocker was found. The author accepted the optional
explanatory clarification and displayed the constant-plus-constant/r
decomposition in the constant-a branch. I checked that expression against the
independent checker's exact decomposition; no mathematical correction was
needed and the PASS verdict is unchanged.

The remaining terminal patterns are `(3,7)`, `(4,5)`, `(4,10)`, and `(6,9)`.
The two alternate presentation types cannot be removed merely because their
semigroups contain smaller generators. The new proof is a closure of one
actual coefficient lane, not an equivalence among the open conjectures.

## Final direct corollary audit

PASS: the proof uses only J independent of y until the last contradiction.
Its three branches therefore show that, in the displayed kappa!=0 normal
form, a nonzero y-independent J has degree congruent to two modulo three.
The values are `12m-1` when a,b both move, `3 deg(b)-1` when a is constant
and b moves, and J=0 when b is constant. Every degree `3s-1`, s>=1, is
realized by a=H=E=d=0 and b=x^s with nonzero kappa. This is a direct
consequence of the audited identities and branch proof; no new finite
extrapolation or classification of every coefficient solution is involved.
