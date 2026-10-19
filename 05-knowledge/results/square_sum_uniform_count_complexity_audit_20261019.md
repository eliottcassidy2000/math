# Independent audit: uniform contact-count complexity

**Verdict: PASS.** No proof-critical gap was found in Sections 1–7 of
`square_sum_uniform_count_contacts_20261019.md`. This audit checks the new
uniform-in-r argument, taking the stated formal-coordinate, signed
recurrence, non-torsion, and orbit-height predecessor theorems as inherited.
It does not certify companion implementation controls added subsequently.

## Proof-critical checks

1. **Integral model and polynomial substitution.** The short-model
   coefficients and both B-coordinate formulas satisfy the curve equation
   identically. Their degrees in h are respectively 4, 8, 12 for
   xi_B, alpha, beta, and 5 or 4 for eta_B. Thus the claimed
   degree-to-weight comparison holds. A monomial of F_n acquires degree
   at most twice its weighted degree; its coefficient norm increases by
   at most H raised to its total degree. The inherited universal norm
   bound therefore gives actual integer arrays of degree O(n²) and
   coefficient bit length O(n² log(n+1)). No bound on arbitrary polynomial
   factors or reduced rational expressions is being assumed.

2. **Multiplication and nonvanishing at the cusp.** The cited
   [MIT Lecture 5, Sections 5.5–5.6](https://ocw.mit.edu/courses/18-783-elliptic-curves-fall-2025/mit18_783_f25_lec05.pdf)
   gives the parity normalization and multiplication coordinates. Its
   definitions yield psi_(2d)=2 psi_d omega_d, so the displayed
   eta-coordinate eta_B G_(2d)/G_d⁴ is correct for odd d. The theorem is
   applied over Q(h), where the curve is nonsingular; specialization at
   h=0 is subsequently only polynomial evaluation.

   At that specialization, the smooth point B_0=(-48,0) has cubic
   derivative 5184, giving G_3(0)=-5184². For every larger odd index,
   precisely the term containing the cubic square vanishes in the
   normalized recurrence; the other term is a product of smaller odd
   entries. Induction proves nonvanishing, without applying a
   nonsingular-curve theorem to the singular fiber. The plus subtraction
   has x-coordinate difference 36 and yields original coordinates (3,4);
   the minus point remains (-1,0). Both inverse-chart denominators are
   nonzero, as are the fixed model-change denominators. Consequently raw
   rational-pair arithmetic needs no cancellation at the cusp.

3. **Uniform positivity cutoff.** The inherited first derivatives make
   the cross-multiplied difference numerators vanish to exactly first
   order, since their denominator constants are nonzero. Removing their
   single factor h by an array shift preserves all coefficient bounds.
   Each resulting numerator and denominator has nonzero integer constant
   term and degree/log-norm bounds O(r²), O(r² log(r+1)). Signing these
   polynomials and choosing b=1+bitlength(M), with the stated lower bound
   b>=4, gives the required positive tail estimate throughout the closed
   cutoff interval. Hence b=O(r² log(r+1)). The elementary estimate
   f(1+h)>=h on this interval validates J_r=b+2 directly.

4. **Fallback budget and output heights.** The inherited r-point signed
   pigeonhole argument uses the positive-z anchor on both ports and
   gives M_r=2r floor(20rW_r)=O(r⁴ log(r+1)). Its guaranteed points have
   distinct spreads, lie in the accepted chart, and differ from the
   anchor. Duplicate filtering therefore preserves termination. Inserting
   M_r into the uniform orbit-height theorem yields
   B=O(r⁸ log²(r+1)(L+log(r+1))) bits for every stored candidate and
   fixed-formula intermediate. Cusp evaluation gives the smaller
   O(r²(L+log(r+1))) bound. Primitive normalization and the rational square
   certificate preserve these orders.

5. **Actual polynomial preprocessing and resource bounds.** Computing
   G_1 through G_(4r+2) by the normalized recurrences is a dependency-safe
   increasing-index computation with O(r) entries. Each has O(r²)
   coefficients of O(r² log(r+1)) bits; all temporary products have the
   same orders. Dense schoolbook convolution therefore costs at most
   O(r⁸ log²(r+1)) per multiplication and O(r⁹ log²(r+1)) overall, with
   O(r⁵ log(r+1)) table storage. The fixed number of rational-pair
   operations per section and cutoff construction obey the same bound.

   Horner partial evaluations obey the same degree-plus-coefficient-height
   estimate as the full evaluation. O(r³) rational operations, each
   bounded by the stated cubic gcd cost, justify the cusp bound
   O(r⁹(L+log(r+1))³). The scan has O(rM_r) operations including duplicate
   checks. Thus O(rM_rB³)=O(r²⁹ log⁷(r+1)(L+log(r+1))³) and working space
   O(rB+r⁵ log(r+1)) are valid conservative bounds.

6. **Quantifiers and scope.** The algorithm is polynomial in r and L,
   not in log r and L. The numerical-height consequence
   C_i<=C raised to O(r⁸ log³(r+1)) follows because C>=5. This is a
   polynomial-in-r exponent for a one-anchor contact construction; it
   does not imply polynomial output size after a growing number of
   successive chain extensions. The draft states these distinctions.

## Targeted independent arithmetic checks

Symbolic expansion independently confirmed both integral B sections and
their common specialized cubic. As a limited recurrence control, odd
normalized constants through index 101 agreed exactly with
(-1)^((d-1)/2) 5184^((d²-1)/4), in particular remaining nonzero. These
checks support the displayed formulas; the induction above supplies the
all-index argument. No predecessor census was repeated.

The initial draft had a cosmetic stray comma in the exponent of equation
(3); removing it has no effect on this verdict. No mathematical correction
is required by this audit.
