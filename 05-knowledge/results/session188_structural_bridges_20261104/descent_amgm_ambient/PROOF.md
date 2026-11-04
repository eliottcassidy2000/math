# Weighted AM-GM extends palindromic-descent ambient separation through order sixteen

**PROVED by the exact certificates retained here:** the common-profile cyclic-child arbitrary-fixed-ambient theorem holds at common child orders thirteen through sixteen. Together with `../descent_profile_cone/PROOF.md`, it holds through sixteen. No inference above sixteen is made.

## 1. The inherited cone and response formula

For fixed n, let `g_(n,d)` be the paired binomial descent generators in `../descent_profile_cone/PROOF.md`, for `0<=d<=floor((n-1)/2)`. Every actual tournament ordered-cover profile is `q=sum_d w_d g_(n,d)` with `w_d>=0`: the coefficient counts permutations with d backward adjacencies, word reversal makes the distribution palindromic, and an ordered k-cover is exactly a word with k-1 cuts containing every backward adjacency. Formal generators need not be realizable by tournaments.

For the actual root masks `a=A1,c=A2`, the complete augmented-core ordered-cover difference is the homogeneous cubic

    D_h(w)=(h/3)L_q[K E0^(h-1)/(1-XYZ)^(h+1)],

with the exact symmetric kernel and coefficient multiplicities specified in that packet. This formula counts all physical two-port interactions. The first-coordinate sign is separately known: `D_1(q)=-P(q)<0` for every nonzero finite nonnegative q.

The raw coefficient-sign test succeeds through twelve but has positive mixed coefficients at larger orders. That does not imply positivity of the cubic on the nonnegative cone. This packet gives a stronger sufficient certificate using its negative pure cubes.

## 2. An exact weighted-cube certificate

Write any homogeneous cubic as

    D(w)=-sum_i A_i w_i^3 + sum_(m in M+) c_m product_i w_i^(m_i)
                               + remaining nonpositive monomials,

where `A_i>=0`, every `c_m>0`, and every exponent vector m has total degree three. Pure positive cubes are not allowed by the certificate. Choose positive rational scales d_i. Weighted AM-GM, with multiplicity m_i, gives

    product_i w_i^(m_i)
       <= (product_i d_i^(m_i))/3
          * sum_i m_i w_i^3/d_i^3.                    (1)

For example, a repeated-index monomial uses the weights 2/3 and 1/3, not three independent variables. The usual arithmetic-geometric mean inequality also covers the three-distinct-index case. All variables are nonnegative, so no sign reversal occurs.

Consequently D is nonpositive on the whole nonnegative orthant if, for each i,

    sum_(m in M+) c_m m_i product_j d_j^(m_j)
          <= 3 A_i d_i^3.                            (2)

This is a finite integer check when the d_i are positive integers. Equation (1) assigns the left-side load in (2) to its negative cube; summing the loads and leaving all other nonpositive monomials in place proves the assertion. No numerical optimization or approximate inequality is needed to verify a certificate.

## 3. Complete exact application at orders thirteen through sixteen

The verifier computes all homogeneous coefficients of every `D_h` for `h=2,...,3n+2` and every n in `{13,14,15,16}`. Coordinates not listed below have every coefficient nonpositive. The only coordinates requiring (2), with generator indices in increasing order, are:

- n=13, h=2: `d=(1,2,4,8,8,16,32)`.
- n=14, h=2: `d=(1,2,4,8,8,16,32)`; h=3: `d=(1,2,4,8,8,16,16)`.
- n=15, h=2: `d=(1,2,4,8,16,16,32,64)`; h=3 and h=4: `d=(1,2,4,8,8,16,32,64)`.
- n=16, h=2: `d=(1,2,4,8,16,16,32,64)`; h=3 and h=4: `d=(1,2,4,8,8,16,32,32)`.

Every scale, complete coefficient list, negative cube, positive monomial, left-side load, right-side capacity and nonnegative slack is retained in the receipts. All arithmetic is integral. The known finite support of each generator permits truncation only above exponent n. After the displayed cover range, the kernel has no support in that box, so no higher coordinate remains unchecked.

Thus every `D_h(q)<=0` for h>=2 on the entire nonnegative palindromic descent cone at these four orders. This statement does not require the mass `sum c_d=n!`, the Stirling cone, or an extremal Hamilton-path estimate.

## 4. Actual named-response consequence

Use the source cyclic order `A1->A2->A3->A1` and the opposite target quotient on the same named child sets, with each source and target child retaining the common complete profile q of common order n. For any fixed rooted exterior Q, give both physical ports the root's incidence to each exterior vertex, as specified in `../equal_small_ambient/PROOF.md`. The augmented core is then a common module, and its Hamilton difference is

    sum_h b_h D_h(q) <= b_1 D_1(q) = -H(Q)P(q)<0,

where every b_h is nonnegative and `b_1=H(Q)>0`. This proves separation by one actual physical query in any fixed ambient, including arbitrary strong exteriors. The new coefficient proof changes none of the physical masks or restoration hypotheses.

The result remains a common-complete-profile cyclic-quotient theorem. It is not an unequal-child result, arbitrary-prime reconstruction, or a claim that positive coefficients imply a collision. The earlier order-thirteen positive mixed coefficient is correctly absorbed here, rather than being relabelled as a false theorem boundary.

## 5. Reproduction and provenance

Run `python verify.py > normal.json` and `python -O verify.py > optimized.json`, then compare. Only the Python standard library is needed. The checker imports no sibling implementation. Its scale lists are explicit constants and every resulting inequality is checked exactly, with checks retained under optimization.

Weighted AM-GM and this sufficient negative-cube allocation criterion are elementary; no novelty claim is made for either. The contribution is the stated complete exact certificate and its consequence via the separately proved response formula, strict scalar identity and physical module lift.
