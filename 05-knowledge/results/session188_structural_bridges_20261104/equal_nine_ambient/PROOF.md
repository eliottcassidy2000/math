# Equal-profile cyclic children through order nine: strict separation in arbitrary ambient tournaments

**PROVED, in the certified range1<=n<=9:** suppose three nonempty named tournaments have the same complete ordered-cover profile and common order n. Reversing their cyclic quotient is detected by the explicit block-uniform two-port query in every fixed ambient quotient, even if the target child interiors are changed in any profile-preserving way. This extends the separately frozen order-seven corollary. It makes no claim beyond order nine.

## 1. An elementary Hamilton bound

Every order-n tournament satisfies

    H(A) <= n! / 2^floor(n/2).                            (1)

Indeed, partition the n! vertex orderings into orbits under independent swaps of positions(1,2),(3,4),...,(2r-1,2r), where r=floor(n/2). Each orbit has2^r elements. The underlying unordered pair at each designated pair of positions is fixed throughout an orbit. At most one ordering in the orbit orients every designated adjacency forward, since each pair has exactly one forward order. Every directed Hamilton ordering must satisfy all those necessary adjacencies. Thus each orbit contains at most one Hamilton path, proving(1).

This is an elementary counting bound, not a claimed new extremal theorem or a claimed best bound. Historical and stronger known bounds were checked before publication; see PROVENANCE.md. The present contribution is its use in an exact, bounded-order response certificate.

## 2. The new continuous certificates at orders eight and nine

Write t^(s)_k=k!S(s,k), with initial coordinate zero. As in the earlier audited cone argument, an actual tournament profile q is a same-parity nonnegative Stirling mixture

    q=t^(n)+sum_(s<n, s=n mod2) w_s t^(s),
    w_s>=0,        q1=1+sum_s w_s.

The top coefficient is one by the classical odd-cycle-cover expansion. Bound(1) therefore puts the lower weights in a simplex of total weight at most

    L8=8!/16-1=2519,
    L9=9!/16-1=22679.                                   (2)

Let D_h(q)=q_h(T+)-q_h(U+) be the complete ordered-cover contrast for the physical core query a=A1,c=A2, with source quotient A1->A2->A3->A1 and target quotient reversed. The exact formula in `../covered_ambient/PROOF.md` is

    D_h(q)=(h/3)L_q[K E0^(h-1)/(1-XYZ)^(h+1)].           (3)

Here K,E0 and L_q have the same definitions as in the earlier certificate.

Set each lower weight w_i=L_n x_i and add a nonnegative slack coordinate so that sum_i x_i=1. Homogenize each term of lower-weight degree p to degree three by multiplying it by(sum_i x_i)^(3-p). The resulting homogeneous cubic equals D_h on the simplex. The complete exact coefficient lists retained by `verify.py` show:

- At n=8, every coefficient is nonpositive for h=2,...,26.
- At n=9, every coefficient is nonpositive for h=2,...,29.

These are53 coordinate polynomials with1,480 retained integer coefficients, including zeros. The checker computes the coefficients directly by symmetric trilinear expansion and the exact denominator series. It does not fit a polynomial to samples. Therefore

    D_h(q)<=0 for every h>=2                             (4)

on the entire continuous domain in(2), not only for integer or realizable weights. Coordinates beyond3n+2 vanish by degree/support. For n<=7, the earlier `../profile_cone_bounded/PROOF.md` certifies the larger domain q1<=n!, so(4) already holds.

For every nonzero finite nonnegative q the independent scalar positive-square theorem gives D1(q)=-P(q)<0. Consequently every actual equal-profile cyclic triple with common order at most nine has complete augmented-cover domination, strict at its first coordinate.

## 3. The precise arbitrary-ambient consequence

Let M=A1 union A2 union A3 be a distinguished module substituted into the root r of an arbitrary fixed tournament Q. All three source children have common complete profile q and size1<=n<=9. The target reverses their cyclic quotient on the same names and retains profile q in each child; their actual interiors need not be unchanged. All exterior vertices and arcs are fixed.

Use x->y. On M put a=1_(A1),c=1_(A2), with a_v=[v->x] and c_v=[v->y]. Give x and y the root's exterior incidences: z->r means z->x,z->y, and r->z means x->z,y->z. Thus M union{x,y} is a common augmented module with the same collapsed quotient Q.

The exact module-run formula has nonnegative coefficients c_h(Q,r) and c1(Q,r)=H(Q)>0. The actual source-minus-target response gap is therefore

    sum_h c_h(Q,r)D_h(q)
      <= H(Q)D1(q)
      =-H(Q)P(q)<0.                                    (5)

The query is an actual named attachment on the whole ambient tournament. Its exterior mask bits are those just specified, not uniformly zero by default. The proof works for arbitrary strong or non-strong Q of any size and for odd or even child orders. No modular nonvanishing is required.

This is the exact same lifting mechanism as `../equal_small_ambient/PROOF.md`; only the certified child-order range is enlarged. The earlier packet's physical count controls remain valid and need not be repackaged as new controls here. The new numerical work is the complete coefficient certificate at n=8,9, not a literal enumeration of all27-name cores or large ambient tournaments.

## 4. Certified boundary and provenance

The simple coefficient-sign test with bound(1) does not certify n=10: already h=2 has positive homogenized coefficients. This is only failure of that sufficient certificate. It is not a counterexample to actual domination, actual separation, or even positivity on the entire constrained domain. No order-ten result is claimed here.

The adjacent-swap inequality, the classical profile cone, homogeneous-simplex nonpositive-coefficient reasoning, and module-run substitution are not claimed as new mathematical methods. The separately replayable application to these complete response coordinates is the result certified here. No historical-priority claim is made.

## 5. Reproduction

Run `python verify.py > normal.json` and `python -O verify.py > optimized.json`, then compare. The standard-library checker is self-contained and imports no producer modules. It retains every coefficient for both new orders, uses arbitrary-size integer arithmetic and explicit exception checks, and binds this proof, its provenance note and source by SHA256. Dependencies for the mathematical corollary are the previously audited exact cover identity, scalar positive-square theorem and small-child ambient substitution proof.
