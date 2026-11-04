# The Hamilton-bounded cone has the domination sign through order seven

**Proved, with an exact finite algebraic certificate:** for every integer n with 1<=n<=7, every profile

    q=t^(n)+sum_(r<n, r=n mod2) w_r t^(r),
    w_r>=0,  1+sum_r w_r<=n!,

has D_h(q)<=0 for every h>=2. Here t^(r)_k=k!S(r,k), and D_h is the exact ordered-cover contrast from `../covered_ambient/PROOF.md`. The weights may be arbitrary real numbers. This is not merely an enumeration of integer profiles, and no all-order assertion is made.

The scalar h=1 sign is already strict for every nonzero nonnegative finite profile by the companion positive-square identity. The order-eleven example in `../profile_cone_hbound/PROOF.md` shows that the present bounded-order statement cannot simply be extended to all n using the stated hypotheses.

## Certificate construction

Write L=n!-1 and list the d-1 lower generator weights. For L>0 put w_i=L x_i and introduce a slack coordinate x_d=1-sum_(i<d) x_i. The Hamilton bound is exactly the standard simplex x_i>=0, sum_i x_i=1. The largest generator has weight one.

The formula

    D_h(q)=(h/3)L_q[K E0^(h-1)/(1-XYZ)^(h+1)]

is a homogeneous cubic in the generator weights before the top weight is set to one. A term containing p lower factors becomes a degree-p monomial in the lower x variables, multiplied by L^p. Homogenize it to degree three by multiplying by (sum_i x_i)^(3-p), which equals one on the simplex. The result is a homogeneous cubic in all d nonnegative simplex coordinates.

The exact certificate stored in the receipt lists every coefficient of this cubic for each n=1,...,7 and h=2,...,3n+2. **Every listed coefficient is nonpositive.** Therefore the cubic, and hence D_h(q), is nonpositive throughout the whole simplex. For n=1, L=0 and there is only the fixed top generator; the same homogeneous calculation applies with a single slack coordinate equal to one.

The checker obtains each original cubic coefficient by the exact symmetric polarization of K E0^(h-1), multiplied by the monomial multiplicity 3!/prod(m_i!). Its denominator expansion uses

    (1-XYZ)^(-h-1)=sum_(k>=0) binom(k+h,h) (XYZ)^k.

All retained generator profiles vanish above n. It is therefore valid to discard numerator monomials with any exponent above n, since all factors have nonnegative exponents. At h>3n+2 every numerator monomial has total degree at least h+2>3n, so all contributions vanish; the verifier also checks the retained numerator has vanished at the boundary. Thus every higher coordinate is covered.

This is an independently replayable integer-polynomial certificate. It neither assumes log-concavity nor estimates a polynomial from sampled values. The receipt retains zeros as well as negative coefficients. The ordinary and optimized runs must agree.

## Scope

The result concerns a sufficient domain for the sign of one complete ordered-cover contrast family. It does not classify realizable tournament profiles or prove named-response reconstruction at these orders. The actual full named response is not replaced by a reduced response to assert a collision. Outside the stated continuous domain and order range, this certificate makes no claim.

Run `python verify.py > receipt.json` and `python -O verify.py > receipt.optimized.json`. The checker uses only the Python standard library and contains its own sparse polynomial multiplication and Stirling recurrence; it imports no producer modules.
