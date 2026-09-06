# The sharp first-order finite-length penalty for signed-root stability

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
For lists of total length at most N, the infimum of the
sharp signed-root stability quotient is

    K3+C/N+o(1/N),
    C=-22-(16/3)sqrt(6)+10sqrt(2)+(40/3)sqrt(3)
     =2.1722010964723645447... .                     (1)

The three-equal-positive-atom family with equal negative dust attains this
asymptotic coefficient. This does not assert that it minimizes the quotient
at any specified finite N, or that a finite-dimensional infimum is attained.

## 1. Domain, inheritance, and the ordered coordinate

For a finite real list r let p_j=sum r_i^j and suppose

    p1=p2=1, E=(1-p4)/2>0.

Let a>=b>=0 be the two largest positive coordinates, padded with zero,
and use

    J=(5-8p3+3p4)/(3(1-p4)),
    c_*=(13-8sqrt(2))/3,
    d2=2-sqrt(2)(a+b),
    R=(J-c_*)/d2,
    K3=4sqrt(3)/[3(1+sqrt(2))(1+sqrt(3))].

The actual normalization makes d2>0: equality would leave exactly two
positive atoms 1/sqrt(2), whose sum is not one. The universal R>=K3 is
**THM-4454**,
`01-canon/theorems/THM-4454-sharp-global-signed-root-duplication-stability.md`.
The second input is the independently audited
[minimizing-sequence classification](long_frontier_sep06_minimizers.md):
R tends to K3 exactly when the three leading positive atoms tend to
1/sqrt(3) and the remaining square mass tends to zero. A reserved theorem
number is not used as a proved dependency.

Define I_N as the infimum of R over this domain with total list length at
most N; zero coordinates may be included or removed. Then

    lim_{N->infinity} N(I_N-K3)=C.                   (2)

The live concept board consists of three-atom concentration, signed tail
first moment, Cauchy saturation, ordered top-two distance, uniform local
Taylor error, and total dimension. The closest proved mechanism is the
classification just cited. Its mixed-sign dust hostile shows that separate
positive and negative tail first moments cannot be inferred from norm
convergence alone. The corrected near miss here is to ignore the ordered
top-two term and treat leading-atom splitting as quadratic. It has a
strictly positive linear cost. Section 3 retains that missing coordinate.

The connection maps norm concentration to a finite-length lower bound by
retaining the signed sum of the whole tail. Square-norm convergence alone
destroys first-moment information; p1=1 and the actual tail count are the
required sidecars. Cauchy--Schwarz preserves the signed first moment but
forgets the detailed tail arrangement; its equality defect records precisely
what must vanish for asymptotic sharpness.

## 2. A length-uniform local expansion

Write z=1/sqrt(3), u=sqrt(2), d0=2-2uz. Near the classified minimizer,
let a>=b>=c>0 be the three largest positive coordinates, and define

    ell=(a+b+c)/3,
    m=sum_{remaining i} r_i^2,
    h=ell-c>=0,
    v=(a-ell)^2+(b-ell)^2+(c-ell)^2.

These symbols refer to actual leading coordinates; no sign pattern is
imposed on the remaining list. If alpha=a-ell and beta=b-ell, then
alpha+beta=h and alpha>=beta>=-h. Consequently

    0<=v<=6h^2,
    3ell^2+v=1-m,
    a+b=2ell+h.                                    (3)

The following expansion is uniform over all finite list lengths and tail
sign patterns as m+h tends to zero:

    R-K3=A*m+B*h+O(m^(3/2)+h^2),                   (4)

where

    A=-2-sqrt(6)/3+2sqrt(2)+(7/3)sqrt(3)
     =4.0533824281458444164...,
    B=K3*u/d0
     =8-5sqrt(2)-4sqrt(3)+3sqrt(6)
     =1.3491981862085498765... .                    (5)

Both coefficients are strictly positive. For A, the rational bounds
sqrt(2)>7/5, sqrt(3)>17/10, sqrt(6)<5/2 already give a positive lower
bound. For B this follows immediately from K3,u,d0>0.

To prove uniformity, the tail bounds are

    |sum_tail r_i^3|<=m^(3/2),
    sum_tail r_i^4<=m^2.

By (3),

    ell=z-(z/2)m+O(m^2+h^2),
    p3=z-(3z/2)m+O(m^(3/2)+h^2),
    p4=1/3-(2/3)m+O(m^2+h^2),
    d2=d0+uz*m-u*h+O(m^2+h^2).                     (6)

For instance, the cubic leading-atom sum is
3ell^3+3ell*v+sum deviations^3; its nonradial terms are O(h^2).
The fourth moment follows in the same way. Every constant in these
estimates is absolute in a fixed neighborhood, independent of length.

The energy and distance denominators have nonzero limits 1/3 and d0.
Taylor expansion of the rational quotient therefore gives (4). Explicitly,

    J0=3-4z,
    J=J0+(10z-4)m+O(m^(3/2)+h^2),
    A=[10z-4-K3*uz]/d0,
    B=K3*u/d0.

Mixed errors mh and m^2 are absorbed into O(m^(3/2)+h^2) in a bounded
neighborhood. In particular there are fixed positive constants and a fixed
neighborhood on which

    R-K3 >= (A/2)m+(B/2)h.                         (7)

The term B*h is the key ordered effect: the top-two sum is 2ell+h,
whereas the moments are symmetric in all three leading atoms.

## 3. The cheap splitting test and the exact upper family

For n negative coordinates and prescribed leading deviations with squared
sum v, the positive mean satisfying both normalizations is

    ell=[1+sqrt(n((n+2)-(n+3)v)/3)]/(n+3),          (8)

provided the radicand and the leading coordinates have the stated signs.
Take the tail to consist of n copies of -(3ell-1)/n. Its signed sum is
1-3ell and its square mass is (3ell-1)^2/n, so (8) gives p1=p2=1.

The first hostile to discarding the ordering coordinate takes either
leading deviation pattern

    (h/2,h/2,-h),      or      (2h,-h,-h),
    h=lambda/n, lambda>0 fixed.

For every fixed lambda these are actual normalized lists for all
sufficiently large n. Their leading roots remain positive and above the
negative dust. Formula (4), or the independent exact derivatives in the
source, gives

    (n+3)(R-K3) -> C+B*lambda.                     (9)

Thus the leading splitting cost cannot be omitted or charged only as h^2.
Neither pattern improves the equal-three coefficient.

For the upper bound set v=h=0 in (8), obtaining three equal atoms

    a_n=[1+sqrt(n(n+2)/3)]/(n+3)

and n equal negative entries -(3a_n-1)/n. For n>=1 the positive root of
3a^2+(3a-1)^2/n=1 lies in (1/3,z); the list is eligible. Its total length
is n+3, rather than n. As n tends to infinity,

    m=(3a_n-1)^2/n=(sqrt(3)-1)^2/n+O(n^-2),
    h=0.

Here the actual tail cubic and fourth moments are O(n^-2) and O(n^-3),
so expansion of this explicit analytic family gives the stronger upper
estimate

    R=K3+A*(sqrt(3)-1)^2/n+O(n^-2)
     =K3+C/(n+3)+O((n+3)^-2).                     (10)

The exact radical identity C=A*(sqrt(3)-1)^2 yields (1). This proves
limsup N(I_N-K3)<=C and supplies an upper family for every N>=4.

## 4. Lower bound for arbitrary approximate minimizers

Choose any sequence N_j tending to infinity and eligible lists of length
at most N_j with

    R_j <= I_{N_j}+N_j^-2.

Such approximate minimizers exist by the definition of an infimum; no
attainment claim is needed. The upper family and THM-4454 imply R_j tends
to K3. The classification forces their three leading positive atoms to
tend to z and their remaining square mass to vanish. Thus h_j tends to
zero and the uniform expansion applies. The upper bound together with
(7) gives

    m_j+h_j=O(1/N_j).

The error in (4), multiplied by N_j, is consequently o(1).

Pad each list with zeros to total length N_j. Its remaining N_j-3
coordinates have signed sum 1-3ell_j, irrespective of their separate
positive and negative first moments. Cauchy--Schwarz gives

    (N_j-3)m_j >= (3ell_j-1)^2.                    (11)

Since ell_j tends to z,

    liminf N_j*m_j >= (sqrt(3)-1)^2.

Now (4), h_j>=0, and A,B>0 imply

    liminf N_j(R_j-K3)
       >= A*(sqrt(3)-1)^2=C.

Together with (10) and the arbitrary subsequence, this proves (2).
The proof does not replace the total length by the number of negative
entries, nor require all tails to be negative. The equal-negative upper
family simply saturates the signed Cauchy inequality asymptotically.

## 5. Exact asymptotic equality conditions

For eligible lists of length at most N, with N tending to infinity,

    N(R-K3) -> C

holds if and only if, using their three leading positive atoms,

    Nm -> (sqrt(3)-1)^2,      Nh -> 0.              (12)

For necessity, the classification and (7) again give m+h=O(1/N).
Then (11), (4), and positivity of A,B force both limits in (12).
For sufficiency, (12) and (3) force ell to z and the local expansion
immediately gives the stated quotient limit.

If S=1-3ell and the padded tail has N-3 labels, its exact Cauchy defect is

    sum_tail (r_i-S/(N-3))^2 = m-S^2/(N-3).

Under (12) this is o(1/N). This is the retained equality sidecar; it does
not say that any finite minimizing tail is exactly constant. In particular
the total length actually used, and the number of negative coordinates,
are both N-o(N): every padded zero or nonnegative tail contributes at
least (S/(N-3))^2 to this defect, while S tends to 1-sqrt(3)<0.

There is also first-moment rigidity at this sharper scale. Cauchy gives

    sum_tail |r_i-S/(N-3)|
       <=sqrt(N-3)*sqrt(m-S^2/(N-3))=o(1).

The reference tail is negative for large N. Therefore the positive tail
first moment tends to zero and the negative tail magnitude tends to
sqrt(3)-1. Root observed this consequence during audit. Plain convergence
R->K3 permits divergent positive and negative dust moments; sharpness at
the first finite-length order rules that behavior out.

## 6. Reproduction and scope

[Source](../../04-computation/long_frontier_sep06_dimension.py),
[output](long_frontier_sep06_dimension.out). The universal lower argument is
the analytic expansion, the inherited classification, and Cauchy--Schwarz.
The source checks exact symbolic identities through two separate derivative
paths: a direct three-atom gradient and the moment/ordered-distance
expansion. It also checks the normalized equal and split families and
their exact first-order derivatives. All algebraic identities are symbolic;
positive constants additionally have exact rational radical bounds.

The finite control universe consists of n=10,100,1000,10000, with each
of the three deviation patterns (0,0,0), (h/2,h/2,-h), (2h,-h,-h), using
h=1/n for the last two. A separate rational interval implementation checks
their actual moments, denominators, and quotient values with outward
rounding. These twelve lists are controls, not an optimization census.
The finite comparison that each named split costs more than its equal
control is not substituted for the asymptotic proof.

```
python3 -B 04-computation/long_frontier_sep06_dimension.py
python3 -B -O 04-computation/long_frontier_sep06_dimension.py
```

Normal and optimized outputs are byte-identical, with 117 explicit gates.

* Source SHA256: `badb2489aeb8e4a27a0e5b61c203323dfb988be693488637056a1f9f3c19e667`.
* Output SHA256: `d2e72609bb6b248bed6d721e75e0a2b0350316b046cd6bed5ad06ddae8973e30`.
* Semantic digest: `15251e617e95cae1f478f4ec83c7ce9b884fde71c9521555a5c78f4c65f21c88`.

The [independent referee](long_frontier_sep06_dimension_audit.md) passed the
complete analytic proof, all length and sign quantifiers, the ordered
splitting term, approximate-infimum argument, and asymptotic equality
conditions, including the final first-moment consequence. The source and
normal, optimized, and frozen outputs and their exact hashes were independently
replayed: PASS. Root also read the complete proof and checked the added
first-moment deduction. All three primary artifacts are frozen.

No finite-N equality theorem or actual Laurent-row transport is inferred
from this abstract normalized moment result.
