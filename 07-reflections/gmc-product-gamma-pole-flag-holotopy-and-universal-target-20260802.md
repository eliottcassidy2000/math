# Product-Gamma pole-flag holotopy and the universal target

**Status:** research reflection.  The finite theorem referred to below is
proved THM-3120.  This reflection and its post-theorem experiments are not
proved dependencies.

## Inheritance pass

- **Closest proved mechanism:** THM-3110 gives the exact 24/25 residual
  alphabet banks, positive degree-five exit, and finite Schur face.
- **Closest positive holotopy:** THM-3115 transports low-degree normalized
  currents through the Young-fibre refinement order, now through degree nine.
- **Canonical near miss:** ordinary numerator coefficients, partial-fraction
  residues, and independent cycle coordinates all have both signs.
- **Least-used sidecar:** pole removal is plethystic virtual-letter subtraction,
  which commutes algebraically with labelled deletion even though it is not a
  positive physical deletion on individual bank atoms.

## Anchor / Niche / Wildcard

The **anchor** is the all-degree one-row response.  THM-3120 replaces a long
coefficient scan by a finite pole-prefix Newton flag.  On the canonical 115
supports, all `8,241` flag coefficients are strictly positive, hence every row
coefficient in every degree `n>=5` is positive.

The **niche** is the mixed pole-prefix/Young-refinement bicomplex.  The flag
uses only prefixes through the numerator degree `d`; longer prefixes are not
part of the certificate and can destroy Hasse positivity.  At `(1,2)`, an
independent exact probe finds all active prefixes `j<=d=2` Hasse-positive,
while the first failure is `j=4` at degree five.  The tempting minimal-active
repair is nevertheless false globally: a subsequent exact 115-support scan
finds its first active failure at `I_2,(a,b)=(1,3),N=5,j=d=5`.  Thus scalar
Newton positivity survives a genuinely nonpositive Hasse current.

The **wildcard** is a discrete B-spline interpretation of reciprocal numerator
division.  The coefficients

```text
Q[r_1,...,r_(j+1)]
```

are generalized Bernstein/Newton coordinates between two already understood
ends: the leading coordinate is the degree-five response, while the constant
Newton coordinate at the other end is the surviving dominant-pole residue.
The open problem is to make the intermediate coordinates into a positive
rooted-current interpolation rather than prove them by coefficient algebra.

## Exact object board

| Object | Representation | Preserved datum | Lost datum | Cheapest next test |
|---|---|---|---|---|
| row response | reduced rational `t^5P/D` | every row coefficient | nonrow partitions | prove flag for all `a,b` |
| pole flag | confluent Newton chain of `Q=x^dP(1/x)` | ordered top-pole geometry | atomwise physical letters | classify active prefixes by chamber |
| suffix decomposition | positive sum `t^(5+k)/(D/q)` | coefficient positivity | normalized row minimum | couple each summand to THM-3115 gaps |
| virtual deletion | alphabet substitution `X -> X-M` | exact deletion commutator | Hasse positivity | test only prefixes `j<=d` |
| cycle current | independent `z_k` coordinates | labelled power-sum data | ordered pole sidecar | retain pole order before cycle collapse |
| cancellation seam | `(a,b)=(2k,3k)`, even `k` | exact local residue parity | raw-envelope top pole | quotient first, then form flag |

## Why the Newton flag is a holotopy

Let the reduced poles be `r_1>=...>=r_E`, let

```text
q_j(t)=product_(ell<=j)(1-r_ell t),
P(t)=sum_(k=0)^d c_k t^k q_(d-k)(t),
Q(x)=x^dP(1/x).
```

Then

```text
Q(x)=sum_(j=0)^d c_(d-j) product_(ell<=j)(x-r_ell).
```

The chain starts at

```text
c_0 = [t^5]F = Phi(h_5),
```

which belongs to THM-3110's globally positive degree-five face.  It ends at
the local top-pole evaluation `c_d=Q(r_1)`, whose positivity is the finite
residue shadow of the positive asymptotic endpoint.  Intermediate `c_k` are
successive rooted divided differences.  They interpolate the low-degree
collision exit to the dominant-pole tail while remembering every repeated
pole seam.  This is more structured than a recurrence barrier and more
faithful than a continuous derivative test: `Q` and its derivatives can be
negative at lower flag nodes even when every divided difference is positive.

The exact algebraic square is

```text
labelled deletion  o  virtual pole subtraction
 = virtual pole subtraction  o  labelled deletion.
```

The square does not itself lie in the positive cone because the pole can be
absent from most atoms.  At `(1,2)`, pole `5` is absent from `21/24` `I_1`
alphabets and `21/25` `I_2` alphabets.  The plausible repair is therefore not
to force atomwise deletion, but to ask whether the *bank-summed active square*
is the boundary of a positive Young-refinement current.

## Enlarged exact evidence beyond the theorem bank

After the canonical THM-3120 run, a separate optimized hostile scan tested

```text
11<=a<=20,
a<b<=min(3a+4,41),
both banks.
```

These are `500` additional support/bank cases and `43,500` additional flag
coefficients, with numerator degrees reaching `137` and `138`.  Every flag
coefficient was strictly positive.  The only denominator cancellations were
the predicted even-parity points `(12,18)`, `(16,24)`, `(20,30)`, once in
each bank.  This is exact integer evidence, but it is deliberately not added
to the finite theorem statement or stored as a canonical transcript.

## Failure boundaries

1. **Ordinary monomials:** `P_1(1,2)=36-108t-72t^2`.
2. **Overlong prefix:** `[t^17](1-5t)(1-4t)F_1=-5,901,696`.
3. **Independent cycles:** the `z_1z_2z_3` coefficient is `-84` in both
   banks at `(1,2)`.
4. **Continuous Peano shortcut:** positivity of derivatives on the node
   interval is false in almost every tested case; the phenomenon is genuinely
   discrete/confluent.
5. **Raw top pole:** the root `2a` cancels on `b=3a/2` exactly when `4|a`.

Each failure removes information the successful flag retains: respectively
pole order, minimal active length, coupling between cycles, node
multiplicity, and reduced-versus-raw pole identity.

## The next decisive experiment

The active-prefix experiment has now selected its second outcome.  Across
`8,241` prefixes and degrees `5..9`, the first active failure is

```text
I_2, (a,b)=(1,3), N=5, j=d=5,
G_(5)=-2,324,160,
maximum positive refinement flow=1,572,480,
negative demand=3,896,640.
```

The coarsest partition already carries unpayable debt.  This is currently an
independent exact incoming result rather than a dependency of this reflection;
its full canonical package should be cited once promoted.  Conceptually it
answers the important question now: the pole flag is a scalar positive
holotopy, but not a prefixwise lift of THM-3115's Hasse-positive holotopy.
The next sidecar must quotient or pair the coarsest debt rather than simply
intersect the two positive cones.

For arbitrary supports, the sharp theorem target is equally finite:

```text
Q_i[r_1,...,r_(j+1)]>0
for 0<=j<=deg(P_i), after exact pole cancellation.             (T)
```

A proof of `(T)` by chamberwise rooted-current transport would connect
THM-3110's degree-five collision face, THM-3115's positive refinement
holotopy, and the all-degree dominant-pole tail in one object.
