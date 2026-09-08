# Independent audit of the complete fixed-DG quartic exclusion

**Status: PASS — complete proof compilation, source and both replays.**

Auditor: root, independently of producer `three_ray_geometry`.
The [primary compilation](planar_jc48_sep08_quartic_closure.md) was
reviewed in its frozen RESERVED state. No case is still pending.

## 1. Accepted conclusions and the exact boundary

All three statements are accepted, with their distinct hypotheses:

1. For every global `H in L2`, `L in L1` with `deg_t H=2`, the
   polynomial `H^2+L` has no polynomial constant-Jacobian mate of any
   degree. The proposed mate need not extend to the surface.
2. On the specified `W=(P1 x P1) minus {z=x^2}`, there is no pair
   of global functions with nonzero constant source Jacobian and one
   coordinate of source fibre degree four. Both coordinates are global
   in this statement; that condition pays the approximate-root gluing.
3. Every nonzero constant output direction of a hypothetical global pair
   avoids `L4`. The already proved degree-at-most-three exclusion is
   used only for a member of that actual output pencil.

This is a theorem about the full stated section spaces on the fixed
surface. There is no reduction of arbitrary planar Keller maps to this
surface, no arbitrary quartic polynomial classification, and no JC(2)
closure. Rational mates in the square-prefix class remain possible.

## 2. Exhaustiveness before applying a case theorem

I checked the full fifteen-parameter quadratic and six-parameter linear
section rows against the proved filtration. The leading coefficient of
`H^2+L` is `N^2`, with nonzero `N` of degree zero through eight.
The general leading-coefficient theorem therefore gives the necessary
exact differential `dx/sqrt(N)` in the actual generated field. A square
`N` uses `C(x)`, not a fictitious connected quadratic cover.

The degree-eight exactness supplier is a proved complete classification,
not a finite-position experiment. Its seventeen finite partition rows
and four labelled position conditions agree with the compilation. I
also generated all partitions independently by descending first part
and remainder, without importing or executing the producer's generator.
This gives all 67 partitions, and its intersection with the literal
classification table is exactly the seventeen displayed rows. The source
itself has a different multiplicity-count enumeration and an independent
Euler-product count. All three agree degree by degree.

The four position-sensitive loci retain their hypotheses: finite distinct
roots, the labelled quadruple/sixfold/triple roots, and the nonzero scalar
and squarefree residual cubic. The C611 example has polynomial
`3k^2+2k+3`, coprime to `k(k-1)`; the C431 example keeps the triple
root at 3 and the simple root at -1. The sparse C5111 residual cubic
has nonzero discriminant. Failed position conditions already fail the
leading necessary differential and are not silently omitted.

The actual infinity order is exactly `8-deg(N)`, from `r^8N(1/r)`.
The same unmarked binary type can consequently require different proof
suppliers. This label is retained throughout; no projective change of
the differential or free swap with infinity is used.

## 3. The eighteen suppliers for the seventeen routes

I checked the statements and current status of every named case supplier.
All eighteen unique route suppliers have proved status, a separate PASS
audit, and exact byte/hash agreement with the inherited 383-artifact
manifest. These checks identify the accepted proof graph; they do not
claim that a matching hash independently reproves each antecedent.

The load-bearing scope checks are:

* Constant `N` uses the original infinity-octuple theorem. A finite
  octuple uses the all-position transport theorem, which explicitly
  includes its already proved `p=0` supplier.
* The finite/infinity `(4,4)` row is covered by the disjoint union
  `p=0` and `p!=0` theorems. The all-finite `(4,4)` partition is
  rejected by nonzero leading residues and is not confused with it.
* The `(5,3)` rows retain all three placements: finite triple,
  finite fivefold, and both finite. Their rational conclusions differ,
  but each gives precisely the polynomial exclusion used here.
* The `(7,1)` rows retain all three placements. The newly audited
  infinity-seven rational classification and both finite-sevenfold
  polynomial theorems have no missing finite-position parameter.
* For finite six/original infinity two, the universal degree-six theorem
  first forces constant boundary value already for a rational mate.
  `M(p)!=0` and `M(p)=0` then exhaust all lower sections. The former
  has the all-position translated M-unit rational exclusion; the latter
  has the complete shared-six polynomial theorem and its retained
  rational exceptions. Neither an active jet nor an unexplained
  constant-boundary assumption is imposed on the initial row.
* The `(5,2,1)`, `(3,3,2)`, `(4,2,1,1)`, `(6,1,1)`, `(5,1,1,1)`
  and `(4,3,1)` routes each have the stated full coefficient and
  boundary-placement coverage. The four positional exactness conditions
  are either paid within the supplier or used only as an already
  necessary restriction.

No case theorem assumes a bound on the degree of a polynomial mate.
Thus every possible leading row gives the actual contradiction required
for statement A, rather than merely a bound on an intermediate statistic.

## 4. The global quartic and output-pencil bridges

I reread the proved `quartic_boundary` theorem and its audit, particularly
Sections 1–4, rather than inferring globality of an approximate root from
its polynomial source formula. The second-chart equation is exactly
`J(Fbar,Gbar)=lambda*r^2`. Its extended finite-pole theorem explicitly
allows a nonzero base-only Jacobian with zeros, and its canonical root
is unique after fixing the leading coefficient. Consequently the
transformed source root agrees with the polynomial root on the whole
second chart. Both `H` and `L=F-H^2` are global, as required by A.
An independent bounded reread by `orthogonal_returns` agrees on this
interface and the necessity that both original coordinates be global.

The source's exact affine fibre-transport check supports uniqueness;
it does not purport to prove the finite-pole supplier by itself. The
globality-of-F-only hostile has a genuine pole in its canonical root
and is correctly retained. Applying A now contradicts the original
polynomial restriction of G, proving B.

For C, a nonzero coefficient row extends to an invertible constant
two-by-two matrix over C; no false assumption about a sum of squares
being nonzero is used. The determinant rescales the Jacobian. The first
coordinate cannot be constant. Degree four is covered by B and lower
degree by the inherited DG corollary with THM-2063, THM-2071 and
THM-2118, each cited with its full slug. At no point is an auxiliary
approximate root substituted into the original output pencil.

## 5. Source, hostile controls and reproduction

The full source has always-active checks and imports no producer from a
case theorem. Its universe, inherited filters, positional conditions,
infinity labels and analytic dependency boundary are explicit. The
binomial residue jets recover the three rational position conditions;
the elliptic operator recovers the sparse fourth condition. The opposite
finite `(4,4)` residues are tested individually.

The actual global rational hostile is checked in both full source charts,
has a nonconstant linear remainder, and has affine principal part
`1/(2x)` in its mate. It therefore blocks the tempting rational
strengthening of A–C without contradicting them. The separate globally
regular quartic with nonglobal root has its exact critical-line test,
so it is not mislabeled as a Keller example.

Independent replays:

    python3 -B 04-computation/planar_jc48_sep08_quartic_closure.py
    python3 -B -O 04-computation/planar_jc48_sep08_quartic_closure.py

Both exit successfully, reproduce **258 gates**, and match the frozen
1,387-byte output exactly. The third partition enumeration and the
eighteen complete supplier pin/status checks also pass.

Frozen producer pins before primary status promotion:

    source bytes 12137
    source SHA256 284fcd8a3f7c2970695ec6fe9d187bc896d07d56c350fd2a6e5523b029079433
    output bytes 1387
    output SHA256 11d0025901ef154426a703e6d501da97fc2e0fc282ff646a476bb9b1349ab253
    primary bytes 20203
    primary SHA256 ae22615103efab8d2224d578ae1426a491f53038735711b2f93e059755ce1b6a

No correction remains. All A/B/C statements are accepted for parent-owned
promotion, with the scope boundaries above preserved.
