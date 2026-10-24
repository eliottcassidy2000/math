# Polynomial weighted bounds from a complete occupancy band

Status: **PROVED**, conditional on the explicitly authenticated native record
interface below. The native twenty-source interface is **FINITE-EXACT** in
`../three_target570_20261024/{THEOREM,WEIGHTS}.md`. These are bounds for a
specified designed chain, not for arbitrary preassigned transition rates.

## 1. Carrier and statement

Retain the original full records, source weights, target-owned currents,
literal pins, normalized capacity units, kernel lifts and converse. Let R be
the existing twenty two-clock sources and E be m other compatible binary
sources with switch cost at most seven. Every option reservoir must be
disjoint from the others and from the immutable remainder. There are no
undeclared cross-source constraints. Charges lie in the original
G=(Z/2)^4 x Z/128, with the authenticated paired-current equivalence; this
does not transfer the harmonic bank to a new tournament order.

At each source i prescribe positive statistical weights a_i(U),a_i(V).
Physical electrical edge weights remain two. Divide out the constant U
weights on E and write w_i=a_i(V)/a_i(U). On R retain both actual weights:

    Z[r,g] = sum_{I subset R, |I|=r, ell(I)=g}
                 product_{i in I} a_i(V) product_{i not in I} a_i(U).

Fix the exact integer V count k and charge q. The visible state is the actual
outside V set J, with r=k-|J| and positive unnormalized mass

    M(J)=product_{i in J} w_i * Z[r,q-ell(J)].                (1)

Zeros are excluded, never filled by invented native witnesses. Conditional
on J, keep and sample the complete actual R fibre with its product weights.

Define, for the complete band 5<=r<=15,

    F[r,g]=Z[r,g]/binom(10,r-5),
    kappa=max F / min F.                                  (2)

All quantities are finite and positive, for any positive reservoir weights.
The following explicit reversible chains exist:

* **Interior:** on all full states with 5<=r<=15, a lazy chain of gap at
  least `1/(2 kappa binom(m+10,2))`. Every event changes at most 22 source
  records and costs at most 54 owner-clock occurrences. Outside weights are
  arbitrary positive numbers; the bound is independent of them.
* **Entire fibre**, if 5<=k<=m+15: put A_r=max_g Z[r,g], B_5=min_g Z[5,g],
  B_15=min_g Z[15,g], and let e_d denote the elementary symmetric polynomial
  in the m indicated actual odds, with e_d=0 when d>m. Set

      H=1+max( sum_{d=1}^5 (A_{5-d}/B_5) e_d(w),
               sum_{d=1}^5 (A_{15+d}/B_15) e_d(w^{-1}) ).   (3)

  A lazy chain on EVERY accepted full state has gap at least

      1 / (16 H kappa binom(m+10,2)).                      (4)

  Its events have the same 22-source / 54-clock budget. This is a
  source-specific bound for arbitrary positive weights. If K^{-1}<=w_i<=K,
  H is bounded by a degree-five polynomial in m, with constants depending
  only on K and the fixed reservoir weights. The gap bound is therefore
  inverse-polynomial of degree at most seven, uniformly in k and q.
* At k<=4 or m+20-k<=4, complete conditional refresh on the full accepted
  fibre costs at most eight sources / 56 clocks and has gap one (or gap
  one-half after lazification). Empty fibres are not states; singleton
  fibres have no nonconstant variance. Thus 22/56 is a uniform event budget
  for the whole collection of nonempty fibres.

No asymptotic claim of sharpness is made. The older 20/45 deterministic
routes remain smaller; the slightly larger budget here permits weighted
exchange at every complete reservoir grade.

For unweighted R the exact extrema in (2) are 67/252 and 16, so
`kappa=4032/67`. Equation (3) becomes

    H=1+9m+4 binom(m,2)+2 binom(m,3)+binom(m,4)+binom(m,5)

when the outside odds are also one. This uses the actual boundary
multiplicities, not a worst-case full-state weight ratio exponential in m.

## 2. The auxiliary exchange model and its coupling

Adjoin TEN AUXILIARY BITS of weight one to the m outside coordinates. They
are a proof representation, not extra native slots or kernel capacity. On
these D=m+10 labelled coordinates take subsets of size k-5 with product
weight, using weights w_i on outside coordinates. A visible set J has
`r-5` auxiliary occupied bits, hence reference mass

    rho(J) proportional to product_{i in J} w_i binom(10,r-5).

Choose an unordered coordinate pair uniformly. If precisely one is
occupied, resample its occupant with probabilities proportional to the two
coordinate weights; otherwise hold. This is a reversible pair heat bath.

Its gap is at least 1/binom(D,2), regardless of the weights. Here is a direct
proof, included so no log-concavity theorem is needed. Two adjacent subsets
differ at precisely coordinates i,j. Couple updates using the same pair and
the same uniform variable for the heat bath. Updating {i,j} coalesces them.
A pair disjoint from {i,j} updates identically. For a pair containing exactly
one discrepancy, one chain is frozen (zero or two occupied coordinates)
and the other either stays or moves the discrepancy. Their exchange
distance remains one. Thus expected distance contracts by exactly
`1-1/binom(D,2)` on adjacent pairs. Concatenate along shortest exchange
paths for arbitrary endpoints. This contracts the Lipschitz seminorm of
the Markov operator by the same factor. Applying it to any nonconstant
eigenfunction proves the asserted spectral bound. The heat-bath operator
is also an average of conditional-expectation projections. A one-state
fixed-size layer needs no spectral estimate.

Tilt the reference measure by F(J) from (2), by Metropolis acceptance
`min(1,F(J')/F(J))`, and then make the chain one-half lazy. The result
lumps under permutation of the auxiliary bits. Write H_0 for its visible
quotient. If f_min<=F<=f_max, then the un-lazified tilted energy is at least
`f_min/E_rho F` times the reference energy, whereas its variance is at most
`f_max/E_rho F` times reference variance. Hence the loss is exactly at most
kappa, NOT kappa squared. Restricting to visible functions cannot decrease
this lower bound. This proves the interior gap gamma>=1/(2 kappa binom(D,2)).

For implementation the nontrivial reference quotient proposals are:

* swap outside i in J with j outside J: `w_j/(w_i+w_j)/binom(D,2)`;
* add outside i: `(r-5) w_i/(1+w_i)/binom(D,2)`;
* remove outside i: `(15-r)/(1+w_i)/binom(D,2)`.

Apply the F acceptance and one-half factor; use the remaining probability
to hold. These are explicit rates on actual visible states. No arbitrary
tiny-rate factor is hidden in the statement.

## 3. Every tilted transition has a legal full-record lift

Given the new outside state, refresh from its complete conditional reservoir
law, including on visible holds. The exact identity in
`../three_target570_20261024/WEIGHTS.md` gives the same lazy spectral gap.
No minimum conditional witness weight is charged a second time.

An outside swap changes at most twenty reservoir sources and two outside
sources: support <=22, clock work <=2*20+7*2=54. A one-bit outside change
changes the reservoir grade by one, so at most nineteen R bits change:
support <=20 and work <=45. A hold refresh costs <=20/40. These estimates
hold for EVERY pair of reservoir representatives, not just selected section
witnesses. Fixed margins, residual charge, shared fixed pins and resource
disjointness hold at both arrivals. The event is simultaneous; no illegal
single-source intermediate is inserted.

## 4. Boundary mass is paid locally, not by a global weight ratio

For each visible boundary state x with r<5, define p(x) by removing the
lexicographically first 5-r selected outside sources, landing at r=5.
For r>15 add the first r-15 unselected sources, landing at r=15. Completeness
of the boundary bands authenticates the destination. The inequalities on k
ensure there are enough sources to change. Set p(c)=c for interior states.

Let L_c={x:p(x)=c}, and h(c)=sum_{x in L_c} M(x)/M(c). For a low-side leaf
with x=c union D, the mass ratio is at most

    (A_{5-|D|}/B_5) product_{i in D} w_i.

Sum over all subsets D of each size, which can only overcount the chosen
parent rule. The high-side bound uses inverse odds. Low and high leaves
attach to different grades, so `1<=h(c)<=H` with H from (3). This is the
source-preserving reason that only a fixed-degree symmetric polynomial
appears. In particular, bounding each entire source product by K^m would
discard the decisive local information.

Normalize pi(J)=M(J)/sum M on ALL visible states, let s=pi(interior), and
pi_C be its conditional interior law. Define a native visible kernel P:

* each leaf x moves to p(x) with probability 1/(4H);
* an interior c moves to its leaf x with probability M(x)/(4H M(c));
* for distinct interior c,d use probability H_0(c,d)/4;
* fill the diagonal.

Detailed balance is immediate. The total leaf probability out of c is
at most (H-1)/(4H)<1/4 and interior probability at most1/4, so P is lazy.
The leaf edges cost <=45 clocks by the original admission calculation;
their reservoir lifts satisfy the same all-representatives estimate.

For completeness, here is the entire variance estimate. For any function f,
write f_p(x)=f(p(x)) and let tilde_pi(c)=sum_{x in L_c} pi(x). Then

    Var_pi f <= 2 sum_x pi(x)(f(x)-f_p(x))^2
                + 2 Var_tilde_pi (f restricted to C).

The first term equals 8H times the leaf-edge energy. Since
`tilde_pi(c)<=H s pi_C(c)`, the second is at most
`2Hs Var_pi_C f <= (8H/gamma) E_interior,P(f,f)`.
As gamma<=1, this proves `gap(P)>=gamma/(8H)`, hence (4).
Conditional full-reservoir refresh again preserves that bound. Both the
stationary measure and every counted event are native; there is no
unpaid three-step jump through a representative.

## 5. Scope of the advance

This supplies an explicit positive projected-chain bound where the previous
packet left the projected bound open. It is uniform under compatible source
extensions with a fixed complete bank and controlled outside odds. The
constants are computable from the actual bank weights; none is called an
independently positive Collatz ROOT score. A different preassigned native
chain needs a lower conductance comparison with these displayed rates.
Arbitrary correlated source laws, extra shared pins, changing harmonic
banks, and strict temporal nonreturn require new proofs.

## 6. Match the reference activity without changing source weights

The ten auxiliary bits may instead have any common positive odds u. Replace
the denominator in (2) by `binom(10,r-5) u^(r-5)`. The gap uses the resulting
kappa(u); the boundary factor H and every native cost stay unchanged. In the
proposal formulas, replace w_i/(1+w_i) by w_i/(u+w_i), and1/(1+w_i) by
u/(u+w_i). Nothing about the actual prescribed law (1) changes. This is
choosing a better reference measure, not reweighting an actual source.

For the explicit heterogeneous control `a_i(U)=1, a_i(V)=1+(source_id mod3)`
on R, exact enumeration gives kappa(1)=27440, whereas kappa(2)=13048/115,
approximately113.46. These are test weights, not claimed ROOT weights. A
single activity choice improves the bound by a factor greater than240.
Global rescaling of all actual V odds by t, together with u->tu, multiplies
all F values by the same t^5; both kappa and H are invariant. The same change
works in the six-bit pinned models with exponent7. This supplies an explicit
normalization rule that retains relative activity instead of discarding it.
