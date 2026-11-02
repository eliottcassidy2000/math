# Independent audit: exact weighted stability for every n>=1500

**Verdict: ACCEPT, including the complete positive-variance equality locus.**
This is a separate third-checkpoint continuation. It uses already accepted
all-order inequalities and no census at order1500 or above.

For R=floor(n/2), the result is

    K_(2R)=(4/R)(1-(3/4)^R-(1/4)^R),
    K_(2R+1)=(4/R)(1-(3/4)^R),             n>=1500.

For V>0 equality consists exactly of ordinal fair pairs, with one singleton
at odd order. Orders9..1499 remain outside this theorem; the threshold1500
is not asserted optimal. The proof does not depend on the earlier10^7
candidate's conclusion: its inputs are the accepted asymptotic baseline,
matching-support theorem, triangle-variance bounds and H1 classification.

## Reviewed dependencies and domains

B=E Delta(T) is retained throughout. For n>=4, B>=tau>=W, the nonlinear bound
is `Delta>=B+(1-exp(-4V))/4`, and the matching-support inequality holds for
arbitrary binary cross orientations. Its coefficient a_m is used only on
nonempty nonperfect retained matchings. The greedy matching selects only
positive-q edges. Its assigned groups have disjoint edge sets, so the
unordered wedge pairs counted inside them are distinct. Both group bounds
`W_group>=wS` and `W_group>=S^2/4` are valid; heavy/light uses do not falsely
add two lower bounds for the same group.

The historical10^7 proof and its exact scalar controls were independently
reviewed. That route is superseded quantitatively, is not a dependency, and
is deliberately excluded from this curated publication packet.

## 1. The below-candidate reduction

It suffices to consider Delta<=kappa V with V>0. The exact fair-pair formula
gives `kappa R<=4` and `kappa<=8/(n-1)<1/187` for n>=1500. The partial
Taylor sum `sum_(j=0)^7 4^j/j!=16319/315>50` gives c0>49/200, so the
nonlinear gap bound excludes V<=1. It then yields

    kappa V>=c0>49/200,
    V>49(n-1)/1600.

The inequality `W>=2V^2/n-V/4`, together with W<=B<=Delta, gives
`V<=n/8+4n/(n-1)<=n/7`. The last inequality is exactly n-1>=224.
Retaining B in the greedy-matching bound gives B<=kappa D; combining with
D<=sqrt(2nW)<=sqrt(2nB) gives B<=2n kappa^2 and D<=16n/(n-1).
The B0 case has W=D0 directly, so no division by zero is hidden.

## 2. Heavy matching and its exponentially small tail

Let eta=1/44 and suppose at least one selected edge is light (weight<=eta).
Delete all L light selected edges. The retained matching has

    r_h<=R-L,
    M_h>49(n-1)/1600-16n/(n-1)-n/88>=n/120.

Subtracting n/120 gives
`f(n)=(577/52800)n-49/1600-16n/(n-1)`. Its derivative is strictly positive,
and its exact value at1500 is9253539/26382400>0. Thus r_h>=4M_h>n/30>0.
Since L>=1, the retained matching is not perfect; the correct coefficient
is a_(r_h), never the full-even b coefficient.

Set epsilon=(3/4)^(n/30). The function n(3/4)^(n/30) decreases for n>=1500:
its logarithmic derivative is less than1/n-1/120<0, since log(4/3)>1/4.
At1500 the required bound is the exact integer comparison
`1500000*3^50<4^50`. Therefore epsilon<=1/(1000n) at every larger real n.
This is an analytic uniform extension, not extrapolation from sample orders.

## 3. The parity-free light-count charge

Using R=floor(n/2), not n/2, is essential to the improved constants. Multiply
the retained-matching inequality by R-L>=r_h and use kappa R<=4. This gives

    L[kappa V-4eta(1-epsilon)]<=4epsilon V+4D.

The term4epsilon V is at most1/1750. Because L>=1 it can be absorbed into
L/1750. The remaining coefficient is strictly greater than

    49/200-1/11-1/1750
      =3/20+271/77000 >3/20.

Consequently `L<(80/3)D`, so D>0. This last strict inference does not hold
by fiat: it follows from the displayed positive coefficient and L>=1.

## 4. The contraction contradiction

On heavy groups, D_h<=W/eta<=44B. On light groups, Cauchy gives
D_l<=2sqrt(LW)<=2sqrt(LB). Hence

    D<=44B+2sqrt(LB).

If B0 this already contradicts D>0. Otherwise substitute the strict light
bound and B<=kappa D, then divide by D>0:

    1<44kappa+2sqrt((80/3)kappa)
      <=4/17+sqrt(320/561)<1.

The final comparison is exact: its remaining positive term is13/17, and
`(13/17)^2-320/561=137/9537>0`. Thus L>=1 is impossible. There is no
uncontrolled square root approximation and no order-dependent limiting step.

## 5. Closure and equality

Every selected edge now has weight w>eta>kappa. If D>0, at least one group
has S>0, so `B>=W>=sum wS>kappa D`, contradicting B<=kappa D. Therefore
D=B0. The uncertainty graph is a matching and the matching theorem gives
Delta>=kappa V, proving the universal inequality. Fair ordinal blocks
supply the reverse upper bound on K_n.

For equality with V>0, B0 means every supported binary completion has zero
strong gap. At this order the only such binary tournaments are transitive;
the order-three exception is irrelevant. Hence H=1, and the already-audited
scalar minimization forces a maximal collection of fair pairs with the
stated singleton parity. The converse is checked by the explicit formula.

## Independent arithmetic controls

`check.py` imports no producer code or data. It separately reconstructs the
Taylor bound, rational kappa/variance endpoints, retained-mass function and
its positive derivative, exact exponential endpoint and negative logarithmic
derivative, absorbed light coefficient, and squared contraction margin.
All12 certificates use Fraction or integer arithmetic and active exceptions.
Normal and Python -O outputs agree byte-for-byte; both processes exited0.

Semantic SHA256:
`06497c79c3f5b0905171f6484197662b809f75932653d415babddd98af2d262a`.
The proof is the preceding inequalities and monotonicity arguments; these
small checks are arithmetic execution controls, not a finite matrix proof.
The reviewed promoted manuscript and dependency snapshots are pinned in
`REVIEWED_SOURCES.json`. No Git or external-state mutation was performed.

## Firewall

None of these arguments uses a termwise lower bound on a Pfaffian containing
one fractional edge. Such a term can cancel (for example3x-2). The accepted
independent-set bounds retain entirely binary principal subsets, and the
robustness lemma uses the complete Gram/Cauchy sum with its actual coefficient
energy. The failed termwise shortcut has no dependency path into this proof.
