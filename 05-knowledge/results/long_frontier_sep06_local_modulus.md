# Sharp local distance and moment moduli at the three-atom limit

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** This connects
[THM-4456's local expansion](../../01-canon/theorems/THM-4456-sharp-finite-length-signed-root-stability-asymptotics.md)
with the incoming [moment/product characterization](continuing3_20260906_stability_near_minimizers.md).
The incoming core rigidity reproof is credited; the claims below concern
an explicit moment-distance remainder and sharp local quotient constants.
No globally sharp stability modulus is asserted.

## 1. Statement and incoming connection

Keep the actual finite-list domain p1=p2=1,E>0 and the quotient R,K3
from THM-4454 and THM-4455. Put z=1/sqrt(3), let a>=b>=c>=0 be the
three largest positive coordinates with zero padding, and define

    Delta=2-2z(a+b+c),
    M=sum_i r_i^2(r_i-z)^2=p4-2z*p3+1/3,
    A=-2-sqrt(6)/3+2sqrt(2)+(7/3)sqrt(3)>0.

For every such list, both Delta and M are positive and

    |M-Delta/3| <= (2/sqrt(3))*Delta^(3/2)+Delta^2.       (1)

In particular Delta<=1/256 implies

    Delta/4 <= M <= (5/12)Delta.                       (2)

Along every sequence with R tending to K3,

    M/Delta -> 1/3,
    liminf (R-K3)/Delta >= A,
    liminf (R-K3)/M >= 3A.                            (3)

Both constants are sharp, attained as limits by the actual three-equal-
positive/equal-negative family of THM-4456. They are approximately
4.053382428 and12.160147284, respectively.

More precisely, let ell=(a+b+c)/3, m be the remaining square mass,
and h=ell-c. Along a minimizing sequence, convergence of either ratio
in (3) to its displayed sharp constant is equivalent to h=o(m).
Equality of only a liminf does not imply this whole-sequence condition.

The source-to-target map retains the ordered three-atom matching and
its signed coordinate errors. The incoming M invariant suppresses
their signs but retains square-weighted proximity to0 and z. The exact
third- and fourth-order remainder below records the lost data. A model
with alternating equal and split upper families is the cheap hostile to
confusing a subsequential optimum with convergence of the full sequence.

## 2. Exact moment-distance remainder

Match the three target entries z to a,b,c, and zeros to every other
coordinate. Write their matched errors as e_i=r_i-z on the first three
slots and e_i=r_i on the tail. Then sum e_i^2=Delta. Direct expansion
gives the exact identity

    M-Delta/3
      =2z(sum_leading e_i^3-sum_tail e_i^3)+sum_all e_i^4. (4)

Consequently its absolute value is at most
2z sum|e_i|^3+sum e_i^4, bounded by2z Delta^(3/2)+Delta^2.
This proves (1) globally, irrespective of length and tail signs.
For Delta<=1/256 the relative error is at most1/(8sqrt3)+1/256<1/12;
the last comparison follows by squaring96/sqrt3<61. This proves (2).

M=0 would force every root to0 or z, with exactly three z entries
by p2=1, contradicting p1=1. Delta=0 has the same contradiction.
Thus every quotient used here is actually defined.

## 3. The sharp local response

THM-4455 gives Delta->0 along every R->K3 sequence. Equation (1)
therefore yields M/Delta->1/3. For the same actual leading coordinates,
THM-4456 proves

    R-K3=A*m+B*h+O(m^(3/2)+h^2),
    B=8-5sqrt2-4sqrt3+3sqrt6>0,
    v=sum_leading(r_i-ell)^2<=6h^2,
    3ell^2=1-m-v.

The exact distance identity is

    Delta=2-2sqrt(1-m-v)=m+O(m^2+h^2).               (5)

Fix any epsilon>0. Subtracting (A-epsilon)Delta from the response
expansion gives epsilon*m+B*h+O(m^(3/2)+h^2), which is nonnegative
in a sufficiently small fixed neighborhood. The constants are independent
of length. This proves the first liminf in (3); the second follows
from M/Delta->1/3. Equivalently, for each epsilon>0 there is a positive
eta such that 0<R-K3<eta implies both local inequalities with constants
A-epsilon and3A-epsilon. This is a sharp local asymptotic statement,
not an explicit optimal global eta or an all-list sharp modulus.

On the equal-three/equal-negative family h=0,m->0. The expansion and
(5) give (R-K3)/Delta->A and (R-K3)/M->3A, proving sharpness.

## 4. Convergence of the sharp ratio and the subsequence boundary

Eventually m>0 on a minimizing sequence: its leading sum tends to
sqrt3>1, so the actual p1 normalization forces nonzero negative tail.
If either ratio in (3) stays bounded, local coercivity and (5) give

    c0(m+h) <= R-K3 <= c1(m+h^2)

for positive fixed constants. Absorb h^2 for small h to obtain h=O(m).
After dividing the expansions by m, uniformly on such a subsequence,

    (R-K3)/Delta=A+B*h/m+o(1),
    (R-K3)/M=3A+3B*h/m+o(1).

Therefore convergence to the respective sharp constant is equivalent
to h/m->0. Conversely h=o(m) directly gives these limits.

If equal-three lists h=0 alternate with the exact normalized split
lists h=1/n from THM-4456, both subsequences still have R->K3. The
liminf ratios are A and3A, but h/m has a positive limit on the split
subsequence. Thus liminf equality alone does not imply h=o(m).
This quantifier distinction was required by the independent initial
analytic review and is retained explicitly.

## 5. Verification and remaining problem

The universal proof consists of the exact identity (4), norm inequalities,
and the audited uniform local expansion. The exact checker verifies the
coordinate remainder, constant identities, distance derivatives and actual
upper-family limit. Its finite controls do not supply a universal quantifier.
The [independent full proof/source audit](long_frontier_sep06_local_modulus_audit.md)
passes without correction. All15 always-active exact gates pass normally
and under optimization, with byte-identical [frozen output](long_frontier_sep06_local_modulus.out).

```bash
python3 -B 04-computation/long_frontier_sep06_local_modulus.py
python3 -B -O 04-computation/long_frontier_sep06_local_modulus.py
```

Raw LF source SHA256:
4757c9d537a3d73ad3093fd94a695fc32eb915264e2b93dcb3c5a3e1e5033366.
Raw LF output SHA256:
4ada2b28cb6d07de820419b5f65c5b2bc6da863eb1a4c5086345713cfd793206.

The incoming entire-product equivalence remains an independent consequence
of signed first-moment normalization. This note improves neither its global
convergence rate nor the global optimal coefficient relating R-K3 to M or
Delta. Those global questions must also account for the one- and two-atom
denominator boundaries.
