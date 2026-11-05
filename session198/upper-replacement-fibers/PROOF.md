# Exact physical seed families for the upper-child replacement

Session198. **PROVED:** a complete physical description of the response-one
fiber behind the C3 upper-child construction. This supplies explicit source
recipes and exact parent readouts, not unrestricted integer-cone solutions.
General odd coverage remains **OPEN**.

Use chronological O/Z/J words and the actual source conventions of the frozen
positive source grammar. All sources in this packet are nontrivial, except the
explicitly named singleton bases used before O blocks. The response-one theorem
deliberately excludes the singleton, whose response C is also one. A source has raw response (D,E,U,d), with positive
coordinates (a,f,e,p)=(D-E-d,U-D+E-d,E,d-E). The nontrivial seed is
(1,0,2,1). Strong means E>0. Write b=(3,1,3,1), the rooted triangle.

## 1. Physical inequalities and balanced normalization

The three letters have the exact positive transfers

    O(a,f,e,p)=(a+e,a+2f+e+2p,0,e+p),
    Z(a,f,e,p)=(2a+f+e+p,f,e+p,0),
    J(a,f,e,p)=(f+p,0,a+e+p,e).

All coordinates are nonnegative; D and d are odd and d>=1. A strong source
ends in J or Z. These formulas therefore prove

    f<=a, p<=e;
    E<=d<=2E, D>=3d/2, D<=U<=2(D-E).             (1)

A balanced source has U=D, equivalently f=p=0. The exact profile identity

    (s . J) . Z = (s . O) . J                    (2)

permits its final Z-run to be moved before J. If the run reaches the seed t,
Z(t)=J(t). Hence every balanced profile has a physical recipe

    r . O^l J, l>=1,                            (3)

where r is strong or is the singleton (the O letters make the latter
nontrivial before J). If r=(D,E,U,d), (3) has profile

    (D+2^l U,D,D+2^l U,D).                       (4)

This is an existence and profile statement, not uniqueness of graph markings.

A nonbalanced strong profile with p=0 has f>0. Trace its final Z-run to the
last weak source. It cannot reach a J without encountering O, because J sets
f=0 and Z preserves f. Thus it has the form

    r . O^l Z^k, l,k>=1,                        (5)

with r strong or singleton. For strong r, the weak source's f is
2^l U-D-d >= D-d >= E >= d/2. Z preserves f and d, so

    f>=e/2 for every nonbalanced strong p=0 profile. (6)

For singleton r, f=2^l-2; the f>0 case has l>=2, and (6) holds as well.
For later use the profile in (5) is

    D'=D+(2^k-1)2^l U-(2^k-2)d,
    E'=d, U'=2^(l+k)U-(2^k-1)d, d'=d.           (7)

## 2. Complete upper-response-one fiber

The reviewed actual replacement identities give

    n=N(t,Oq)=3D_q+12U_q-2d_q,
    J0(t,Jq)=2n+C(q),
    C(q)=9D_q-3E_q-9U_q+4d_q=e_q-9f_q-5p_q.  (8)

**Theorem (source profiles).** The response profiles of actual nontrivial
sources q with C(q)=1 are exactly the union of:

1. q=(2^j+1,1,2^j+1,1), j>=1, with physical word
   t . O^(j-1) J. The parent is n=15*2^j+13.
2. q=s . O O J J, where s is an actual strong source satisfying
   4U_s-5D_s=1. Here

       q=(11D_s+2,5D_s+1,12D_s+2,6D_s+1),
       n=165D_s+28.                              (9)

3. q=s . O Z J, where s is an actual strong source satisfying
   D_s+2U_s-6d_s=1. Here

       q=(4U_s+4d_s+1,5d_s+1,4U_s+5d_s+1,6d_s+1),
       n=60(U_s+d_s)+13.                         (10)

Every displayed source gives an actual C3 graph with H=2n+1. The lists need
not be disjoint or unique as graph recipes.

**Proof of necessity.** C=1 excludes weak q. By (6), a nonbalanced Z profile
has C=e-9f<=-7e/2<0. If q is balanced, C=e=d, so d=1. Formula (4), and D_r=1
forcing singleton/transitive r, give precisely case1.

Otherwise q ends in J with a strong predecessor r, because p_q=E_r>0.
The equation C(Jr)=1 is exactly

    D_r-6E_r=1.                                 (11)

If r is balanced, use (3)-(4). Equation (11) becomes 2^l U_s-5D_s=1.
The strong-source bounds D_s<=U_s<2D_s force l=2; singleton s would require
2^l=6 and is impossible. This is case2.

If r is nonbalanced, it cannot be a J profile. A strong predecessor v would
make D_r-6E_r=U_v-5(D_v-E_v)<=-3(D_v-E_v)<0. A weak predecessor would make
r balanced. Thus r has form (5), with k,l>=1. For strong s, (7) and d<=2D/3
give

    D_r-6E_r
      =D+(2^k-1)2^l U-(2^k+4)d.

If l>=2 this is at least D[(10/3)2^k-17/3]>=D>=3.
If l=1 and k>=2 it is at least D[(4/3)2^k-11/3]>=5D/3>=5.
Only l=k=1 survives, giving D+2U-6d=1 and case3.
For singleton s, (7) gives D_r=2+(2^l-1)(2^k-1) and E_r=1.
Equation (11) would require the product of two numbers of the form 2^m-1
to equal5. One factor must be5, which has no such form. So there is no
singleton exception. This proves necessity. Conversely, direct substitution
in O/Z/J proves all profiles, readouts and C=1 identities in cases1-3.

## 3. An exact bridge from the lower C2 fiber to case3

Case3's seed condition has a further complete physical description:

    s strong, D_s+2U_s-6d_s=1
    iff s has the profile of y . J J,
        where y is actual balanced and D_y=4d_y-1. (12)

The converse constructs s physically. Necessity preserves a genuine source
recipe, not just an algebraic predecessor.

To prove it, a Z profile has d=E and D>=2d. Its expression D+2U-6d is at least
3(D-2d)>=3, since D and d are odd. Hence s must end in J. Write s=Jx.
Its equation is 3(U_x-D_x)-E_x=1, which rules out weak x. For a Z profile x,
(6) would give E_x=3f_x-1<=2f_x, hence f_x<=1. The f_x=0 case fails; f_x=1
would give E_x=d_x=2, contrary to odd deletion count. Thus x=Jy, with y strong,
and the equation becomes

    D_y=4E_y-1.                                  (13)

A nonbalanced J profile y is impossible: a strong predecessor z would require
U_z-3(D_z-E_z)=-1, contrary to U_z<=2(D_z-E_z) and D_z-E_z>=2; a weak
predecessor makes y balanced. A nonbalanced Z profile has form (5). For strong
base z, (7) gives

    D_y-4E_y
      >=D[(4/3)2^k-7/3]>=D/3>0,

using l>=1, whereas (13) requires -1. For singleton base, the same expression
is (2^l-1)(2^k-1)-2. Equality -1 forces l=k=1, the balanced triangle.
Therefore y is balanced, proving (12).

Write y=(4h-1,h,4h-1,h). Then the case3 source is concretely

    q=y . J J O Z J,
    n=1320h-347,      J0(t,Jq)=2640h-693.         (14)

The seed y lies on the negative C2 replacement fiber:

    3D_y-4E_y-2U_y=-1.

Thus (14) physically transports a known lower-response-one witness to a
specified upper-response-one family. It does not claim that its new parent
readout equals the old parent, or that all integers h are allowed.

## 4. Exact arithmetic correction, rather than residue coverage

Both seed equations in the preceding description have uniform actual solutions.
For any source with deletion sidecar d, direct substitution gives

    (2U-3D)(s . O Z) = (2U-3D)(s)-2d,
    (4U-5D)(s . O O Z Z) = (4U-5D)(s)-2d.        (15)

Both words preserve d and finish at a strong source. Starting at t and using
only O and Z preserves the genuine unit deletion sidecar d=1: deleting the
root leaves a transitive tournament. Consequently:

**Exact correction theorem.** For every odd integer r, positive or negative,
there is an actual strong source with d=E=1 and 2U-3D=r. There is also such a
source with 4U-5D=r. These are exact equalities over the integers, not modular
or relaxed-profile assertions.

For the first construction, choose a power u=2^(l+1), l>=0, with 2u-3>=r+2.
The actual prefix t . O^l has D=d=1 and U=u. Set k=(2u-3-r)/2>=1 and append
(OZ)^k. Formula (15) proves the claim. Its exact final cover count is

    U=[(3u-1)4^k+1]/3,       D=(2U-r)/3.         (16)

For the second construction choose such a power u with 4u-5>=r+2, set
k=(4u-5-r)/2>=1, and use t . O^l (OOZZ)^k. Then

    U=[(5u-1)16^k+1]/5,      D=(4U-r)/5.         (17)

The displayed quotients are integers because they came from the physical
recurrences U'=4U-1 and U'=16U-3. Choosing the smallest allowed power gives
O(|r|+1) vertices; the positive Hamilton and cover counts can be exponential
in |r|. No bound comparable to |r| for those counts is claimed.

Taking r=1 in (17) supplies case2 of the main theorem uniformly. Taking r=-1
in (16), then applying OJ, supplies a balanced y with D_y=4d_y-1; this gives
case3 through (12). Thus both nontrivial branches have explicit infinite
primitive physical families, without relying on a finite source search.
The minimal-power choices in these two pump constructions return the seeds
(23,1,29,1) for case2 and (19,5,19,5) for the balanced negative seed y in
case3. Their upper-child parents are respectively3823 and6253. These are not
the smallest physical seeds in the complete branches: the literal checks
also include parents2503 and973.

This mechanism controls a signed response, not an independently specified
positive count n. In particular the exponent k fixes the large U in (16) or
(17); it is not a free second coordinate after r has been prescribed.

## 5. Reproduction and scope

Run python3 verify.py and python3 -O verify.py. The verifier checks the cones,
the balanced and residual formulas, all actual source profiles through
U<=5000, every upper-response-one profile under that bound, both exact seed
conditions, the necessity/sufficiency of (12) under the same complete bound,
the two arithmetic-correction laws on a basis and exact target tests for all
odd r from -101 through101, and actual parent/child graphs for three seed types using subset
Hamilton-path and cycle counts. It also checks all symbolic profile/readout
identities on a coordinate basis and sample parameter values. These finite
controls supplement the all-order proofs above.

The positive source model and Python/proof dependencies are unchanged copies
of the reviewed target43 replacement packet. No numerical saturation claim is
made. The remaining arithmetic task is still physical: determine which D occur
in case2, and which h occur in balanced sources (4h-1,h,4h-1,h). The theorem
provides exact construction branches and a complete description of this one
response fiber; it does not establish cofinite C1/C2/C3 coverage.
