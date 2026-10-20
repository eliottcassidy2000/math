# Smooth common-root noncancellation for outward edits across a one-step gap

20 October 2026. Status: **PROVED**, with a cited Minkowski input and independent mathematical and exact-computational audit. See [the independent audit](audit/OUTWARD_GAP_ONE_AUDIT.md). This packet treats a specified mixed-direction disjoint family. The replacement for the external-prime argument is a valuation-and-norm contradiction using exact legality.

## 1. The actual family and proposed content theorem

Let p>=3, 1<=r<p, gcd(p,r)=1, q=p+r, and start from the nonconstant lower zero-intercept mechanical word. Source ones have positions s_j=j+ceil(jr/p)-1 modulo q and residue labels a_j=-jr modulo p.

Choose consecutive source ones separated by gap one, with residues u+r followed by u. Move the first one left and the second one right. The actual local change is

    0110 -> 1010 -> 1001,

or the commuting order 0110 -> 0101 -> 1001. The complete parameter range is

    max(0,p-2r)<=u<min(r,p-r).                       (1.1)

Both affected edges are disjoint and legal in the source. All nonempty cases here have q>=5. The residue inequalities alone would also admit p=2,r=1,u=0, but at q=3 both edges share the same target zero; those two swaps cannot both be executed. That boundary is excluded explicitly, not counted as a disjoint path.

For the resulting word w put D=2^q-3^p, B=B(w), and g=gcd(|D|,B). The proposed theorem is

    g<=22^(4 p^(2/3)).                              (1.2)

This is a full integer-content bound including prime powers and both clock signs. It is not a norm bound or an all-length integer exclusion by itself.

## 2. Legality and a smooth root gate

A gap-one consecutive pair has residues u+r,u with 0<=u<p-r. A source left move is legal exactly at residue at least p-r; a source right move is legal exactly below r. Applying these criteria gives u+r>=p-r and u<r, proving (1.1). Conversely those inequalities make the two outer positions zeros. At q>=5 the two affected pairs are disjoint, so the swaps commute and both orders remain actual.

Every legal parameter satisfies the three strict inequalities

    u<p-r,  r>p/3,  r+u<=p-1.                       (2.1)

The second follows because if r<=p/2, nonemptiness requires p-2r<r; if r>p/2 it is immediate.

For any h>1 dividing |D|, the canonical unit root satisfies

    z^p=2, z^r=3/2,

and z-1 is also a unit. The [moving-anchor carry identity](../collatz_two_swap_norm_barrier_20261020/symbolic/TWO_SWAP_SPARSE_FORMS_AND_BARRIER.md) gives carry vanishing iff

    1+(z-1)(z^u-z^(u+r)/2)=0.

Using the clock reduces this exactly to

    z^u(z-1)=-4.                                   (2.2)

Thus if z is reconstructed modulo the full |D|, then

    g=gcd(|D|,4+z^u(z-1)).                          (2.3)

No prime outside 2 and 3 occurs on the right of (2.2). The previously used single external-prime valuation witness therefore does not apply. The following argument retains the extra arithmetic and geometric information instead.

## 3. The exponent lattice and cleared polynomials

Use the determinant-p lattice

    L={(b,c,m) in Z^3: m-br-cu is divisible by p}.

For a vector set a=(m-br-cu)/p and H=max(|b|,|c|,|m|). Exact legality gives |a|<=H because 1+r+u<=p. Negate the vector if needed so c>=0. The modular equations imply

    z^m(z-1)^c=(-1)^c 2^(a-b+2c)3^b.                (3.1)

Writing e=a-b+2c, clear only powers of 2 and 3:

    F_v(X)=U X^max(m,0)(X-1)^c - V X^max(-m,0),
    U=2^max(-e,0)3^max(-b,0),
    V=(-1)^c2^max(e,0)3^max(b,0).                   (3.2)

The whole actual g divides F_v(z).

### Zero, one, constants, and degree

No nonzero legal relation produces the zero polynomial. Neither 0 nor 1 is a root. At 0 a possible cancellation requires m=0 and equality of the rational side of (3.1) with (-1)^c. Unique factorization forces b=0,a=-2c, whence 0=m=c(u-2p). If c>0 this contradicts u<p-r; if c=0 it is the zero vector. Other m cases have a single nonzero constant term. At 1, c>0 leaves the second term, while c=0 would force a=b=m=0.

A possible leading cancellation when m<0 and -m=c requires U=V, hence b=0,a=-2c and u=2p-1 (as well as an appropriate even sign). Again this is outside the legal range. Therefore F_v is constant exactly when m=c=0, and has degree at most |m|+c<=2H. Two independent relation vectors cannot both produce constants.

## 4. Coprimality without an external prime

**Lemma.** Polynomials (3.2) from any two independent vectors of L are coprime over Q.

A nonzero constant is harmless, so suppose a common algebraic root x exists. Section 3 shows x!=0,1. Both rational identities (3.1) hold at x.

If the projected pairs (m_i,c_i) are dependent, a nonzero integer combination kills both m and c. Independence of the full rows leaves nonzero (a,b) with ap+br=0. Thus (a,b)=(rt,-pt), t!=0, and the same multiplicative combination would give 1=2^(qt)3^(-pt), impossible. Hence only independent projections require a new argument.

Let K=Q(x), of degree d. At any finite place above a rational prime ell, normalize the valuation by nu(ell)=1. Independence of the two (m,c) rows implies that nu(x) and nu(x-1) are uniquely determined by the two equations. In particular they have the same values at every place over a given ell.

For ell outside {2,3}, both rational right sides have valuation zero, so

    nu(x)=nu(x-1)=0.

At 2 write these common values as alpha_2,beta_2; at 3 as alpha_3,beta_3. The row identities give

    m_i alpha_2+c_i beta_2=a_i-b_i+2c_i,
    m_i alpha_3+c_i beta_3=b_i.

Multiply the first equation by p and the second by q=p+r, then use m_i=a_i p+b_i r+c_i u. Since the two projected rows are independent,

    p alpha_2+q alpha_3=1,
    p beta_2 +q beta_3 =2p-u.                       (4.1)

### Only two valuation patterns are possible

If alpha_2<0, then beta_2=alpha_2 by the valuation triangle law. The first equation in (4.1) forces alpha_3>0, hence beta_3=0. The second equation would have a negative left side and positive right side, a contradiction. The case alpha_3<0 is symmetric.

If both alpha values are positive, both beta values are zero, again impossible. Both alpha values zero contradict the first equation. Therefore exactly one alpha is positive and the other is zero. The beta at the positive alpha is zero, and the remaining beta is nonnegative. All finite valuations of x and x-1 are consequently nonnegative, so both are algebraic integers.

Since every place over ell has the same normalized beta_ell, the norm satisfies

    |Norm_K/Q(x-1)|=2^(d beta_2)3^(d beta_3).         (4.2)

Indeed the exponent at ell is sum_{P|ell} f_P e_P beta_ell=d beta_ell. No other prime occurs.

### Eliminating the edit relation

Define

    M_1=2/x^p, M_2=3/(2x^r), M_3=-4/[x^u(x-1)].

The two row relations are M_1^a_i M_2^b_i M_3^c_i=1. Eliminating c produces

    M_1^A M_2^B=1, (A,B)!=(0,0).                   (4.3)

The combination is nonzero because the original full rows are independent; alternatively its m-coordinate is the nonzero determinant of the projected rows.

### Case I: the positive x-valuation is above 2

Then

    alpha_2=1/p, alpha_3=0,
    beta_2=0, beta_3=(2p-u)/q>1,

where the strict inequality uses u<p-r. At a place over 3, M_1 has valuation zero and M_2 has valuation one. Equation (4.3) forces B=0 and A!=0. Thus M_1 is a root of unity. Every conjugate of x therefore has absolute value 2^(1/p)<=sqrt(2). Hence

    |Norm(x-1)| <= (1+sqrt(2))^d < 3^d.

But (4.2) says |Norm(x-1)|=3^(d beta_3)>3^d, a contradiction.

### Case II: the positive x-valuation is above 3

Then

    alpha_2=0, alpha_3=1/q,
    beta_3=0, beta_2=2-u/p>q/p>4/3.

The two strict inequalities use u<p-r and r>p/3. At a place over 2, M_1 has valuation one and M_2 has valuation minus one. Equation (4.3) forces A=B!=0, so M_1 M_2=3/x^q is a root of unity. Every conjugate of x has absolute value 3^(1/q). Since q>=3,

    |Norm(x-1)| <= (1+3^(1/q))^d
                 <=(1+cuberoot(3))^d
                 <(5/2)^d<2^(4d/3).

The last strict comparison is exact since (5/2)^3=125/8<16. Equation (4.2), however, gives |Norm(x-1)|=2^(d beta_2)>2^(4d/3). This contradiction finishes the coprimality proof.

The proof uses legality essentially; it is not a general smooth-S-unit noncancellation claim.

## 5. Height and full-content resultant bound

From (3.2),

    ||F_v||_1 <= 2^(max(-e,0)+c)3^max(-b,0)
                 +2^max(e,0)3^max(b,0).

The right side is convex on |a|,|b|<=H, 0<=c<=H. At the eight normalized vertices its H=1 values in lexicographic order are

    4,10,7,5,7,22,4,14.

At height H>=1, each vertex is A^H+B^H<=(A+B)^H, so ||F_v||_1<=22^H. Degree is at most 2H.

The [audited Minkowski cube argument](../collatz_two_swap_overlap_lattice_content_20261020/OVERLAPPING_TWO_ONE_CONTENT.md), Section 5, gives two independent vectors with H_1 H_2<=p^(2/3). The coprimality lemma makes their integer resultant nonzero. Its integral Bezout identity, including the one-constant case, shows that full g divides it. Sylvester-Hadamard yields

    g<=|Res(F_1,F_2)|<=22^(4H_1H_2)<=22^(4p^(2/3)),

proving (1.2). This final resultant still need not recover actual content exactly.

## 6. Exact controls and boundaries

`verify_outward_gap_one.py` checks all 1,101 coprime shapes 2<=p<=60, 1<=r<p. The p=2,q=3 false-disjoint boundary contributes no path. It compares the residue universe with direct physical source edges, executes both commuting orders, verifies two carry algorithms and the whole-gcd iff, and constructs short relation polynomials and exact Sylvester resultants.

The universe contains 7,367 unordered edit pairs and 14,734 ordered actual paths, including 1,100 original-cut crossings. The maximum observed actual content is 2785 at p=59,r=30,u=17. No zero resultant or integral output was observed. This is a finite control, not an all-length exclusion.

An explicit illegal-parameter hostile shows why no unrestricted smooth-root lemma may be substituted. Take p=3,r=1,u=9, outside (1.1), and relations (b,c,m)=(1,0,1) and (0,2,0). They give

    F_1=2X-3,
    F_2=4(X-1)^2-1=(2X-1)(2X-3),

so the resultant is zero. The common root x=3/2 has a negative 2-adic valuation. In (4.1), 2p-u is now negative, so the contradiction eliminating that valuation is unavailable. The legal range is a mathematical hypothesis, not an optimization or cosmetic filter.

Reproduce the finite controls from this directory:

    python3 verify_outward_gap_one.py > receipt.replay.json
    python3 -O verify_outward_gap_one.py > receipt.optimized.replay.json
    cmp receipt.json receipt.replay.json
    cmp receipt.json receipt.optimized.replay.json

A separate closure could combine this content theorem with a signed finite-minimum argument and the already audited absolute-gap certificate. That closure is not claimed here. Neither this family nor its proof covers arbitrary disjoint pairs.
