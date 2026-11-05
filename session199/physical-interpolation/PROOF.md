# A sharp local interpolation barrier and exact one-edge rebasing

Session199. **PROVED, pending independent audit:** local unit-step interpolation
inside a sufficiently long retained-root block is impossible: its unit-step
graph is a matching. In contrast every even Hamilton-path difference is
realized by a single edge flip between two strong C1 tournaments, and an
explicit digit family gives exact target-preserving rebases. **OPEN:** a
construction at an independently prescribed parent count, and universal odd
coverage. The positive families below do not close those quantifiers.

## 1. Inherited physical setting

Use `session199/source-cut-rearrangement/PROOF.md`, published at
`7812db73cdac9526ba5e9aa7849a4d2d5a78cb5d`, including its independent audit.
A nontrivial actual source has (D,E,U,d), with d>=1 and U>=2d, and

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).

These are actual one-vertex moves. The source seed is tau=(1,0,2,1).
Write F=D-E+U for its C1 closure and G=D+2U for its balanced C1 closure.
The inherited paper proves that every fixed-partner classified C1/C2/C3
readout, including after any fixed O/Z/J suffix, has integer coefficients

    q=alpha D+beta E+gamma U+eta d,
    alpha>=1, gamma>=1, tau_q=2alpha+beta+gamma>0.

Moreover tau_q is even: it is exactly the difference of q on tau.OZ and
tau.ZO, two actual tournaments with odd Hamilton-path counts. Hence tau_q>=2.
The letter tau for the edge and the scalar tau_q are distinct.

## 2. No unit-step paths inside a long block

Fix an actual prefix source, an integer k>=5, and one such final readout q.
The vertices of the *word cube* are the 2^k O/Z blocks of length k. Join two
words if they differ at one letter and their final path counts differ by2.

**Theorem 1.** This graph is a matching: every word has degree at most one.
The conclusion also permits any fixed later O/Z/J suffix, because it can be
transported into q by the inherited strict-cone theorem.

**Proof.** Let d be the deletion count of the input; it stays constant in
this block. The digit formula makes the scalar a quadratic bit function with
mixed coefficient `-alpha d w_ij`, where
`w_ij=2^(j-i-1)+[j=i+1]`. Define m_i to be the scalar change from O to Z at
position i, with all other bits fixed. This is independent of the bit at i.
For adjacent positions i,i+1, write Delta for the positive difference
`q(...OZ...)-q(...ZO...)` with their other bits held fixed. Their mixed
coefficient is -2alpha d, so at any word

    m_(i+1)-m_i = Delta+2alpha d(epsilon_(i+1)-epsilon_i)
                 >= Delta-2alpha d.                    (1)

We show that the last quantity is at least6d. Put m=k-i-2, the number of
letters after the adjacent pair. For i>=1 the actual state immediately
before that pair satisfies

    U_i+E_i >= 2^i(U_0-d)+2d >= (2^i+2)d.                (2)

Indeed U_i+E_i=2U_(i-1) regardless of its last O/Z letter, and iterating the
smallest possible cover update gives
`U_(i-1)>=2^(i-1)(U_0-d)+d`.

If i,m>=1, the exact exchange formula from the inherited paper gives

    Delta-2alpha d >= [alpha2^i+gamma2^m]d >=6d,

because i+m=k-2>=3. If i=0, m>=3, it gives
`Delta-2alpha d>=gamma2^m d>=8d`. If m=0, i>=3, it gives

    Delta-2alpha d >= [alpha(2^i-2)+tau_q]d >=8d.

Thus the consecutive bit marginals at any word are strictly increasing with
gap at least6. At most one can belong to {-2,+2}. Such marginals are exactly
the incident unit-step edges of that word, proving degree at most one. QED.

The threshold5 is sharp. For tau, the readout F, and k=4, the words
`ZOOO`, `ZZOO`, `ZZZO` have counts27,25,27 and form a three-vertex path.
The first and third differ at two letters. Hence the graph at k=4 need not
be a matching.

This is a limitation of **local one-letter interpolation in a fixed block**.
It does not bound the size of an interval in the image, does not forbid using
several blocks, and does not imply a gap in the full tournament spectrum.
Distinct components can have adjacent scalar values.

## 3. Sharp adjacent-exchange scales

For tau and the F readout, among all adjacent OZ-to-ZO exchanges inside
k-letter blocks, the minimum positive decrement is

    2 at k=2, 4 at k=3,
    2+2^floor((k-2)/2)+2^ceil((k-2)/2) at k>=4.          (3)

For i,m>=1 the exact lower bound is `2^i+2+2^m`, attained with prefix Z^i
and suffix O^m. The two end cases give `2+2^(k-2)`. Balancing i,m yields (3),
with the end cases tying at k=4. Thus even the smallest ordered-exchange
step grows on the square-root scale in block cardinality. A greatest-common-
divisor argument about the steps cannot turn this into interpolation.

## 4. Every even difference is a literal single-edge difference

Let g be the reviewed one-apex function. Every U>=2 gives an actual source

    w(U)=(g(U), U mod2, U,1).

It is obtained from tau using O at n->2n and Z at n->2n-1. The inverse
reduces n to n/2 if even and (n+1)/2 if odd, terminating at2. If U is odd
and at least3, its final Z makes it strong with E=d=1.

**Theorem 2.** Every even integer delta is the difference H(T_Z)-H(T_O)
of two strong C1 tournaments on the same vertex set differing at exactly one
edge. Infinitely many such pairs exist for each delta.

Choose any B=2^m, m>=1, so large that U=B+1+delta>=3. This U is odd. Start
with the same actual source w(U), append respectively O or Z, then append
m further O letters and take the C1 closure J. These two graphs differ only
in the edge between the first newborn and the old marked root. Every later
edge and vertex is identical. Both have exactly one Hamilton cycle. Their
counts are exactly

    H(T_O)=g(U)+2BU,
    H(T_Z)=g(U)+2BU+(U-1-B)=H(T_O)+delta.                (4)

This proves the assertion. B is arbitrarily large, and the count tends to
infinity, so the pairs are unbounded. Their order is
`m+ceil(log2 U)+3`. Nothing prescribes the first count independently of delta.
For example delta=2 gives pairs25/27 at B=2,61/63 at B=4,185/187 at B=8.

There is no conflict with Theorem1: these are isolated unit-step edges when
the full block is long, not a connected path through a long odd interval.

## 5. A parameterized exact scalar rebase

For a nonempty O/Z suffix W, let its digit data be A,p,B=2^m as in the
inherited paper. A single O-to-Z toggle immediately before W gives the full
profile difference

    (s.Z.W)-(s.O.W)
       =(U-E-(A+p)d, 0, -Bd, 0).                       (5)

This follows by applying the exact linear block transfer to the initial
difference `(U-E,d,-d,0)`. In particular its G difference is

    U-E-(A+p+2B)d.                                     (6)

Choose the actual unit-deletion strong source w(U) with

    U=A+p+2B+1.

A+p is even, so U is an allowed odd integer. Equations(5)-(6) give two
physical sources with exactly equal G and full profile difference

    (2B,0,-B,0).                                      (7)

Their balanced OJ closures therefore have the same Hamilton-path count and
different endpoint counts, separated by2B. This is an explicit scalar-
preserving rebase, not a claim that every prescribed scalar has such a rebase.

There is also a universal six-letter identity. For every actual source s,

    G(s.OJJOOO)=G(s.OJJOZO)=33D(s)+68U(s).              (8)

Indeed after OJ put h=D+2U and e=D; the profile is(h,e,h,e). After JO it is
(2h-e,0,4h,h), satisfying U-E=4d. The common final suffix is O, so(6) vanishes.
The two final profiles differ by(4h,0,-2h,0), and the common G is34h-e.
No full-profile equality, arbitrary-suffix equality, or graph isomorphism is
asserted in(8).

## 6. Exact two-block probe and remaining scope

The accompanying finite probe fixes the word shape `B_0 J B_1` from tau,
with block lengths(k,l), and uses G. All 2^(k+l) bit assignments are tested,
with an edge exactly for a one-bit flip of scalar difference2. For k=l=5,
a component gives the three consecutive values6679,6681,6683. Equal lengths
6 through9, and the listed unequal lengths, have only isolated unit edges.
These are complete finite word-cube statements, not all-source censuses or
an all-length bound on the two-block case.

The single-block matching theorem explains why naive swap paths did not
fill the gaps. The one-edge construction and(7)-(8) supply genuine physical
survivors, but their parent count is not freely selectable. A positive
coverage proof still needs target-dependent global rebasing, midpoint
transport, or another exact-level argument. No all-odd coverage theorem is
claimed.

## 7. Why complementing digits cannot convert the signed readout

The auxiliary signed response `2U-D` has positive mixed digit coefficients
`+d w_ij`, whereas every fixed classified physical Hamilton readout has
negative mixed coefficients `-alpha d w_ij`, alpha>0. This remains true
for a fixed later source suffix by the inherited strict-cone theorem.

**Proposition 3.** On a cube of k>=3 independent retained-root bits, no
permutation and independent complementation of the bits converts `2U-D`
into such a physical fixed-gate scalar identity, even if one permits an
arbitrary additive affine-linear correction to the digit expression.

Under a bit complementation the mixed coefficient on ij is multiplied by
sigma_i sigma_j, with sigma_i in{-1,+1}; permutation just relabels edges.
The product of coefficient signs on any triangle is therefore invariant.
It is+1 for the signed response and-1 for the physical response. Every
coefficient is nonzero. This proves the obstruction. Affine-linear terms
have no mixed coefficients and cannot change it.

This is a switching/frustration invariant of the **digit interaction graph**,
not a tournament invariant. It does not cover input-dependent piecewise
changes, extra root births controlled by input bits, duplicated variables,
or a source construction changing the cube itself. In particular it does
not refute signed-response coverage or a future genuinely structural
converter. It only rules out the tempting bare complement/reversal repair.
