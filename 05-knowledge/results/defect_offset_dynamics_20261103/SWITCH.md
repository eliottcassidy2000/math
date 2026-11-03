# A reusable joint-incidence switch invisible to the entire golden bank

**PROVED all-size source identities and traced implementation.** The finite
modulo-seven minimality certificate is **FINITE-EXACT**. This is not a claim
of closure under every native complement or a completed ternary constructor.

## S1. The two sources and their preserved observer

Let A_m have block sizes (m,1,1,1,2,2) over Q, and B_m have sizes
(m,1,1,1,3,1), m>=1. Keep the natural vertex labels; both have order m+7.
They differ in exactly one arc: with u=m+5 and v=m+2, A_m has u->v and
B_m has v->u. The vertex u has moved from the first position of block5 to
the last position of block4. Its other incidences are unchanged because
quotient vertices4 and5 differ against the remaining core only at vertex3.

At the separately named processing root p=m (quotient vertex1), both sources
have the **same integer** golden observer

    (H,X,Y)=(24*2^m+23,12*2^m,12*2^m+23).                (S1)

This follows from W6/W8: their right wings are (12,1) and (12,3), while the
left wing is (2^m+1,0). Its zero second coordinate hides the right-wing bit.
Thus all five Q4 matrix presentations, and every word in that bank, have
identical integer responses from the two marked sources. Their F49/charge
triples agree as well; no phase identification is being guessed from H alone.

## S2. A second selector sees the precise defect change

At root r=v (quotient vertex3), write P=2^m. The exact (h,s,q,F) profiles are

    A_m: (24P+23, P+1, P+1, 29P+29),
    B_m: (24P+23, P+1, P+1, 25P+25).                     (S2)

The all-size run coefficients for h,s,q,F are checked symbolically in
`wing_switch.py`; for example F has coefficients (87,29) for A and(75,25)
for B in the basis 1,2^m-2. Therefore their neutral defects are respectively
5P+6 and P+2. B_m is the literal word 11^3,00^m,01 in the old two-in-port
chart. In particular **B_1 is a defect-four source and B_2 is a defect-six
source**. The switch reaches three-in-port roots of defect16 and26 while
preserving the separate processing observer exactly.

This supplies an actual controlled enlargement of the defect families.
It does not make the defect change invisible to the joint (p,r) observer.
The missing datum is precisely an incidence that one selector could not see.

## S3. One legal birth exposes the hidden bit

Clone p by adding a vertex just after it in its transitive block: all old
in-neighbours of p and p point to the new vertex, which beats every old
out-neighbour of p. This changes the left block length from1 to2. The counts
are now respectively

    H(A_m plus clone)=72*2^m+1,
    H(B_m plus clone)=72*2^m+5.                           (S3)

Both targets are strong. So the entire integer golden observer, even with
all of its future Q4 responses, is insufficient once this actual birth is
allowed. Given m, one additional binary choice is necessary and sufficient
to distinguish these two sources and answer this query. At m=1 the pair is
H71 with outputs145 and149; at m=2 it is H119 with outputs289 and293.

More generally allow the left block length l>=1. Retain the right-wing bit
epsilon, choosing (12,1+2epsilon), and the left response (x,b)=wing(m,l).
The full chart is closed under either left block growth and under toggling
the same right-edge incidence, with

    H=24x+2(1+2epsilon)b-1.

The discrepancy between the two sources after t left-length growths from
l=1 is exactly4(2^t-1). The bit is hidden on the boundary l=1 and observable
after leaving that boundary. Parameters (m,l,epsilon), with their explicit
vertex-address transport, reconstruct the whole source and every joint query.
Modulo7, this chart has an exact minimal18-state carrier for Hamilton responses
under the two growths and the toggle: all18 future signatures are distinct.
The finite upper bound is (m mod3,l mod3,epsilon). It is not a claim about
arbitrary complement/re-marking closure.

## S4. Actual complement realization and the helper clock

In A_m the following Hamilton path contains the switched arc u->v:

    (m+3,u,v,m+4,m+6,0,1,...,m-1,m,m+1).                  (S4)

Apply the inherited traced-edge compiler: append one temporary vertex,
perform its two specified legal native path complements, and restrict to
the original named vertices. The resulting old relation is exactly B_m.
The two legal paths and all names are stored in `wing_switch.json`.
The peak order is m+8; the retained target order is m+7. Reverse the recorded
two moves to implement the opposite switch with its corresponding helper
incidence. The same construction works for every l after expanding block1
in S4, with peak order m+l+7.

**Restriction of the helper is explicit.** This is a traced source operation,
not a claim that two native moves at the original order reverse one edge.
Indeed a native complement reverses (n-1)(n-2)/2 edges. At n congruent to1 or2
modulo4 this number is even, so any sequence at fixed n preserves edge parity
and cannot make this single-edge switch. The helper/trace capability cannot
be silently dropped from a claimed closed carrier.

The pure native-plus-birth escapes in ESCAPES.md need no helper deletion.
Combining those with this parameterized joint switch is a concrete next target;
the present packet does not assume the two operation languages coincide.
