# The infinite entry adapter: joint recurrence, integral guards, and phase

**PROVED with exact finite certificates in a proved response space.**
Identifiers JR1--JR4 are scoped to this file. Let D_m and X_m=N^3(D_m,w)
be [JE5](ENTRY.md), m>=1. Name the original C3 vertices a,b,c, the old
mark w, the m-vertex transitive module Z, and the three new vertices A,B,C.
The seven vertices outside Z have a fixed relation at every m.

Bulk insertion is a literal birth: add the next Z vertex beating exactly
a,b,c, and losing to every other existing vertex. All old names and arcs
remain. In the saved canonical bit rows, the old A,B,C indices move up
one and the newborn takes their previous first index; the receipt retains
this role-to-physical-index permutation. It is a coordinate change, not
permission to rename or discard an old physical port.

## JR1. Forty-eight actual joint responses have a three-state recurrence

For each pair (a,w),(b,w),(c,w), retain the sixteen fields in the inherited
endpoint convention:

    H, P(2), [E(1),S(1),q,E(2),S(2)] at each port,
    W2_aa,W2_aw,W2_wa,W2_ww.

W2 uses the end of bin0 and start of bin1; its repeated-port entries are
zero. Empty bins are allowed unless an endpoint requirement forbids them.
For all m>=1 the combined 48-vector is

    y_m = C (1,2^m,3^m)^T,                                (10)

where the integer 48-by-3 matrix C is exported in
`response_recurrence.json` (the first three coefficient columns; the
remaining six are zero). Thus

    y_(m+3)=6y_(m+2)-11y_(m+1)+6y_m.                      (11)

In particular

    H=18*2^m+15,
    E_a(2)=76*3^m+96*2^m-9,
    E_b(2)=76*3^m+96*2^m-27,
    E_c(2)=96*2^m-9.                                     (12)

**Why nine samples prove an all-size statement.** Contract maximal Z
runs in any counted ordered path cover. Its exterior skeleton has at
most seven exterior vertices and two bins, hence at most nine run slots.
For r runs the m ordered module labels can be distributed in exactly
r! S(m,r) ways; within each run their order is forced. All endpoint
constraints concern exterior vertices and remain literal on the skeleton.
Thus every field has a finite expansion sum_(r=0)^9 B_r r! S(m,r).
For m>=1 this is an integer linear combination of 1^m,...,9^m.
Deletion and one-bin fields have a smaller bound and belong to the same
space. The m=0 boundary may carry a separate 0^m contribution.

The integer subset-DP endpoint counter supplies all 48 values at m=1..9,
on actual graphs of orders8..16. Exact rational elimination of the
nine-by-nine power matrix gives C and zero coefficients at4,...,9.
Its determinant is nonzero (a Vandermonde determinant times the product
of its nonzero bases), so no other sequence in the proved space fits.
An independent gap counter on the actual order27 graph, m=20, checks
all48 values without supplying a fitting sample. The full certificate
and field ordering are saved, including the exceptional m=0 table.

This uses the familiar maximal-run substitution law, recovered in
[Session188, skeleton PROOF](../session188_structural_bridges_20261104/skeleton/PROOF.md).
The new content is its application to this fixed, named joint interface
and the explicit integral response matrix, not a claim to invent the
Stirling transform. This is a certified finite determination in a proved
finite-dimensional space, not an unconstrained recurrence fit.

Every quantity in (12) is divisible by3 for m>=1. Consequently wrapping
X_m at any of these three restored unique-start marks, with any positive
prefix clock, still gives a multiple of3. This particular interface alone
does not supply prime or7p lanes with p!=3. A successful return must alter
this interface or its incidences; the entry adapter by itself is insufficient.

## JR2. Integer observation lattice and its exact operation scope

For arbitrary boundary vector u=(u_1,u_2,u_3), the rows of C generate

    Z x 6Z x 38Z.                                       (13)

Indeed every second column entry is divisible by6 and every third by38.
Three exact row combinations attain the generators:

    S_a(1) = u_1,
    E_a(1)-5 S_a(1) = 6u_2,
    P(2)-3E_a(2)-16E_a(1)+107S_a(1) = 38u_3.              (14)

These are checked as integer row identities in the complete C. Its Smith
factors are (1,2,114), since gcd(6,38)=2 and lcm(6,38)=114. The prime19
is an actual observation divisor here: it comes from38, not a numerical
analogy with an unrelated nineteen-object construction.

Under bulk insertion in Z, u maps to diag(1,2,3)u. The lattice (13) is
stable under this update. Therefore for arbitrary modular boundary weights
the instantaneous observation kernel is already closed under all future
bulk insertions. For modulus M, the exact number of distinguishable states is

    M * [M/gcd(M,6)] * [M/gcd(M,38)]
      = M * [M/gcd(M,2)] * [M/gcd(M,114)].                (15)

This is the smallest state carrier for **these linear weights, these48
observations, and forward bulk insertion only**: equal observations remain
equal, and different instantaneous observations already distinguish states.
Over Q the minimal linear dimension is three, by (14) and the three
independent exponential modes. On the actual family a single unbounded
integer m reconstructs everything; a linear-dimension bound does not
assert a nonlinear memory lower bound.

## JR3. A free coordinate at19, a size-supplied bit at7

Actual sources impose u_1=1,u_2=2^m,u_3=3^m. Modulo19 the response loses
the third coordinate for arbitrary weights, but on this realized orbit

    u_3=u_2^13 mod19,                                    (16)

because 2^13=3 mod19. The observable E_a(1)-5 recovers6u_2; 2 has order18
mod19. Thus the actual bulk trajectory has exactly18 distinguishable
states, and its otherwise invisible coordinate is recovered at no extra
state cost. This inverse requires the realization guard; it is false for
arbitrary boundary weights.

Modulo7 the full response instead has exactly six states, with

    3^m=(-1)^m (2^m)^2 mod7.                             (17)

The rooted golden triples (H,(S(2)-E(1))/2,(S(2)+E(1))/2) at all four
tracked marks a,b,c,w use only1 and2^m, as their rows in C show. Both
the inherited F49/charge codec and F343 codec therefore have period three
along this family. They lose the parity needed for the full joint response.
For m=1 and4 the actual graphs have identical charts at all four marks,
but E_a(2) is411 and7683. Birth a vertex beaten only by a. The resulting
Hamilton counts are461 and7985, which are6 and5 modulo7. `repair.json`
and `response_recurrence.json` retain the actual graphs, fields and birth
cuts; this is a legal one-command distinction, not a formal weight change.

The examples have different orders,8 and11. They do **not** refute a
carrier that already retains size. In this seven-exterior family the
known size gives m=n-7, so the missing parity bit is free from that
retained clock. Without size, phase-only compression must keep this bit.
After an arbitrary complement, the seven-exterior realization guard may
fail; size alone then licenses no such recovery. The size clock, the
quadratic Frobenius, and the cubic Frobenius remain distinct operations.

## JR4. An older kernel becomes visible after two actual complements

Let z be the first vertex of Z. In D_m its sole in-neighbor is w. The
only old-to-new arc of X_m is w->A, and C is a universal sink of
{B,C} with C->A. These are exactly the hypotheses of the inherited
[Session103 bridge--sink star lemma](../session103_words_necklaces_20261029/bridge/BRIDGE_SINK_STAR_LEMMA.md),
with R=D_m, Q={A,B,C}, its a=z, s=q=w and r=A. Hence its ordered
two-nonempty-path matrix satisfies

    M2(X_m) [(e_z-e_w) wedge e_A] = 0.                    (18)

Apply the literal RP2 repair reversing C->A, obtaining Y_m. Then

    [M2(Y_m) ((e_z-e_w) wedge e_A)]_(z,B) = 1.             (19)

For the positive term with endpoints(z,A), the first path is the singleton
z. The other is B,C, followed by the unique Hamilton path of D_m-z ending
at w, then A. That old path is the remaining transitive Z block followed
by a,b,c,w. The term with endpoints(A,z) is impossible because z is
already the first start. Both terms using endpoints w and A in different
paths are impossible: after repair A has sole in-neighbor w. This proves
(19) for every m>=1. The inherited proof supplies (18) at every size.

The same explicit two-complement guards apply using the cycle
(w,Z,a,b,c,w). `star.py` checks them at m=1,2,4,9,20 and independently
enumerates the complete star action by permutations and cuts at m=1.
Thus a fixed integral kernel direction really becomes observable under
actual complements at every family size. It concerns weighted port
responses, not two distinct tournaments proved to share an entire carrier.
The three-state bulk observation result does not authorize dropping such
directions when the operation language includes these complements.
