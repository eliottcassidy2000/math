# Paid triangle repair, absorption, and its exact clock boundary

**PROVED, scoped all-size operation.** Local identifiers RP1--RP3 belong
to this file. This operation restores entry for an infinite family but
is not the missing ternary return.

## RP1. Keep the old mark and break the new triangle

Start from any nontrivial strong (D,w), with
(h,s,q,F)=(H(D),start_w(D),H(D-w),E_w(D,2)). Form the inherited N^3 dial
X on new a,b,c, then reverse just c->a. Call the result Y and keep w.
The three new vertices now form the transitive chain a->b->c; w->a is
the only old-to-new arc. Then

    (h_Y,s_Y,q_Y,F_Y)=(h+4F, q, q, 8F-q).                  (6)

**Proof.** A Hamilton path with all new vertices before all old vertices
contributes h. In every other path, its only old-to-new crossing is w->a.
The old vertices split into an ordered two-bin cover with the first path
ending at w. The new block after w contains a and any subset of {b,c};
the other new vertices form the initial block. All four choices are legal
and exhaustive. They contribute 4F.

Deleting w leaves TT3 -> (D-w), so q_Y=q. A path starting at w must go
through a,b,c before D-w; hence s_Y=q. For two-bin covers ending their
first bin at w, distribute a,b,c among the two bins in their chain order,
ahead of their respective old paths. There are 8F possibilities, except
when the first old path is the singleton w and its new bin is exactly
{a}; the edge a->w is absent. These q configurations are the only losses.
Thus F_Y=8F-q. Strongness follows from the old strong component and
w->a->b->c->D. QED.

This is literally the inherited zero-ready wrapper Z_2(D,w), with a as
its outer root and b,c as its prefix, now viewed as a component at old w.
The scalar h+4F is inherited from the wrapper law. The new use is a paid
entry repair and a complete retained-mark recurrence; it does not create
a fresh realization of an already counted wrapper.

At q=1 the operation returns to the same unique-start/transitive-deletion
class indefinitely. Its one-apex word transformation is u -> 100u, where
bit1 means w beats the corresponding transitive-chain vertex. Starting
with C3 gives

    (h,F)=(3,3) -> (15,23) -> (107,183) -> (839,1463) -> ... .

For general q it preserves the deletion tournament up to a transitive
prefix and makes its entry budget exactly q. Iterating gives

    F_j=(F_0-q/7)8^j+q/7,
    h_j=h_0 + 4(7F_0-q)(8^j-1)/49 + 4qj/7.                (7)

These rational-looking expressions are integers by the recurrence.
The scalar h_j satisfies the recurrence with characteristic
(X-1)^2(X-8). The visible linear term records the fixed deletion budget.

For the actual original component at name31,
(h,s,q,F)=(2195,1,53,42153), the repaired/absorbed component has
(170807,53,53,337171) at that same name. It is the known full H170807
graph under an explicit change of component boundary; it is not ready
with budget one at31. The original outer root39 remains a different mark.

## RP2. Uniform legal native guards

Choose a directed Hamilton cycle (w,u_1,...,u_m,w) of D; m>=2. The word

    P=(b,c,a,u_m,w,u_1,...,u_(m-1))

is a Hamilton path of X and contains c->a. Birth a temporary helper z
that beats only c. Use the actual Hamilton-complement words

    P followed by z,
    P with z inserted between c and a.

The first word is legal since its last vertex is not c. In the first
off-path complement, the unprotected arcs z->c and a->z reverse; hence
the second word is legal. Two complements and deletion of z reverse
only c->a among permanent vertices. This is the inherited traced-arc
gadget, here with an explicit uniform path. The neighbors of w agree in
both words. In any zero-ready wrapper, prepend its outer root and its
ordered prefix to both words to keep the outer-root frame as well.

All old names and every old-old arc survive the completed repair.
There are three permanent births and two native complements, with one
temporary helper. The final component is strong. No invariant H or
strongness is claimed inside the helper program. Other ports retain their
physical names and actual incidences, but their adjacent path frames need
not be unchanged. A supplied Hamilton cycle suffices; existence is the
standard strong-tournament cycle theorem.

## RP3. Why this repair cannot itself supply either ternary branch

At any input prefix clock M=2^s, s>=1, the wrapped count is h+MF.
At an arbitrary output clock M'=2^t, t>=1, the repaired count is

    h+4F+M'(8F-q).

Equality with 3(h+MF)+delta, delta in {2,4}, would require

    (8F-q)M'=F(3M-4)+2h+delta.                            (8)

Always 1<=q<=h<=F; the middle inequality follows by canonical vertex
insertion and the last by cutting a Hamilton path immediately after w.
Strongness gives h>=3. The right side of (8) is strictly smaller than
(8F-q)M/2, because their difference is

    F(M+4)-Mq/2-2h-delta
      >= F(M/2+2)-delta >0.

For M>=4 it is strictly larger than (8F-q)M/4, since the difference is
F(M-4)+Mq/4+2h+delta>0. There is no power of two strictly between M/4 and
M/2. For M=2 the upper bound already gives M'<1. Thus (8) never holds.
This excludes the specified repair with any two prefix clocks, not other
source rewrites, mark switches, or compositions of several different repairs.

The inverse arithmetic step is guarded:

    F=(F_Y+q)/8,   h=h_Y-4F,   F_Y=-q mod8.

Integer divisibility is necessary; an actual inverse also needs the
named a,b,c incidences and a strong recovered D. Counts alone do not
authorize deleting three vertices or resetting a component.

Modulo7, (h,F,q) -> (h+4F,F-q,q) is unipotent and its seventh power is
the identity. Its orbit has exact period7 unless F=q=0 modulo7, when
it is fixed. This is an arithmetic recurrence clock. The retained F49
Frobenius has order dividing2 and the F343 Frobenius order dividing3;
neither is identified with this seven-step clock. Actual chart changes
are retained separately in the named-return certificates.
