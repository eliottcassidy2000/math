# Three path-bin polynomials and a two-source seam

**PROVED, for the stated operation language.** For a nonempty rooted
tournament (T,r), define P_T(m) as the number of ordered m-tuples of
possibly empty directed paths partitioning V(T). For m>=1 define E_T(m)
and S_T(m) by requiring the first path to be nonempty and respectively
end or start at r. The root is a physical vertex.

P is an integer-valued polynomial of degree n, and E,S have degree at most
n-1. For example expand over nonempty path covers and their assignments
to the available bins; fixing the root bin leaves m-1 available labels.
They have integral coefficients in the Newton basis binom(m,j).
P(0)=0. The *formal* values E(0),S(0) need not vanish: for a singleton
P(m)=m, E(m)=S(m)=1. Rooted zero/negative values are polynomial evaluations,
not counts of covers into a nonexistent distinguished bin.

## C1. Polarize the reflection construction

Take disjoint nonempty rooted tournaments (X,a),(Y,b), retain X, take a
converse copy Y*, and add c. All arcs between copies point Y* -> X. Set
a->c, c->X-a, c->b*, Y*-b*->c. Call this actual tournament Z=X diamond Y,
rooted at c. Then

    P_Z(m) = m[E_X(m+1)E_Y(m+1)+E_X(m+1)P_Y(m)
                 +P_X(m)E_Y(m+1)
                 +(P_X(m)-S_X(m))(P_Y(m)-S_Y(m))],
    E_Z(m) = E_X(m)P_Y(m)+P_X(m-1)(P_Y(m)-S_Y(m)),
    S_Z(m) = E_Y(m)P_X(m)+P_Y(m-1)(P_X(m)-S_X(m)).       (4)

Proof: select the root bin, then split at c. Every ordinary bin traverses
Y* before X. In the root bin, distinguish whether X-before and Y*-after
are nonempty. If both are, they respectively end at a and start at b*;
each child uses m+1 slots, giving the first term. If exactly one is
nonempty there is only the corresponding marked restriction. If neither
is nonempty, exclude a as the start of X-after and b* as the end of
Y*-before. These are precisely the four disjoint terms. For a bin ending
at c, either its X part ends at a or is empty and its Y* part does not end
at b*. This gives E_Z; reversal and swapping X,Y gives S_Z.

If supplied Hamilton paths of X and Y end at a and b, their concatenation
through c, using the reversed Y path, closes by a Y*->X arc. Hence Z is
strong. No additional strongness decision is hidden in the formula.
Moreover E_Z(1)=E_X(1)P_Y(1)>0 and S_Z(1)=E_Y(1)P_X(1)>0, so repeated
seams retain the required endpoint guard.

This is literally an actual complement line: put A=K_p(X), B=K_q(Y),
form A -> c -> B*, and use the path p,c,reverse(q*). Its actual complement
is Z. For X=Y and equal frames this is Session158's blue reflection line.
For different inputs self-converseness is not asserted.

Set h=P(1),s=S(1),F=E(2). Then

    H(X diamond Y)=F_X F_Y+h_X F_Y+F_X h_Y+(h_X-s_X)(h_Y-s_Y).

In coordinates u=h+F,v=h-s,w=h this is u_Xu_Y+v_Xv_Y-w_Xw_Y.
The coordinate change is integral unimodular. This is a literal polarized
quadratic response, not an assertion that every integral Lorentz vector
is realized by a tournament. Its diagonal recovers the inherited sandwich.

## C2. Births and the graded carrier

Adjoin x with x->r and T-r->x, then mark x. Its exact polynomials are

    P_new(m)=m[S(m+1)+P(m)-E(m)],
    E_new(m)=P(m)-E(m),
    S_new(m)=S(m)+P(m-1).                             (5)

Delete x in its bin: a nonempty suffix must start at r; an empty suffix
requires that its predecessor is not r. Endpoint restrictions give the
other formulas. The opposite birth follows by exchanging E and S before
and after (5). Converse itself exchanges E,S and preserves P.

Equations (2) in FUNNELS, (4), (5), and converse close on three polynomials
of unbounded degree. They compose as polynomial circuits with shifts,
products and sums. They do not close on three scalar values. Incoming
GRADED_OBSERVER already proves that both polarities force unbounded
linear memory over all strong sources; this construction supplies a
specific bin-polynomial quotient for sole-new-root observation.

For an ordinal sum X->Y, P multiplies. With the root retained in X,

    E_(X->Y)(m)=E_X(m)P_Y(m-1),
    S_(X->Y)(m)=S_X(m)P_Y(m).

With the root retained in Y, the formulas are

    E_(X->Y)(m)=P_X(m)E_Y(m),
    S_(X->Y)(m)=P_X(m-1)S_Y(m).

Each bin traverses X before Y. Ending in X forbids Y in the distinguished
bin, and starting in Y forbids X there; that is exactly the m-1 shift.
These identities also preserve the actual selected root.

## C3. Retain the inherited Q4 and F49 phase exactly

On the actual source T at r the golden observer, defined through its
converse rooted profile, is

    v=(h,x,y)=(P_T(1), [S_T(2)-E_T(1)]/2,
                       [S_T(2)+E_T(1)]/2).             (6)

The inherited portal parity makes these integers. All five Q4 presentation
matrices therefore act on the extracted v with their original orientation
convention. The field coordinate remains

    alpha=4(h+2x-y)+y theta in F49,
    theta^2=4theta+2,    g=h-y mod7.

The inherited q2 command multiplies alpha by theta and preserves g;
Frobenius sends theta to 4-theta. A new seam or funnel is evaluated over
the integers first, then mapped by (6). Its action is not silently declared
to be multiplication by one fixed field element. In particular a closed
polynomial calculus does not make the old 31/95/97 phase controls vanish.

## C4. Preserved information and remaining operations

The carrier preserves all future H responses under these specified rooted
joins, births, funnels, ordinal sums and converse, by substitution in the
proved laws.
It forgets individual path supports, source isomorphism, and arbitrary
old-root re-marking. Retain literal source constructions, root names, and
supplied path frames when using the metagraph realization in C1. No claim
of closure under arbitrary actual complements or minimal general strong-
atom attachment data is made. The joint correction in FUNNELS shows where
even a modest operation enlargement asks for more information.
