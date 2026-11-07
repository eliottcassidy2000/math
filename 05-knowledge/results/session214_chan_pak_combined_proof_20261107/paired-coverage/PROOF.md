# The exact common language of the four staggered paired carriers

**Status: PROVED candidate, awaiting independent audit.** For each of the four
published physical carrier pairs at bases 677, 759, 1547 and 821, the common
O/Z/J words that preserve equality of the two adjusted bases are exactly
the words over O,Z. Once the first common J is used, no common suffix can
restore equality. This includes arbitrarily long, target-dependent suffixes.
It does not concern different words on the two sources, changing their fixed
start labels, replacing a source, new carriers, or new graph operations.
The Chan--Pak all-odd spectrum and F(8,137) remain OPEN.

## 1. Inheritance and notation

The eight explicit rooted graphs and their two-residue interpretation are
in `dependencies/PAIRED_CARRIERS.md` (Session203). The all-input physical
nine-coordinate transfer is in `dependencies/ROOTED_CLOSURE.md`. Both were
read from current main at commit `4c63fd58b566f3e006eb87eaedcf48824741cc84`;
the frozen copies have Git blob hashes `ec6f992a72efd909f0af5e05cbabaf26e6aa2563`
and `4360809d176a820160e2d201c464bcba920d103f`, respectively.

Use coordinates (h,s,e,L,R,d,C,b,T), with their path-count meanings from
the dependency. In particular B=h+L-s, and for residue r=1 or 2 the fixed
start s=r+3j gives adjusted base

    X=B+2j=h+L-(s+2r)/3.

Subtract source B from source A in each named pair; write delta for this
difference and D=delta X. The j values stay fixed under all three moves.
Hence D=delta(h+L) for pairs 759 and 1547, while
D=delta(h+L)+1 for pairs 677 and 821.

The initial difference vectors, in the displayed nine-coordinate order, are

    677:  (-10,-1,0, 9,-12,0, -88,12,-57),
    759:  (  0, 2,0, 0,  0,0,-168, 0,-84),
    1547: (  0, 2,0, 0,  0,0, -24, 0,-12),
    821:  (  2,-1,0,-3,  0,0,   8, 0,  3).

All four have D=0 initially. The inherited paired-carrier theorem proves
D=0 after every O/Z word. The new assertion is the converse even after
arbitrary later attempts to repair a J.

## 2. A useful nonnegative cover block

Put A=C-2L and W=T-b. The four cover coordinates (L,A,W,b) satisfy

    O: (L, 2L+2A, 3W+b, 2b),
    Z: (L+b, 2A+2W, 3W, 2b),
    J: (L+A, 2W+2b, W+b, 2L+A).             (1)

These identities follow by substitution in the inherited nine-coordinate
law. Their matrices have nonnegative coefficients. Both signs of the
coordinate orthant are therefore preserved. If L is strictly positive
(respectively negative) and the other coordinates have the same weak sign,
its sign remains strict forever, since every new L contains old L.

For completeness, actual L,A,W,b are nonnegative: C>=2L by the disjoint
choices of the labeled bin ending at the root, and T>=b by making the
designated root-starting path the singleton root. These interpretations
are not needed to treat differences, which can have either sign. The block
is a small change of the prior positive cover block, not a new cover law.

For the ordinary block put

    a=h-e-d, f=R-h+e-d, p=d-e.

Its algebraic maps in order (a,f,e,p) are

    O: (a+e, a+2f+e+2p, 0, e+p),
    Z: (2a+f+e+p, f, e+p, 0),
    J: (f+p, 0, a+e+p, e).                 (2)

They too have nonnegative coefficients. No claim that this cone contains
every rooted tournament is made or needed. We use it only for a particular
difference vector after its first J. Here h=a+2e+p; a weakly nonpositive
vector with h<0 retains h<0: O keeps h, and Z/J add
R-e=a+f+e+2p<=0.

## 3. The two purely cover-separated pairs

For pair 759 or 1547 the ordinary coordinates h,e,R,d have zero difference
and remain equal forever. Initially their cover difference is
(L,A,W,b)=(0,-2t,-t,0), with t=84 or 12.
Every O/Z prefix retains L=b=0, A<0 and W<0 by (1).
Its first J creates L=A<0 with all cover differences nonpositive.
Every later suffix retains delta L<0. Therefore D=delta L<0 after
every word containing J. No resynchronization is possible.

## 4. Pair 677

Along every O/Z prefix write u=delta h and v=delta R. Direct use of the
transfer maintains

    v<u<0, delta e=delta d=0,
    delta L=-1-u>0, delta b=-v>0,
    delta C<0, delta T<0.                   (3)

At the empty prefix these are the displayed initial values. O fixes u
and doubles v; Z sends (u,v) to (u+v,2v), preserving v<u<0.
For the remaining strict signs, O scales C,T by 2,3; Z sends
C to 2C+2T-2L and T to 3T-b, with all quantities interpreted as
differences. This proves (3) inductively.

Immediately after the first J, the ordinary difference (a,f,e,p) is

    (v-u,0,u,0),

which is nonpositive and has h=u+v<0. By (2), every later suffix has
delta h<0. The cover difference after that J is

    (delta C-delta L, 2 delta T, delta T, delta C),

strictly negative in all four entries by (3). Thus delta L<0 forever.
The integral quantities then give
D=delta h+delta L+1<=-1<0 at every later stage.

## 5. Pair 821

Every O/Z prefix preserves

    delta h=2, delta L=-3,
    delta e=delta R=delta d=delta b=0,
    delta C>0, delta T>0.

The last signs follow because O scales C,T, while Z sends them to
2C+2T+6 and 3T. After the first J the ordinary difference
(h,e,R,d) equals (2,2,2,2), twice the actual singleton profile.
The ordinary transfer laws hold also at a singleton; the only singleton
exception in the nine-coordinate law is the unused start coordinate of Z.
Consequently after any common suffix this ordinary difference is twice an
actual tournament profile. In particular delta h>=2.

The cover difference after the first J is

    (delta C+3, 2 delta T, delta T, delta C),

strictly positive. Equation (1) preserves delta L>0. Therefore
D=delta h+delta L+1>0 forever.

Sections 3--5 exhaust the four pairs and prove the exact accepted language.

## 6. Consequence and loss ledger

For these fixed eight sources and fixed residue labels, the richer common
alphabet O/Z/J supplies no new synchronized base at all. Thus the previously
proved four O/Z arithmetic images and their sparse upper bound remain the
entire output of this common-word paired strategy. This is stronger than
noting a one-step J failure, and differs from Session203's criterion for
synchronization under *every* word: here each individual candidate word
containing J fails, including words chosen after seeing its intermediate
state.

Source: an ordered pair of actual rooted graphs with fixed starts and j.
Map: append the same literal word to both and apply the four-vertex compiler.
Preserved: physical graphs, starts, j bounds, and the signed separation above.
Lost if only X,e,R+b,d are retained: the cover difference which becomes
visible at J. The required sidecar is the signed four-cover block.
Escape routes not excluded: asynchronous words, new source pairs, multi-entry
operations, and replacing either source rather than extending both identically.

`verify.py` reconstructs the eight seed profiles by subset path counting,
checks the difference table, checks (1)--(2) symbolically, tests every common
word through a declared length, and directly recounts small literal children.
These controls do not replace the invariant proof above. No numerical
coverage interval is promoted to an all-order occupancy claim.
