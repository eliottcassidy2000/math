# A second alternating round and the missing cover state in the old cap

**FINITE-EXACT:** an actual alternating23→73→221→667→2003 path with one
original root and one original port surviving every return. **PROVED:**
the exact recursive identification and cover response of the older
growing-width cap, and the failure of its repeated one-letter funnel to
produce an unlimited+2 ready return. Chan–Pak remains **OPEN**.

## 1. A direct+2 bridge across the667 ready obstruction

Use outgoing-bit rows and zero-based labels. The inherited (+4,+2)
construction's final H221 graph has physical names

    (P,x,y,a,b,r,u,v).

Here P,a,b,r came from the original H23 source; x,y were born on its+4
edge and u,v on the following+2 edge. The original two tail vertices
were the deletions already specified in the inherited instructions.

Prepare the2-ready221 graph

    B221=(34,108,105,113,103,64,129,63)

by retaining old labels(0,1,3,5,6,7), mapping them to(4,2,6,7,0,1), and
inserting the two missing vertices. Its labels are then

    (u,v,x,f,P,g,a,r).

Its marked root is still physical r, with sole in-neighbour physical a.
The count-preserving preparation trace is221→23→99→221. Perform two
standard births at r and its newborn; this is the lawful+4 step to667.

Prepare the unique strong2-ready667 class

    B667=(6,196,200,211,231,207,129,257,127)

by retaining labels(2,3,4,5,6,7,8) of that ten-vertex child and mapping
them to(1,2,3,6,7,8,5). Insert the two missing target vertices. Its names
are now

    (e,x,f,P,h,c,g,a,r).

The preparation trace is667→51→197→667. Physical r and a again survive
as the marked root and its unique in-neighbour. The return deletes the
most recent newborn along with two earlier vertices. It does not merely
rename a scalar667 witness.

There is no1-ready tournament of value667, even allowing weak sources or
converse; see [TARGET_RETURN §4](TARGET_RETURN.md). Nevertheless a different
literal operation gives3·667+2 with a1-ready target. Define

    B2003=(514,284,505,497,481,451,387,259,1,510).

Retain B667 labels(1,2,4,5,6,7,8), and map them to B2003 labels
(4,5,2,3,6,0,9). Delete old0,3, and insert the missing target labels1,7,8
with their indicated cuts. The exact trace is

    667 → 33 → 321 → 997 → 2003.

The final physical names are

    (a,j,h,c,x,f,g,l,m,r).

The marked root r and its sole in-neighbour a are still the same two
vertices from the original H23 graph. Its defect has changed from2 to1.
All three replacements preserve every comparison among retained vertices.
All completed sources in the displayed alternating path are strong;
the intermediate deletion cards need not be strong.

Thus a finite actual second alternating round exists:

    23 →[+4]73 →[+2]221
       →[same-count2-ready preparation]221 →[+4]667
       →[same-count2-ready preparation]667 →[direct+2]2003.

The first failed implication in the naive ready induction was that+2 at667
must use a1-ready source. This direct graph replacement bypasses that
specific interface while preserving the original root and port. Its
corrected scope is one explicit finite bridge, not closure of all such
sources.

### The cheapest uniform lift fails

For each of the9 vertices of B667 and each of the10 vertices of B2003,
compute the complete positive transitive-module run polynomial. None of
the90 pairs has target coefficient equal to three times source coefficient
at every run number j>=2. Therefore replacing one source vertex and one
target vertex by a common TT_k never makes this particular bridge satisfy
H_target(k)=3H_source(k)+2 for every k>=1. This includes pairs not compatible
with the retained-vertex map, so it is a conservative obstruction.

No inference is made about two moving modules, parameter-dependent cuts,
different667/2003 targets, or other growing cores. Exact equality at the
first count is not enough to certify an infinite source-return family.

## 2. Recovering the old cap in the new six-coordinate grammar

The older [CAP_WIDTH theorem](../../compensation_returns_20261103/CAP_WIDTH.md)
constructs actual D_j with H(D_j)=3^j, end_(x0)=1, and feedback vertex
number ceil(j/2). Its displayed order is

    x0,y1,x1,y2,x2,...,yj,xj,

with exactly the backward arcs xi→x_(i−1). Let C_j be its converse, rooted
at x0. Delete x0 from C_j. The first strong core is C_(j−1), on names
x1,y2,x2,...,yj,xj, followed by the singleton y1. The root x0 beats just
x1 in that core and loses to y1. Hence C_j is EXACTLY the one-port funnel
with tail word(1) applied to C_(j−1), beginning with the singleton C_0.
This identifies actual named graphs, not just the number3^j.

Apply the six-state compiler from [TARGET_RETURN §3](TARGET_RETURN.md).
For the one-letter word(1),

    h'=h+2s2,       e'=h,       s2'=h+2s2,
    p2'=4p2−2e2+6s3,
    e2'=2p2−e2,
    s3'=2p2+3s3.                               (1)

It follows that h_j=s2_j=3^j. The omitted cover coordinate e2_j satisfies

    e2_(j+3)=6e2_(j+2)+3e2_(j+1)+12e2_j,
    (e2_0,e2_1,e2_2)=(1,3,21).                  (2)

Indeed the last three coordinates in(1) have matrix

    [4 −2 6]
    [2 −1 0]
    [2  0 3]

and characteristic polynomial X^3−6X²−3X−12. Its first e2 values are
1,3,21,147,981,6579. Appending a vertex beaten only by the unique-start
root creates a1-ready tournament with count

    H_ready(j)=3^j+e2_j−1
              =1,5,29,173,1061,6821,... .       (3)

The j0 output is the weak H1 boundary; all j>=1 outputs are strong. The
parent's growing feedback width is inherited exactly from CAP_WIDTH.
Constant endpoint mass therefore coexists with both growing structural
width and a growing path-cover coordinate.

### No seed makes this fixed funnel program an unlimited+2 ready return

This failure is not confined to the singleton seed. Start from ANY fixed
nonempty unique-start rooted seed and repeatedly apply the same word(1)
funnel before taking its ready newborn. After the first step, s2=h, so
h_j is a positive constant times3^j. The e2 sequence still satisfies(2).
An eventual exact ready recurrence R_(j+1)=3R_j+2 would force
e2_j=c·3^j eventually. Substitution into(2) gives

    (27−54−9−12)c = −48c = 0.

But e2_j counts nonempty rooted path covers and is positive. Contradiction.
Thus no seed of this fixed funnel program supplies an unlimited exact+2
ready return. Other funnel words, changing word programs, and other
preparations are outside this negative theorem.

The missing datum was e2, not the endpoint count or H alone. The strongest
survivor is the exact six-state construction and the complete recursive
family(1)–(3), which now supplies a faithful source for a compensation
attempt. The outstanding obligation is to control its actual baseline.

## 3. Controls and scope

`python -B finite_bridge.py` checks every retained arc, actual appended
cut, root/port name, complete count trace, readiness and strongness, plus
all90 fixed-module pairs. `python -B cap_cover_state.py` independently
checks the cap/funnel graph identity through depth12, all six coordinates
by literal subset counts through depth6, the recurrence, ready newborns,
and the separate growing-defect birth orbit of TARGET_RETURN §5.

The infinite statements use the graph decomposition and exact transfer
identities; the finite controls are not extrapolated into an infinite
target-return claim. Neither construction proves Chan–Pak.
