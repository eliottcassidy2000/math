# A rooted bridge that lets windows and the actual matrix readers interact

**Proved physical construction.** The existing bilateral window method
constructs the initial source; the exact Wave-20 matrices provide
different final readers. Their amplitudes must match the full rooted
profile, not just its Hamiltonian count.

## Rooting a window source at the apex

Let C be the old puncture with contact mask U. In the bilateral
notation (h,f,P,Q,S,T,K), f counts ordered complementary path pairs,
including empty parts, P requires a path to start in U, S is the
whole-path end in U, and K is the count with just the marked core and
apex. Insert transitive wings before and after C. Use positive binary
addresses x,y, with the right contact word the reversed complement of y.
Let m=2^(bitlength y), and rev(y) the reversal at that width.
At the literal apex the rooted orbit coordinates are

    W = m(f x+P+h) − h rev(y),           d = h.            (1)

For W, the path from the apex is either a singleton or begins with
a contacted first vertex on a nonempty complementary part. Under an
ordinal join of old punctures A before B, the P count is
P_A f_B+h_A P_B (split according as the root path has an A part).
The left wing has P=x and the right has P=m−1−rev(y). This proves W.

For d, the same sum has an alternating sign. For every nonempty B
the signed free complementary cover sum is zero:
sum over U subset B of (−1)^|U| H(B[U])H(B−U)=0. Use the
[odd-cycle-packing expansion](https://arxiv.org/html/2412.10572v3#S3.SS4).
For any packing on B, each odd cycle and each unused singleton is
assigned wholly to one of the two parts. Each such odd object has
signed choice 1−1=0. Thus every root path that meets A has zero
signed sum over the unrestricted B slots. Only root paths with
empty A part survive, and give h_A d_B. For the terminal transitive
wing its d_B is 1 minus the last contact bit: any earlier contacted
first vertex leaves a nonempty free suffix with cancelling parity.
The last bit is 0 by positivity of y, which proves d=h.

The first-coordinate row of any actual reader product can now be
applied to (H,W,d). This turns a scalar witness into a state for
both constructions rather than assuming its source port survives.

## Two actual reader words agree except for their d coefficient

The single block J OOOZ has first row (9,9,−1).
The construction-order pair (J O) then (J OZ) has row (9,9,−3).
Thus on every rooted source they give difference 2d, exactly.
For any transitive puncture whose apex sees the last vertex as an
in-neighbor, d=1, so this supplies a two-odd window with no restriction
on its previous source length. This does not assert that every such
initial source occurs within the capped initial run of the strict orbit.

## A three-window requiring both changing the source and changing its reader

For independently chosen powers u,v>=2, there are actual tournaments
of counts

    3456uv+459u+2430v+{319,321,323}.                        (2)

Take x=u,y=v so g(u)=u+1, g(v)=v+1 and all terminal parity terms vanish.
The source count is f uv+(Q+h)u+(P+h)v+K, and by (1) its W has
only 2fuv+2(P+h)v−h. If a reader has row (a,b,c), its coefficients become

    ((a+2b)f, a(Q+h), (a+2b)(P+h), aK+(c−b)h).

The actual three certificates (core rows, contact, reader) are

    319 (2,4,9,3,15,31)    46    J OOOZ
    321 (0,1,3,7,15,31,63) 89    J O ; J OZ
    323 (0,1,3,7,15,31,63) 89    J OOOZ

The first core is nontransitive, with h=5; the latter is transitive,
with h=1. The lower two consecutive repairs come from different
mechanisms: the nontransitive puncture changes both W and H, and the
last two use the identical source but switch the reader's d term.
At u=v=2 this gives 19921,19923,19925. It is a literal merged
construction, not a union of unrelated residue projections.
One tests prescribed N by removing each offset, enumerating the
possible u powers and exactly testing the quotient v=(N−delta−459u)/
(3456u+2430). These pairs still realize only O((log N)^2) values.

## A free-height window can survive a reader, but not every offset does

In the inherited [slope-908 four-point bank](../chanpak_absolute_window_reanchoring_20261114/SLOPE_908.md),
use its four certified punctures in the displayed order. Their
source H values are 908u+9g(u)+752,754,756,758. At the common
chosen apex d=9. Formula (1) gives W=1344u+1099 for the first
three and W=1344u+1105 for the fourth. A single J therefore
gives exactly the three consecutive odd outputs

    2252u+9g(u)+{1842,1844,1846},                    u>=1,

and a fourth offset 1854. The +8 hole is a precise example
of why an unrooted exact window is not automatically an exact
matrix-orbit window. It is not a nonexistence result for other
choices of the fourth source or reader.

Replay with python3 rooted_window.py, including literal graph
construction, rooted subset counts and final Hamiltonian counts.
