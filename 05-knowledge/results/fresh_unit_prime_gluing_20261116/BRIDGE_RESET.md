# Width-four window: exact obstacles to a reusable bounded deposit

This note uses the reverse-bridge convention of Session radix_three_shields:
join a nonempty left tournament U to a nonempty right tournament V by
all forward arcs except y→x, where x is the retained left port and y is
the incoming right port. Write H,E,O,g for the old count, the weighted
even/odd root-start path sums, and H(U−x). Empty complements count as 1.
For V write h=H(V), h_-=H(V−y), K=P_V(2)/2, and E_s,O_s and E_e,O_e
for weighted root-start and root-end sums. Both optional linear slots
in P_V(2) may be empty, but neither can be empty simultaneously for
nonempty V.

## Root change that recovers the modular freedom

The scalar output, independent of the chosen next port, is
`H'=hH+2(O_e E+E_e O)`. Put X=E+O and D=E−O. Keeping the old port,
the continuation is

```
X' = P_V(2) X − (E_s+O_s)g,
D' = (E_s−O_s)g,
g' = hg.
```

Keeping instead the right port y, when V−y is nonempty, gives

```
X' = (E_s+O_s)H + P_(V−y)(2) X,
D' = (E_s−O_s)H,
g' = h_- H.
```

Indeed, a run starting at y either stays in V, with left complement
factor H, or crosses immediately to x. In the latter case there
is a run in U starting at x, followed by an optional run in V−y.
The complementary full forward join factors. The weighted
parity totals of the optional V−y run coincide: the ordered
two-slot sum has zero alternating difference by the standard
odd-cycle packing cancellation. The switched join has disappeared
from the complement, because y has already been consumed. When
V is a singleton, the optional run is empty, so it instead gives
X'=H+X, D'=−H−D and g'=H.

For V a directed triangle, its profiles are E_s=E_e=1,
O_s=O_e=2, and P_(V−y)(2)=4, so literal right reroot gives

```
(H',X',D',g')=(3H+3X+D, 3H+4X, −H, H).
```

Modulo 3 this rotates (H,D)↦(D,−H) and preserves X. For instance,
the transitive edge rooted at its first vertex has (H,X,D)=(1,2,0),
so take instead a tournament with H and D both units to get a
perpetual unit pair. Such an input is the transitive 3-vertex
tournament rooted at its last vertex, (H,X,D)=(1,1,−1).
This escape is physical, but it does not give a low independent digit.

## A nonnegative transfer and a no-return theorem

Take positive coordinates (H,E,B,g), where B=O−g≥0: B excludes the
singleton root-start path. In these coordinates both port policies
have entrywise nonnegative linear transfers. For retention:

```
E' = K(E+B)+(K−O_s)g,
B' = K(E+B)+(K−E_s−h)g,
g' = hg.
```

The K for either first-slot parity includes all paths, not just
those starting at y; K≥O_s and K≥E_s+h since the even empty
first slot contributes h outside the rooted paths. For right
reroot with V−y nonempty, put K_-=P_(V−y)(2)/2:

```
E' = E_s H+K_-(E+B+g),
B' = (O_s−h_-)H+K_-(E+B+g),
g' = h_- H.
```

O_s≥h_- by its singleton root term. For singleton right reroot
the exact formulas are E'=B+g, B'=E, g'=H. The scalar formula
is nonnegative in all these coordinates and has coefficient
h≥1 on H.

**Consequence.** Apply a common word of such bridges to two
inputs that differ only by a positive scalar increment in H.
Their difference remains coordinatewise nonnegative and has
strictly positive H at every step. The first right reroot
transfers it into strictly positive g; subsequent steps
retain that strict positivity in g, whether they retain the
port or use a right reroot. So no word containing a right
reroot can transport even a two-value H-window into another
window with identical continuation E,O,g (or E,B,g). This
statement is structural and applies to arbitrary right counts
and arbitrarily long words. It relies on using the same word
on both inputs and on reading the full continuation.

## Fixed-port 3-adic obstruction at arbitrary macro length

Suppose the root is retained. Count the right modules of count three; other right modules are allowed as intermediate steps.
A right count-three tournament consists of one directed
triangle and singleton strong components. Its P_V(2) and
start total E_s+O_s are multiples of 3, independently of
the port. Consequently after t count-three bridges,
interspersed with any number of other fixed-root bridges,

```
3^t | X,g,       3^(t−1) | D   (t≥1).
```

Proof by induction on bridge letters: a count-three bridge
raises simultaneous divisibility of X,g by one, while D'
uses old g; any other bridge preserves the simultaneous
X,g divisibility and does not worsen D. On such a state,
any fixed-port word with scalar multiplication coefficient
M (product of its right counts) outputs H'≡MH modulo 3^(t−1).
In particular two alternative words with M=3 give a
difference divisible by 2·3^(t−1), since both physical
counts are odd. There is no fresh +2 at height t≥2, even if
the alternative words have unequal lengths and grow
unboundedly in order. No equality of their ending profiles
is needed.

Together these results exclude the common-macro strategy:
for each digit choice let a word act on *all* old settings,
insist on identical full continuation, and let the pure
transport coefficient be 3 per stage. A right reroot can
break the 3-adic trap but cannot synchronize the old settings;
keeping every port cannot deposit a unit at depth. A
candidate escape must bring in a new ingredient, e.g. words
chosen according to old settings with cancellation across
the resulting g offset, a coarser proven future interface,
or a different physical splice. This is not an obstruction
to arbitrary ports, bridges, or Chan–Pak.

## No matched one-step right-reroot count-three companions

In fact two right count-three choices using right reroot
cannot give different scalars on the same input while matching
all continuation coordinates. First, any physical X≥H by
splitting each Hamilton path just before its port; the suffix
and complementary prefix inject into the paths counted by X.
If g' matches, the two h_- agree (1 for a triangle root,
3 otherwise). With fixed h_-, the possible P_(V−y)(2)
are in a geometric sequence of ratio 2. Each start sum
a=E_s+O_s satisfies 0<a≤P_(V−y)(2) by assigning the part
after y and the complementary part to two slots of V−y.
If p_1<p_2, then (p_2−p_1)X≥p_1 H>(a_1−a_2)H,
so X' cannot match. Thus p and a match, and D' forces the
same signed start difference as well.

For completeness, a right root in a triangle with u/v
singletons before/after has a=3·2^v and end total 3·2^u,
with end difference −1 if u=0 and 0 otherwise. Equal p fixes
u+v, so a fixes u and the entire end profile. For a singleton
root let u/v count other singleton components before/after.
Before the triangle one has a=12·2^v and end total 3·2^u,
with end sign −3 when u=0 and 0 otherwise; after the triangle
one has a=3·2^v and end total 12·2^u with sign zero. Within a
type, fixed p,a fix u,v. For mixed types with matching a and
fixed p the after-triangle choice has two fewer preceding
singletons than the before-triangle choice, so the totals
coincide and the before-triangle sign is zero. Thus even
mixed types give the same end profile and the same H'.

## Mixed one-step policies: exact +4 but never +2

Allow two count-three companions on the same physical input,
one retaining the old port (F) and one rerooting to the incoming
port (R). There can now be a genuine synchronized choice, but
not a primitive +2. Use p_F=P_V(2) for F and p_R=P_(V−y)(2)
for R; let a denote start total, c end total, d signed end
difference. Equal X' requires

    (p_F−p_R) X = a_F g+a_R H.

Equal g' requires H=g for an R singleton port or H=3g for an
R triangle port. Put q=p_R/p_F<1 and use X≥H, |D|≤X. The
difference in H' is (c_R−c_F)X−(d_R−d_F)D. Every c is
divisible by three, and the only d not divisible by three
is −1 at a triangle port with no preceding singletons.

* If H=g, R uses a singleton port. A gap not divisible by 3
  needs d_F=−1, so a_F/p_F=1/4. Powers of two in the slot
  counts force q≤1/2, and 1−q≤1/4+q forces q=1/2. Then
  X/g=1/2+a_R/p_R. The latter ratio is a nonpositive power
  of 2. To get X/g≥1 it must be 1 or 1/2; the larger gives
  X/g=3/2, impossible for odd g. So X=g, D=−g, and R is the
  singleton-before-triangle case with exactly one singleton
  preceding its port. Hence (c_R,d_R)=(6,0), while
  (c_F,d_F)=(3,−1). The gap is precisely 4g.
* If H=3g, R uses a triangle port and q≤2/3, with next
  smaller possibilities at most 1/3. If d_F=−1 and d_R=0,
  the start ratios are a_F/p_F=1/4 and a_R/p_R≤3/8.
  Therefore 3(1−q)≤1/4+9q/8 forces q=2/3 and equality,
  i.e. X=3g and exactly one singleton preceding the R port.
  The absolute gap is at least 2X≥6g. If d_R=−1 and d_F
  is divisible by 3, then c_R=3. When c_F≥6 we have d_F=0,
  so again the gap is at least 2X≥6g. The only remaining
  possibility for c_F=3 is the singleton-before-triangle
  case with no prefix, having d_F=−3. Here
  a_F/p_F=1/2, so (1−q)X/g=1/2+(9/4)q. Then
  q≥10/21 forces q=2/3 and X=6g. Since D has the parity
  of X, the remaining gap −2D is divisible by 4, not ±2.

Thus no matched single-step pair of count-three companions
using at least one new entrance port can deposit +2. The +4
case is sharp: for U a singleton, F = C3→TT2 at the triangle
port gives H'=5, and R = TT2→C3 at the second singleton gives
H'=9. Both end with (X',D',g')=(36,0,3). This gives a small
synchronized signed move, not a repeatable window.
