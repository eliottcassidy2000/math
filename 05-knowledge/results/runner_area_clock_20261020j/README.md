# Runner count, polygon symmetry, and the area clock

**Status:** `PROVED` elementary clock and threshold reductions; `FINITE-EXACT`
literal controls; general LRC and the internal LRC(14) proof obligation remain `OPEN`. This packet changes
coordinates and tests their exact guards. It supplies no new common lonely
time for an arbitrary speed set.

Incoming literature qualification: [Allikvere, arXiv:2609.02604v2, Theorem 1.1](https://arxiv.org/html/2609.02604v2)
states the fourteen- and fifteen-total-runner cases. This is a `CITED`
external preprint result, not an independently reproduced computation in
this packet. The repository's internal LRC(14) status must not be read as
an unqualified claim about the current literature. General LRC remains
the target here.

## 1. Inheritance and the retained object

The closest existing mechanism is the AP witness orbit in
[THM-403 / cyclotomic-witness-orbit-of-the-AP](../../../01-canon/theorems/THM-403-cyclotomic-witness-orbit-of-the-AP.md).
For speeds $0,1,\ldots,n-1$, the critical times are precisely $a/n$ with
$\gcd(a,n)=1$. Section 6 supplies a short proof of completeness in the full
continuous clock, not just at those rational samples. The current anchor is
the [LRC(14) proof map](../../../00-navigation/LRC14-PROOF-MAP.md).

The old [regular-polygon reflection](../../../07-reflections/regular-polygon-shared-stage-lrc-twin-goldbach-s517.md)
is provenance for the geometric idea, not a proof that arbitrary runner
trajectories attain a regular configuration. The corrected near miss is the
difference between every polygon vertex and the primitive AP clock points:
at $n=6$ the only AP witnesses are $1/6,5/6$. The safe arithmetic content of
[THM-404 / doubling-rigidity-dichotomy-of-the-witness-orbit](../../../01-canon/theorems/THM-404-doubling-rigidity-dichotomy-of-the-witness-orbit.md)
is that doubling permutes the unit clock exactly for odd $n$, with cycles
the cosets of $\langle2\rangle$. Permutation need not mean a single cycle:
modulo 7 the cycles are $\{1,2,4\}$ and $\{3,6,5\}$. Its July 21 correction
also separates a bad period-14 group object in characteristic 2 from
ramification of $\mathbb Q(\zeta_{14})$.

The live concepts here are literal runner histories, the loneliness
threshold, boundary arclength, swept area, central angle, and the shared
clock. The least-used sidecars are the marked observer, the origin/focus,
and each point's radius. Their retention determines which maps below are
equivalences.

Kepler's second law supplies the physical motivation: equal time increments
sweep equal areas about the attracting focus; see
[NASA, Kepler's second law](https://pwg.gsfc.nasa.gov/stargaze/Skepl2nd.htm).
The polygon statements below are kinematic area identities. Polygon motion
with corners is not asserted to be an inverse-square gravitational orbit.

Write $d(x,y)=\min(\{x-y\},1-\{x-y\})$ for circular distance, with one lap
normalized to 1. For runner histories $x_i(t)$ and a literal observer $k$,
the target is

$$\min_{i\ne k}d(x_i(t),x_k(t))\ge 1/n. \tag{1}$$

## 2. Exactly when equal area is equal perimeter

Let $P$ be a strictly convex polygon, oriented counterclockwise, and let
$O$ lie strictly inside it. Mark the same initial boundary point for both
clocks. Let $s$ be boundary arclength, $L$ the perimeter, $A$ the area, and
$h_j>0$ the perpendicular distance from $O$ to side $j$. The swept-area
increment along that side is

$$dA_O=\tfrac12 h_j\,ds. \tag{2}$$

Indeed this is the area of the thin triangle with base $ds$ and height
$h_j$. For normalized arclength $x=s/L$ and normalized swept area $a(x)$,

$$a'(x)=\frac{Lh_j}{2A}. \tag{3}$$

Consequently $a(x)=x$ on the entire marked lap **iff all side distances
$h_j$ are equal**. The common distance is then $2A/L$. This is precisely
the condition that the polygon have an incircle centered at $O$.
Every centered regular polygon satisfies it, as do nonregular polygons
with such an incircle. In those cases constant area speed and constant
perimeter speed are literally the same normalized speed.

For a general polygon set
$\alpha=\min_j Lh_j/(2A)$ and $\beta=\max_j Lh_j/(2A)$. Applying the
derivative bounds to each of the two complementary arcs gives

$$\alpha d(x,y)\le d(a(x),a(y))\le\beta d(x,y). \tag{4}$$

Thus area clearance $1/n$ implies perimeter clearance at least
$1/(\beta n)$; the exact critical threshold need not survive. For the
centered $4\times2$ rectangle the four derivatives are
$3/4,3/2,3/4,3/2$, and $a(1/3)=1/4$. This is an exact metric loss, despite
both clocks being continuous, increasing, and degree one.

## 3. A sharp runner-count/symmetry theorem

Let $F$ be an orientation-preserving circle homeomorphism and use a
degree-one increasing lift to $\mathbb R$. For every integer $n\ge2$,

$$\begin{split}
&d(x,y)\ge1/n\quad\Longleftrightarrow\quad
  d(F(x),F(y))\ge1/n\quad\text{for all }x,y\\
&\hspace{35mm}\Longleftrightarrow\quad
  F(x+1/n)=F(x)+1/n\quad\text{for all }x.
\end{split} \tag{5}$$

For necessity, apply the preserved predicate to the $n$ points
$x+j/n$. In their preserved cyclic order every image gap is at least
$1/n$, while the gaps sum to 1. Every gap is therefore exactly $1/n$.
For $n=2$, the image pair must remain antipodal, giving the same conclusion.
For sufficiency, for fixed $x$ the safe arc is
$[x+1/n,x+1-1/n]$. Monotonicity and the displayed equivariance map this
arc exactly onto the safe arc based at $F(x)$.

Now let $F_m$ send normalized arclength on a centered regular $m$-gon to
the normalized central polar angle, where $m\ge3$. Then

$$F_m\text{ preserves the critical }1/n\text{ predicate }
\quad\Longleftrightarrow\quad n\mid m. \tag{6}$$

Rotational symmetry proves the forward equivariance when $n\mid m$.
For the converse let $\rho$ be the apothem and $r(x)$ the radius at the
boundary point. On every side,

$$F_m'(x)=\frac{L\rho}{2\pi r(x)^2}. \tag{7}$$

The derivative extends continuously across the vertices. Its minima are
exactly the $m$ equally spaced vertices, and its maxima are exactly the
$m$ side midpoints. If (5) holds, the derivative has period $1/n$, so this
shift preserves the vertex set. Hence $1/n$ is an integral multiple of
$1/m$, which says $n\mid m$.

In particular, **for $n\ge3$, $n$ runners on a regular $n$-gon have exactly the same
critical loneliness predicate in normalized perimeter, swept area, and
central angle, at every same time and for every observer**. More generally
the polygon can have any multiple of $n$ sides, with at least three sides;
for two runners a square works. The central angular
histories usually have variable speeds: this theorem is an exact
reformulation of the perimeter problem, not an independent existence proof.

Radially projecting all points to one circle of radius $R$ gives the
equivalent chord test

$$|R e^{2\pi iF_m(x)}-R e^{2\pi iF_m(y)}|
       \ge2R\sin(\pi/n). \tag{8}$$

The literal polygon chords are different. In the square with vertices
$(\pm1,\pm1)$, adjacent vertices and adjacent side midpoints both have
perimeter and angular separation $1/4$, but their squared chord lengths
are respectively 4 and 2. Their radii are the necessary sidecar. A
mismatched count also fails: on this square the points at arclengths
$0,1/3$ are perimeter-safe for three runners but their central angle is
less than $2\pi/3$.

### The focus of a genuine ellipse

For a noncircular ellipse about one focus, let $\mu$ be normalized swept
area and $F(\mu)$ the normalized true angle. If its semiaxes are $a,b$,
the same area derivative gives $F'(\mu)=ab/r(\mu)^2$. There is exactly one
closest and one farthest boundary point from that focus per lap, so
$F'$ has a unique maximum and a unique minimum. It has no period $1/n$
for any $n\ge2$. By (5), the area-to-angle map of a noncircular Kepler
ellipse does **not** preserve every critical runner predicate. The circle
is the positive boundary case. This identifies a clock obstruction; it
does not disprove loneliness for any particular family of elliptic motions.

## 4. Common time changes survive; individual ones can fail sharply

For any configuration curve $X(t)$ and any common surjective time map
$\psi$, the two sets of configurations $\{X(t)\}$ and $\{X(\psi(t))\}$
are equal. Thus existence of (1), and the supremum clearance, are
unchanged. A degree-one increasing periodic clock maps a full period onto
a full period. All runners must receive the same clock. This is a map of
the entire source configuration, not separate freedom to place each runner.

Here is a literal hostile with the original start and average speeds
retained. Keep observer 0 and the first runner $x_1(t)=t$. The constant
speed second runner $x_2(t)=2t$ has maximum common clearance $1/3$.
Let $G$ be the piecewise-linear degree-one lift with knots

$$ (0,0),\quad(1/3,1/6),\quad(2/3,5/6),\quad(1,1),$$

and replace only the second history by $\widetilde x_2(t)=G(2t)$.
Equivalently its private clock is $\psi_2(t)=G(2t)/2$. This clock is
strictly increasing, fixes 0, and satisfies $\psi_2(t+1)=\psi_2(t)+1$;
the second runner still has average speed 2. Nevertheless

$$\max_{0\le t\le1}\min(d(t,0),d(G(2t),0))=3/10<1/3, \tag{9}$$

with exactly the maximizers $3/10,7/10$. To see this directly on the first
half-period: for $0\le t\le1/6$ the minimum is $t$; for
$1/6\le t\le1/3$, $G(2t)=4t-1/2$, and comparison with the tent function
$d(t,0)$ gives a maximum $3/10$ at $t=3/10$; for
$1/3\le t\le1/2$, $G(2t)=t+1/2$, giving minimum $1/2-t\le1/6$.
The second half is its reflection. This is a counterexample to the
individual-clock transfer, not to constant-speed LRC.

## 5. An intrinsic orientation is still too small an observer

Orient a pair by whether its positive circular difference lies in
$(0,1/2)$ or $(1/2,1)$, retaining difference 0 or $1/2$ as a tie.
The configurations $(0,1/10,1/5)$ and $(0,7/20,2/5)$ have the same three
strict orientations, but observer 0 has clearance respectively below and
above $1/3$. Thus even the honest half-turn tournament does not determine
(1). Gap sizes and the observer must remain. For even regular polygons
antipodal pairs are ties; a cosmetic orientation of them is not licensed.

## 6. The surviving general LRC obligation

For the AP, $n$ points $0,t,\ldots,(n-1)t$ satisfying critical observer
clearance have **all** pairwise distances at least $1/n$: every difference
is another speed from $1,\ldots,n-1$. All cyclic gaps must then be exactly
$1/n$, giving a regular configuration and $t=a/n$ with $\gcd(a,n)=1$.
Conversely each such clock gives precisely the full regular polygon.
This proves the exact AP maximum and complete witness set.

An arbitrary speed set has no such difference-closure, and LRC only asks
for a star of separations from one observer. Requiring a complete regular
configuration is stronger. The unresolved target remains, for each marked
observer, a single literal time satisfying that observer's inequalities
for the original speeds; different observers may use different times.
The polygon equivalence carries this obligation faithfully; Kepler's law
does not supply that simultaneous time.

The [Conway tail certificate](../conway_tail_barrier_20261013c.md) gives its
own source-indexed $1/20$ result. Its recurrence, finite boundary, and the
number 1024 are not runner speeds or a substitute for the LRC threshold
$1/n$. The useful shared pattern is preserving the exact clock and
checking the whole boundary, not identifying the constants.

## 7. Exact implementation and reproduction

After the local proof was frozen, incoming
[polygon areal transport](../polygon_areal_transport_20261021/REPORT.md)
supplied the same threshold criterion with a broader rational-threshold
and radial-curve formulation. The matching result is shared, not a
separate increment in LRC coverage. That packet's LCM symmetry criterion
also shows how to preserve half-turn ties together with1/n: use a polygon
whose side count is divisible by lcm(2,n).

[check.py](check.py) uses rational polygon coordinates with metric
$dx^2+Ddy^2$, exact cosine inequalities, and literal rational
piecewise-linear histories. Its lower-envelope reader partitions at all
phase breakpoints and half-integer crossings, then checks all affine
intersection points. This proves each supplied finite maximum exactly;
it does not enumerate arbitrary speeds or arbitrary nonlinear clocks.
For a general history with a plateau, the reported `witnesses` are
maximizing certificate points, not an interval encoding of all maximizers.

The finite universe comprises all 60-by-60 pairs for regular triangles,
squares and hexagons at thresholds $1/2,1/3,1/4,1/6$; rectangle distortion;
the sharp clock hostile; six common-clock speed sets; complete AP clocks
through 16 runners; orientation and strict type hostiles. The universal
statements (2)–(8) are proved above, not inferred from this census.

```powershell
python -B -X utf8 05-knowledge/results/runner_area_clock_20261020j/check.py
python -B -O -X utf8 05-knowledge/results/runner_area_clock_20261020j/check.py
```

The immutable reader recomputes and compares [receipt.json](receipt.json).
Only explicit `--write` replaces that receipt. Normal and optimized replay
both pass **25,317** always-active checks with saved LF SHA-256
`a219f5e29a66055e87e8de098e76c8a32112499ea5886d99aaf566521a76477e`.
No imported literature
headline or unresolved LRC claim is tested or promoted by this packet.

**Independent audit:** `family_algebra` read the complete proofs and public
APIs and independently reproduced the optimized output and saved hash:
**PASS**. The root review also passed. The final prose explicitly restricts
the regular $n$-gon sentence to $n\ge3$ and retains a separately marked
observer; these clarification edits do not change the computation.

## Incoming Kepler work: closed endpoints and the rotation sidecar

The [boundary and colouring supplement](boundary_and_colouring.md) repairs
THM-4737 / HYP-9350 at three exact interfaces: one-sided endpoint rates,
conditional rational-witness existence, and the additive group action
needed to recover a runner time from a graph colouring. Its explicit
distance graph has chi_f=4 while1/kappa=5. An LCM family also proves that
finite witness denominators are unbounded at each fixed runner count with
at least two nonobservers. These are scoped proofs, not general LRC
coverage. The supplement has its own exact checker and receipt and does
not alter this packet's original25,317-check universe.
