# Complete owner-line dictionaries are count-safe at every height

**Status: PROVED, with exact finite head and independent mathematical
audit.** This closes a class with arbitrarily many primitive
directions by retaining its affine owner lines. No universal count bound,
network inequality, physical entry, or LRC(14) proof is claimed.

Let `w=(a,b,c)` be primitive, sorted, distinct, positive, odd, and nonzero
modulo three. Set

```text
L={C in Z^3 : C dot w=0},
K={C in L tensor R : 14|C_i|<3(w_j+w_k), i=1,2,3},
Omega={C in K intersect L : C_i!=0 mod3 for all i},
tau(C)=aC_1=bC_2=cC_3 mod3 in {1,2}.
```

Suppose `Omega` spans the plane and the points of each owner color are
collinear. Then

```text
N=|Omega| < 2c/11.                                     (1)
```

Thus the count gate in
[THM-4422 / projection-deficit-and-Beatty-row-reduction](../../01-canon/theorems/THM-4422-lrc14-projection-deficit-and-beatty-row-reduction.md)
already proves all three network projections strictly below `6/77` on this class.
The hypothesis concerns two parallel affine lines, one per owner color,
and permits arbitrarily many directions through the origin.

The [independent audit](overnight2_20260906_lrc_owner_lines_audit.md)
accepted the lattice saturation proof, reconstructed all seven chord bounds
using a different affine section, and rebuilt all 678 finite-head triples
by solving a different coordinate. Its 3,880 exact gates passed. There is
no equality case in (1).

## Inheritance and the changed object

The closest proved mechanism is
[the exact midpoint payment and affine residual](overnight_20260906_lrc_midpoint_deficit.md).
It exposed the complete `(1,m,2m+3)` family, whose entire positive deficit
is affine on each owner line. Counting or optimizing second differences
cannot see this baseline. The canonical excluded `A2` dictionary
`(19,23,29)` has a noncollinear triangle in each owner color; it is already
closed by the [full-cap theorem](overnight_20260906_lrc_cap_carriers.md).
The incoming `(41,47,49)` hostile has an affine AP with coefficients
`(1,1,2)`, not an additive `(1,1,1)` circuit.

The corrected near miss is treating a count of origin-rays as lattice
rank, or treating all three-ray configurations as A2. The recovered
operation is to inspect the **index of a transverse lattice basis** after
convex saturation. The least-used sidecar is the primitive invisible
direction between same-owner points, rather than a live origin-ray.

| Board object | New operation | Preserved information / boundary |
|---|---|---|
| Full convex carrier dictionary | Insert a lattice point at transverse height one | Completeness is essential; selected subsets do not suffice |
| Owner color | Keep the index-three admitted lattice and its kernel | Primitive invisible steps have size three in the admitted lattice |
| Affine deficit | Replace its zero curvature by the whole supporting line | A line's offset now forces a determinant exactly equal to `w` |
| Raw projection target | Pay by the sufficient count gate | No fixed average or physical-mass substitute is used |
| Other live lanes | Compare saturation and cost before quotienting | This is a methodological comparison, not an implication to Smith or moments |

## 1. The two owner lines force a saturated transverse basis

Define the admitted rank-two lattice

```text
Gamma={C in L : aC_1=bC_2=cC_3 mod3}.
```

Reduction modulo three maps `L` onto its two-dimensional kernel in
`F_3^3`, since `w` is primitive. The equal-coordinate owner condition is
a one-dimensional subspace of that kernel. Therefore

```text
[L:Gamma]=3,
phi:Gamma -> F_3, phi(C)=aC_1 mod3,
ker(phi)=3L,       [Gamma:3L]=3.                      (2)
```

In particular, `Omega=(K intersect Gamma) minus 3L`. These lattice facts
are also proved in the cap companion; no reserved theorem is a dependency.

Let `A` be one owner color. Central symmetry gives the other as `-A`.
Write `n=|A|=N/2`. Plane span forces `n>=2` and ensures that the affine
line containing `A` differs from its negative and does not pass through
zero. Let `r` be a primitive vector of `L` parallel to that line.

**First claim: `r` is not in `Gamma`.** If it were, there are two cases.
If `phi(r)=0`, then `r in 3L`, contradicting its primitivity in `L`.
If `phi(r)!=0`, two points of `A` differ by a multiple of `3r`; the
segment between them contains an intermediate point of the other
nonzero owner color. Convexity makes that point live on the same affine
line, whereas all other-color points lie on the distinct negative line.
This is impossible. Thus `r notin Gamma`.

Since `Gamma` has index three in `L`, its primitive vector parallel to
`r` is `v=3r`. It belongs to `3L`. The full set `A` is a consecutive
interval of the affine lattice `u+Z v`: every intermediate point has the
same owner and is in the convex open set. In particular, some consecutive
pair is `u,u+v`.

Extend `v` to a basis `(v,z)` of `Gamma`, and call its second integer
coordinate the transverse height. The functional `phi` vanishes on `v`
and is nonzero on `z`, so every height not divisible by three is owner-live.
The line of `A` has a nonzero integer height `k`; its negative has height
`-k`. Negating the choice if necessary, assume `k>0`.

If `k>=2`, the parallelogram with vertices

```text
u, u+v, -u, -u-v
```

is contained in `K`. At transverse height one its horizontal section
has length exactly one in units of `v`. A closed interval of length one
contains an integer, so this section contains a point of `Gamma`.
Its height one is not divisible by three: it is owner-live, but belongs
to neither of the two permitted lines at heights `+-k`. This contradiction
proves `k=1`.

Consequently `(u,v)` is a basis of `Gamma`. Since an integral basis of
`L=w^perp intersect Z^3` has cross product `+-w`, equation (2) gives

```text
u cross (3r)=+-3w,       u cross r=+-w.               (3)
```

This is the key new conclusion. The affine line cannot float at an
arbitrary lattice distance from the origin. Convex completeness and owner
color force the minimum transverse distance. The proof uses the actual
full dictionary, not a selected collinear subset.

It also follows that **every carrier is primitive**, and that the number
of unoriented origin-directions is exactly `N/2`. Indeed every point `C`
on one owner line has `C cross r=+-w`; a common divisor of its coordinates
would divide the coordinates of primitive `w`. Two points on that line
cannot be scalar multiples, since their cross products with `r` are the
same nonzero vector. The negative line supplies just the opposite rays.
Thus `N=4` is exactly the two-ray intersection, `N=6` the three-ray
intersection, and larger owner-line dictionaries go beyond either bound
on the number of rays.

## 2. The general primitive-direction bound

Let `M=max_i |r_i|`. Because `A` has `n` consecutive points with step `3r`,
the strict roof in a coordinate attaining `M` gives

```text
3(n-1)M < 3(w_j+w_k)/7 < 6c/7,
N=2n < 4c/(7M)+2.                                   (4)
```

For `M>=4`, this implies `N<c/7+2`.

It remains to classify `M<=3`. An invisible relation `r notin Gamma`
must have some coordinate zero modulo three: if all coordinates were
nonzero modulo three, their three nonzero weighted residues would sum
to zero, hence all be equal, putting `r` in `Gamma`.

If a coordinate of `r` is zero as an integer, the other two coefficients
have the same parity because the speeds are odd. Primitivity, distinct
speeds, and ternary units exclude all possibilities of magnitude at most
three: `(1,1)` forces equal speeds, `(2,2)` is not primitive, and `(1,3)`
forces a speed divisible by three. Two zero coordinates cannot annihilate
positive speeds. Thus every coordinate is nonzero and one has magnitude
three. The same parity condition forces the other magnitudes to be one
and two or three and two; the latter is impossible modulo three. Hence

```text
M<=3  ==>  {|r_1|,|r_2|,|r_3|}={1,2,3}.              (5)
```

## 3. Seven elementary chord certificates close the small directions

Negate `r` so `r_1>0`. Equation (3) allows choosing either owner line so
`u cross r=w`. Normalize by `c`, writing `b/c=beta`, `a/c=alpha`.
The relation `r dot w=0` gives

```text
alpha=-(r_2 beta+r_3)/r_1,
u(t)=(-beta/r_3,alpha/r_3,0)+t*r.                     (6)
```

This parameterizes the normalized line with cross product `(alpha,beta,1)`.
Each strict roof is an open interval `L_i(beta)<t<U_i(beta)`.
Write `ell` for the length of their common interval; if it is empty there
are no carriers. For any `i,j`, `ell<=U_i-L_j`. The table below lists all
signed vectors in (5) compatible with `0<alpha<beta<1`, enlarging to the
closed ratio endpoints solely for an upper bound.

| r | Allowed closed beta interval | Chosen `U_i-L_j` | Resulting bound on ell |
|---|---|---|---|
| `(1,-3,2)` | `[2/3,1]` | `U_3-L_1=beta/7` | `1/7` |
| `(1,3,-2)` | `[1/2,2/3]` | `U_1-L_2=2/21` | `2/21` |
| `(2,-3,1)` | `[1/3,1]` | `U_3-L_1=beta/7` | `1/7` |
| `(2,3,-1)` | `[1/5,1/3]` | `U_1-L_2=1/21` | `1/21` |
| `(3,-2,1)` | `[1/2,1]` | `U_3-L_1=2*beta/21` | `2/21` |
| `(3,1,-2)` | `[1/2,1]` | `U_1-L_2=2/21` | `2/21` |
| `(3,2,-1)` | `[1/5,1/2]` | `U_1-L_2=1/21` | `1/21` |

Every entry follows by substituting (6) into
`14|u_i(t)|<3(w_j+w_k)` with `c=1` and reversing an interval when
`r_i<0`. The verifier independently generates all signed permutations,
their admissible ratio intervals, and the affine roof differences; it
does not infer these bounds from sampled speeds.

The normalized parameter spacing between same-owner carriers is `3/c`.
Thus `3(n-1)/c<ell<=1/7`, yielding

```text
N<2c/21+2 <= c/7+2.                                 (7)
```

Combining (4)--(7) gives the uniform `N<c/7+2`. For every `c>=53`,

```text
c/7+2 < 2c/11.                                      (8)
```

This closes the all-height tail using line geometry, not a height fit.

## 4. Complete small head and controls

The independent-coordinate verifier enumerates all **678** primitive
sorted distinct positive odd ternary-unit speed triples with `c<=51`.
There are **139** plane-spanning dictionaries whose owner colors are
collinear: 137 have four points and two have six. Every one satisfies
`N<2c/11`; the largest `N/c` in this finite head is `4/25`, first at
`(17,23,25)`. This is not a claim of the sharp all-height ratio.

Carriers are reconstructed by scanning the first and third coordinates
inside their strict integer roofs and solving the middle coordinate
exactly. There is no direction, cap, or carrier-count prefilter before the
full set is built. The owner-line test then checks plane span and affine
collinearity of the complete owner class. Every selected row also checks
consecutive spacing, invisibility, (3), and the strict tail estimate.

The named controls are `(17,23,25)`, `(41,47,49)`, the flat ten- and
thirty-carrier examples `(1,137,277)`, `(1,499,1001)`, and the excluded
owner triangles `(19,23,29)`, `(23,29,37)`. The latter two prevent treating
all small multi-ray configurations as owner lines. At `(41,47,49)` the
primitive invisible direction is `(1,-4,3)`, outside the special small
direction table, while the flat family uses `(3,2,-1)`.

## 5. Consequence and continuation boundary

The map from a complete two-owner-line dictionary to `(u,r)` preserves
every point through its exact interval of multiples of `3r`. It turns
convex completeness into a determinant equality and then into a strict
length bound. Discarding fullness loses the forced height-one point;
discarding the owner color loses the correct primitive step and index.

Combining (1) with the previously audited cap, one/two-ray, and payment
closures, any remaining possible network failure has at least three
origin-directions, contains an affine collinear triple, has **noncollinear
points in each owner color**, and has positive unpaid count surplus.
The general case requires a genuinely two-dimensional owner polygon.
The next decisive object is its lattice width and boundary count, with the
canonical owner triangle retained as a hostile to an unqualified count
bound. LRC(14) and physical entry remain open.

For that next object there is an elementary exact count identity. Let
`P=conv(A)` be two-dimensional, let its area be normalized so a fundamental
parallelogram of `3L` has area one, and let `B(P)` count **all boundary
points in the affine owner coset**. Then

```text
N=2 Area_(3L)(P)+B(P)+2.                             (9)
```

Indeed `P` is contained in the open convex region `K`, so fullness gives
`A=P intersect (u+3L)`. Translate a vertex to zero and triangulate using
every lattice point. An empty triangle is primitive: otherwise a nonzero
lattice representative `x=alpha*v+beta*z` in its half-open fundamental
parallelogram lies in the triangle when `alpha+beta<=1`, or `v+z-x` does
when `alpha+beta>1`. Both contradict emptiness. Every small triangle
therefore has normalized area `1/2`. With `n=|A|`, Euler and edge counting
give `T=2n-B(P)-2`, proving (9). In the first two ambient coordinates,
the covolume of `3L` is `9c`, since `w` is primitive.

This is an elementary polygon-count identity, not a novelty claim or a
bound on the target. An independent review accepted its coset and boundary
normalization. The canonical A2 owner triangle has area `1/2`, boundary
count three and `N=6`. Formula (9) locates the missing boundary contribution
that an area-only continuation would lose.

Reproduce with

```text
python3 -B 04-computation/overnight2_20260906_lrc_owner_lines.py
python3 -B -O 04-computation/overnight2_20260906_lrc_owner_lines.py
```

Companions: [source](../../04-computation/overnight2_20260906_lrc_owner_lines.py)
and [frozen output](overnight2_20260906_lrc_owner_lines.out).
The 1,127 explicit gates are finite evidence for the small head and
symbolic seven-case identities. The all-height saturation and chord proofs
above are separate from that finite enumeration. No external priority or
Lean formalization claim is made.
