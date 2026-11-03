# Session 145: when pseudoforests do and do not carry fence geometry

**Status:** PROVED elementary reductions and exact construction verification; CITED classical polygonal isoperimetry and Friedman's best-known construction. No global fence or product optimum is claimed. The reductions are new to this recovered repository lane; historical novelty is not claimed.

## 1. The two problems are different

The primary pages were read on 2026-11-03:

- [Maximum Fence Area](https://erich-friedman.github.io/packing/fence/): original objects are unit straight segments; all their ends are supported on another segment; crossings are forbidden; each bounded field has area at most one; maximize the **sum** of field areas. A supported endpoint may be in the interior of its supporting segment.
- [Maximum Product of Areas](https://erich-friedman.github.io/packing/product/): unit edges may intersect only at common endpoints; maximize the **product** of the areas of the bounded polygons. There is no area-at-most-one condition. We require at least one bounded polygon, as the page's n=3 value does; interpreting an empty product as one would change that case.

Both pages explicitly describe best-known configurations. Displayed decimals are not certificates of global optimality.

Throughout the fence results, distinct original fences meet in at most one point; collinear overlaps are excluded, as in the inherited THM-4505/4509 model.

For fences, first split every original segment at every junction. The resulting plane graph's edges need not have length one. For the product problem, the original unit edges are already the graph edges. This difference is load-bearing.

The recovered general fence frontier is THM-4509, `01-canon/theorems/THM-4509-fence-corner-lemma-and-density-bound.md`: the asymptotic density is at most 0.5224525, or 0.5168084 when every field has a convex outer boundary. Its per-fence note records exploratory 0.50116-class LP values as EMPIRICAL, not universal theorems. The conjectured density one half remains OPEN.

## 2. A sharp fence pseudoforest obstruction

**Theorem F.** If the split plane graph of an admissible fence configuration is a pseudoforest (each connected component has at most one cycle), its total field area is at most n/4. Equality holds exactly for non-nested, mutually disjoint unit squares, with no other fences.

Every graph vertex has degree at least two: a fence endpoint is supported by another fence, and a point interior to a fence has its two opposite rays. A connected pseudoforest with minimum degree at least two is a cycle: the degree sum is at least 2V, while E<=V, so E=V and every degree is two. In particular no T junction survives. Every split edge is consequently an original unit fence, although consecutive fences may be collinear.

For a component cycle of k edges, k>=3. If k=3, its enclosed disk is an equilateral triangle of area sqrt(3)/4<3/4, and its field has no larger area. If k>=4, the area cap gives at most 1<=k/4. Each cycle corresponds to one bounded field; if cycles are nested, the inner regions are subtracted from the containing field. Summing proves the bound. Equality forces k=4 in every component, field area one, and no nesting. A unit equilateral quadrilateral has area at most one, with equality only for a square. Conversely disjoint unit squares attain n/4.

**Hostile control.** Two adjacent unit squares use seven original fences and enclose area two, which exceeds 7/4. Their graph is a theta graph with three paths of lengths 3,3,1 between the endpoints of the common edge, and has cycle rank two. The theorem therefore does not extend to unrestricted fence configurations.

## 3. Exact product maximum for all endpoint-only cactus graphs

A cactus is a graph whose edges belong to at most one simple cycle. It includes every pseudoforest, and also allows cycles sharing vertices and bridges between cycles. Write

    a(k) = k cot(pi/k)/4,  k>=3.

**Theorem P.** Among product configurations with n>=3 unit edges, at least one bounded face, and a cactus underlying graph, the exact maximum is

    Q_n = max over 1<=r<=floor(n/3) of
          a(q)^(r-s) a(q+1)^s,
    where q=floor(n/r) and s=n-rq.

The same maximum holds if the graph is required to be a pseudoforest. It is attained by disjoint regular polygons with the displayed balanced cycle lengths.

**Classical input.** The polygonal isoperimetric inequality says a simple k-gon of perimeter k has area at most a(k), with equality only for the regular polygon. This is classical, not a new theorem here; see Beniamin Bogosel, [A Geometric proof for the Polygonal Isoperimetric Inequality](https://beniamin-bogosel.github.io/pdfs/DiscreteIsopShorter_Rev.pdf), author-hosted revision dated 2024-05-16, inspected 2026-11-03. Degenerate collinear vertices and nonconvex polygons do not improve the bound.

**Reduction.** Each cycle in an endpoint-only cactus has k original unit edges. Every bounded face is enclosed by one cycle, possibly with nested cycle interiors removed. The cycle edge sets are disjoint. Therefore the product is at most the product of a(k_i), with sum k_i<=n. Conversely disjoint regular polygons attain any such product. Bridges, pendant edges, and nesting can be removed for this upper bound; no claim is made that this preserves a specific given embedding. Since a(x) is strictly increasing for x>=3, any unused edge budget can be added to a cycle to improve the bound. Thus the maximizing cycle sizes sum to n.

**Balancing.** Let f(x)=log a(x), t=pi/x, so 0<t<=pi/3. Direct differentiation gives

    f'(x) = [1+2t/sin(2t)]/x > 0,
    f''(x) = -[1+4t/sin(2t)-4t^2 cos(2t)/sin(2t)^2]/x^2 < 0.

For the second sign, multiply the bracket by sin(2t)^2. If cos(2t)<=0 all terms are positive. Otherwise use sin(u)>u cos(u) for 0<u<pi/2 (its difference has derivative u sin(u)>0), with u=2t. The numerator sin(2t)^2+4t sin(2t)-4t^2 cos(2t) is then positive. Thus f is strictly concave. Replacing cycle sizes p,q with p+1,q-1 whenever q>=p+2 strictly increases the product. For fixed r, all sizes therefore differ by at most one, proving the formula.

The exact formula is the theorem; numerical tables in the checker are certified finite evaluations, not a replacement for this all-n proof.

## 4. The first displayed product improvement has an exact structural explanation

At n=12 the balanced possibilities are

    r=1: a(12)=6+3 sqrt(3),
    r=2: a(6)^2=27/4,
    r=3: a(4)^3=1,
    r=4: a(3)^4=9/256.

Hence Q_12=6+3 sqrt(3)=11.1961524... . Friedman's n=12 image and formula are a non-cactus construction with product

    (89+26 sqrt(13))/16 = 11.4215208... .

Here is a standalone exact realization, reconstructed from the inspected image, not a new record. Put

    R=sqrt(13), D=sqrt((R-1)/6),
    A=(R+1)/4, B=3D/2, H=D(R+7)/4.

The upper field has boundary vertices

    (-1,0), (-A,B), (-1/2,H), (1/2,H), (A,B), (1,0),

and the lower field is its reflection across the x-axis. Include the ten outer edges, and split the common length-two base into the two unit edges from (-1,0) to (0,0) to (1,0). All twelve original edges have length one; intersections occur only at common endpoints. The upper polygon has five unit sides and one length-two side. The latter is not one of the original unit edges: it comprises two collinear edges.

Direct squared-distance checks use R^2=13 and D^2=(R-1)/6. The two polygons are strictly convex and lie on opposite sides of their common base. A shoelace calculation gives each area

    H(A+1/2)+B/2 = D(23+5R)/8.

Its square is (89+26R)/16. The connected graph has E=12, V=11, cycle rank two. It is the theta graph with path lengths 5,5,2, not a cactus.

For an elementary strict comparison, R>18/5 and sqrt(3)<7/4 give

    26R-7-48 sqrt(3) > 26(18/5)-7-48(7/4)=13/5>0.

Thus the displayed product strictly exceeds the maximum of every unit-edge cactus configuration, including every pseudoforest. This proves a genuine gain from sharing a path between cycles. It does not prove the displayed construction globally optimal. Its fields have area greater than one, so it is not an admissible fence configuration.

## 5. Why the cactus extension is unsafe for fences

Theorem P extends to cactus graphs because endpoint-only intersections preserve unit cycle edges. Theorem F has not been extended to general fence cactus graphs.

An exact counterexample to the *proposed unit-cycle proof*, not to the unproved bound, is this valid five-fence arrangement. One original fence is [(0,0),(1,0)]. Two further unit fences form the upper isosceles triangle on base [(0,0),(1/3,0)] with apex (1/6,sqrt(35)/6). Two others form the triangle on base [(2/3,0),(1,0)] with apex (5/6,sqrt(35)/6). Every original fence endpoint is supported, and each field has area sqrt(35)/36<1.

The split graph is a cactus with two triangles joined by a bridge, but each triangular cycle has edge lengths 1,1,1/3, not 1,1,1. It has seven split edges, not five. The two T junctions are exactly the lost coordinate. Applying a unit-cycle partition to this graph would use the wrong edge budget.

**Safe corollary.** The fence n/4 bound does extend to cactus configurations with endpoint-only intersections (no through-fences): each cycle edge is then an original unit edge, the area cap bounds k>=4 cycles by k/4, and k=3 cycles contribute sqrt(3)/4. Bridges only consume unused budget. General T-junction cactus fences remain outside this corollary.

## 6. Connection contract and stopping boundary

Source: an abstract pseudoforest / cactus. Target: a unit-segment plane configuration. The map retains graph incidence and cycle rank, but does not supply coordinates, segment lengths, an embedding without crossings, supported endpoints, or field areas. For fences it must also retain which split graph pieces belong to each original unit fence and their order/lengths. A valid enumeration of abstract pseudoforests cannot by itself count fence configurations or certify the product objective.

The useful transport is the opposite direction: after checking actual geometry, the graph class gives the exact restrictions above. No relationship to 1729, solfeggio frequencies, knots, or unrestricted Ramsey numbers is inferred from a shared integer or from graph terminology alone.

Cheapest next target: optimize product on prescribed theta graphs with path lengths p,q,r, retaining the common path's shape, rather than abstract cycle counts alone. The n=12 witness is the (5,5,2) seed. A fence continuation needs the per-original-fence bin-length constraints already identified by THM-4509; the abstract cactus extension is not established here.
