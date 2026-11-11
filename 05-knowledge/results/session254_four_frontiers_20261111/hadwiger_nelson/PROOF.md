# Explicit algebraic soft amplifiers at fixed planar ports

**Session254, 2026-11-11. PROVED:** the all-size induced unit-distance construction and its exact conditional colouring counts. **FINITE-EXACT:** the attached rational-parameter and independent count controls. **OPEN:** any improvement to the unrestricted chromatic number of the plane. No novelty/priority claim relative to the literature is made.

## 1. Inheritance and the changed operation

Source pins are repository-relative, at inherited commit `7b7c5965cc1e`:

* `05-knowledge/results/boundary_palette_compression_20261011.md`, Sections 1–3: complete boundary equality relations, geometric sidecars and circle saturation.
* `05-knowledge/results/boundary_context_weights_20261011d.md`, Section 2: conditional integer weights and their product under private-interior gluing.
* `05-knowledge/results/existence_count_amplifiers_20261011e.md`, Sections 2–3: serial/parallel algebra, positive-support obstruction and the two-common-neighbour geometric cap.
* `05-knowledge/results/unit_equality_quantized_guards_20261011e.md`, Sections 2–3: absolute integer quanta, fixed-universe deletion versus new-vertex amplification.
* `05-knowledge/results/unit_triangle_product_boundary_20261013c.md`, Sections 2–3: exact incidence and palette-specific preservation as separate obligations.

The closest mechanism is private-interior multiplication. The hostile is `K_{2,r}`: excellent equality bias abstractly, but only two common neighbours can be placed between distinct unit-separated planar ports. The corrected near miss is extrapolating this cap to longer paths, or treating a probability tending to zero as an exact zero. The underused sidecar is **continuous geometric freedom of an entire path with fixed ports**, rather than freedom of one common neighbour.

Board: fixed metric ports / path-length parity / exact induced incidence / integer transfer channels / positive support. Anchor: realizable palette operations. Niche: an explicit quadratic arch family. Wildcard: positive-integer product realizability. We replace length-two paths by length-four paths and keep all cross distances.

## 2. All-size induced unit-distance theorem

**Theorem A.** For every integer r>=1 there is an explicit finite algebraic planar point set whose full unit-distance graph consists precisely of r internally vertex-disjoint length-four paths between two ports at distance **two**. It has 3r+2 vertices, 4r edges, and chromatic number exactly two. Both equal and unequal port assignments are genuine geometric boundary states when k>=3.

**Construction.** Put u=(0,0), v=(2,0). For j=1,...,r define

    t_j = 1/2 - 1/[20(j+1)],
    c_j = (1-t_j^2)/(1+t_j^2),
    s_j = 2t_j/(1+t_j^2),
    h_j = s_j + sqrt(1-(1-c_j)^2),
    x_j = (c_j,s_j), y_j = (1,h_j), z_j = (2-c_j,s_j).

The jth path is u,x_j,y_j,z_j,v. Every coordinate is real algebraic, with x_j,z_j rational and y_j quadratic over Q (possibly degree one). No fixed uniform field-degree bound for the whole r-path configuration is claimed.

**Proof of complete incidence.** All c_j are distinct and satisfy 3/5<c_j<2/3. Their positive sine coordinates obey c_j^2+s_j^2=1, so ux_j and z_jv have length one. The definition of h_j gives |x_j-y_j|=|z_j-y_j|=1.

The low vertices x_j,z_j lie in the open rectangle

    (3/5,7/5) x (sqrt(5)/3,4/5).

Its width is 4/5 and its height is less than 1/10 (since sqrt(5)/3>7/10), so its diameter is less than sqrt(65)/10<1. There are no edges among any low vertices. Distinct x_j have x coordinate <2/3, distinct z_j have x coordinate >4/3, and all their parameters are different, so these vertices never coincide.

For all j, 1<h_j<2: the lower bound follows from s_j>7/10 and sqrt(1-(1-c_j)^2)>sqrt(21)/5>9/10; the upper bound follows from both summands being less than one. As a function of c in (3/5,2/3),

    h(c)=sqrt(1-c^2)+sqrt(1-(1-c)^2)

is strictly decreasing. Indeed its derivative is

    -c/sqrt(1-c^2) + (1-c)/sqrt(1-(1-c)^2) < 0,

because a/sqrt(1-a^2) is strictly increasing on (0,1) and c>1-c. Thus the y_j are distinct, and their pairwise distances are strictly less than one. Each y_j has distance greater than one from both u and v.

The two intersections of the unit circles centered at u and y_j are exactly x_j and y_j-x_j: the latter has x coordinate 1-c_j<2/5, and both are on both circles. They are distinct. Every x_i has x coordinate >3/5, so |x_i-y_j|=1 forces x_i=x_j and hence i=j. Similarly, the two intersections of the unit circles centered at v and y_j are z_j and v+y_j-z_j. The latter has x coordinate 1+c_j>8/5, whereas every z_i has x coordinate <7/5. Hence |z_i-y_j|=1 iff i=j.

Finally,

    |u-z_j|^2 = |v-x_j|^2 = 5-4c_j > 7/3 > 1.

No unwanted edge to either port exists. The ports themselves have distance two, so there is no port edge. This exhausts all unordered pairs and proves the exact induced graph. A two-colouring gives u,v,y_j colour 0 and x_j,z_j colour 1; an edge exists, so the chromatic number is exactly two. QED.

**Unit-separated auxiliary variant.** Setting v=(1,0), z_j=(1-c_j,s_j), y_j=(1/2,s_j+sqrt(1-(c_j-1/2)^2)), and choosing distinct c_j in (1/2,2/3) instead realizes the r paths plus the unavoidable port edge. Its path-only colouring channels require deliberately omitting that actual unit edge, and equal-port assignments are not proper colourings of the full point-set graph or plane. This variant is an edge-deleted subgraph model only. The distance-two theorem above does not have that limitation.

**Boundary cases.** r=0 gives two isolated ports. The strict parameter interval is load-bearing. Unit-distance graphs allow edge crossings; no noncrossing drawing requirement is imposed.

## 3. Exact palette amplification that is genuinely realizable

Let F_r be the **full induced unit-distance graph** in Theorem A, with its distance-two ports. No unit edge is omitted.

For k>=2 let E_r(k) count proper k-colourings of F_r with fixed port colours (0,0), and D_r(k) count those with fixed port colours (0,1).

**Theorem B.** For every k>=2 and r>=1,

    E_r(k) = [((k-1)^4+(k-1))/k]^r
           = [(k-1)(k^2-3k+3)]^r,
    D_r(k) = [((k-1)^4-1)/k]^r
           = [(k-2)(k^2-2k+2)]^r.

**Proof.** The single-edge palette matrix is M=J-I. On the constant vector its eigenvalue is k-1; on its orthogonal complement the eigenvalue is -1. Thus M^4 has diagonal e=((k-1)^4+k-1)/k and off-diagonal d=((k-1)^4-1)/k. These entries count the colourings of a length-four path conditional on its endpoint colours, by summing over its three internal colours. With port colours fixed, the r private interiors are independent, giving e^r,d^r. These are the inherited transfer formulas applied to the newly authenticated geometric family. QED.

For k>=3, e=d+1 and d>0. Consequently

    D_r/E_r=(d/(d+1))^r -> 0,
    Pr[different ports] = (k-1)d^r / [(d+1)^r+(k-1)d^r] -> 0,

but D_r>0 for every finite r. At palette six the channels are exactly

    (E_r,D_r)=(105^r,104^r).

Thus **unbounded equality bias really is possible at fixed distance-two planar ports**, despite the two-common-neighbour cap on the inherited short-path bank. The correct obstruction is not geometric impossibility of all amplification: it is the unchanged positive support and growing absolute count. The full induced graph has kE_r(k)+k(k-1)D_r(k) colourings, and remains two-chromatic for all r. Its two-colouring forces equal ports, but at k>=3 every unequal-port count remains positive. At k=2, d=0 and e=1; this is exactly the even-path parity exception already predicted by the inherited algebra.

The improvement over the inherited construction is the all-r geometric realization with full cross-edge exclusion. The positivity argument, matrix eigenvalues, and conditional multiplication are inherited elementary mechanisms, not claimed as new general algebra.

## 4. Exact bridge and what it does not transfer

**To positive-integer constructive arithmetic (usable alongside Chan–Pak):** the map sends r labelled independent path interiors, after conditioning on a fixed port assignment, to their Cartesian product of colouring sets. It preserves cardinality and supplies a literal bijection, restriction/gluing. The attainable channel semigroup is {(e^r,d^r):r>=1}; at palette six it is {(105^r,104^r)}. Geometry, integrality, and multiplication are all paid here. It does not realize arbitrary positive integers, addition of channels, or arbitrary graphs with a specified statistic. There is no supplied map from these colouring sets to the Chan–Pak native counting statistic. Export the Cartesian-product mechanism and this explicit realizability witness, not a Chan–Pak conclusion.

**To LRC:** normalized rarity is insufficient for eliminating a positive integer fibre; the precise shared lemma is N>=1 implies N/Z>=1/Z. In our construction Z=E_r+(k-1)D_r grows. An LRC use must provide its own finite source universe and absolute denominator bound; no map from these vertices to runners or owner intervals is supplied. This is a proof-obligation transfer, not a problem reduction.

**To planar Jacobian:** no geometric map is asserted. Euclidean unit-distance coordinates and polynomial maps with constant Jacobian are different objects. A finite fibre product has a count-multiplication mechanism only after its actual compatible fibres are identified. Nothing here certifies Jacobian properness or injectivity.

The native next target is now sharply separated: add **controlled interactions between different path interiors** while preserving one equal-port colouring and deleting every unequal-port colouring. The current exact arch family deliberately forbids those interactions, hence cannot meet that target. Any useful modification must alter this incidence theorem, not merely increase r. A second concrete search target is a flexible two-port gadget with e>d and smaller d than the four-path value, to amplify faster while still acknowledging the support barrier.

## 5. Focused checks and reproduction

Run from repository root:

    python3 05-knowledge/results/session254_four_frontiers_20261111/hadwiger_nelson/check.py
    python3 -O 05-knowledge/results/session254_four_frontiers_20261111/hadwiger_nelson/check.py

Both outputs agree with `output.json`. The finite universe is 100 explicit rational arch parameters; palettes 2 through 8 for a single path, with independent literal internal-colour enumeration versus matrix multiplication; and two parallel paths at palettes 2,3,4 by literal enumeration. The rational checks authenticate the strict hypotheses used in the analytic all-pair incidence proof. An additional 10,000 ordered cross-incidence checks use the exact rational norm A^2-B^2 rad_j of squared-distance minus one, independently excluding x_i-y_j and (by reflection) z_i-y_j cross edges for i!=j; low-low distances are also checked rationally. High-high and port-high exclusions use the proved strict height bounds, rather than a claimed full symbolic distance census. The all-r realization and every cross-edge exclusion are proved in Section 2. The k=2 boundary and r=2 positive-support examples are hostile controls against overextending the k>=3 conclusion. No inherited historical audit campaign was conducted.
