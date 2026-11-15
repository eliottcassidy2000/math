# One reversed ordinal seam: a sharp gap obstruction and a root-switch construction

Here h(A) counts directed Hamilton paths, including h(empty)=1. For x in A, let S_x sum h(A\P) over directed paths P starting at x, including the singleton path. Its signed version d_x uses sign (-1)^(|P|-1). Define T_y and e_y dually for paths ending at y. Put E_x=(S_x-d_x)/2, the weighted sum over even rooted paths.

## Signed endpoint identity

Exactly d_x is the number of Hamilton paths of A ending at x; dually e_y is the number starting at y.

An algebraic proof works in the commutative ring of variables z_v with z_v²=0. Let B=A−x, M its adjacency matrix, Z=diag(z), u indicate outneighbors of x in B, and 1 denote the all-ones vector. The ordinary path series (empty path included) is F_B=1+1ᵀZ(I−MZ)⁻¹1. Squarefree coefficients count exactly paths on the selected set. Reversing paths shows that the same series uses Mᵀ. With G=I−MᵀZ−Z, one has ZG⁻¹=Z(I−MᵀZ)⁻¹: in the expanded walks, a diagonal step repeats the same variable and dies, including a diagonal step at the first vertex because the prefactor Z includes that vertex. Since M+Mᵀ=11ᵀ−I, one has I+MZ=G+11ᵀZ. The rank-one inverse formula, applied to 1, now gives
```text
F_B Z(I+MZ)⁻¹1 = Z(I−MᵀZ)⁻¹1.
```
All inverses here exist as finite nilpotent series or inverses of scalars with constant term one. The rooted signed sum is the full monomial coefficient of (1−uᵀZ(I+MZ)⁻¹1) F_B. By the displayed identity it is h(B) minus those Hamilton paths of B ending at an outneighbor of x. The remainder can append x, which proves the assertion, also when B is empty.

## Single-seam gap theorem

Given the ordinal join A→B, replace x→y by y→x for x in A, y in B. The new count is
```text
h(A) h(B) + S_x T_y − d_x e_y.
```
Indeed, apply the positive odd-cycle packing expansion (the Grinberg–Stanley/Rédei–Berge identity; repository THM-466). Old cycle packings give the product h(A)h(B). A crossing cycle uses the unique backward arc, and decomposes into a nonempty path rooted at x in A followed by a nonempty path ending at y in B. The parity of the two path sizes is opposite for an odd cycle; its weight 2 turns the parity filter into 1−(-1)^(|P|+|Q|). All residual packings are independently internal to the two sides. This gives the formula.

Cutting an ordinary Hamilton path of A before x shows S_x≥h(A), and symmetrically T_y≥h(B). If x has any outneighbor in A, then 2d_x≤h(A): delete x from each Hamilton path ending there and reinsert it just before the first outneighbor in the resulting ordering. The preceding vertex, if any, is an inneighbor by definition. The map is injective because deleting x recovers the original remaining order. By Rédei parity this gives d_x≤(h(A)−1)/2. A sink x instead satisfies S_x=d_x=h(A). Similarly e_y≤(h(B)−1)/2 unless y is a source, in which case T_y=e_y=h(B).

Consequently a reversal changes nothing when x is a sink and y a source, and in every other case the increment is at least half the old count. If both local endpoint exceptions fail, the increment is at least h(A)h(B)−(h(A)−1)(h(B)−1)/4, more than three quarters of the old count. In particular one seam cannot act as an additive, sublinear repair to an arbitrarily large fixed ordinal baseline.

For singleton B, the increment is exactly 2E_x and 0<2E_x implies h(A)≤4E_x−1. The constant 4 is asymptotically sharp: use the n-vertex near-transitive tournament with last→first as its only backward arc, and root x at the penultimate vertex. Its count is h=2^(n−2)+1; its even rooted paths consist of x→last, or begin x→last→first and then take an odd subset of the n−3 other sites. Every nonempty complement in this computation is transitive, so E_x=1+2^(n−4) for n≥4 and h=4E_x−3. Taking the dual copy on the right shows the three-quarters coefficient for two active endpoints is also asymptotically sharp: in these examples S_x=h+1 and d_x=(h−1)/2. A source-side root whose ordinary count is h cannot produce h+2 by this surgery unless h≤3; h=1 and h=3 both occur.

## A bounded-increment exact inventory

Decompose A=P→C→R by the strong component C containing x, with P and R optional. A rooted path cannot enter P and can use an arbitrary optional path in R after visiting C. Thus S_A=h(P)S_C Q(R), where Q(R) counts ordered two-slot covers including empty paths. The signed identity gives d_A=0 for nonempty R, because a Hamilton path cannot end in C; when R is empty d_A=h(P)d_C. It follows that
```text
E_A=h(P)E_C                        if R is empty,
2E_A=h(P)S_C Q(R)                  if R is nonempty.
```
In the latter case S_C≥h(C), Q(R)≥2h(R), and Q(R)≥2^|R|, so h(A)≤E_A and |R|≤log₂(2E_A). In the former case for E_A>0 the root component is nontrivial, and Camion's Hamilton cycle yields from x a rooted path of every prefix size; all complement counts are positive. Thus |C|≤2E_C+1. In the latter case nontrivial |C| is bounded logarithmically by h(C)≤2E_A and the strong path-count floor. Prefix singletons have no effect, while every nonsingleton prefix component contributes an odd factor ≥3, and has bounded order by the same floor. This gives a finite exhaustive inverse for fixed E_A, even allowing arbitrarily large A.

For E_A=1 the possible old counts are 1,3. For E_A=2 they are 1,5. To see the only potentially nontrivial second case, with R empty note h(P) is an odd factor of 2 and hence 1. Then C has order at most five. The triangle has E_C=1 and the unique strong order-four class has root E_C∈{2,3,4} and H=5 (replayed below). At order five, E_C=2 would require exactly one rooted 2-path x→o and one rooted 4-path; the three-vertex complement of {x,o} must be transitive, say p₀→p₁→p₂. If o beats p₀ then at least two rooted 4-paths exist, so x and o cannot reach p₀ via any of the remaining vertices. This contradicts strength.

## A nonmirror exact root switch

The obstruction above applies when *adding* a seam to the old ordinal baseline. Comparing two seams can instead yield small differences. For any nonempty tournament D, construct A on x₀→D→x→z with only the extreme arc changed to z→x₀. Let h=h(D), q=Q(D). Then h(A)=h+2q, and the outdegree-one roots z and x satisfy
```text
S_z=h(A),   S_x=h(A)+h.
```
The count h(A) comes from the seam theorem applied to x₀→(D→x) on the left and the singleton z on the right. The optional signed weight on the nonempty D→x is zero because the apex root x₀ is not a possible Hamilton endpoint there (equivalently d_{x₀}=0 in the enlarged source); and its two-slot cover is 2q. For the root equalities, S_w−h(A) counts the covers obtained by cutting before w that fail the joining test. If w has a unique outneighbor u, such a failure can occur only if the rooted path is the singleton w, and then it is exactly a Hamilton path of A−w ending at u. For z this count is zero since its sole outneighbor x₀ is a source in A−z. For x it is h(D), since the only exceptional edge in A−x leaves z and therefore cannot be used in a path ending there.

Choose an arbitrary nonempty tournament B (count b, two-slot count p), and form the ordinal A→B→y with y singleton. Reverse either z→y or x→y, but not both. For the two roots, their total path weight across B becomes S_z p or S_x p, while the signed weight is zero since the root cannot be a Hamilton endpoint of A→B. The two exact counts are
```text
N₀=(h+2q)(b+p),      N₁=(h+2q)(b+p)+hp.
```
Both tournaments are strong: A is strong as x₀→Hamilton path through D→x→z→x₀ is a spanning cycle, every A vertex can reach every latter vertex via B, and from every B vertex the remaining arc to y returns to A. Thus this gives a physical two-arm target selector: the unswitched arm is multiplicative, but the switched arm adds hp while reusing both operands. With D transitive of m≥1 vertices and B singleton, it gives the strong consecutive odd pair (3·2^(m+1)+3, 3·2^(m+1)+5) on m+5 vertices, at a target-selected power-of-two scale. With general D it yields the line N₁=5h+6q for singleton B. None of these structured arithmetic conditions places a repair at every prescribed target.

Run `python verify_seam.py` for direct directed-path DPs, including every tournament through order 5 for the endpoint identity and rooted E controls, deterministic literal seam reversals, sharpness tests, and the two-arm construction.

### Exact target inversion and a limit of degree-one roots

The switched output is symmetric in the two operand profiles after
expansion:
```text
N₁ = 2(h+q)(b+p) − hb.
```
This gives an exact target-first decoder for any chosen finite bank of
un-padded operand profiles (h,q₀),(b,p₀) with freely available ordinal
singleton padding. Factor (N+hb)/2 and test divisor pairs (u,v). They
are implemented by that bank exactly when
((u-h)/q₀) and ((v-b)/p₀) are both nonnegative integer powers of two,
meaning 1, 2, 4, … . For if the ratios are 2^s and 2^t, append s and t
ordinary singleton components to the bank operands before building
the strong two-arm construction. Conversely any such padding produces
that divisor pair. Arithmetic tests without both literal bank profiles
supply no witness. This decoder is strictly conditional on the bank
and does not make its represented targets dense.

Within strong A of order at least four, using two distinct roots of
outdegree one cannot give a richer unit-difference baseline. Orient
their mutual arc as x→z; then z has unique outneighbor u in
D=A−{x,z}. In A−z, x is a sink, so S_z=h(A). In A−x, a Hamilton path
ending at z corresponds precisely to a path of D not ending at u.
Hence S_x−S_z=h(D)−end_u(D). Strength requires u to reach every vertex
of D (a path from u that detours through x,z returns to u before
entering D again), so u is not a sink in D of size at least two.
By the endpoint injection this difference can be one only when
h(D)=1. There is then no directed triangle in D, by the positive
odd-cycle expansion, so D is transitive, with u its source by
reachability. The resulting A is the near-transitive backward-arc
construction already displayed. Thus a more flexible unit repair
using two roots of the same first operand needs at least one root
of outdegree two or more.

## Higher-degree root switches, signed repair, and a ternary growth scale

One can escape the degree-one rigidity with a module inversion. Use
the order-five template with outgoing row masks [24,25,3,4,12]
(orientation code 108, low bits for lexicographic unordered pairs,
1 meaning the smaller vertex wins). Replace vertex 4 by a tournament
D, keeping uniform arcs to it. Number the remaining vertices 0–3
as in the template. Vertices 2 and 3 give the exact identities
```text
        H(D)   P₂(D)
h(A)      11      6
S₂        13      6
S₃        12      6
```
All higher columns are zero. To compute the matrix, enumerate the
orders of all non-module vertices. For h(A) put them in one path;
for S_i put them in an ordered pair of paths, first starting at i,
second allowed empty. In every internal gap between old vertices
a,b, assign either no module visit, allowed iff a→b, or one visit,
allowed iff a→4→b. At each extreme, a visit is allowed by the
corresponding arc, except that a rooted path cannot begin with
the module; an empty second path may freely have zero or one
visit. Multiplying the gap polynomials counts compressed patterns by
total visits j. For an actual D these visits are replaced by an
ordered cover of D by j nonempty paths, giving P_j(D). Conversely
compress every actual cover, so the prescription is exact.
The code below replays the finite polynomial enumeration and
independently tests the inflated literal graph.

Write q=Q(D)=P₂(D)+2h and h=H(D). The table says
h(A)=6q−h and (S₃,S₂)=(6q,6q+h). Add a buffer singleton
then a closer, with ordinal arcs except for one backward arc
from the closer to root 3 or 2. Its two outputs are
```text
18q−h, 18q+h.
```
Thus both signed repairs are now physical from the **same**
absolute baseline and actual source profile; for a transitive
D of m vertices this gives adjacent odds 18·2^m ±1. Both output
tournaments are strong because the template is strong, substitution
by a nonempty D preserves strong connectivity, and the buffer and
closer return to the core along the chosen edge. This does not show
that a prescribed N has an appropriate literal (h,q). It moves the
problem to matching a signed difference and the same q.

There is also a ternary-scale switch, not just a different numerical
coefficient of binary inflation. On order six take orientation
code 2330 and replace vertex 5 by D. Roots 4 and 3 have the following
coefficients, computed by the identical rules:
```text
        H(D)   P₂(D)   P₃(D)
h(A)      43     54      16
S₄        74     66      16
S₃        75     66      16
```
Adding the same buffer and closer gives
191h+186P₂(D)+48P₃(D) and that count plus 2h. In particular for
D transitive, exact surjection counts give the consecutive odds
48·3^m+42·2^m−37 and 48·3^m+42·2^m−35, on m+7 vertices.
This establishes arithmetic scale diversity in a single-root
repair without giving positive-density coverage.
