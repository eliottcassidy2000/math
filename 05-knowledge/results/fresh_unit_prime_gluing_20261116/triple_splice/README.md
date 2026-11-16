# An exact tripling splice with a fresh +2

**Result.** On an unbounded family of rooted tournaments there are two physical extensions preserving the source induced and keeping the same root. They implement
$$
(H,E,O,g)\longmapsto 3(H,E,O,g)+(d,0,0,0),\qquad d=0,2.
$$
They differ by one arc. The two children have identical rooted continuation coordinates $(E,O,g)$, not merely the same total response or the same residue. This is one exact digit choice, not a reset or an inductive ternary generator.

## Construction and certificate

Use the outgoing bit masks below, with vertex 0 the low bit. The root is vertex 4. In both children the source occupies vertices 0 through 5; the alternatives reverse just arc 0–6. Both children are strong.

| Tournament | Outgoing bit masks | $(H,E,O,g)$ |
|---|---|---|
| Source $C$ | 30, 60, 8, 16, 36, 13 | (17, 9, 12, 3) |
| Child $A$ | 158, 252, 8, 16, 36, 13, 189, 60 | (51, 27, 36, 9) |
| Child $B$ | 222, 252, 8, 16, 36, 13, 188, 60 | (53, 27, 36, 9) |

The [bridge protocol](../BRIDGE_RESET.md) defines $E,O$ as weighted even/odd root-start path sums, with each path weighted by the Hamilton count of its complement. Equivalently they count two ordered paths covering the vertices, the first starting at the root and the second possibly empty. The coordinate $g$ is the Hamilton count after deleting the root. These are the continuation coordinates for retaining the left root. Allowing an arbitrary change of root can read the scalar difference too, so no larger continuation claim is made.

Inflate vertices 2 and 3 to transitive blocks of sizes $s,t\ge1$ in all three tournaments. Put $X=2^s, Y=2^t$. The source profile is
$$
H=3(2X-1)(Y-1)+8,\quad E=3(2X-1)Y/2,\quad O=E+3,\quad g=3.
$$
For example, at $s=t=2$ the source has 71 paths, and the two children have 213 and 215 and common continuation (126,135,9). Fixing $t=1$ already gives an unbounded family.

**Exact finite verification.** Compress each maximal run of an inflated module in either path to a module symbol. Every nonterminal run of module 3 must be followed by vertex 4, so across two paths there are at most three such runs. Every nonterminal run of module 2 is followed by module 3, giving a cap of five. The verifier exhausts valid words within these caps, recording ordinary-vertex masks and both run counts. A module in one run has one filling; in two runs it has $X-2$ or $Y-2$ fillings by ordered nonempty assignment to the two runs. A module split between the first and second path contributes signed parity sum $-1-(-1)^s$ or its size-$t$ analogue. A module wholly in the first path contributes its size parity sign times the unsigned count; one wholly in the second contributes just the unsigned count.

For single paths the run pairs $(1,1),(2,1),(1,2),(2,2)$ have coefficients $(17,6,9,6)$ on $C$, $(51,18,27,18)$ on $A$, and $(53,18,27,18)$ on $B$. After deleting the root the only run pair is $(1,1)$, with respective coefficients 3, 9, 9. The two-path verifier prints all unsigned and signed coefficients, including how the two signed groups cancel the change at the fresh arc. Only terms with one or two total runs per module survive. For fixed parities, each profile polynomial is therefore multi-affine in $X,Y$. Matching at a two-by-two formal grid for each parity pair establishes the identities for every size. Independent subset DP checks six actual inflations and verifies the literal witnesses.

Run **python3 check_splice.py**. The calculation requires no external packages.

## Arbitrary modules: exact digit cost

The same coefficient proof works with any two nonempty substituted tournaments $M,N$. Write $u=H(M)$, $v=H(N)$ and $R=P_M(2)$, $T=P_N(2)$ for the two optional ordered path-slot counts. Then the source has
$$
(H,E,O,g)=\bigl(3(2R-u)(T-v)+8uv,\;3(2R-u)T/2,\;3(2R-u)T/2+3uv,\;3uv\bigr).
$$
The children have three times these coordinates, with respective scalar additions 0 and $2uv$.

Indeed for one run in $M$ the filling count is $u$; for two it is $R-2u$. When the two runs are split between the paths, the even-minus-odd filling count is $-(1+(-1)^{|M|})u$. To see this, include the two empty-side terms: the alternating sum over vertex subsets of the products of the two Hamilton counts is zero. In the [positive odd-cycle packing formula](https://arxiv.org/abs/2307.05569), assign each odd cycle or singleton of a full packing to one side or the other. Flipping its side flips the sign and not its weight. Remove the two empty-side terms to get the displayed value. Thus the filling table is precisely the earlier one with $X=R/u$ and all terms multiplied by $u$, and likewise for $N$. The independent DP also checks examples with nontransitive substituted modules.

This gives a limitation and an exact relation. The fresh difference is +2 precisely for $uv=1$, which forces both modules transitive. The difference cannot be +4 in this family, since $u,v$ are odd. Further, no +2 child can already be a source in this same generalized family with its rooted data: its $g=9$ would require $u'v'=3$, forcing the next source's count to be divisible by 3, while the +2 child count is 2 modulo 3. This is a failure of return even with arbitrary substitutions into the two marked sites. It does not rule out a different next construction.

## Why one reverse arc was not enough

For a word of the [single-reverse-edge bridges](../BRIDGE_RESET.md) whose right counts multiply to three, starting at $H\ge3$, the final count is either $3H$ or at least $3H+6$, regardless of the two allowed root policies or the choice of word per input. To prove it, let a right module have count $h$ and entry $y$. Its scalar surplus over $hH$ is $2(O_e E+E_e O)$. If $y$ has an inneighbor, both $O_e$ and $E_e$ are at least 1, from the singleton and length-two paths, and the surplus is at least $2(E+O)\ge2H$. Otherwise $y$ is a source, $E_e=0$ and $O_e=h$, giving a multiple of $2h$. There is one letter of count 3, the others have count 1. Any positive surplus before the 3 is multiplied by at least 3 and costs at least 6 at the end, and a positive surplus at the 3 likewise costs at least 6.

If the 3 letter has zero surplus, then $y$ was a source and the old root was a sink. Its deletion count is the old $H$. After retaining the old root, an even rooted path can take one forward arc to a right vertex other than $y$; its complement has at least the old $H$ paths via the forward ordinal paths. After switching to the right port use $y$ followed by the old root, with the same bound. Therefore immediately after the 3 letter the even rooted coordinate is at least the old $H$. If there is a next count-1 letter, its surplus is at least twice that old $H$ in either entry case. This rules out both +2 and +4 in the bridge-word architecture. The new splice adds two vertices with several coordinated contacts, even though its two settings differ by just one arc.

## Search reduction for the missing digit

Suppose an induced strong $q$-path source with $q\ge11$ is contained in a $3q+2$ or $3q+4$ host. The exterior is necessarily transitive, and at least one exterior vertex has mixed contact to the core.

If all contacts are uniform, compress the core to one quotient vertex. If it lies on no quotient triangle, it is a singleton strong component of the quotient, and the host count is a multiple of $q$. Otherwise choose a quotient directed triangle containing that vertex. On these two exterior vertices together with the core, the three quotient paths with a contiguous core give $3q$ paths. Splitting any old Hamilton path at any internal cut and taking the two outside vertices consecutively between the pieces gives another $(|C|-1)q$. Thus the count is at least $5q$, also impossible.

Now for every exterior vertex $x$ put $d_x=(H(C+x)-q)/2$. The positive packing formula gives the lower bound $qH(\mathrm{exterior})+2\sum_x d_x$: independent old-core and exterior packs are disjoint from those using a crossing cycle with exactly one exterior vertex, and the latter single-contact packs are disjoint by that vertex. The [triangle budget](https://github.com/eliottcassidy2000/math/blob/c3cf550db9ba849e8d1dfb9cf2c6eefc61b550e5/05-knowledge/results/triangular_backbone_chanpak_20261116/CORE.md) gives $d_x\ge |C|-2\ge3$ at any mixed contact. An exterior containing a triangle would therefore give at least $3q+6$ paths, a contradiction. This leaves a fixed ordinal exterior and an explicit local contact budget, but does not bound all possible ordinal lengths. A viable +4 sibling and a reset for further digits remain open.

## Follow-through: a signed carry, not just a fork

The [carry note](CARRY.md) gives two strong two-stage routes that meet at exactly the same profile $9P+(6,0,0,0)$: take +2 and then zero, or take zero and then deposit +6. It also gives a further strong zero step. The all-order literal continuation scan in that note distinguishes actual second-step scalar hits from the missing synchronized fresh digit.
