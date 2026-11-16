# One splice for all three exact ternary digits

**Result.** There is an unbounded family of strong sources with a single physical step realizing all three counts $3H,3H+2,3H+4$. The three outputs share the same tripled rooted continuation $(E,O,g)$. The middle setting uses one old-edge reversal as well as the new vertex; the outer settings retain the source induced.

Here $E,O$ count ordered path covers in which the first path starts at the root and has even or odd vertex count, respectively, and the second path may be empty. The coordinate $g$ counts Hamilton paths after deleting the root. This is the [fixed-root continuation](../BRIDGE_RESET.md), not a claim about arbitrary new roots.

## Physical templates

The following rows encode outgoing-neighbor masks, vertex 0 in the low bit. Root at vertex 7 and mark vertex 6 for substitution.

| Setting | Row masks | Base $(H,E,O,g)$ |
|---|---|---|
| Source | 118,124,248,113,224,192,128,11 | (163,84,89,5) |
| 0 | 374,380,504,369,224,448,128,11,208 | (489,252,267,15) |
| 2 | 118,372,504,115,224,192,128,11,249 | (491,252,267,15) |
| 4 | 118,124,504,369,224,192,128,11,243 | (493,252,267,15) |

Each setting adds vertex 8, with its specified contacts. Setting 2 also reverses $1\to3$ to $3\to1$ in the old graph. None changes the marked slot or its only outward neighbor, which is the root. Replacing vertex 6 by a transitive block of length $s\ge1$ gives
$$
q=84\cdot2^s-5,\qquad
(H',E',O',g')=(3q+d,\ 126\cdot2^s,\ 126\cdot2^s+15,\ 15),\qquad d\in\{0,2,4\}.
$$
Deleting vertex 6 gives the $s=0$ case: $79\to237,239,241$, all with continuation $(126,141,15)$. The source's continuation there is $(42,47,5)$. The certificate checks this case literally; no run argument for a nonempty module is assumed at the empty boundary.

## Arbitrary marked module and finite proof

Let $M$ be nonempty, $u=H(M)$ and $R=P_M(2)$ the number of two ordered path slots covering $M$, each allowed to be empty individually. Write $P=(H,E,O,g)$. Then
$$
P(C[M])=(84R-5u,\ 42R,\ 42R+5u,\ 5u),\qquad
P(A_d[M])=3P(C[M])+(du,0,0,0).
$$
In particular the increment is the absolute $2$ or $4$ when $u=1$.

Compress each maximal run of the marked module to one symbol in a path. The marked slot can have at most two runs in a Hamilton path, because the only nonterminal exit consumes the root. In an ordered cover whose first path starts at the root, a run in either path must be terminal: the sole exit is already used at the beginning of the first path. Thus at most one run occurs in each path. For deletion of the root, any marked run must be terminal.

Enumeration of the legal fixed-vertex and run words gives:
* Hamilton coefficients at one and two runs: $(163,84)$ for the source; $(489,252),(491,252),(493,252)$ for the children.
* Root deletion has just one run, with coefficient 5 in the source, 15 in the children.
* Rooted two-path patterns $(0,1),(1,0),(1,1)$ have unsigned coefficients $(89,84,84)$ and signed coefficients $(-5,0,0)$ in the source. The sign is parity of the fixed vertices in the first path, even positive. Each child has precisely triple those tables.

One run fills in $u$ ways, two runs in $R-2u$. The unsigned source sum is therefore $84R+5u$. The first signed coefficient is for the module entirely in the second path, giving $-5u$ with no further parity factor; the other signed coefficients vanish. These give the two rooted counts and prove the identities. Both cases with a marked module and the empty boundary are checked by independent subset DP in **node check_ternary.js**, which also verifies the literal arc-change claim and several nontransitive substitutions.

## A real reset for zero, but not a fresh second digit

For ordinal sums, both $H$ and $P(2)$ multiply. Every path crossing ordinal components crosses in their order, and an ordered cover restricts to an ordered cover on each component with empty slots allowed. A transitive pair has $(H,P(2))=(1,4)$, and a directed triangle has $(3,12)$ by its eight subset splits. If the module has a consecutive pair of ordinal singleton factors, replace that pair with a directed triangle. The new module has $(u',R')=(3u,3R)$. Evaluating it in the source template gives **the exact whole profile of the zero child**, not just its scalar. This is a re-encoding; the zero child need not equal the new source graph.

Starting with a transitive block of length $s$, this allows $\lfloor s/2\rfloor$ successive zero returns within the same family, each consuming two singleton factors. After $k$ such returns the current module has $u=3^k$ and $R=3^k2^s$, so its three new counts are $3^{k+1}q+d\,3^k$. Those are carried digits, not a fresh unit at that later stage.

The nonzero children cannot return to this source family with their full profiles, even via different nontransitive modules. Matching $E,g$ forces $R'=3R,u'=3u$, in which case the source formula necessarily gives the zero-child scalar. Equivalently the source equation is $H=2E-g$; the three outputs have exact defects $H-2E+g=du$. Thus the obstacle to iteration is now explicit: a new ready source must carry the nonzero defect without forcing digit spacing to grow. The full one-step deposit is proved; indefinite fresh deposits and Chan–Pak are not.


## Why the reversed contact matters at the literal seeds

Flipping $1\to3$ before adding the new vertex is not by itself a neutral preparation. Its intermediate profile is $(88R-7u,44R,44R+5u,5u)$, from the same run table. The new contacts compensate for the change. As a separate check, enumerate all contacts of one new vertex to the unchanged old graph. At $q=79$ there is exactly one induced single-vertex extension giving 239, but its continuation is $(127,142,15)$ instead of $(126,141,15)$. At $q=163$ there is no induced single-vertex extension even at the scalar 491. This explains why the one-edge recontact matters for this step, without making a no-go claim about larger exteriors.

At the literal seeds $q=79,163$, the exact second-round boundary is stronger than the same-template obstruction. See [SECOND_ROUND.md](SECOND_ROUND.md): even without an induced retention constraint, only the three carried offsets $0,6,12$ realize the demanded tripled continuation within the nine-count window.

