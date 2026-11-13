# Mirrored bridges: four exact full-count and blue-path laws

Let $T$ be nonempty and let $a$ be a named vertex. Write $H=H(T)$,
$E=E_a(T)$ for the number of Hamilton paths ending at $a$,
$M=H(T-a)$ (with $H(\varnothing)=1$), and
$$
 L=\sum_{a\in X\subseteq V(T)}
       S_a(T[X])H(T[V\setminus X]).
$$
Here $S_a$ counts paths starting at $a$. In other words, $L$ counts
two ordered disjoint paths that cover $T$, the second starting at $a$
and the first allowed to be empty. No arc is required at their seam.
Cutting a full path immediately before $a$ gives $L\ge H$, and
deleting $a$ from a path ending there gives $E\le M$. In addition
$L\ge M$ by taking the second path to be the singleton $a$.
For a nonsink a term $X=\{a,v\}$ with $a\to v$ is positive by
Rédei, so $L>E$; for a sink, $L=M=E=H$.

Start with $T\to T^*$, optionally with a middle singleton $c$ and
$T\to c\to T^*$. Reflect vertices to their mates (fix $c$); this is
an anti-involution. A **bridge flip** reverses $a\to a^*$.
A **middle flip** reverses both $a\to c$ and $c\to a^*$.
Both preserve the anti-involution. Below $F$ counts Hamilton paths
fixed by reflecting the labels and reversing the path order.

| Construction | Total | Fixed paths $F$ |
|---|---:|---:|
| even, bridge flip | $H^2+L^2-E^2$ | $H+L-E$ |
| odd, bridge flip only | $H^2+2L^2$ | $H$ |
| odd, middle flip only | $(H-E)^2+2HL+L^2$ | $H-E+L$ |
| odd, both flips | $(H-E)^2+2HL+3L^2-2M(E+L)$ | $H-E+L$ |

**Even.** The only lost paths of the ordinal sum use the original bridge,
giving $E^2$. A new path has chronology
$L_0R_0a^*aL_1R_1$, with the left path split just before $a$
and the right path split dually just after $a^*$. More precisely, the
four blocks are the left prefix, the right prefix ending at $a^*$, the
left suffix starting at $a$, and the right suffix. Every surviving
forward cross seam is automatic since $a^*$ precedes $a$. This gives
a Cartesian product of two $L$-sets. Reflection transposes the product;
the diagonal contributes $L$ fixed paths. Likewise the old product
contributes $H-E$ fixed paths.

**Odd, bridge only.** An old path places $c$ at the left/right seam,
so none of them are lost. Each new even path has exactly two possible
slots for $c$, the two forward side changes (including the first or
last position if the corresponding outer block is empty). Reflection
exchanges the slots; no additional fixed paths.

**Odd, middle only.** For a nonreturning path with $c$ at the ordinary
left/right seam, its left endpoint cannot be $a$ and its right start
cannot be $a^*$, giving $(H-E)^2. With $c$ inside the left side,
it is immediately before $a$, and the seam between the two left paths
is free: this gives $LH$, including an initial $c$. Inside the right
side gives the reflected $HL$, including a terminal $c$. A returning
path uses $a^*ca$, and the two unconstrained half splits give $L^2$.
The seam and return terms contribute diagonals $H-E$ and $L$.

**Both flips.** From the preceding law lose $ME$ paths in each of the
two inside cases: the left suffix at $a$ (or its reflected counterpart)
is a singleton and the opposite full path starts at $a^*$ (or ends
at $a$). The ordinary seam term and $a^*ca$ term survive. New paths
use the direct backward bridge. They have two slots for $c$, except
the first slot is blocked if the right prefix is the singleton $a^*$,
and the second is blocked if the left suffix is the singleton $a$.
So the new term is $2L^2-2ML$, paired off under reflection.

For the first three constructions, if $|T|\ge2$ the result is strong
iff $a$ reaches every vertex of $T$. In each case all entries to the
left from the other side land at $a$; conversely the reachability from
$a$ and its dual on the right, together with an available forward
crossing, close a spanning return. For $|T|=1$ the even case is not
strong; both odd single-flip cases are directed triangles. The count
formulas include these boundaries.

## What the laws do and do not cover

For the even bridge, $L\equiv E\pmod 2$ by Rédei on the output.
Consequently every output is $1\bmod 4$. The odd bridge is an
$x^2+2y^2$ norm with $y\ge x>0$, hence is $1$ or $3\bmod 8$,
and misses even some admissible residues as absolute targets. For
example $11,17,19,25$ are impossible in this family: if $H\ge3$
the output is at least 27, and if $H=1$ the input is transitive
(a directed triangle's three paths each extend by insertion).
For a transitive input $L=2^{d^+(a)}$; the resulting counts are
$1+2\cdot 4^{d^+(a)}$. A directed triangle input has
$(H,E,L)=(3,1,3)$, giving 17 in the even bridge, 27 in the odd bridge,
and 31 in the middle-only construction. Both single-bridge families
miss 23. The separate paired-site census establishes that 23 is actually
absent in the entire self-converse sector, although ordinary tournaments
realize it.

Appending $k\ge1$ sinks to $T$ leaves $H$ fixed, sets $E=0$,
and sends $L\mapsto2^kL$, since each sink can go to either half and
there is one permissible increasing order within each half. Prepending
sources leaves $H,E,L$ unchanged, as every new source must go in
the first half. Thus the single-bridge outputs have literal dyadic rays,
but each such ray is sparse. The port $M$ is needed once the two
sorts of flip interact; the term $-2M(E+L)$ is an exact obstruction
to treating them as independent digits.
