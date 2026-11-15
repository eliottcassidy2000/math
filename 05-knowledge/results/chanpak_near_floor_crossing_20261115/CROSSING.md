# The Chan–Pak crossing lemma is false, even for an actual unit cut

[HYP-9732](../../hypotheses/HYP-9732-crossing-lemma-a-double-crossing-cut-has-defect-at-least-three-so-unit-cuts-are-monotone.md)
proposed that a tournament cut $X$ admitting a Hamiltonian path with two
$X\to X^c$ steps forces $H(T+x)-2H(T)\ge3$, where the new vertex $x$
beats exactly $X$. This fails for strong tournaments, even if the defect
is exactly $-1$. The finite observations through order 8 in THM-5192(C)
are unaffected. No minimal-order claim is made.

## Literal strong unit witness

Vertices are 0 through 10, with each outgoing row encoded by the sum of
$2^j$ over vertices $j$ beaten by its vertex. Use

    rows = (1982,1976,2042,432,2016,1984,1931,1792,1536,1032,8)
    X = {3,5,7,8,9,10}  (mask 1960)

A crossing witness is $(0,2,1,3,4,5,6,7,8,9,10)$. The two forward
color crossings are $3\to4$ and $5\to6$. Direct traversal gives

    number of X-to-complement steps   0     1     2
    Hamiltonian paths                194   465   156

and $H(T)=815$, $H(T+x)=1629=2\cdot815-1$. Exactly 37 of the enlarged
tournament's Hamiltonian paths have $x$ between two endpoints whose
junction would fail on deleting $x$. In the insertion identity the
defect is therefore $E+R-Z=156+37-194=-1$. Thus **unit cuts need not
be monotone**. This is not a counterexample to sign freezing itself,
which requires the existence of additional grandchildren.

## An exact two-shelf family

The witness is in a two-parameter family. The seven-vertex base $0,\ldots,6$
has every arc increasing except the reversals $2\to1,6\to0,6\to1,6\to3$.
Follow it by a transitive shelf $A$ of size $r\ge1$ and a second
transitive shelf $B$ of size $s\ge1$. All interblock arcs go forward
except that every member of $B$ beats 3. Set $X=\{3,5\}\cup A\cup B$.
Write $U=2^r,V=2^s$. Then
$$
\begin{aligned}
H(T)&=68UV-68U-19V+75,\\
H(T+x)&=133UV-133U-31V+157,\\
H(T+x)-2H(T)&=(7-3U)V+3U+7.
\end{aligned}
$$

To prove the formulas, every exit from $B$ goes through vertex 3, so
there can be at most two maximal consecutive runs from $B$. Every
exit from $A$ goes into $B$; to return to $A$ must therefore use
$B\to3$. There can be at most two runs from $A$ as well. These facts
remain true after the insertion since $x$ beats both shelves. Collapse
the individual runs to shelf symbols while retaining the sequence of
other vertices. All interblock adjacencies are independent of which
shelf vertices each run contains, and inside each run the order is
unique. A single run has one assignment, and two runs have $2^r-2$
or $2^s-2$ assignments by choosing a proper nonempty subset for the
first run. Hence both counts are affine separately in $U$ and $V$
for all positive shelf sizes. Evaluate by literal DFS at the four
interpolation vertices:

    r   s   H     Hplus
    1   1   173   361
    1   2   407   831
    2   1   309   627
    2   2   815   1629

The four values determine the displayed identities uniquely. A
second way to recover their coefficients is to traverse with two
vertices in each shelf, count paths by their number of shelf runs,
and divide by two for each shelf appearing in two runs. The
exact verifier confirms equality of these two computations.

For every $(r,s)$ the old tournament is strong: vertex 3 reaches
4, 5 and 6, vertex 6 reaches 0 and 0 reaches 2 which reaches 1; layer
$B$ returns to 3. The early vertices reach 3 either directly or via 6, and $A$
does so via $B$. The enlarged tournament is strong as well since
$x$ has both an incoming and outgoing arc to the old strong
component. For every $(r,s)$ the old path $(0,2,1,3,4,5,6,A,B)$,
with each shelf in increasing order, has the two required crossings.

For $r=1$ the defect is $V+13\ge15$. For $r=2$ it is $19-5V$;
the only unit value is $-1$, at $s=2$. For $r\ge3$ it is at
most $21-3U\le-3$, and tends to $-\infty$ as either shelf grows
with the other appropriately fixed. Thus no constant lower bound,
even a negative one, holds for crossing cuts in strong tournaments.
The same marked interface changes the sign of its defect under inflation.
When $U,V\to\infty$, the defect is not negligible:
$\delta/H(T)\to-3/68$.

The insertion identity itself splits into three bilinear polynomials:
$$
 E=12UV-12U+12,\quad
 Z=17UV-17U-8V+22,\quad
 R=2UV-2U-V+17.
$$
In the old tournament every arc from a shelf to the exterior of the
shelves can only go to 3, which is marked. So the only vertices that
can originate forward color crossings are 3 and 5, and the histogram
has no entry above 2. The shelf slot expansion applies to each
crossing class and also to the reverse-junction class separately;
its defining adjacency has endpoints in fixed or uniformly joined
blocks. Evaluate at the same four interpolation vertices, as the
verifier checks. The actual deficiency is at the leading $UV$
coefficient: the reverse covers add 2 to the excess-crossing
coefficient 12, but the missing-path coefficient is 17.

## Independent direct replay

Run python3 verify_crossing.py (or python3 -O verify_crossing.py). This is a
literal permutation DFS following allowed arcs, different from the
subset endpoint dynamic program and cut-table evaluation used to
find the witnesses. It checks the tournament conditions, the arc
count histogram and reverse-junction count, both witness defects,
the four interpolation vertices, and the shelf-slot expansion.
