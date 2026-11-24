# A literal Chan–Pak seam in a pentagon observer tournament

Use the pentagon labels in the order $(0,1,4,3,2)$. On source vertices $0,\ldots,4$ use the ascending transitive chain; on observer vertices $5,\ldots,9$ use its reverse. A source beats an observer exactly when their pentagon labels agree or are adjacent. This is a strong blue tournament. Its Hamilton-path height is 6461; reversing the single source arc $(3,4)$, or its reflected observer arc $(8,9)$, gives height 6463. Neither reversal preserves the original blue presentation.

The subset recurrence is $f_{\{v\},v}=1$ and $f_{S\cup\{w\},w}\mathrel{+}=f_{S,v}$ for every $v\to w$ and $w\notin S$. Summing the full-set vector gives the height. The replay enumerates the 12 label orders modulo the dihedral action (fix initial label 0 and canonicalize reflection), and all 45 single flips in each. The only exact difference-two flips are the two above. Both changed tournaments are strong.

## One-vertex lifts

Adjoin an apex beating the core vertices indicated by mask $w$. For $w\notin\{0,1023\}$ both sides of the seam stay strong and the exact same core arc is flipped. The complete exact-difference-two list is:

| $w$ | lower core extension | flipped core extension |
|---:|---:|---:|
|60|24445|24447|
|102|29895|29897|
|155|48479|48481|
|226|29455|29453|
|259|52083|52085|
|482|32669|32667|
|510|16895|16893|
|511|13681|13679|
|614|34321|34319|
|673|43245|43247|
|678|38217|38219|

For a split into sets $(A,B)$ before and after the apex, multiply the number of paths through $A$ ending outside the mask by the number through $B$ starting inside the mask. Use value one for empty $A$ or $B$ in its respective contribution. Sum over all ordered splits. This gives the complete 1022-mask replay, not an independent random sample.

## Two limits on immediate renewal

1. Replace any marked vertex by a chain of $k$ twins, all carrying the old external incident orientations. For a skeleton with $t$ separated marked runs, the $k$ chain vertices realize it in $F_t(k)=t!S(k,t)$ ways: assign each chain vertex to one of the $t$ nonempty ordered runs, whose internal order is forced. A subset DP counts the skeletons by $t$. In nine sites, the $t\ge2$ nonzero coefficient differences from the physical seam all have one sign and the $t=2$ coefficient has magnitude at least four. Thus for $k\ge2$ their differences are never $\pm2$. At site 8 the signed difference is $2-15F_2-2F_3+12F_4$, giving $-28,-100,8$ for $k=2,3,4$. At $k=5$ we have $F_3\ge3F_2$ and $F_4\ge F_3$; the recurrence $F_t(k+1)=t(F_t(k)+F_{t-1}(k))$ preserves these inequalities, so the difference is greater than two for all $k\ge5$.

2. Put a transitive tail of length $k\ge1$ after all core sites and before a new apex, and let the apex beat precisely the core mask $W$. Define $U$ to count splits whose after-apex part is nonempty and starts in $W$, with no condition on the before-apex path; define $V$ similarly but require the before-apex path, when nonempty, to end outside $W$. Then
$$H(k)=H_{\mathrm{core}}+V+(2^k-1)U.$$
If the apex is last, count $H_{\mathrm{core}}$. Otherwise there is a nonempty after-apex core path. Each nonempty selection of tails before the apex removes the condition on the before-apex core path, and the complement follows the after-apex path. If all tails are after the apex exactly the $V$ channel applies. Thus a uniform exact gap requires the two integer slopes $U$ to agree and the intercept difference to be $\pm2$. The complete search over the 12 pentagon observer orders, all 45 core flips each, and all 1023 nonempty common masks finds no such pair in this one-tail construction. The single compact marked Petersen witness with outgoing rows $(966,940,888,753,483,193,386,772,552,112)$ has height 8265, no single-arc difference two, and no affine fork in this same one-tail test. These are scoped obstructions, not an obstruction to other lifts or to Chan–Pak.

## A renewable fork after one preparatory apex

The one-tail obstruction above is evaded by a single shared preparatory apex.
1. Start from the ten-site rows $(254,636,440,912,832,24,44,114,227,485)$ and retain the marked source arc $(3,4)$.
2. Adjoin a vertex 10 with outgoing mask 194: it beats exactly vertices $1,6,7$ of the old core.
3. Place a chain of length $k\ge1$ after all 11 core vertices. Place a final apex behind the chain, so all chain vertices beat it. Its outgoing mask into the eleven sites is 1240: it beats exactly vertices $3,4,6,7,10$.

The exact tail parameters are $(U,K)=(245121,-74075)$ and $(245121,-74073)$, where $K=H_\mathrm{core}+V-U$. Thus the single marked core arc gives the strong-tournament fork
$$245121\cdot2^k-74075\quad\longleftrightarrow\quad245121\cdot2^k-74073,$$
on $k+12$ vertices. The last apex connects the tail back to the core through a positive contact, so both sides are strong. The first three pairs are $(416167,416169)$, $(906409,906411)$, and $(1886893,1886895)$; direct subset DP checks all three. The lower point is prime for $k=1$ and $k=6$ by trial division, while the upper rail is always divisible by 3.

More generally the growing slot can be any nonempty tournament $A$, carrying uniform arc orientations to all core vertices and the final apex. Set
$$Q_2(A)=\sum_{S\subseteq V(A)}H(A[S])H(A[V(A)\setminus S]),\qquad H(\varnothing)=1.$$
The exact heights are $245121Q_2(A)-74075H(A)$ and $245121Q_2(A)-74073H(A)$, giving gap $2H(A)$ with the same physical reversal. This follows by weighting the nonempty-before slot choices by $Q_2(A)-H(A)$ instead of $2^k-1$, and multiplying the other two channels by $H(A)$. The directed triangle checks the nontransitive module case directly.
