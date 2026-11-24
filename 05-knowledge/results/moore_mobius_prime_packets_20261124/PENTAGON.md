# Pentagon arrays: an exact completion test and a conditional residual spectrum

Take a pentagon $C$ in a Moore graph, let $t=d-2>0$, and let $X_i$ be its $t$-element sets of exterior neighbors. Index $C$ in cyclic order, so the graph on $X=\bigcup_i X_i$ is a lift of the **pentagram**. Let $Y$ be the remaining $t(t-1)$ vertices. Every $y\in Y$ has one contact in every $X_i$, by its unique two-step route to $c_i$, giving a 5-column row array. The induced graph $B$ on $Y$ is $(t-3)$-regular when $Y$ is nonempty.

## Not just a minimum case: exactly why completion can fail

For any pair of columns $i,j$, the rows realize all $t(t-1)$ ordered pairs except one forbidden matching. If the columns are adjacent in the pentagram, the forbidden matching consists of already adjacent exterior vertices. If their pentagram distance is two, it consists of the vertices matched by the length-two exterior route. For a pair not in these matchings its unique common neighbor cannot be central or exterior, so there is a unique row of $Y$. For the forbidden matching the short route already exists; no row can realize it.

An extension by $t$ rows to an $\operatorname{OA}(t^2,5,t,2)$ exists **if and only if the pentagram monodromy is the identity**. Necessity: each added row must use the forbidden pairs on all five successive pentagram links. It must close; in addition each value of a fixed column has to occur exactly once among the new rows (each forbidden pair is to occur once). Thus all $t$ sheets are fixed. Conversely if the monodromy is the identity, label each exterior 5-cycle by a symbol. The uncovered pairs are exactly the equal-symbol ones for every column pair, so add the $t$ constant rows. In fact the maximum number of extra rows that can be added without repeating an ordered pair is exactly the number $c_1$ of fixed sheets: each fixed sheet gives one consistent row and two distinct fixed sheets use different values in every column. The minimum $9^t$ is therefore an exact array-completability test, rather than only a small numerical height; the decoded deficit $t-c_1$ is the number of rows still unobtainable after every permissible addition.

In the minimum case each existing row has five distinct symbols (because the equal-symbol pairs are forbidden). For a row $a$, its neighbors in $Y$, in column $i$, take each symbol outside
$\{a_i,a_{i-2},a_{i+2}\}$ exactly once. For a fixed exterior vertex $z\in X_i$, the routes to it from row $a$ already using zero or one exterior edges occur at exactly those three symbols: the direct contact and the two matched contacts via columns $i\pm2$. All other $z$'s must have one common neighbor with $a$ in $Y$. This is stronger than only forbidding an equal coordinate between adjacent rows.

## The rank-and-trace prediction at degree 57

Write $N$ for the $Y$-by-$X$ incidence matrix of the row contacts, and $A_X$ for the adjacency of the exterior graph. In the minimum case, $N^{\mathsf T}N$ has diagonal blocks $(t-1)I$ and all off-diagonal blocks $J-I$. On the five copies of the symbol-centered space the Gram is $tI_5-J_5$ tensored with identity. On that centered sector it is invertible for $t=55$; the column-constant sector maps onto one constant direction. Thus $N$ has rank $5(55-1)+1=271$. In the Hoffman–Singleton calibration $t=5$, the centered Gram is singular instead, so the degree-57 rank calculation must not be extrapolated directly from it.

The exact Moore block equations give
$$
 BN=J-N(I+A_X),\qquad B^2+B=(d-1)I+J-NN^{\mathsf T}.
$$
On the image of centered column vectors, $-I-A_X$ controls the action. The five-column cycle has eigenvalues $2$ and $(-1\pm\sqrt5)/2$, the latter twice each, so $B$ has $-3$ with multiplicity 54, and the two conjugate roots $(-1\pm\sqrt5)/2$ with multiplicity 108 each. The constant mode is 52. On the orthogonal complement of $\operatorname{im}N$, dimension $2970-271=2699$, the second identity leaves only 7 and $-8$. Trace zero determines their multiplicities: the contribution from the image is $52-162-108=-218$, so the remaining multiplicities are 1454 and 1245. This is an exact necessary residual spectrum for a **minimum pentagon**, not an assumed feature of every degree-57 candidate. The array and row-neighbor rule can be checked without computing any tournament paths; a failure would reject this minimum branch, not the whole Moore problem.
