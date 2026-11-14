# Optimal domination at order 36, and a uniform lower bound

This strengthens the previously exhibited 12-vertex dominating set in the fixed cubic lift of the three slot roles. On its 36 vertices the **ordinary domination number is exactly 10**. The very same graph and optimal dominating set still admit seven different role-label projections, hence seven different frame weights. The domination statistic is not sufficient to select the overlap, even when optimal.

Write the lift $G_r$ with points $(i,j)$, lines $(d,t)$, $i,d\in\{0,1,2\}$, and second coordinates in $\mathbb Z/r\mathbb Z$, with incidence $j=t+di$. For $3\mid r$, select the points $(0,j)$ with $j\not\equiv2\pmod3$ and all lines with $t\equiv2\pmod3$. Every unselected line meets a selected point at $i=0$. Every point with $i=1,2$ has a selected line because the three $j-di$ exhaust the residues mod 3; a point with $i=0$ not selected also meets a selected line. Thus $\gamma(G_r)\le 5r/3$. In particular the six-column graph has a dominating set of size 10. Nine centers on 36 vertices would have disjoint closed neighborhoods, so with $x$ point centers and $y$ line centers one needs $x+3y=3x+y=18$, impossible in integers.

## An integer slack obstruction at every order

In fact
$$
                       \gamma(G_r)\ \ge\ \left\lceil\frac{8r}{5}\right\rceil.
$$
Use the shift operator $zf(j)=f(j-1)$ on the cyclic columns. In the six point and line coordinates the closed-neighborhood operator is
$$
 M(z)=\begin{pmatrix}I_3&B(z)\\ B(z^{-1})^{\!T}&I_3\end{pmatrix},
 \qquad B(z)_{i,d}=z^{id}.
$$
Let $f$ be the binary indicator of a dominating set. Its integer coverage slack $s=Mf-\mathbf1$ is nonnegative entrywise. Expansion gives
$$
 \det M=-4(z^3+z^{-3}),\qquad
 \operatorname{adj}(M)\mathbf1=-2\mathbf1.
$$
The second equality on the constant vector uses $M(1)\mathbf1=4\mathbf1$ and $\det M(1)=-8$. Multiplying $Mf=\mathbf1+s$ by the adjugate and reducing modulo 4 gives a local congruence
$$
           (\operatorname{adj}M)s\equiv 2\mathbf1\pmod4.
$$
For just the four rows $i=0,2$ of the points and $d=0,2$ of the lines, replace each Laurent coefficient of the adjugate by a residue of minimal absolute magnitude. In each of the six columns the sum of absolute coefficients in these four rows is exactly 20, as verified by the integer Laurent expansion in the replay. After this replacement the congruence still holds, so the absolute value at each of the $4r$ output coordinates is at least 2. On the other hand its summed absolute value is at most $20\sum s$ by the triangle inequality (cyclic reduction can only improve this upper bound). Thus $8r\le20\sum s$. Cubicity gives $\sum s=4|D|-6r$, proving the bound.

For $r=3,6,9,12$ the rounded lower bound agrees with the constructive upper bound. The exact domination numbers at the corresponding graph orders $18,36,54,72$ are $5,10,15,20$. For later multiples of three the proof only gives the stated interval. There is no claim about the general cubic domination conjectures from this particular family.

Run **python replay_domination.py** for the determinant, coefficient norms and constructions, with no nonstandard dependencies. The geometry and role projection are in [the previous cubic comparison](../cubic/BALANCE_CUBIC_OVERLAP.md).
