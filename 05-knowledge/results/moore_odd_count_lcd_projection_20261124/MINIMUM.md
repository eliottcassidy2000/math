# Minimum codewords identify actual neighborhoods

Use an odd-degree diameter-two Moore graph of degree $d\ge3$, $n=d^2+1$ and adjacency $A$, and work over $\mathbb F_2$. Write $U=\ker A$, $V=\operatorname{im}A=U^\perp$. The odd vectors in $V$ satisfy $Ay=y+\mathbf1$ over $\mathbb F_2$, by $A(A+I)=J$ and $A\mathbf1=\mathbf1$. The even vectors in $V$ are exactly $\ker(A+I)$: one direction follows from the parity identity, and for the other, $y=Ay$ puts $y$ in the image and $Jy=(A^2+A)y=0$.

## Odd band: no extra first-height impostors

**The smallest odd weight in $V$ is $d$, and its words are exactly the $n$ open neighborhoods.** For an odd word of support $S$ and weight $m$, let $t_z=|N(z)\cap S|$ over the integers. Outside $S$ these counts are odd so at least one, while inside they are even. Thus $dm\ge n-m$, ruling out $m<d$. Neighborhood rows attain $d$.

For equality $m=d$, the total excess over the exterior baseline is
$$
 \sum_{z\notin S}(t_z-1)+\sum_{z\in S}t_z=d-1.
$$
Count pairs of $S$ via common neighbors, adding back the adjacent pairs (which have none):
$$
 \binom d2=\sum_{z\notin S}\binom{t_z}2+
      \sum_{z\in S}\left(\binom{t_z}2+\frac{t_z}{2}\right).
$$
Every term is at most $d/2$ times its corresponding excess, and their sum already attains $d(d-1)/2$. Thus every positive interior $t_z$ would have to equal $d$, impossible because $d$ is odd and the interior counts are even. Every exterior count exceeding one must equal $d$; the budget forces exactly one. Therefore $S$ is a whole open neighborhood. Distinct such neighborhoods cannot coincide because two vertices have at most one common neighbor.

At degree 57 this requires **exactly 3,250 words of weight 57**, and none at smaller odd weights; since $\mathbf1\in V$ and $n$ is even, the complementary band at weight 3,193 also has exactly 3,250 words. Under the triangle-cell readout in the main note, the first odd-exponent band is $3^{57}$, but the value itself has been known constructible without Moore. What the degree-57 construction must supply is the certified *multiplicity and arrangement*, not a witness for that odd integer. At degree 3 there are exactly 10 weight-three words, and at degree 7 exactly 50 weight-seven words, also checked by full enumeration.

## The two even codes have distinct equality geometries

For a support $S$ of size $m>0$, the pair count is $\sum_z\binom{t_z}2=\binom m2-e(S)$.

- In $\ker A$, all $t_z$ are even, so the left side is at least $dm/2$. Therefore $m\ge d+1$. At equality $S$ is independent, and all $t_z$ are 0 or 2.
- In $\ker(A+I)$, $t_z$ is odd on $S$ and even off it. Therefore the left side is at least $(dm-m)/2$, and $e(S)\ge m/2$. Again $m\ge d+1$. At equality the induced graph on $S$ is a perfect matching, with $t_z=0,2$ outside. This is different from the [previous maximum cap](../moore_chanpak_counts_20261124/degree57-ports.md), whose support is a union of cycles.

For a hypothetical weight-58 word of $\ker(A+I)$ at degree 57 take $x=1_S$, $y=1_Y$ with $Y$ the outside vertices having two neighbors in $S$. Then $Ax=x+2y$ over the integers and the Moore identity gives $Ay=27x+29\mathbf1-2y$. The partition $(S,Y,V(G)\setminus(S\cup Y))$ is therefore equitable, with sizes $(58,1624,1568)$ and matrix
$$
 \begin{pmatrix}
 1&56&0\\
 2&27&28\\
 0&29&28
 \end{pmatrix}.
$$
The calculation is conditional on attaining the code bound. The known degree-seven case is a useful hostile test: its open null code has 750 words of weight 8, but its closed null code has **none** of weight 8 or 10; its minimum is 12, attained 1,050 times. At degree three both lower bounds are attained at weight four (5 and 15 words, respectively). Thus there is no reason to hard-code existence of a closed weight-58 word into the degree-57 search. The complete small code weight distributions are in **census.out**.

**Practical use.** For the degree-57 LCD search, the odd first-band count is a mandatory rejection test, while the two even minimum patterns are separate conditional branch tests. This parallels the earlier Golay contact work in looking at a whole weight spectrum rather than a lone height, but the code here is complementary dual, not self-dual. No code of the required size satisfying the squared-degree criterion is produced.
