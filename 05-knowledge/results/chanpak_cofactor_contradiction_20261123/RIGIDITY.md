# Rigidity of the physical both-sign fork over the factor seven

The existing fork has central ordinary rows $(16,9,67,21,38,15,59)$ and apex mask 68. Its marked lower and upper arcs are $(2,4)$ and $(0,2)$, respectively. With the inherited right address $y=2^{s+1}-1$, it gives center $3h$ for $h=280\cdot2^s+7$. We ask whether either arbitrary tagged ordinal grafts, or alternative fixed marked arcs with two unrestricted binary flanks, give a genuinely new free parameter. For this center the answer to each precise question is no.

## Arbitrary ordinal grafts with the two designated edges

Write the seven tagged-apex profile coordinates as $(h,q,c,t,e,r,k)$: ordinary height, two-slot count, positive start and end reads, signed start and end reads, and apex height. For a join of ordinary blocks $X\to Y$ under the same apex, the apex-height read is
$$k_{XY}=h_Yk_X+h_X(k_Y-h_Y)+c_Xt_Y+e_Xr_Y.$$
For the start coordinates the corresponding reads are $c_{XY}=q_Yc_X+h_Xc_Y$, $e_{XY}=h_Xe_Y$. Here $k_X\ge h_X$ by the positive odd-cycle packing formula, $|e_X|\le c_X$, and $t_X\ge0$. Empty blocks have $h=q=k=1$ and all other reads zero.

Let $C$ be the center and let $L,U$ be after its designated lower and upper reversals. Direct subset DP gives
$$
\begin{aligned}
C&=(81,1140,420,873,10,45,431),\\
L-C&=(-14,-148,0,-95,-2,-13,0),\\
U-C&=(28,224,0,115,8,27,-6).
\end{aligned}
$$
In particular the combination $U-C+4(L-C)=(-28,-368,0,-265,0,-25,-6)$ kills the signed right input; that is what controls arbitrary grafts.

Take *any* ordinary tournaments $A,B$, any apex contacts to them, and put $A\to C\to B$ with the corresponding two marked edits in $C$. From the join identities, the combined difference of the three apex heights is
$$
-6h_Ah_B
-28\{h_B(k_A-h_A)+h_A(k_B-h_B)\}
-(265c_A+25e_A)h_B-368c_A t_B.
$$
Every term is nonpositive, and the first is at most $-6$. To have both desired unit exits this expression must equal $2+4(-2)=-6$. Equality forces $h_A=h_B=1$, $k_A=k_B=1$, and $c_A=0$ (because ($265c_A+25e_A\ge240c_A$)). Thus both grafts are transitive, and A has no apex out-neighbor since each such vertex contributes a positive one-vertex starter path to $c_A$. A is only ordinal source padding and never changes the height. The lower difference now becomes $-2r_B$, so $r_B=1$. Conversely these conditions give both desired exits.

A transitive B with $k_B=1$ has tag word $I^aO^b$ (no $O$ before $I$, since that would be an apex triangle). Here I means a vertex that beats the apex, O one beaten by it. In such a transitive B the signed end read $r_B$ is one exactly when its first vertex is I: the contribution of each possible endpoint cancels by the choice of each earlier vertex, except at the first one. Therefore $a\ge1$. The final O's are only ordinal sinks, whereas $I^a$ is precisely the old address $2^a-1$. **The full arbitrary ordinal-graft result collapses to the original dyadic center values $420\cdot2^a+21$, with no extra width parameter.**

## All alternative fixed arcs, with two arbitrary binary flanks

Allow arbitrary nonnegative binary addresses $(x,y)$ to the left and right of the same center. Among all 28 edges in the seven-vertex core or between that core and its apex, a pair of opposite physical unit exits exists **if and only if** $x=0$ and $y=2^a-1$ for some $a\ge1$. The eligible signs are again realized by the two designated edges. This is a finite-certified classification of an unbounded address space, not a box scan.

For any candidate edge the height difference has the exact form
$$
a(f(x)+f(y))+qxy+tx+cy+r[x\text{ odd}]+e[y\text{ odd}]+k,
$$
with the inherited binary read $0\le f(v)\le v$. If $q\ne0$ and both addresses are positive, $|\Delta|=2$ forces
$$
\min(x,y)\le \left\lfloor\frac{2|a|+|c|+|t|+|e|+|r|+|k|+2}{|q|}\right\rfloor.
$$
This follows after dividing the required equality by $xy$; the same covering by fixed rows or columns includes both axes automatically. On a fixed row, the remaining slope $u=q x+c$ gives a finite upper bound on $y$ whenever $u,u+a$ have the same strict sign, using $0\le f(y)\le y$. The transpose handles columns. The only unbounded rows in this test are
$$
(x,\text{edge})=(0,02),(0,24),(0,36),(0,46),(1,12),(1,25);
$$
there are no unbounded columns for this center. Pair each exceptional row with each possible second arc. Once parity is fixed the two equations are linear in $y,f(y)$. A nonzero determinant gives at most one exact candidate. The only compatible rank-one pair is the old $(0,2),(2,4)$ at $x=0$, where the equations reduce to $f(y)=0$ and odd $y$. The apex edits (the $q=0$ cases) have $a=0$, so any pair of them is instead a two-equation linear system in $x,y$ after specifying both parities. The independent subset-DP and exact elimination scripts recover all coefficients and check the finitely many residual candidates.

    python3 sevenfork-all-arc.py

This second theorem allows the pair of exits to depend on the address, but not edits touching a variable flank. The first allows completely general tagged ordinal grafts, but retains the designated pair of unit edges. Neither rules out a different core, a nonordinal attachment, or surgery on the variable vertices.
