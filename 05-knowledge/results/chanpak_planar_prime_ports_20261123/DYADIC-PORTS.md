# Exact dyadic port inverses and the factor-seven lane

Build a tournament from a transitive chain and an apex. A contact is $O$ if the apex beats that chain vertex, and $I$ if it loses to it. Every word beginning $O$ and ending $I$ gives a strong tournament. For any such word $w$, odd-cycle packing, or directly assignment of chain vertices to the two sides of the apex, gives
$$
H(w)=1+\#\operatorname{adj}_w(OI)+
\sum_{i<j:w_i=O,w_j=I}2^{j-i-1}.
$$
An odd cycle chooses these two endpoints and an even subset of the intervening vertices. When the two endpoints are adjacent, the unique empty interior subset has weight two rather than one; that corrects the geometric sum by one.

## First chart: unique positive inverse, but no factor seven

For $O^aI^b$, $a\ge2,\ b\ge1$, put $x=2^a,y=2^b$. Geometric summation gives $H=(x-1)(y-1)+2$. Flipping the final $O$ moves to $O^{a-1}I^{b+1}$, so its signed gap is $\delta=x/2-y$. Thus
$$
\det\frac{\partial(H,\delta)}{\partial(x,y)}
=-\bigl(x/2+y-3/2\bigr)<0
$$
on the physical quadrant. The exact inverse on the positive branch needs only a single square test: with $w=x-1$,
$$
w^2-(2\delta+1)w-2(H-2)=0.
$$
Take its unique positive root and require it to be integral, $x=w+1$ and $y=x/2-\delta$ to be powers of two, and the two widths to obey their bounds. The discriminant is $(2\delta+1)^2+8(H-2)$. This is an actual global positive inverse on the chart, not just local formal invertibility.

But **seven never divides the height in this chart**, even with $a=1$ allowed. Each of $2^a-1,2^b-1$ belongs to $\{0,1,3\}$ modulo seven; their products never equal $5=-2\bmod7$. A nonsingular real inverse therefore need not touch even a single point in the target factor-seven lane.

## Three-run repair: exact paired-port inverse and a prime/7p fork

Take $W=OI^aOI^bOI^c$ with all $a,b,c\ge1$, and write $x=2^a,y=2^b,z=2^c$. Geometric summation gives
$$
H=4xyz+(2y+1)(z-x)+1-y.
$$
For clarity, the elementary two-out-contact read is $H(OI^rOI^s)=2^{r+s+1}+2^s-2^r+1$: its two adjacencies plus the constant give $3+(2^r-1)+(2^{r+1}+1)(2^s-1)$. Therefore flipping the middle or last $O$ in $W$ to $I$, respectively, gives the literal heights
$$
H_1=4xyz+z-2xy+1,\qquad H_2=4xyz+2yz-x+1.
$$
For fixed $y$ the signed gaps $g_i=H_i-H$ give a genuine constant-Jacobian affine map:
$$
(g_1,g_2)=(x-2yz+y,\ 2yx-z+y),\qquad J=4y^2-1.
$$
Putting $A=g_1-y,B=g_2-y$ gives exact inverse
$$
x=\frac{2yB-A}{4y^2-1},\qquad z=\frac{B-2yA}{4y^2-1}.
$$
Require integral powers of two $x,z\ge2$ and, if the center was prescribed, verify its displayed formula. If $y$ is unknown, enumerate the Mersenne divisors $2y-1=2^{b+1}-1$ of $g_2-g_1$, since
$$
g_2-g_1=(2y-1)(x+z),\qquad
g_1+g_2=(2y+1)(x-z)+2y.
$$
This is a finite exact inverse for arbitrary two-gap inputs to this family. The difference is at least 12, so these ports cannot be simultaneous opposite unit exits.

Unlike the two-block chart, factor seven occurs. Since the residues of $2^a,2^b,2^c$ have period three, there are only 27 cases to substitute:
$$
\begin{array}{c|c|c}
(a,b,c)\bmod3& (g_1,g_2)\bmod7 & J\bmod7\\
\hline
(2,0,0)&(3,1)&3\\
(0,1,0)&(6,5)&1
\end{array}
$$
These are exactly the two classes with $7\mid H$. Thus at every physical center on the seven lane the paired port inverse is nonsingular modulo seven, and both literal reversals leave that lane, with fully determined residues. This is a theorem about distinct physical endpoints on this particular chart, not about all one-arc reversals.

For $(a,b,c)=(5,3,6)$, the same 18-vertex strong tournament has center $66073=7\cdot9439$ and the two literal one-arc neighbors $65089,66529$, both prime; $9439$ is prime as well. This is an actual prime–$7p$–prime fork (not unit-spacing). For a separate exact example obtained symmetrically, fix $b=1$ and $a=c$; the center is $2^{2a+3}-1$. Seven divides it just when $3\mid a$; with $2a+3=3r$, $r$ odd and $r\ge3$, factor
$$
2^{3r}-1=(2^r-1)(2^{2r}+2^r+1).
$$
For $r\ge5$ both factors exceed seven and seven divides exactly one, so the cofactor is composite. For $r=3$ the cofactor is $73$, the unique prime on this symmetric subfamily.

The earlier marked-cut Jacobian and two-apex packets are local contact-space results. Here the parameters are actual dyadic widths, the map is inverted globally on its explicitly specified discrete chart, and the separate arithmetic exclusions and prime witnesses are literal. Neither chart proves coverage of arbitrary primes or of the $7p$ lane.
