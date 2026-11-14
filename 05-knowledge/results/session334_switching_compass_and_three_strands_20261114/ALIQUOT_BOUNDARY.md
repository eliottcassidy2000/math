# Aliquot boundary of the smooth target lattice

In the offset-one chart write the **increment**, not the path count, as $m=H-1=A B$, with $A=2^a$, $B=3^b$, $a\ge3,\ b\ge1$. The prime-power rectangle of divisors factors the divisor sum, giving the exact absolute identity
$$
 s(m)=2m-A-\frac{3B-1}{2}.
$$
The two costs of moving from the bulk $2m$ to the aliquot successor depend on the two axes separately. This is a second, arithmetic occurrence of the same **anchored sum of two boundary prices** as the compass chamber $C-\alpha X-\beta L/X$; the two formulas have different geometry. In the compass the moving powers have fixed product, so the boundary charges trade off. In the divisor rectangle both powers can independently grow. We do not identify these moves on actual tournaments.

## A necessary resonance for staying in the rail

Almost every exponent pair in this chart leaves the $2,3$-smooth monoid after the aliquot step. More precisely, if $s(m)$ is $2,3$-smooth, it is necessary that at least one of
$$
 \boxed{a=v_2((3^{b+1}-1)/2)}
 \qquad\text{or}\qquad
 \boxed{b=v_3(2^{a+1}-1)}
$$
hold. This is a boundary-cancellation test, not a statement that every resonance works.

To prove it, in
$$
 s(m)=A(2B-1)-(3B-1)/2,\qquad
 2s(m)=(4A-3)B-(2A-1),
$$
the factors $2B-1$ and $4A-3$ are units at 2 and 3 respectively. With unequal valuations at the two subtraction operands, the valuation of the difference is their minimum. If neither resonance holds, $v_2(s(m))\le a$ and $v_3(s(m))\le b$; smoothness would force $s(m)\le AB$. But
$$
 s(m)-AB=(A-\tfrac32)B-A+\tfrac12>0
$$
already for $A=8,B=3$, and increases in both axes. Contradiction.

The first resonance is impossible for even $b$; for odd $b$ it is $a=1+v_2(b+1)$. The second is impossible for even $a$; for odd $a$ it is $b=1+v_3((a+1)/2)$. These follow by factoring the odd-exponent geometric sum and repeatedly doubling at 2 or cubing at 3. Among $3\le a\le U,\ 1\le b\le V$ there are at most $U+V$ candidates. In fact each resonance chooses one $a$ for a given $b$ or vice versa. Recasting the test in target coordinates, almost every $1+2^a3^b$ leaves the offset-one smooth rail under $H\mapsto1+s(H-1)$. No assertion about leaving the entire tournament spectrum follows.

Both an early success and a nontrivial control lie on these resonance faces: $24\mapsto36$, $216\mapsto384$, whereas $72\mapsto123=3\cdot41$. Thus the first macro rail point $217=1+36\cdot2\cdot3$ has aliquot-increment neighbor $385=1+2^7\cdot3$, already furnished by the low ternary sector in Session333. It is the aliquot sum of **216**, not of 217, that realizes this precise comparison. Persistence away from the resonance faces cannot be obtained just by repeating the same smooth rail rule.

## Historical reading

The old [Vieta aliquot carrier](../../../07-reflections/vieta-perfect-aliquot-carriers-s644.md) isolates the same issue of comparing different observers: on $\binom{14}{2}=91$, the square-discriminant observer sees $27$, the aliquot observer sees $21$. No nonrealization theorem about 21 follows. The useful transfer is to retain which observer is being applied, and to look for exact cancellation where the two terms share a valuation. The [older summand closure](../../../07-reflections/summand-graph-fermat-zeckendorf.md) likewise preserves the actual parent-pair guard rather than asserting a number can be generated because one of its incoming arcs is available.

## The transitivity comparison, with the parent pair retained

For distinct positive summands, the projected directed summand graph has edge $x\to z$ iff $x<z$ and $z\ne2x$. So it is a transitive tournament **minus directed doubling chains**, rather than itself transitive: $3\to4\to6$ is legal but $3\to6$ is not. For $x\ge2$, the omitted relation $x\to2x$ is nevertheless forced in the transitive closure by $x\to x+1\to2x$. Thus on a finite initial segment of size at least two exactly two transitive tournament completions are possible: orient $1,2$ either way, and all other pairs have forced order. This shows concretely what is lost by taking the closure: all the forbidden equal-parent events except the first disappear.

The unordered parent pairs of a new node $z$ are parametrized by
$$
 (u,z-u),\quad 1\le u\le \lfloor(z-1)/2\rfloor .
$$
Color by $u\pmod3$. Advancing $z$ by 6 grows each color strand by exactly one parent pair, on top of the translated older pairs $ (u,z-u)\mapsto(u,z+6-u)$. This is a three-lane growth law, not the assertion that the residue of a parent pair identifies an exact child independent of its row. In fact unordered *residue pairs alone* give at most two pair types per target residue, since $u\mapsto z-u$ on three residue classes has one fixed class and exchanges the other two. Labelling the smaller parent is therefore a genuine coordinate choice. The marked-fibre counterpart is to name the three runs before tracking the assignment of a site; forgetting the name preserves the total but loses the charge used in the three-strand theorem.
