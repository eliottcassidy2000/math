# The 8–9 square-sum endpoints really do carry a centered Collatz branch

The distinction is between a literal raw Collatz edge (really labeled by the two successive odd integers) and a whole square-sum subpath carrying that edge after a phase-dependent change of coordinates. Here the former has a complete obstruction, but the latter produces arbitrarily long exact walks. All positive constructions below are on explicit vertex sets; except for the first unscaled block, the sets are sparse, not intervals.

## 1. Raw square-sum edges are isolated on every accelerated Collatz trajectory

Let $x$ be positive odd, let $j=v_2(3x+1)$, and let $x'=(3x+1)/2^j$ be its next odd state. Set $A=2^j+3$. Then
\[
 x+x'=m^2 \quad\Longleftrightarrow\quad A x+1=2^j m^2 .
\]
No odd guard $j$ can satisfy this. Indeed $2^j m^2\equiv1\pmod A$, so $m$ is a unit modulo $A$. If $j$ is odd the Jacobi symbol $(2/A)$ would have to be $+1$. But $A\equiv5\pmod8$ at $j=1$ and $A\equiv3\pmod8$ at all larger odd $j$, giving $(2/A)=-1$ in both cases.

Every even guard $j$ has infinitely many exact examples. Choose a positive even $m\equiv 2^{-j/2}\pmod A$ and put
\[
 x=\frac{2^j m^2-1}{A},\qquad x'=\frac{3m^2+1}{A}.
\]
Both quotients are integral, positive, and odd, and $3x+1=2^j x'$ proves the guard is exact. Since $m$ is even and $A\equiv3\pmod4$, every such edge has $x\equiv1\pmod4$ and $x'\equiv3\pmod4$. Conversely any square edge between odd states has even $m$, so these residue conditions apply to all raw square-sum Collatz edges. Consequently no accelerated Collatz trajectory contains two consecutive square-sum edges, even with changing guards. Replacing a two-step branch by two direct square edges cannot represent the $9/8$ scale law.

## 2. Why the first Pythagorean block is uniquely relevant

Use the endpoint formula of [THM-4510](../../../01-canon/theorems/THM-4510-sum-graphs-as-reflection-orbits-pythagorean-zigzags.md) for the primitive Pythagorean triple $s^2+t^2=u^2$, $s<t$, and its three-reflection chain at order $t^2-1$. If $t$ is odd, the endpoint ratio or its inverse is $2$ or $1/2$, which has no positive power of 3. If $t$ is even, its endpoints are $s^2$ and $t^2/2$. Suppose their ratio in one orientation is $3^a/2^b$ with $a,b>0$. In the orientation $2s^2/t^2$, primitivity gives $s=3^m$, $t=2^n$, $m,n>0$. In the other orientation it forces the odd leg $s=1$, which cannot be a leg of a nontrivial integer right triangle.

In the viable orientation,
\[
 (u-3^m)(u+3^m)=2^{2n}.
\]
The two factors have gcd 2, and therefore are 2 and $2^{2n-1}$. It follows that
\[
 3^m=2^{2n-2}-1=(2^{n-1}-1)(2^{n-1}+1).
\]
For $n=1$ this is impossible. For $n\ge2$ the two odd factors are coprime, so the smaller equals 1 and $(m,n)=(1,2)$. Thus the $3$–$4$–$5$ triple is uniquely the primitive Pythagorean block of this family with endpoint ratio a nontrivial Collatz slope. The ratio is exactly $9/8$, with exponents $a=2,b=3$. This classification is only for the stated three-reflection family; it does not assert uniqueness among all square-sum paths.

## 3. Two exact guards at each joint of a path of 8–9 blocks

The two accelerated guard words $(1,2)$ and $(2,1)$ have affine maps $(9x+c)/8$, with $c=5$ and $c=7$, respectively. Their exact starting residue classes are $x\equiv-c\pmod{16}$. In the centered coordinate $y=x+c$, each map is just $y\mapsto9y/8$. For $y=8\lambda$ with $\lambda$ even, the two actual trajectories are
\[
 8\lambda-5\ \to\ 12\lambda-7\ \to\ 9\lambda-5,\qquad
 8\lambda-7\ \to\ 6\lambda-5\ \to\ 9\lambda-7 .
\]
Checking $3x+1$ against $2^j x'$ and oddness verifies each stated guard, not merely total division by 8.

Fix the square-sum path
\[
 (8,1,15,10,6,3,13,12,4,5,11,14,2,7,9)
\]
on $[15]$. There is a literal middle-out identity: its middle vertex is $12$, exactly the centered intermediate state between $8$ and $9$ for guards $(1,2)$. At any scale $\lambda$, the two guards are displayed by equal seven-edge halves $8\lambda\leadsto12\lambda\leadsto9\lambda$, using phase centers $5,7,5$. For guards $(2,1)$ the corresponding internal checkpoint is $6$, at path index 4 with zero-based indexing, giving a $4/10$ split.

The alternating vertices beginning at 8 are $(8,15,6,13,4,11,2,9)$. A two-edge zigzag is a rotation by $7$ modulo 16, because the first reflection uses target 9 or 25 and the second target 16. Explicitly a low alternate vertex $z\le8$ goes to $z+7$, a high one $z\ge10$ goes to $z-9$, and 9 is the terminating alternative. Seven rotations send $8$ to $9$. This is the reflection mechanism realizing the scale ratio in this particular path; a fixed composition of reflections still has slope $\pm1$ on the number line and is not the global Collatz map.

For any desired number $k$ of blocks, take
\[
 \lambda_i=2^{6k-3i}3^{2i},\qquad 0\le i<k.
\]
Each $\lambda_i$ is even, and $\lambda_{i+1}/\lambda_i=9/8$. Adjacent sets $\lambda_i[15]$ and $\lambda_{i+1}[15]$ meet only at $9\lambda_i=8\lambda_{i+1}$: an equality requires $8a=9b$, hence $a=9,b=8$ since both indices are at most 15. Sets two or more indices apart are disjoint because the analogous equality requires an index at least $9^2>15$. Concatenating all scaled paths at their common endpoints therefore gives a genuine spanning path on their union, with $14k+1$ vertices. At the joints, subtracting $c$ recovers the odd states reached by exact consecutive copies of the same two-guard branch.

The edge sums alternate by block between squares and twice squares, because $\lambda_i$ alternates between those two square classes. This alternation is necessary for this specific one-block homothetic scheme: if both neighboring copies had only square sums, both scales would be rational squares, which is incompatible with a ratio of $9/8$. Translating every vertex by $-c$ changes the square/twice-square tests to centered sums $x+y+2c$; it does not give raw square sums.

## 4. Every individual guard checkpoint on a pure-square path

A two-phase coordinate change removes the affine constants from the individual transitions. For guard 1 take incoming center 5 and outgoing center 7; for guard 2 take incoming center 7 and outgoing center 5. In each case,
\[
 y=x+c_{\rm in},\qquad y'=x'+c_{\rm out},\qquad
 y'=\frac3{2^j}y.
\]
Thus on an alternating sequence the center is determined by the phase. After four guards the net scale is $81/64=(9/8)^2$, now a rational square. The following literal pure-square paths contain all centered checkpoints in order, with exactly four edges between successive checkpoints:

| Guard word | Centered checkpoints at scale 1 | Pure-square spanning path on its displayed support |
|---|---|---|
| $(1,2,1,2)$ | $64,96,72,108,81$ | $(64,17,8,73,96,4,21,28,72,49,51,13,108,36,45,19,81)$ |
| $(2,1,2,1)$ | $64,48,72,54,81$ | $(64,17,8,1,48,33,3,97,72,49,15,10,54,27,9,40,81)$ |

Take one row and repeat it at scale $\mu_i=2^{6(k-i)}3^{4i}$ for $0\le i<k$. All scales are even perfect squares and $\mu_{i+1}/\mu_i=81/64$. The four guards at every checkpoint are exact: the first and third, second and fourth form the two branches in §3 with parameters $8\mu_i$ and $9\mu_i$, both even. All adjacent sums of the displayed paths are squares and stay squares at these scales. The adjacent supports meet only at $81\mu_i=64\mu_{i+1}$ because $64a=81b$ and all labels are below $2\cdot81$. Nonadjacent supports are disjoint because an index would be at least $81^2$. So again the concatenated path is spanning on the union, with $16k+1$ vertices, and it displays every accelerated step of $4k$ consecutive steps at uniform four-edge distance. Extra edges between support vertices are allowed: this is a spanning path, not a claim that the induced graph is itself a path.

There is a shorter pure-square representative if one only needs the two-step checkpoints: $(64,36,28,72,9,40,81)$, with sums $(100,64,100,81,49,121)$. It has $6k+1$ vertices after gluing. The scalar $9/8$ relation therefore lifts in three different senses: one scaled first-onset block per two-step branch if twice squares are allowed; pure squares with a checkpoint per two steps; or pure squares with every individual guard checkpoint after alternating the affine centers. None of these claims Collatz convergence, a square-sum result on a full large interval, or transport to ordinary tournament Hamilton counts.

## The middle on the general even-leg reflection path

Let \(a=s^2,b=t^2\) in the primitive Pythagorean construction of [the earlier reflection-orbit note](../procgen_sumgraph_20260926_reflection_orbits.md), with \(s<t\) and \(t\) even. On \([b-1]\) its three targets are \(a,b,a+b\). Start at \(e_0=b/2\). Alternate a reflection at \(a\) or \(a+b\) (their available domains are disjoint) with a reflection at \(b\). With path indices counted from zero, the two sheets satisfy
\[
   e_j\equiv b/2-a j,\qquad o_j\equiv a-e_j\pmod b .
\]
Here there are \(b/2\) even-indexed labels and \(b/2-1\) odd-indexed labels, with endpoints \(b/2,a\). Reversing the path pairs the \(e\)-labels by sum \(a+b/2\bmod b\) and the \(o\)-labels by sum \(b/2\bmod b\); this follows by substituting \(b/2-1-j\) and \(b/2-2-j\) respectively. Since \(a\equiv1\bmod4\), the exact central label is \(3b/4\), with neighbors \(b/4\) and \(a+b/4\bmod b\) as positive representatives. At \(Q_{15}\) the center is 12 and the neighbors are 4 and 13. The two sheets really matter: the extra square-sum edge \(1\!-\!3\) maps under this path mirror to \(7\!-\!5\), whose sum is 12. It is not a full-graph automorphism.

This exact endpoint mechanism is distinct from the [golden-field clock](../five_eight_nine_transfer_20261005.md), where 9 is a field size and 8 a unit period, and the [cubic square-sum contact chains](../square_sum_cubic_contact_chains_20261019.md), where \(9/8\) bounds **normalized cubed-height growth** and is not an endpoint ratio. The new transport above uses the Pythagorean path, not either unrelated occurrence. Its middle-out lesson for the tournament selector is to keep both ports and their interaction, not to use the midpoint as a replacement for the two readings.

## One common pure-square path for both gap-preserving branches

The two paths above can in fact share **one and the same** pure-square carrier. Use the following ordered sparse support:
\[
 (64,105,16,48,73,96,100,21,28,\ \mathbf{72},\
                 9,27,54,46,3,61,108,88,81).
\]
Every consecutive sum is a square. The middle vertex 72 divides the path into two equal nine-edge halves, and the left half contains **both** interior ports 48 and 96 while the right half contains both ports 54 and 108. Thus the single carrier displays the two different chronological checkpoint subsequences
\[
 (64,48,72,54,81),\qquad (64,96,72,108,81).
\]
For an even square scale \(\mu\), these are not merely a pair of formal words: with alternating phase centers \(7,5,7,5,7\) and \(5,7,5,7,5\), respectively, they recover the ordinary integer trajectories of the input pair \(64\mu-7,64\mu-5\), down to each exact division guard. After two guards the pair is \(72\mu-7,72\mu-5\); after four it is \(81\mu-7,81\mu-5\). The gap of two is preserved at the shared centered checkpoints, not generally at the individual intervening odd states.

Repeat the *same* path at scales \(\mu_i=2^{6(k-i)}3^{4i}\), \(0\le i<k\). The sparse supports glue at \(81\mu_i=64\mu_{i+1}\), and only there: all labels are at most 108, so the overlap calculation of the preceding section applies exactly. This gives a single pure-square spanning path on \(18k+1\) labels that represents **both** trajectories of the positive Collatz unit pair through \(4k\) consecutive guards. The first four-step block preserves the same common support while retaining distinct interior ports. This is the exact square-sum version of the middle-out organizing idea, not a map from its vertices to an apex contact in a tournament.

A slightly longer carrier gives an exact positional mirror between the two guard words:
\[
 (64,105,16,33,48,52,69,100,96,73,27,9,\mathbf{72},
 97,3,46,54,10,39,61,108,36,45,19,81).
\]
The seven ordered ports \(64,48,96,72,54,108,81\) are exactly four path edges apart at indices \(0,4,8,12,16,20,24\). Path reversal exchanges the **positions** of the two branch checkpoints in reverse order. In particular the two lifted guards use alternating subpath lengths \(4,8,4,8\) and \(8,4,8,4\), not one scalar length pretending to encode both. It glues at the same square scales, on \(24k+1\) labels. The mirror is on the marked path positions; no claim is made that it is an affine reflection of the numerical labels or a symmetry of the whole graph.
