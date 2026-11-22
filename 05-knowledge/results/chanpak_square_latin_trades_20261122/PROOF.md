# Five square sums and a binomial family of Hamilton paths

**PROVED.** At infinitely many orders $N$, five specified square sums support at least $\exp(cN^{2/5})$ distinct undirected Hamilton paths on $1,\ldots,N$, for an absolute $c>0$. Here the five sums depend on $N$. The proof also counts exactly the paths with prescribed endpoints and a prescribed matching on these orders. The count is a binomial sum with a single coprimality test.

## Starting point

The 3–4–5 triple gives a path from 8 to 9 on 1 through 15 with sums 9, 16, 25. Its inherited generalization uses a primitive triple to produce a path by composing two modular reflections. The scaled version classifies what happens when the resulting translation has several orbits. We instead take **two different primitive triples at the same even leg**, keep the common reflection and exchange the two other reflections separately on each alternating cycle. This brings three formerly separate pieces of information together. They are the two Berggren addresses, the cyclic addition table at the common leg, and the number of connected components after the exchanges. [THM-4510, reflection orbits and Pythagorean zigzags](../../../01-canon/theorems/THM-4510-sum-graphs-as-reflection-orbits-pythagorean-zigzags.md) · [Scaled three-square chains](../chanpak_berggren_square_bridge_20261117/squares/PROOF.md) · [Romik, *The dynamics of Pythagorean triples*](https://arxiv.org/abs/math/0406512)

## Theorem 1. An exact count with a prescribed matching

Let $B$ be even and let $0<A<C<B$ be units modulo $B$. Write

    g = gcd(B,C−A),    K = g/2,    M = B/g,    D = (C−A)/g.

In the graph on $[B]=\{1,\ldots,B\}$ admitting only the five sums $B,A,A+B,C,C+B$, require endpoints $B,B/2$ and require every edge of sum $B$. The exact number of undirected Hamilton paths with these properties is

    Σ binom(K,ℓ),     over 0 ≤ ℓ ≤ K with gcd(A+2ℓD,M)=1.

Every choice on layer $\ell$, including a failed choice, makes one path on $B/d$ vertices and $(d-1)/2$ disjoint cycles on $2B/d$ vertices each, where $d=\gcd(A+2\ell D,M)$. Removing vertex $B$ shortens just the path by one. Connected choices remain distinct after this deletion.

**Proof.** Identify the number $B$ with zero in $\mathbb Z/B$. With loops omitted, $t(x)=-x$ supplies exactly the required matching. It has two fixed points $0,B/2$. Each of $r_A(x)=A-x$ and $r_C(x)=C-x$ is a fixed-point-free involution. Its edges have precisely the allowed ordinary sums of the corresponding residue. This includes the edge from $B$ to $A$, of ordinary sum $B+A$, and its analogue for $C$.

A path with the specified matching and endpoints uses exactly one edge of $r_A\cup r_C$ at every vertex. This union consists of even alternating cycles. The product of the two reflections translates by $C-A$, whose orbits are the individual residue classes modulo $g$. A reflection pairs the classes $i,A-i$. These classes are distinct because $g$ is even and $A$ odd. Thus there are $K$ alternating cycles, each with exactly two choices of perfect matching. Choose the $C$-matching on $\ell$ cycles and the $A$-matching on the rest. Denote the resulting involution by $R$.

Modulo $g$, the permutation $F=tR$ subtracts $A$, so it visits every class once in each $g$ steps. In $\mathbb Z/B$, the step subtracts $C$ on the $2\ell$ chosen classes and $A$ otherwise. Its first return to each class is therefore translation by

    −[(g−2ℓ)A+2ℓC] = −g(A+2ℓD).

Consequently $F$ has exactly $d$ orbits, each of size $B/d$. The union of the loopless edges of $t,R$ has degree two except at $0,B/2$, where it has degree one. It consists of a path between these vertices and alternating cycles. On the path, $tR$ visits its vertices in a single orbit. In path order it goes up by two on one parity and comes back by two on the other, turning at the endpoints. On each cycle it has two orbits of equal size. This proves the claimed lengths and count. On removing $B$ from a connected choice, the other endpoint besides $B/2$ is precisely the former neighbor of $B$, either $A$ or $C$. Reattaching $B$ there recovers the choice. ∎

**Latin-square interpretation.** In the cyclic Latin square $L(x,y)=x+y\pmod B$, exchange symbols $A,C$ in rows whose residues modulo $g$ lie in a chosen set of pairs $\{i,A-i\}$. Each changed row still uses every symbol once. The two affected rows in any column belong to the same class modulo $g$, so both exchange their entries together and all column constraints survive. Closure under $i\mapsto A-i$ makes the new table symmetric. The new positions of $A$ give exactly $R$. Thus a genuine two-symbol Latin trade preserves row, column and symmetry constraints for **every** choice. The extra coprimality test decides which choices preserve path connectivity. Paths using additional edges or other endpoints can occur in the full graph.

## A concrete collision of Berggren addresses

Both $(41,840,841)$ and $(559,840,1009)$ are primitive. Put $B=705600$, $A=1681$, $C=312481$. Then $g=8400$, $K=4200$, $M=84$, $D=37$. The sums are $1681,312481,705600,707281,1018081$, all squares. Every one of the $\binom{4200}{2100}>2^{4187}$ central choices gives a Hamilton path at both orders 705600 and 705599. At the larger order all paths have endpoints 352800 and 705600. Every individual one-cycle switch instead gives $d=\gcd(1755,84)=3$, one path and one cycle. The arithmetic change is invisible to a count of matching incidences at individual vertices.

## Theorem 2. Infinitely many five-square palettes

For all sufficiently large odd positive integers $k$, put

    u = 2k+1,    v = 3k+1,    s = v⁴−u⁴,    r = s+u²v²,
    b = 2uvrs,   B = b².

On each of $[B]$ and $[B-1]$ there are at least $\exp(cB^{2/5})$ distinct Hamilton paths using only five square sums. On $[B]$ they all have the same two endpoints and contain every edge of sum $B$.

**Proof.** The four integers $u,v,r,s$ are pairwise coprime. Indeed $\gcd(u,v)=1$ by their definitions, $s=v^4-u^4$ is coprime to $uv$, and $r=s+u^2v^2$ is coprime both to $uv$ and to $s$. Only $v$ is even. The two Euclid parameter pairs $(us,vr)$ and $(vs,ur)$ give primitive triples with common leg $b$. Their other legs are the absolute differences of the respective parameter squares. To see these positive legs are strictly smaller than $b$, note $v/u\to3/2$ and $r/s\to101/65$. The ratios of the larger to smaller parameter in the first and second pairs therefore approach $303/130$ and $202/195$, respectively. Both limits lie strictly between 1 and $1+\sqrt2$, the condition for the difference of squares to be smaller than twice the product.

Call their squared legs $A,C$ in increasing order. Both are odd units modulo $B$, and $B+A,B+C$ are the squared hypotenuses. The difference of squared legs is

    C−A = (v⁴−u⁴)(r⁴−s⁴).

It is positive. Since $r\equiv s\pmod{u^2v^2}$, this gives $g\ge s u^2v^2$, as $s u^2v^2$ divides both $B$ and the difference. Also $8\mid g$ because both odd legs have square 1 modulo 8, and $b$ is divisible by 4. Thus the central choice $\ell=K/2=g/4$ is integral and its translation parameter is $(A+C)/2$.

We show that $(A+C)/2$ is coprime even to $B$. It is odd. Expanding gives

    A+C = (u⁴+v⁴)(r⁴+s⁴)−4u²v²r²s².

For an odd prime $p\mid uv$ we have $r\equiv s\pmod p$, so the right hand side is a unit modulo $p$. For $p\mid s$ it is a unit since $v^4\equiv u^4\pmod p$. For $p\mid r$, vanishing would require $v^4=-u^4\pmod p$. Substituting in $r=v^4-u^4+u^2v^2\equiv0$ gives $(v/u)^2=2\pmod p$, hence $4=-1\pmod p$. But 2 is a nonsquare modulo 5. This exhausts the primes of $B$.

By Theorem 1 all central choices yield distinct paths. Their number is $\binom K{K/2}\ge 2^K/(K+1)$. Finally $u,v\asymp k$, $r,s\asymp k^4$, $B\asymp k^{20}$ and $K\gg k^8\gg B^{2/5}$. These estimates prove the bound on an unbounded sequence of distinct orders. ∎

## Exponent, dimension and the tournament boundary

The permutation theorem only uses two units and a cyclic group. Squares enter in making each lifted ordinary sum admissible. Requiring the same additive closure $a^e+b^e=c^e$ for a palette of positive $e$-th powers gives the classical Fermat obstruction for every integer $e>2$. The projective curve remains one-dimensional as $e$ changes. Its degree and genus change. Passing to a larger cyclic group or allowing more target sums are different operations. The independent complex Jacobian dimension boundary comes from retaining a noninjective constant-Jacobian map under identity padding. Neither that map nor its padding supplies an operation on the Hamilton path count. [Wiles, *Modular elliptic curves and Fermat’s last theorem*](https://doi.org/10.2307/2118559) · [Fermat curve geometry](https://mcs.unife.it/alex.massarenti/files/FLT.pdf) · [Isabelle verification of the Jacobian counterexample and its padding](https://isa-afp.org/entries/Jacobian_Counterexample.html)

For tournaments, two existing facts delimit the application. A rooted tournament construction already uses disjoint inverse chambers to decode ternary branch words. And the direct proposal $H(T)=2p(Q_n)+1$, for $p$ the unrestricted square-sum path count, fails at $Q_{23},Q_{25}$. Here the actual shared operation is exchange in a fixed interface. The count of a *successful* exchange requires the first-return permutation. Its reduction to a single integer is proved above because the quotient is a cyclic translation. This argument gives a large family of paths **inside a graph**, with no implication that these cardinalities occur as Hamilton counts **of tournaments**. [Rooted tournament tree](../chanpak_berggren_square_bridge_20261117/root/PROOF.md) · [Session 211, scalar obstruction and path-cover ports](../../../session211/square-sum-port-bridge/PROOF.md)

The accompanying check compares the component formula with direct graph traversal on every even modulus through 32 for every pair of admissible units and every switch choice. It also traverses three specified choices of the five-square example and checks the parameter identities for all odd $3\le k\le 101$. Run `python3 05-knowledge/results/chanpak_square_latin_trades_20261122/check.py` from the repository root.
