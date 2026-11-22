# Five square sums and a binomial family of Hamilton paths

**PROVED.** At infinitely many orders $N$, five specified square sums support at least $\exp(cN^{5/8})$ distinct undirected Hamilton paths on $1,\ldots,N$, for an absolute $c>0$. Here the five sums depend on $N$. The proof also counts exactly the paths with prescribed endpoints and a prescribed matching on these orders. The count is a binomial sum with a single coprimality test.

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

Take a primitive Pythagorean triple $u^2+v^2=s^2$ with $u$ odd and $v$ even. Put $r=u+v$, $b=2uvrs$, and $B=b^2$. Define the two other legs by Euclid's parameter pairs $(us,vr)$ and $(vs,ur)$. Whenever both other legs are positive and smaller than $b$, the five corresponding squares support at least
$\binom{2uvrs^2}{uvrs^2}$ distinct Hamilton paths on each of $[B]$ and $[B-1]$. On $[B]$ they have common endpoints $B,B/2$ and contain every edge of sum $B$.

In particular the original **3–4–5 triangle itself** gives the order-705600 example above. An infinite sequence of primitive triples approaching the shape 3–4–5 gives at least $\exp(cN^{5/8})$ paths at each of the two successive orders $N=B-1,B$.

**Proof.** The four integers $u,v,r,s$ are pairwise coprime. Only $v$ is even. The only new coprimality to check is that of $r,s$. An odd prime dividing both would give $v=-u$ and $0=s^2=2u^2$ modulo that prime, which is impossible. Thus both parameter pairs give primitive triples with common even leg $b$. Write $A,C$ for their squared other legs in increasing order. In the asserted chamber $0<A<C<B$, and $B+A,B+C$ are squared hypotenuses. Using $r>s$, the difference of squares factors as

    C−A = |v⁴−u⁴|(r⁴−s⁴)
        = |v−u|·r·s²·4uv(u²+uv+v²).

Here $u$ and $v$ are distinct by parity. We claim that the switch parameters are exactly

    g = 4uvrs²,       M = uvr,       D = |v−u|(u²+uv+v²).

Indeed $B=gM$ and $C-A=gD$. Each prime of $uv$ divides neither factor of $D$. At a prime $p\mid r$ we have $v=-u\pmod p$, so the two factors reduce up to sign to $2u$ and $u^2$. Since $r$ is odd, these are also units. Hence $\gcd(M,D)=1$.

Now $K=g/2=2Ms^2$, so the central binomial layer has $\ell=Ms^2$. Its translation parameter $A+2\ell D$ reduces to $A$ modulo $M$, which is a unit by primitivity. All the central choices give the claimed distinct paths by Theorem 1. More generally every layer with $\ell\equiv0\pmod M$ is connected.

For $u=3,v=4,s=5,r=7$ the pairs are $(15,28)$ and $(20,21)$, so the two other legs are 559 and 41. To obtain distinct unbounded examples, for even positive $k$ use the Euclid pair $m=2k+1,n=k$ and set $u=m^2-n^2,v=2mn,s=m^2+n^2$. Then $v/u\to4/3$ and $r/s\to7/5$. The ratios of larger to smaller parameter for $(us,vr)$ and $(vs,ur)$ approach $28/15$ and $21/20$, respectively. Both limits are strictly between 1 and $1+\sqrt2$, the condition making the difference of squares positive and smaller than twice the product. Thus the chamber holds eventually. Here $u,v,r,s\asymp k^2$, so $B\asymp k^{16}$ and $K\asymp k^{10}\asymp B^{5/8}$. Finally $\binom K{K/2}\ge 2^K/(K+1)$. ∎

## Theorem 3. Connected components of the graph of trades

Use the $K$ alternating cycles of Theorem 1 as named coordinates. The vertices of a cube are choices of matchings, and a single edge changes the chosen matching on one cycle. Keep only vertices whose chosen graph is Hamiltonian.

For each maximal consecutive run $[i,j]$ of layers satisfying the coprimality condition, the resulting graph has the following components. If $i=j$, each of the $\binom K i$ choices is isolated. If $i<j$, it has a single component of size $\sum_{\ell=i}^j\binom K\ell$. Indeed all vertices of layer $i$ are connected by exchanging one element at a time through layer $i+1$. Every vertex in the run can descend to layer $i$ by deletions. A missing whole layer separates consecutive runs because one move changes weight by one.

At $b=840$ the forbidden layers are precisely $\ell\equiv1\pmod3$ or $\ell\equiv5\pmod7$. For example every one of the $\binom{4200}{2106}$ paths at layer 2106 is isolated. Its lower neighbors fail at 7 and its upper neighbors fail at 3.

There are also infinitely many orders with $\exp(cN^{5/8})$ **isolated Hamiltonian states** under these moves. In the construction of Theorem 2 take $n=42j+22$, $m=2n+21$ for positive $j$. They are coprime of opposite parity and $m/n\to2$, so the same limiting shape and growth estimates hold. Modulo 21 the new Euclid parameters $(m,n)$ equal $(2,1)$, giving $(u,v,s,r)=(3,4,5,7)\pmod{21}$. In particular $21\mid M$. Both possible squared legs are $1\bmod21$. Eventually $v>u$, giving $D\equiv37\pmod{21}$ by Theorem 2. Therefore the obstructions 3 and 7 at the indicated neighboring layers persist.

To ensure success of the layer itself at all the other primes, choose by CRT an integer $0\le q<M$ with $q\equiv6\pmod{21}$ and $q\equiv0\pmod p$ for every prime $p\mid M$ other than 3 and 7. Put $\ell=K/2+q=Ms^2+q$. Then $A+2\ell D$ is a unit at 3 and 7 by the original residue check, and at every other prime dividing $M$ it equals the unit $A$ modulo that prime. Layer $\ell$ is good and both neighbors are forbidden. Since $\ell$ is between $K/2$ and $2K/3$ for sufficiently large $j$, it contains at least $2^{\lfloor K/3\rfloor}$ distinct choices: with $h=\min(\ell,K-\ell)$ one has $\binom K\ell=\binom Kh\ge(K/h)^h\ge2^h$. This proves the claim. Isolation here is relative to these single alternating-cycle exchanges with the common matching and endpoints fixed. More general path transformations can introduce other edges.

## Exponent, dimension and the tournament boundary

The permutation theorem only uses two units and a cyclic group. Squares enter in making each lifted ordinary sum admissible. Requiring the same additive closure $a^e+b^e=c^e$ for a palette of positive $e$-th powers gives the classical Fermat obstruction for every integer $e>2$. The projective curve remains one-dimensional as $e$ changes. Its degree and genus change. Passing to a larger cyclic group or allowing more target sums are different operations. The independent complex Jacobian dimension boundary comes from retaining a noninjective constant-Jacobian map under identity padding. Neither that map nor its padding supplies an operation on the Hamilton path count. [Wiles, *Modular elliptic curves and Fermat’s last theorem*](https://doi.org/10.2307/2118559) · [Fermat curve geometry](https://mcs.unife.it/alex.massarenti/files/FLT.pdf) · [Isabelle verification of the Jacobian counterexample and its padding](https://isa-afp.org/entries/Jacobian_Counterexample.html)

For tournaments, two existing facts delimit the application. A rooted tournament construction already uses disjoint inverse chambers to decode ternary branch words. And the direct proposal $H(T)=2p(Q_n)+1$, for $p$ the unrestricted square-sum path count, fails at $Q_{23},Q_{25}$. Here the actual shared operation is exchange in a fixed interface. The count of a *successful* exchange requires the first-return permutation. Its reduction to a single integer is proved above because the quotient is a cyclic translation. This argument gives a large family of paths **inside a graph**, with no implication that these cardinalities occur as Hamilton counts **of tournaments**. [Rooted tournament tree](../chanpak_berggren_square_bridge_20261117/root/PROOF.md) · [Session 211, scalar obstruction and path-cover ports](../../../session211/square-sum-port-bridge/PROOF.md)

The accompanying check compares the component formula with direct graph traversal on every even modulus through 32 for every pair of admissible units and every switch choice. It also traverses four specified choices of the five-square example and checks the parameter identities for all even $4\le k\le 102$. Run `python3 05-knowledge/results/chanpak_square_latin_trades_20261122/check.py` from the repository root. The same check verifies the isolated-layer congruences at the first 30 parameter pairs with $j=10,\ldots,39$.
