# Exact self-converse gaps and an automorphism sieve on count fibers

## An all-order low-count classification

The twelve holes among positive odd counts below 75 in the **self-converse**
sector are
$$
7,19,21,23,35,39,43,47,53,55,59,63.
$$
This is an all-order statement for that target range, not just a finite
order observation. It does not exclude the ordinary 19, 23, or 63.

To see why a finite certificate suffices, write a self-converse tournament
as the unique ordinal sum $L\to C\to L^*$, with $C$ strong or empty.
Then $H=H(L)^2H(C)$. By [Busch's sharp strong bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3)
a strong tournament of order at least 9 has at least 75 paths.
Thus below 75 the middle factor is found by the complete paired-site
census through order 8, with the empty factor equal to 1. The only
possible positive values of $H(L)$ below $\sqrt{75}$ are 1, 3, 5:
they are odd, and the ordinary exclusion of 7 is known. All three are
attained by the singleton, triangle, and strong order-4 tournament.
Multiply the paired-site count set by 1, 9, and 25 and truncate at 75.
This gives precisely the displayed list. In fact the two nonidentity
multipliers add no value not already attained at order at most 8.
The companion replay constructs each pair-site assignment and computes
its Hamiltonian paths by endpoint subset DP, with no quotient
classification needed to decide support. It also checks an ordinary
23-path example with outgoing masks $(48,1,3,7,14,30)$.
The base 75 is a strict cutoff; this is not a claim about 75 itself.

The contrast with the dyadic construction is exact: every listed hole
has a strong self-converse representative in *each* of its finite dyadic
residue neighborhoods, without being an exact count. The original
7–21–63 chain thus has three, not two, holes in this sector.

## Primes in joint coefficients

Let $\rho$ be a fixed anti-involution and $F$ the number of paths fixed
by reversal composed with $\rho$. The centralizer $C_{\operatorname{Aut}T}(\rho)$
acts freely on both the $H$ paths and the $F$ fixed paths. Hence in the
joint decorated-site Burnside coefficient indexed by $(H,F)$, site
permutations with order not dividing $\gcd(H,F)$ contribute zero. This
follows because the order of the permutation equals the order of its
canonical lift that preserves tail and head labels. $F$ is odd by
pairing against nonfixed paths, so this sieve never needs a convention
for $F=0$. For ordinary $H$-stratified Burnside, allowable permutation
orders must divide $H$. This is a coefficient restriction, not a
claim that the surviving terms can be counted using only $H,F$.

In fact $F=H$ occurs only for transitive tournaments. If a self-converse
strong core has at least three vertices, take a Hamilton cycle. Cutting
it at any of its $n$ seams gives $n$ distinct Hamilton paths.
If one is fixed, $\rho$ acts as a reflection on the cycle, and at most
two seams are fixed (a reflection of an $n$-gon). Thus some path
is not fixed. The ordinal decomposition gives
$F=H(L)F(C)$ and $H=H(L)^2H(C)$, so equality requires
$H(L)=1$ and a transitive middle. A tournament with $H=1$ is
transitive: each of the three paths of any directed triangle could
otherwise be extended by successive insertions. This also shows that
every prime-count self-converse chart has trivial *marked*
stabilizer: $0<F<H$ and $\gcd(H,F)=1$.

An additional filter does not require converse symmetry. If a
prime $p$ divides $|\operatorname{Aut}T|$, take an order-$p$
automorphism with $r\ge1$ nontrivial orbits. Each such orbit is
a circulant tournament on $p$ vertices. Its $(p-1)/2$ allowed
differences each give $p$ Hamilton paths by walking in arithmetic
progression, and these paths are distinct. Thus
$$
p\mid H,\qquad H\ge [p(p-1)/2]^r.
$$
To justify multiplying the bounds without any module assumption,
two arbitrary disjoint directed paths can be merged while preserving
both internal orders: at each step put the loser of the two current
last vertices in the final free position; both candidates at the
next step beat it. The result is recoverable by restriction to each
vertex set, so each tuple gives a distinct merged path. Merge one
path per orbit plus one per fixed singleton. This recovers the
inherited rigidity of prime counts $p\ge5$, and filters the large
prime factor from the symmetry palette of $H=7p$ as soon as $p>15$.
Neither restriction constructs the remaining asymmetric target cores.

## One pointed cycle recurrence for both A000568 and the paired sector

Give a permutation type $\lambda$ with odd cycle lengths a weight
$a^{\ell(\lambda)}b^{o(\lambda)}\mathbf x^{\mathbf m}/z_\lambda$;
here $o(\lambda)$ counts unordered-pair orbits, as in
[the paired-site note](PAIRED_SITES.md). Let $P_n^{a,b}(\mathbf x)$
be the sum at size $n$ and $P_0=1$. Mark the cycle containing
a distinguished label. Deleting a marked $d$-cycle removes
$(d-1)/2$ pair orbits inside it and $\gcd(d,e)$ per remaining
$e$-cycle. In the partition weight the marking cancels the
multiplicity and length of that cycle: $d(m_d+1)$. Thus
$$
 nP_n^{a,b}(\mathbf x)=
 a\sum_{\substack{d\le n\\d\ {\rm odd}}}
 x_d b^{(d-1)/2}
 P_{n-d}^{a,b}\bigl((b^{\gcd(d,e)}x_e)_{e\ {\rm odd}}\bigr).
$$
At all $x_e=1$, $P_n^{1,2}$ gives A000568 and
$P_m^{2^\epsilon,4}$ gives the self-converse count at order
$2m+\epsilon$. Consequently the merged node count at size $n$
is $\big(P_n^{1,2}(\mathbf 1)+P_{\lfloor n/2\rfloor}^{2^{n\bmod2},4}(\mathbf 1)\big)/2$.
Restricting $d$ to an allowed prime palette computes the corresponding
truncation without changing the recurrence. This is an exact shared
recursive algebra, not a bounded-state scalar recurrence: future
steps apply different gcd substitutions to the full cycle inventory.

To retain Hamilton counts across SCC composition, use formal
Dirichlet weight $h^{-s}$ so products of counts add no information
loss. The ordinary series $\mathcal T(z;s)$ and strong series
$\mathcal S(z;s)$ obey $\mathcal T=(1-\mathcal S)^{-1}$
by unique ordinal SCC words. The self-converse series satisfies
$$
 \mathcal B(z;s)=\mathcal T(z^2;2s)(1+\mathcal S_{\rm sc}(z;s)).
$$
This is a coefficientwise exact projection and recovers the count
in the merged sector by $(\mathcal T+\mathcal B)/2$. It isolates
the remaining exact spectrum question in the strong cores, rather
than mistaking an SCC identity for a method to generate those cores.
