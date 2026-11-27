# Guarded machines for Moore, Pollock, Collatz and Cullen

**Terminology update:** The requested comparison was Woodall versus Cullen, not Wilson versus Cullen. See the [signed Cullen–Woodall continuation](../cullen_woodall_signed_collatz_clocks_20261127/README.md). The Wilson calculation below is retained as a separate observation.

A common operation across the four problems is not merely "run a machine." It is to write one more datum while preserving a previously checked constraint. The possibilities differ: Moore routing disallows a common abelian clock, Pollock needs a bound on actual indices, and Cullen divisibility admits a uniquely writable digit at every precision even when the corresponding Collatz power clock stalls.

## 1. Moore: a terminating branch machine, not a halting conjecture

For an undirected diameter-two Moore graph of degree \(d\ge2\), root at a vertex. The \(d\) parents each have \(b=d-1\) leaves. Leaves in one branch have no mutual edges, and each pair of different branches has a perfect matching \(P_{ij}\) (with \(P_{ji}=P_{ij}^{-1}\)). Relabel the slot sets so that \(P_{0j}=I\). The machine's tape records the matchings, current permutation candidate, and backtrack stack. Its legal transition writes one new pair of mutually inverse matchings. For every pair \(i,j\), and every source slot \(x\), reject if any two of the already specified routes
\[
P_{ij}(x),\qquad P_{hj}(P_{ih}(x))\quad(h\ne i,j)
\]
coincide.

At completion, there are exactly \(d-1=b\) routes, so injectivity is sufficient for full reach. The other root, parent and same-branch pairs have the required routes by construction. Thus walking this finite search tree accepts exactly when a Moore graph exists at that degree; exhausting it is a proof of nonexistence. The included code lazily enumerates permutations and stops only at a witness, exhaustion, or an explicit optional budget. It finds degree 2 and 3 directly and independently checks the Hoffman–Singleton bank extracted from the pentagon/pentagram construction. The bounded degree-7 search is not an impossibility test. Degree 57 is a finite decision problem, although nothing here makes its full search practical.

**Why one common cyclic clock cannot program all the routes.** Suppose all normalized \(P_{ij}\) commute. The generated group is transitive, since the routes for one pair move any one slot to all slots. A faithful abelian transitive permutation group is regular: any element fixing one point fixes all other points by commutation. Thus identify the slots with an abelian group \(G\) and the transports with translations \(s_{ij}=-s_{ji}\). Define \(R_i=\sum_{h\ne i}s_{ih}\) and \(\sigma=\sum_{g\in G}g\). Summing the route translations, which enumerate \(G\) exactly once, gives
\[
\sigma=s_{ij}+\sum_{h\ne i,j}(s_{ih}+s_{hj})=R_i-R_j-s_{ij}.
\]
So \(s_{ih}+s_{hj}=R_i-R_j-2\sigma=R_i-R_j\), the same for all intermediates (negation permutes \(G\), so \(2\sigma=0\)). If \(d\ge4\) this repeats a destination. It rules out even an **abelian transitive subgroup centralizing the whole bank**: its centralizer is a regular abelian group when the subgroup is abelian and transitive. The obstruction applies in particular to a common full-length slot cycle. A genuine completion needs noncommuting transport somewhere; the recovered Hoffman–Singleton bank generates all 720 slot permutations. This supplements the [mixed exact-routing result](../mixed_moore_four_state_tiles_20261127/README.md), which permits reciprocal and missing pairs but still keeps diagonal and off-diagonal constraints distinct.

## 2. Pollock: a six-track Boolean tape

For a height bound \(X\), the tape interval is \(0,\ldots,X\). Cell \(n\) carries bits \(b_0,\ldots,b_5\), where \(b_j(n)\) records a representation by exactly \(j\) positive tetrahedra so far. Start with \(b_0(0)=1\) and everything else zero. Generate tetrahedral atoms by the Pascal recurrence
\[
(f,g,h)=(1,3,3),\qquad (f,g,h)\gets(f+g,g+h,h+1);
\]
the first register is \(1,4,10,20,\ldots\). For each \(f\le X\), sweep \(j=1,\ldots,5\) in that order:
\[
b_j(n)\gets b_j(n)\ \vee\ b_{j-1}(n-f).
\]
This is a Turing-machine program on a finite segment. At the local-rule level, a scratch track copies \(b_{j-1}\), shifts right by one cell per tick with zero boundary for \(f\) ticks, then writes the OR to \(b_j\). Serial sweeps with a controller implement the same rule on a one-head tape; the bitset replay performs those shifts in bulk. The ascending sweep is critical: it allows repeated use of the new atom, since \(b_{j-1}\) has already been updated. The alternative descending sweep would test a different problem, distinct tetrahedra. This is an exact decision procedure for each input and bound; it does not assume the conjecture or use modular searches as substitutes for positive indices.

A second machine advances \(X\) and halts if it sees a zero across all six tracks at any positive \(n\le X\). It halts exactly if Pollock has a counterexample. Thus "no halt" here is precisely the open universal statement, not a proof of Turing undecidability. Replay through \(10^6\) matches the repository's earlier frontier, with 241 four-term defects and no five-term defect. The [two-residual inverse](../pollock_chanpak_five_active_and_two_residual_20261127/README.md) gives constructive witnesses in its selected shells without relying on the tape scan.

## 3. Two odd Collatz steps have two distinct clocks

For the halved Collatz map, two successive odd inputs are legal exactly on \(x\equiv-1\bmod4\). Their return is \(x\mapsto(9x+5)/4\), so the centered coordinate \(y=x+1\) is multiplied by \(9/4\). At least \(r\) successive blocks are allowed exactly when \(4^r\mid x+1\); for *maximal* length \(r\) additionally require \(4^{r+1}\nmid x+1\). This follows by losing two factors of 2 at each block. All finite guard depths intersect in the \(2\)-adic state \(x=-1\), not in the positive integers. At an odd prime \(p\ne3\), the same quotient is a unit and its orbit period is an independent multiplicative clock. Neither this orbit nor its period enforces the integer parity guard.

## 4. Cullen digits can be written even when the Collatz clock stalls

Let \(C_a(m)=ma^m+1\), with base \(a\) a positive integer prime to an odd prime \(p\); \(a=2\) is the [Cullen sequence](https://faculty.nps.edu/pstanica/research/cullen.pdf). Let \(d=\operatorname{ord}_p(a)\) and fix a clock phase \(m\equiv j\bmod d\). For \(r\ge0\), adding \(hdp^r\) gives
\[
C_a(m+hdp^r)\equiv C_a(m)+hdp^r a^m\pmod{p^{r+1}}
\]
for nonnegative \(m,h\). Indeed \(a^{dp^r}\equiv1\bmod p^{r+1}\) by iterating the binomial theorem. The relevant digit slope \(da^j\bmod p\) is a unit. Consequently within any phase, any target mod \(p^r\) has exactly one index class mod \(dp^r\). The next digit is obtained by subtraction and multiplication by a fixed inverse, with no branching. Taking target zero indefinitely gives compatible local divisibility certificates but no nonnegative integer index where \(C_a(m)=0\). The clock has not controlled the ordinary size of the representative. This is a concrete counterpart to Pollock's local versus bounded-positive distinction, not a claim about how often Cullen values themselves are prime.

**Wilson versus Cullen, through the Collatz pair.** Put
\[
D_a=\frac{(p-1)a^{p-1}+1}{p}\pmod p,\qquad
q_p(a)=\frac{a^{p-1}-1}{p}\pmod p,
\]
so \(D_a=1-q_p(a)\). In particular \(p^2\mid C_2(p-1)\) asks for \(q_p(2)=1\), not the base-2 Wieferich condition \(q_p(2)=0\). Wilson tests a checksum of **all** bases, not a selected base: Lerch’s formula gives
\[
\sum_{a=1}^{p-1}D_a\equiv -1-W_p\pmod p,\qquad
W_p=\frac{(p-1)!+1}{p}.
\]
For completeness \(\prod_{a=1}^{p-1}a^{p-1}\equiv1+p\sum q_p(a)\equiv1+pW_p\bmod p^2\) proves the needed [Lerch relation](https://arxiv.org/abs/1110.3113).

The two-odd clock is a **difference** of two of those digits:
\[
q_p(9/4)=2(q_p(3)-q_p(2))=2(D_2-D_3)\pmod p
\]
using modular inverses for \(4\). Its order mod \(p^2\) equals its order mod \(p\) exactly when \(D_2=D_3\): \(p-1\) is a prime-to-\(p\) multiple of the mod-\(p\) order. For instance at \(p=23\), the two digits are both 7, and the actual orders of \(9/4\) modulo \(23,23^2,23^3\) are \(11,11,253\). Thus the Collatz power clock stalls for one level, while each individual Cullen digit controller still has a unit slope and writes every next digit. Among \(5\le p<10000\), 23 is the only stall at the first lift in the replay, not a statement about larger primes or integer Collatz trajectories.

Run the JavaScript replay with Node. It uses ordinary integer/BigInt arithmetic, not primality assumptions for any produced Cullen value.
