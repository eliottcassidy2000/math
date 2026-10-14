# Tournament path counts, Fermat quotients, and the missing action marker

**PROVED:** the marked-path orbit interpretation and its symmetry boundary; the elementary quotient identities below. **INHERITED PROVED / CITED:** the exclusions $H\ne7,21$ and the strong-tournament lower bound used to recover them. **FINITE-EXACT:** the declared small graph and arithmetic controls. Realization of every other odd Hamiltonian-path count remains **OPEN**. No aliquot-period or Collatz conclusion is asserted.

Artifacts: [exact checker](../../04-computation/experiments/tournament_path_defects_20261014e.py), [JSON](tournament_path_defects_20261014e.json).

## 1. Inheritance and the objects that must stay distinct

The closest proved mechanism is the endpoint-reversal family in [arithmetic seams, sections 4–5](arithmetic_seams_20260921_tournaments.md): its Hamiltonian paths are explicitly coded by binary subsets, and its Fermat-number subsequence is already proved there. The canonical hostile is $H=63$: divisibility by seven does not force a strong component with exactly seven paths. The corrected near miss is the rejected $7\cdot3^j$ gap tower; see [THM-344, the two order-eight $H=63$ classes](../../01-canon/theorems/THM-344-n8-h63-complete-omega.md) and the 2026-09-17 arithmetic-braids correction in [MISTAKES](../../01-canon/MISTAKES.md). The least-used sidecar is the distinction between an action on **marked path addresses** and an automorphism of the underlying tournament.

The working board is **ordered vertex / path address / orbit size / strong component / prime quotient / residual term**. The anchor is the owner's seven-and-21 comparison; the niche is an exact Fermat-quotient path count; the wildcard recovers the independent Burnside formula for tournament *isomorphism classes*.

Three different quantities occur:

| Object | Observable | Required marker |
|---|---|---|
| One tournament $T$ | $H(T)$, its directed Hamiltonian-path count | Its actual oriented arcs |
| Marked paths in a specified family | Number of orbits under an address action | Distinguished edge, interior order, action |
| All tournaments on $p$ vertices | $\mathcal T(p)$, the number of isomorphism classes | Relabelling action on the whole class of labelled graphs |

Aliquot dynamics supplies a functional directed graph, usually with missing pairwise comparisons. It does not supply a tournament orientation. A cycle's length is not any of the three counts above.

## 2. What the forbidden values seven and 21 actually say

[THM-1370, all-order omission of seven and 21](../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md) proves $H(T)\notin\{7,21\}$ for every finite tournament. Its last consequence had retained the stronger phrase “EXACTLY two” odd holes, despite its own explicit completeness conjecture. The authorized correction now says **at least** these two holes; it uses the already recorded September17 correction lineage, not a new conjecture resolution.

The mechanism can be recovered compactly. Strong components have a unique transitive order and every Hamiltonian path concatenates their internal paths, so

\[
H(T)=\prod_i H(T_i).
\tag{1}
\]

Busch's [published strong-tournament bound](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3) gives $H(S)\ge5^{(n-1)/3}$ for a strong $n$-vertex tournament, with the sharper residue-class minima. The weaker displayed bound already excludes seven for $n\ge5$ and 21 for $n\ge7$. The remaining small orders have neither value; the checker independently exhausts all labelled tournaments through order six. In a product of positive odd integers, a value seven requires a factor seven, and a value 21 requires a factor seven or 21. Equation (1) completes the exclusion.

This is a strong-component and small-floor obstruction, not a prime-power congruence. In particular, the explicit eight-row tournament

\[
(238,220,216,16,225,206,8,72)
\]

has $H=63$, where bit $j$ of row $i$ means $i\to j$. Adding one vertex that beats every old vertex preserves $H=63$ and makes the graph non-strong. This directly refutes the old gap tower. The supplied row is replayed by subset dynamic programming; no new classification of order-eight tournaments is claimed.

## 3. A lossless path address and the Fermat numbers

For $n\ge3$, define $E_n$ on $0,\ldots,n-1$ by $i\to j$ when $i<j$, except replace $0\to n-1$ by $n-1\to0$. This is an intrinsic binary relation once the total order and distinguished extreme pair are specified. There are no ties.

Exactly one Hamiltonian path avoids the reversed edge: $(0,1,\ldots,n-1)$. Every other one has the unique form

\[
P_S=(S\text{ in increasing order}),\ n-1,\ 0,
       \ (\{1,\ldots,n-2\}\setminus S\text{ in increasing order}).
\tag{2}
\]

Both increasing runs are forced. Conversely every subset $S$ gives a valid path. Thus

\[
H(E_n)=1+2^{n-2}.
\tag{3}
\]

This proof and family are inherited. The API now makes the marker executable: `Address(order, subset)` uses an exact bit mask, or `None` for the exceptional increasing path. `encode` and `decode` are inverse; `transport` first authenticates the entire original adjacency tuple as $E_n$. An address for a different graph is not silently reused.

At $n=2^j+2$, equation (3) is the Fermat number $F_j=2^{2^j}+1$. Orders $3,4,6,10,18,34$ yield $3,5,17,257,65537,4294967297$. The last value equals $641\cdot6700417$. The graph construction and address bijection do not change at this composite value. A prime value of $2^m+1$ forces $m$ to be a power of two, by factoring $x^d+1$ when $m$ has an odd divisor $d>1$; the converse is false.

## 4. A Wieferich criterion for a genuine path-orbit count

Let $p$ be an odd prime and use $E_{p+2}$. Cyclically rotate the $p$ interior **subset coordinates** in (2), and fix the exceptional increasing path. This is a well-defined $C_p$ action on actual Hamiltonian paths because each rotated subset is decoded again using (2).

Its fixed paths are precisely the exceptional path and the empty/full subsets. Every other subset has orbit size $p$: a nonidentity rotation generates $C_p$, so an invariant subset is empty or full. Therefore the number of free path orbits is

\[
\mathcal N_p=\frac{H(E_{p+2})-3}{p}
 =\frac{2^p-2}{p}=2q_p(2),\qquad
q_p(2)=\frac{2^{p-1}-1}{p}.
\tag{4}
\]

The classical Wieferich condition is $2^{p-1}\equiv1\pmod{p^2}$; see [Wieferich's original paper](https://doi.org/10.1515/crll.1909.136.293) for its historical source. Equation (4), proved here without importing any theorem about its frequency, gives

\[
p\mid\mathcal N_p\quad\Longleftrightarrow\quad p\mid q_p(2)
\quad\Longleftrightarrow\quad p\text{ is Wieferich to base }2.
\tag{5}
\]

The checker verifies this condition for the supplied primes 1093 and 3511 and negative controls $3,5,7,17,257,65537$. This is not a claim that no other Wieferich primes exist. It counts a quotient of paths after removing three fixed addresses; $H(E_{p+2})$ itself is congruent to **three**, not zero, modulo $p$.

The primality hypothesis matters. At four interior coordinates, the nonconstant subset mask $0101$ has rotation period two. Dividing all nonconstant addresses by four would not count its orbits correctly. The usual primitive-necklace Möbius inversion would retain those shorter periods; the prime case is exactly where only the two constant subsets need removal.

### Why this is not a hidden tournament automorphism

For every $n\ge4$, $E_n$ has trivial automorphism group. Every directed cycle uses the sole backward edge, so the only directed Hamiltonian cycle is
$0\to1\to\cdots\to n-1\to0$. Any automorphism is a rotation of this cycle. The maximum-score vertices are exactly the adjacent pair $\{0,1\}$; no nontrivial rotation of a cycle of length at least four preserves this pair. Thus only the identity remains. The boundary $E_3=C_3$ has three automorphisms.

For a concrete failure, at $p=3$ the address $S=\{1\}$ yields $(1,4,0,2,3)$. Literal vertex rotation $1\mapsto2\mapsto3\mapsto1$ gives $(2,4,0,3,1)$, which is not a path because $3\not\to1$. Rotating the address and **re-sorting** instead gives the valid $(2,4,0,1,3)$. The missing operation is the ordered decoder.

In contrast, a genuine automorphism group acts freely on directed Hamiltonian paths: fixing an ordered path fixes every vertex. Consequently its order divides $H(T)$. The action in (4) has three fixed paths and is a different action. An arithmetic divisibility of its free-orbit count does not establish a symmetry of the graph or explain the forbidden spectrum values.

### An exact Fermat/Wieferich interface

The companion order-lifting lane supplied a useful quotient comparison. It actually holds for every Fermat **number** $P=2^{2^j}+1$, with primality unnecessary for this algebraic sentence. Put $d=2^{j+1}$; then $d\mid P-1$ and

\[
2^d=(P-1)^2=1+P(P-2).
\]

Raising to $(P-1)/d$ and reducing modulo $P^2$ shows

\[
P\mid2^{P-1}-1,\qquad
\frac{2^{P-1}-1}{P}\equiv2^{-j}\pmod P.
\tag{6}
\]

Here $d$ is invertible modulo the odd number $P$, and $((P-1)/d)(P-2)\equiv2/d=2^{-j}$. Thus when $P$ is prime, the actual path-orbit count (4) is $2^{1-j}\pmod P$, always nonzero. Every Fermat prime is non-Wieferich. At the composite $F_5$, (6) still holds, but the prime-length orbit proof is unavailable. Modular exponentiation checks this boundary without constructing the enormous integer $2^{F_5-1}$.

## 5. A second, inherited arithmetic formula counts different objects

[THM-307, Burnside–QR decomposition](../../01-canon/theorems/THM-307-burnside-qr-decomposition.md) already gives a genuine quotient identity for $\mathcal T(p)$, the number of **unlabelled tournaments on $p$ vertices**. Its residual term must be retained.

For an odd-part partition $\lambda=(\lambda_1,\ldots,\lambda_t)$ of $n$, set

\[
z_\lambda=\prod_a a^{m_a}m_a!,\qquad
e(\lambda)=\sum_i\frac{\lambda_i-1}{2}+\sum_{i<j}\gcd(\lambda_i,\lambda_j).
\]

A vertex permutation with an even cycle fixes no tournament: a half-turn in that cycle reverses a pair. With all cycles odd, its unordered-pair orbits are counted by $e(\lambda)$ and none reverses an orientation. Thus Burnside's counting argument gives

\[
n!\,\mathcal T(n)=\sum_{\lambda\text{ all odd}}\frac{n!}{z_\lambda}2^{e(\lambda)}.
\tag{7}
\]

For odd prime $p$, write

\[
\chi=(2/p),\quad f_p=\frac{2^{(p-1)/2}-\chi}{p},\quad
W_p=\frac{(p-1)!+1}{p},\quad
A_p=\frac1p\sum_{\lambda\ne(1^p),(p)}\frac{p!}{z_\lambda}2^{e(\lambda)}.
\]

Every remaining $z_\lambda$ is prime to $p$, so $A_p$ is an integer. The identity contribution is $2^{p(p-1)/2}\equiv\chi\pmod{p^2}$. The $p$-cycle contribution is $(p-1)!2^{(p-1)/2}\equiv-\chi+p(\chi W_p-f_p)\pmod{p^2}$. Divide their combined sum by $p$ and use $(p-1)!\equiv-1\pmod p$ to obtain the inherited formula

\[
\mathcal T(p)\equiv f_p-\chi W_p-A_p\pmod p,
\qquad q_p(2)\equiv2\chi f_p\pmod p.
\tag{8}
\]

Wieferich cancels the Euler-quotient term but does not cancel the Wilson or remaining-partition terms. At $p=5$, $f_p=1$, $\chi=-1$, $W_p\equiv0$, and $A_p\equiv4$ modulo five. Dropping $A_p$ predicts one, whereas $\mathcal T(5)=12\equiv2$. At seven vertices, $\mathcal T(7)=456$ and $H(E_7)=33$: the two observations are visibly different.

No priority claim is made for (7) or (8); the new contribution here is the explicit type separation and independently checked marked-path interface. The omitted Burnside sum is the same kind of missing-coordinate failure as dropping a cycle length or an endpoint phase, but no cross-domain implication follows from that methodological similarity alone.

## 6. Scope, next obligation, and reproduction

The source-to-target map in (2) is lossless with its order and distinguished edge. The quotient to an unmarked orbit loses the exact subset and requires an orbit representative plus rotation phase to reconstruct the path; at fixed addresses the phase is nonunique. The map to the integer $H$ loses the entire orientation. Formula (8) uses a different source, namely all labelled tournaments under relabelling. None of these maps preserves an aliquot successor, a Collatz native guard, or a ROOT certificate.

The next useful question needs an independently specified arithmetic pairwise observable that selects oriented arcs and whose path statistic preserves a desired arithmetic predicate. The present construction already preserves the Wieferich congruence **for its marked path-orbit count**, but does not turn that congruence into an obstruction to a graph's Hamiltonian-path spectrum or an arithmetic convergence theorem.

The checker exhausts all 33,867 labelled tournaments on one through six vertices; it tests the endpoint family through order12 by subset dynamic programming, and literally enumerates all Hamiltonian permutations through order8. It checks all binary subset orbits for prime lengths $3,5,7,11,13$, automorphisms through order7, the explicit $H=63$ graph and its non-strong padding, six Fermat-number values, and the supplied prime congruences. The Burnside partition formula is independently checked by direct unordered-pair orbits of every vertex permutation through order6, and its prime residual is tested through prime23. Malformed Boolean, float-like/noninteger, wrong-source, composite-prime and out-of-cap inputs are rejected where relevant. The graph API cap is64 vertices, literal path-count cap12, prime-count API cap100000, and partition-count cap23; these are implementation limits, not theorem limits.

```text
python -B -X utf8 04-computation/experiments/tournament_path_defects_20261014e.py
python -B -O -X utf8 04-computation/experiments/tournament_path_defects_20261014e.py
```

Both runs give **8,368 always-active checks** and identical LF-normalized JSON,
SHA256 `d98b193e9c9de6db73d59129d320434bfb65afc47b6d04e96fbea16aa3d97687`.
