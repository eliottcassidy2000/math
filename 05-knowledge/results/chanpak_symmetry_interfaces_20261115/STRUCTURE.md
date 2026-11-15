# Signed run transfer for tournament substitution

**Results.** A finite modular substitution has an exact signed transfer that counts not just Hamilton paths but all ordered path covers fixed by relabeling followed by reversal. On a self-converse tournament this gives the number of Hamilton-path presentations fixed by converse, by a centralizer division. A triple inflation $C_3[S,S,S]$ has a superexponentially vanishing proportion of fixed presentations as the size of its self-converse child $S$ grows. A different symmetric construction $C_3[M,M^{\mathrm{op}},1]$ gives a quadratic norm law: for every choice of $M$ of size at least two, its Hamiltonian count is exactly 3 modulo 8, an obstruction to using that natural composition to select general Chan–Pak targets.

The ordinary run-content factorization is inherited from THM-3121. The additional state is the signed cover profile. The group conventions use the existing free-action and converse-presentation dictionary (LEM-003, THM-781 and THM-4741), not a new claim that quotient modules classify prime tournaments.

## 1. Ordered and signed covers

Write $O_T(k)$ for the number of ordered $k$-tuples of nonempty vertex-disjoint directed paths covering $T$. The paths themselves are directed vertex lists, so $O_T(k)=k!\,pc_T(k)$, with $O_T(1)=H(T)$. Suppose that $\tau$ is an anti-involution: it reverses every arc and $\tau^2=1$. Define $Z_{T,\tau}(k)$ to count the ordered covers satisfying
$$
 P_{k+1-j}=\tau(\operatorname{rev} P_j) \quad(1\le j\le k).
$$
An anti-involution has at most one fixed vertex, since any arc between two fixed vertices would have to reverse. A self-converse tournament has such an involution: its ordinary automorphism group has odd order (an involutive vertex permutation would reverse the arc across a transposed pair), and its group of automorphisms and anti-automorphisms has twice that order. Cauchy's theorem supplies an involution outside the odd subgroup. At order one we just use vertex reversal with identity relabeling.

## 2. Exact transfer over a quotient and a rational signed kernel

Suppose $T=Q[T_1,\ldots,T_d]$, with the displayed nonempty modules permuted by $\tau$, which induces the quotient anti-involution $\sigma $. A canonical modular partition has the required preservation since both isomorphism and converse preserve modules, inclusion, and intersections (hence strong modules and their maximal proper blocks at a nontrivial node). The cover transfer itself needs only the displayed permutation property, not the canonical partition.

Take any ordered tuple $w_1,\ldots,w_k$ of nonempty directed quotient walks. Each walk is a list of visited modules with no consecutive equal letters. Let $c_i$ be its aggregate content, counting all visits to module $i$ across the tuple. Discard any tuple with $c_i=0$ or $c_i>|T_i|$. Read occurrences of each letter in order through the tuple. Into these occurrences insert the components of any ordered $c_i$-cover of $T_i$. Uniform cross-arcs make all concatenations legal, and collapsing maximal consecutive runs in each path recovers the quotient walks and the factor covers. This is the ordinary THM-3121 kernel, with ordered instead of unordered factors.

For the signed $Z$ transfer, retain only quotient tuples satisfying
$$
 w_{k+1-j}=\sigma (\operatorname{rev} w_j).
$$
For each such tuple, its number of lifts is
$$
 \prod_{i<\sigma (i)} O_{T_i}(c_i)
 \prod_{i=\sigma (i)} Z_{T_i,\tau|_{T_i}}(c_i).
$$
Indeed occurrences of $i$ correspond in reverse order to occurrences of $\sigma (i)$. Choose the cover arbitrarily on exactly one representative of each swapped pair; the reverse-label action forces the mate. On a fixed module its ordered cover must itself satisfy the same signed condition. Summation proves the all-order signed transfer.

The signed quotient kernel is rational for fixed $k,Q,\sigma $. Put $y_i=x_i x_{\sigma (i)}$, and let $w_i(y)$ denote the content-generating function of nonempty quotient walks ending at $i$; it is the $i$-th coordinate of the standard quotient-walk resolvent in THM-3121. Write $W_Q(y)=\sum_i w_i(y)$. A symmetric ordered quotient tuple with $k=2r+\eta$ is given by its first $r$ arbitrary walks, their forced mates, and, if $\eta =1$, one symmetric central walk. The series for the central walk is
$$
 \operatorname{Pal}_{Q,\sigma }(x)=
 \sum_{i\to\sigma (i)} w_i(y)+
 \sum_{c=\sigma (c)}x_c\left(1+\sum_{i\to c}w_i(y)\right).
$$
An even-length central walk is a nonempty first half ending at $i$, followed by its image in reverse order, so the joining arc is $i\to\sigma (i)$. In the odd case the central letter $c$ is fixed; the incoming arc $i\to c$ automatically implies its outgoing reflected arc. Thus the symmetric quotient tuple series is $W_Q(y)^r\operatorname{Pal}_{Q,\sigma }(x)^\eta $. Coefficients are multiplied by the factor-profile weights above; no additional factorials occur.

## 3. Divide by the centralizer, not by the full automorphism group

Let $G=\operatorname{Aut}(T)$ and $C=C_G(\tau)$. Ordinary automorphisms act freely on ordered covers, since a vertex list in every assigned path is fixed pointwise if the cover is fixed. Relabeling by $\tau$ and reversing all lists and the list of paths induces an involution on the ordinary $G$-orbits. The number of those orbits fixed by the involution is
$$
 z_{T,\tau}(k)=Z_{T,\tau}(k)/|C|.
$$
To see this, count pairs consisting of a cover and an anti-automorphism fixing that cover under the signed action. In a stable ordinary orbit of size $|G|$, each cover has exactly one such map (existence by stability, uniqueness by freeness); unstable orbits contribute none. A fixing anti-automorphism must square to the identity by freeness. Sylow conjugacy in the group of order $2|G|$ shows all such involutions are conjugate by $G$: if a conjugating element is anti, multiply it by the source involution to make it ordinary without changing the conjugation. There are $|G|/|C|$ involutions, each fixing $Z_{T,\tau}(k)$ covers. Equating counts proves the quotient, in particular $0\le z(k)\le O_T(k)/|G|$. The $k=1$ specialization is the number of converse-fixed Hamilton-path presentations, not the total number of labeled antisymmetric paths.

When the given modular partition is preserved by all of $G$, let $G_Q$ be the subgroup of quotient automorphisms preserving the isomorphism types of its modules. The usual exact sequence has kernel $\prod_i\operatorname{Aut}(T_i)$: a type-preserving quotient permutation lifts by arbitrary isomorphisms between corresponding modules. For the centralizer the exact order is
$$
 |C|=|C_{G_Q}(\sigma)|\,
 \prod_{i<\sigma (i)}|\operatorname{Aut}(T_i)|\,
 \prod_{i=\sigma (i)}|C_{\operatorname{Aut}(T_i)}(\tau|_{T_i})|.
$$
Projection must commute with $\sigma $. Conversely, for any such quotient element, choose a lift on one representative of each swapped pair and force the map on its mate by commuting with $\tau$; in the sole possible fixed module take any internal map commuting with the restricted involution. Count independent choices. This formula keeps the symmetry of the quotient and that of the children separate.

## 4. Repeated cyclic triples: total and mirror-fixed counts separate

Use $A\to B\to C\to A$. For $T=C_3[S,S,S]$, with $S$ self-converse of order $m$, choose $\sigma $ to swap A and B and fix C. Put $\tau$ on each component by an internal anti-involution of $S$ combined with this block permutation. Write
$$
 a_j=O_S(j)/|\operatorname{Aut}S|,\qquad
 z_j=Z_{S,\tau|_C}(j)/|C_{\operatorname{Aut}S}(\tau|_C)|,
$$
with zeros at $j=0$ and beyond $m$. Any proper module of this triple lies in one of the three factors: if it meets two factors, every vertex in the third distinguishes those two and thus belongs to it; then any excluded vertex distinguishes vertices in two of the other factors. Consequently the natural three factors are canonical. The ordinary quotient group is $C_3$, and the centralizer of the reflection in this group is trivial.

A quotient Hamilton word with three positive contents has counts differing by at most one. Equal contents have three possible starts; otherwise a specified content triple has one start. This recovers THM-3121's ordinary cyclic diagonal and gives the total presentation count
$$
 q(T)=\sum_{j=1}^m(a_j^3+a_{j+1}a_j^2+a_{j+1}^2a_j).
$$
A symmetric quotient word is unique at every sufficiently long length once its center is fixed: the center letter is C in odd length, the center arc is A→B in even length, and it extends uniquely in the directed cycle. For paired content $j$, the three possibilities for the C count are $j-1,j,j+1$; content zero is simply discarded. The signed formula and centralizer cancellation therefore give
$$
 p(T)=\sum_{j=1}^m a_j(z_{j-1}+z_j+z_{j+1}).
$$
Since $0\le z_j\le a_j$, Cauchy's inequality on adjacent products and Hölder give
$$
 p(T)\le3\sum_j a_j^2\le3m^{1/3}\left(\sum_j a_j^3\right)^{2/3}
 \le 3m^{1/3}q(T)^{2/3}.
$$
Moreover $a_m=m!/|\operatorname{Aut}S|\ge 2^{v_2(m!)}$ since the denominator is an odd divisor of $m!$. The elementary identity $v_2(m!)=m-s_2(m)\ge m-1-\lfloor\log_2m\rfloor$ then yields
$$
 p(T)/q(T) \le 3m^{1/3}\,2^{-m+1+\lfloor\log_2m\rfloor}.
$$
The theorem applies in particular to ternary towers of increasing depth. It does not say that every self-converse composition has few fixed views; it is a quantified conclusion for equal cyclic triples and does not require any assumption on the internals of $S$ beyond being self-converse. The finite control $S=TT_3$ already has $H(T)=2721$, $q(T)=907$, $p(T)=31$.

An elementary automorphism bound improves the exponential estimate to a factorial one. Every tournament on $m\ge3$ vertices has a vertex whose out- and indegrees $a,b$ are both at least $(m-2)/4$: otherwise the sets of vertices respectively below those thresholds would each have size less than $m/2$, by their internal average out- or indegree, and cannot cover the vertices. For such a vertex $m^2\le9ab$: check $m=3,4,5,6$ using $a,b$ positive integers, and for $m\ge7$ use $ab\ge(m-2)(3m-2)/16$. Its stabilizer injects into the product of the ordinary automorphism groups on its in- and out-neighborhoods. Induct from $m=1,2$ to obtain
$$
 |\operatorname{Aut}(S)| \le m\;\frac{3^{a-1}}a\;\frac{3^{b-1}}b
 \le \frac{3^{m-1}}m.
$$
Since $q(T)^{1/3}\ge a_m=m!/|\operatorname{Aut}(S)|$, the stronger bound is
$$
 \frac{p(T)}{q(T)} \le \frac{3^m m^{-2/3}}{m!}.
$$
In particular the decay is superexponential. For the ternary tower $C_3[C_3,C_3,C_3]$, direct replay gives $H=3159$, ordinary group order 81, labeled fixed path count 33, centralizer order 3, and hence $(q,p)=(39,11)$.

The same scarcity holds for any **fixed strongly connected self-converse quotient** $Q$ of order $d\ge3$ and repeated self-converse child $S$. Choose compatible quotient and child anti-involutions. The quotient has a Hamiltonian cycle, so repeating its directed order $m$ times gives at least $(m!)^d$ ordinary lifts: each child is traversed as $m$ singleton runs that can be assigned freely. A fixed quotient word has at most $\lceil d/2\rceil$ independently chosen factor covers (one per swapped pair, plus at most one fixed factor). Any ordered cover count of $S$ is at most $m!2^{m-1}$, by encoding the concatenated vertex permutation and its cut set. There are at most $dm\cdot d^{dm}$ possible quotient words. Finally converting the labeled fixed/ordinary ratio to the fixed/total presentation ratio multiplies it by at most $|\operatorname{Aut}(T)|\le 3^{dm}$. Hence
$$
 \frac{p(Q[S,\ldots,S])}{q(Q[S,\ldots,S])}
 \le (dm)(3d)^{dm}2^{m\lceil d/2\rceil}(m!)^{-\lfloor d/2\rfloor}.
$$
This proof does not require the chosen repeated-module partition to be canonical. In the convention of THM-4741 the degree at such a self-converse class in the fundamental black graph is $(q-p)/2$, so the fixed quotient inflation asymptotically approaches the full $q/2$ degree scale. It explains why increasing a symmetric module's size, while keeping a strong quotient fixed, need not produce more mirror-fixed views as a proportion. No corresponding assertion for general, possibly transitive quotients is made.

## 5. The mirror norm and a Chan–Pak boundary

Let $M$ be any tournament of order at least two, $g=|\operatorname{Aut}M|$, $h=H(M)$, and $u=O_M(2)$. Put $T=C_3[M,M^{\mathrm{op}},1]$. Matching labels of the two mirrored modules and fixing the singleton is an anti-involution. The three modules are again canonical, and the singleton prevents every nontrivial quotient rotation. Thus $|\operatorname{Aut}T|=g^2$ and its anti-centralizer has order $g$.

Because the singleton is visited once, the ordinary run-content terms are $3h^2,2hu,u^2$. The signed kernel forces a reflected pair with either one or two runs, giving $Z_{T,\tau}(1)=h+u$. Therefore
$$
 H(T)=(h+u)^2+2h^2,\qquad
 p(T)=(h+u)/g,\qquad
 q(T)=p(T)^2+2(h/g)^2.
$$
The paired run gives a real arithmetic restriction on this entire all-size construction, not merely a similarity in names to a quadratic norm. Here $h$ is odd by Rédei, $u$ is even by exchanging the two paths, and $g$ is odd. Hence $H(T)\equiv q(T)\equiv3\pmod 8$ for all choices of the mirrored core. For a transitive $m$-vertex core $h=g=1,u=2^m-2$ and we get the effective infinite family $(2^m-1)^2+2$, starting $11,51,227,963,3971,15875,64519,261123\ldots$. These are not assertions of primality or spectrum completeness. In particular changing the mirrored core alone cannot construct a prescribed odd target outside the residue 3 modulo 8: at least one symmetry condition must be relaxed. The one-vertex exceptional core would make all three factors equal, so its extra quotient symmetry must not be divided as though the singleton were unique.

**Replay.** Run `python verify_structure.py` from this directory. It compares every ordered and signed cover coordinate against independent literal permutation-and-cut counting for ten substitutions (66 coordinates per profile), checks both ordinary automorphism and anti-centralizer orders by degree-preserving permutation search, and includes a mirrored modular norm, a fixed middle module, a quotient involution with no fixed vertex, and a prime five-vertex quotient.

## 6. A fixed mirror with a bounded seed palette hits every odd dyadic residue

First, for any center $C$ there is an exact residue law independent of the mirrored core:
$$
 H(C_3[M,M^{\mathrm{op}},C])\equiv 3H(C)+O_C(2)=H(C_3[1,1,C])\pmod 8.
$$
In the cyclic diagonal, center runs $k\ge4$ contribute an $O_C(k)$ divisible by $8$ (ordered covers have the factor $k!$). Center run $3$ has even $O_C(3)$, and both mirrored-factor run counts are at least $2$, so its terms are divisible by $8$ as well. At center run $1$, the bracket is $3h^2+2hu+u^2\equiv3\pmod8$ with $h=H(M)$ odd and $u=O_M(2)$ even. At center run $2$, the bracket is $h^2+2hu+3u^2+2uw+w^2\equiv1\pmod4$, where $w=O_M(3)$ is even. Multiply by the even $O_C(2)$. This proves the law. The self-converse choices $1,TT_2,TT_3,C_3$ already give all four odd residues modulo 8 for every mirrored core.

A stronger exact statement keeps **one** mirrored core $M=TT_2$. There are sixteen bounded self-converse seed choices $D$ (one is empty) such that the actual tournaments
$$
 T_{\ell,D}=C_3[TT_2,TT_2,TT_\ell\to D\to TT_\ell],\qquad \ell\ge1
$$
hit every odd residue modulo every power of two, even with arbitrarily large $\ell$. Here $\to$ denotes an ordinal stack, with all arcs from the earlier block to the later. All these tournaments are self-converse: the two $TT_2$ modules are exchanged by mirrored reversal, and the center is a self-converse stack.

For a tournament $D$, let
$$
 A_D(j)=\sum_{f:V(D)\to[j]}\;\prod_{i=1}^j H(D[f^{-1}(i)]),
 \qquad H(\varnothing)=1 .
$$
Then $A_D(1)=H(D)$, and the exact nonempty-color counts are
$$
 O_D(2)=A_D(2)-2A_D(1),\qquad
 O_D(3)=A_D(3)-3A_D(2)+3A_D(1).
$$
This is inclusion-exclusion on unused colors, followed by independently choosing the Hamiltonian path of each nonempty color class. Color-class counts multiply under stacks, so $A$ multiplies under stacks. Also $A_{TT_\ell}(j)=j^\ell$.

In the cyclic kernel with both external modules $TT_2$, the central run-count coefficients at $1,2,3$ are exactly $11,17,4$, and larger counts do not occur. Indeed $O_{TT_2}(1)=1,O_{TT_2}(2)=2$; using the cyclic contents with differences at most one gives $3+2\cdot2+2^2=11$, $1+2\cdot2+3\cdot2^2=17$, and $2^2=4$. Inclusion-exclusion and the equal transitive wings therefore give the exact count
$$
 H(T_{\ell,D})=4A_D(3)9^\ell+5A_D(2)4^\ell-11H(D).
$$

For any $r\ge5$ and $2\ell\ge r$ the $4^\ell$ term vanishes modulo $2^r$. The powers of $9$ fill the subgroup $1+8\mathbb{Z}$ modulo $2^{r-2}$: for positive $k$, $v_2(9^k-1)=3+v_2(k)$, by factoring a square minus one (its plus-one factor has valuation 1) and factoring an odd power (the geometric factor is odd). Thus the order is the size $2^{r-5}$ of that subgroup. Also $A_D(3)$ is odd, since every term is odd and their number is the odd integer $3^{|D|}$. Consequently this seed's outputs run through precisely all lifts of
$$
 4A_D(3)-11H(D)\pmod {32}.
$$
Adding a period to $\ell$ allows the lower bound $2\ell\ge r$ or any other lower bound.

The seeds cover all sixteen odd residues modulo 32, as follows.

- For $0\le t\le7$ take a stack of $t$ directed triangles. Each triangle has $H=3$ and $A(3)=33$ (three monochromatic assignments each add two to the $27$ all-unit contributions). Thus $A_D(3)\equiv1\pmod8$ and $H(D)=3^t$. Powers of 3 modulo 32 have order 8 and run through all residues that are 1 or 3 modulo 8; the map $u\mapsto4-11u$ sends those eight values to the same eight residue classes modulo 32.
- Use a nearly transitive four-vertex tournament $S_4$ as the central seed; it is obtained by reversing the end arc of $TT_4$, has $H=5$ by literal path enumeration, and is self-converse by vertex-order reversal. Put a stack of $t$ triangles on each side of $S_4$, $0\le t\le3$. Then $A_D(3)$ is constant modulo 8 because $33^{2t}\equiv1$, and $H(D)=5\cdot9^t$. The modular outputs form all four residues congruent to 5 modulo 8: $9^t$ runs through $1+8\mathbb Z$ modulo 32, $A_D(3)$ is odd, and the coefficient of $H(D)$ is an odd unit.
- Repeat the same paired stacks with the circulant five-vertex tournament whose arcs have differences $1,2$ modulo 5. Negation is an anti-involution. Starting at zero, the three four-step strings without a repeated vertex are $1111,2112,2222$, so it has $H=15$. Its corresponding residues are all four congruent to 7 modulo 8.

Reduction from $r=5$ handles the smaller powers. This is an exact mathematical dyadic targeting result, not absolute targeting: since the seed palette is finite and the actual counts in each sequence grow as a positive constant times $9^\ell$, the union produces at most $O(\log X)$ distinct counts up to $X$ and has zero density. Thus an all-modulus dyadic selector and severe global sparsity can coexist in one physically realized symmetry grammar.

## 7. Finite arithmetic states at every modulus

The ordinary profile has $k!\mid O_T(k)$. For the signed profile put $r=\lfloor k/2\rfloor$. The group of permutations commuting with reversal of the $k$ index positions contains the group of order $2^r r!$ obtained by permuting the $r$ mirrored pairs and independently exchanging the two members of each pair, fixing the central index if present. It acts freely on ordered covers (the component lists are disjoint, nonempty, and hence distinct) and preserves the signed condition. Thus
$$
 2^{\lfloor k/2\rfloor}\lfloor k/2\rfloor!\mid Z_{T,\tau}(k).
$$
Modulo any fixed integer, both profiles truncate at finite indices independent of tournament size. With finitely many quotient templates (and compatible involutions in the signed case), the exact substitution rules then give a finite residue-state computation up the decomposition tree: for each modulus and each retained output coordinate it uses only finitely many quotient-word tuples and child-profile residues. This finite-state statement concerns the labeled counts; in odd-modulus calculations one must not divide away automorphism orders that are nonunits. It does not imply a finite unweighted-state description of the absolute spectrum.
