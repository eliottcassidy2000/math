# A strong 90-vertex witness for the exponent-61 lower endpoint

Let $M=2^{60}-1$ and $N=E_{61}-1=M(2^{61}+1)=2658455991569831744654692615953842175$. This is a literal strong tournament of count $N$, not just an ordinal factorization or a modular hit. Its section-1 witness has 90 vertices, improving this repository's 124-vertex canonical three-port witness by 34 vertices. No minimal-order claim is made.

## 1. Two simultaneous multiplicative coordinates

For a tournament $U$, let $Z_U(q)$ count ordered covers by $q$ paths with empty slots allowed. Thus $Z_U(1)=H(U)$, where the empty tournament is assigned count 1. If $A\to B$ is an ordinal sum, independently split each slot into its part in $A$ and its part in $B$; this gives
\[
 Z_{A\to B}(q)=Z_A(q)Z_B(q).
\]
Let the outer quotient $D=(8,5,17,6,11)$ be specified by outgoing masks in zero-based order. Mark vertex 0. Every return between two consecutive visits to the marked module needs a crossing from an outneighbor of 0 to an inneighbor of 0. Vertex 0's only outneighbor is 3, so at most one return is possible, giving at most two visits. There are 13 ordinary Hamilton skeletons and six with the mark appearing twice, checked by enumerating the permutations of the four unmarked vertices and their gaps for the mark. Replacing the mark by $U$ therefore gives
\[
 H(D[0\mapsto U])=13H(U)+6[Z_U(2)-2H(U)]=H(U)+6Z_U(2).
\]

Use the following six strong blocks \(A_i\), followed by 33 singleton blocks, all in ordinal sum. Each row certifies \(H_i=H(A_i)\) and \(Z_{A_i}(2)=2^{v_i}u_i\).

| \(H_i\) | \(u_i\) | \(v_i\) | Outgoing masks |
|---:|---:|---:|---|
| 1321 | 775 | 5 | (256,289,67,23,295,141,59,95,236) |
| 1655 | 3641 | 3 | (258,12,393,305,295,71,31,123,226) |
| 8723 | 1321 | 7 | (768,269,73,561,551,391,59,351,604,230) |
| 1359 | 1599 | 4 | (258,288,283,403,99,141,303,87,176) |
| 155 | 427 | 3 | (2,40,35,21,39,201,31,95) |
| 287 | 151 | 5 | (10,20,25,82,161,15,55,111) |

The exact products are \(\prod H_i=M\), \(\prod u_i=M/3\), and \(\sum v_i=27\). Hence \(H(U)=M\) and \(6Z_U(2)=2^{61}M\). The sum of the six orders is 53, and \(53+33+4=90\). The quotient is strong. A quotient path lifts to reach distinct blocks, and a directed cycle through the mark lifts to connect vertices within \(U\); hence the full tournament is strong despite the internal ordinal decomposition of \(U\).

The two-coordinate composition is how the count was found. For each directed quotient with \(c_1\) one-run and \(c_2>0\) two-run skeletons, write \(c_1-2c_2=d\), \(c_2=2^e f\) with \(f\) odd. Seek signatures of small blocks whose coordinate products are \((M/d,M/f)\) and whose accumulated valuation is at most \(61-e\); insert enough transitive singletons to fill the remainder. This is divisor-state exact cover in two coordinates, not a numerical search over \(2^{121}\) path counts. The same method also gives [91- and 93-vertex cyclic-quotient families](STRONG-93.md).

## 2. All component orders are genuine classes

The quotient is prime: across the 25 nontrivial proper subsets, every subset has an exterior distinguisher. The inflated marked block \(U\) is the unique maximal proper nontrivial module. Indeed any module meeting both \(U\) and its exterior projects to a nonsingleton module in \(D\). Its projection must be all of \(D\); but if it omitted a vertex of \(U\), that vertex sees both directions to the exterior of \(U\), a contradiction. A module disjoint from \(U\) also projects to a module in \(D\) and therefore is a singleton. So any isomorphism of inflations restricts to \(U\). Its SCC condensation is totally ordered; all six nonsingleton strong factors have distinct Hamilton counts. Thus ordering these six types and the 33 indistinguishable singletons gives exactly
\[
 \frac{39!}{33!}=2349088560
\]
distinct strong 90-vertex classes within this displayed family. This is a lower bound on the full 90-vertex fiber, not a claim of completeness.

## 3. Literal arc check and its finite reduction

The companion code gives exact subset-endpoint counts on blocks of order at most 10, checks the quotient grammar, and covers all \(\binom{90}2=4005\) reversals of the displayed representative. Its counts by type are six between quotient exterior vertices, 344 contacts with \(U\), 209 internal to a strong block of \(U\), and 3446 between distinct blocks of \(U\). None produces \(N+2\). In the exact class-wide extension, the intra-\(U\) proof does not require any arc scan by order:

- If an arc lies within a strong block \(A_i\), its change upstairs is \((M/H_i)\Delta H_i+6(Z_U(2)/Z_i(2))\Delta Z_i(2)\). Its odd divisor \(M/\operatorname{lcm}(H_i,u_i)>1\) rules out \(+2\).
- If an arc crosses an ordinal cut \(A\to B\), orient it first as \(u\to v\). In any induced tournament including both endpoints the number of lost paths is the product of counts in \(A\) ending \(u\) and in \(B\) starting \(v\). New paths concatenate two independent ordered two-slot covers across \(v\to u\), with one slot starting at \(u\) on the \(A\) side and one ending at \(v\) on the \(B\) side. Each relevant one-path read injects into its two-slot counterpart by putting the specified endpoint in a singleton slot and its complement in the other. Thus every induced Hamilton count is nondecreasing. In particular the two-slot increase is at least twice the one-slot increase and is divisible by four (use Rédei on the unordered nonempty bipartitions and on the two terms with an empty slot). Therefore \(H+6Z(2)\) cannot rise by two.
- An arc between exterior quotient vertices induces the same folded slot formula for every order, since all weak slot totals of \(U\) are products independent of summand order. The six exact tests therefore apply to the whole family.

For a contact between a vertex \(u\) of \(U\) and an exterior vertex \(a\), enumerate the exterior order and the gaps where \(U\) is visited. The only nonuniform contact is at \(a\). Its incident boundary needs either a designated visit starting or ending at \(u\), or the opposite constraint excluding \(u\). If both adjacent symbols are visits to \(U\), one constraint requires \(u\) and the other then excludes it automatically. All other incidences follow the quotient.

Let \(S_U(u;j)\) denote ordered weak \(j\)-slot covers with a designated slot starting at \(u\), and similarly \(F_U(u;j)\) for ending at \(u\). For \(U=A\to C\to B\) and \(u\in C\), the exact endpoint formulas are
\[
 S_U(u;j)=Z_A(j-1)S_C(u;j)Z_B(j),\qquad
 F_U(u;j)=Z_A(j)F_C(u;j)Z_B(j-1).
\]
Excluding empty slots by inclusion–exclusion gives the \(r\)-run rooted fill
\(\sum_{j=1}^r(-1)^{r-j}\binom{r-1}{j-1}S_U(u;j)\), and the same expression for ending reads. The gap enumeration therefore depends only on which nonsingleton blocks lie before \(C\), and the count of singletons on that side, rather than on the order within either side. The companion code compresses the entire contact test to those finite configurations. No claim is made that \(N+2\) has no witness outside these classes.

There are exactly $53\cdot32\cdot34+64\cdot33=59776$ contact configurations per exterior vertex: a nonsingleton-root contact chooses one of 53 local vertices, a prefix subset of the other five nontrivial factors, and a split of all 33 singletons; a singleton-root contact chooses a prefix subset of six nontrivial factors and a split of 32 remaining singletons. The exact contact search finds no change $+2$ in all $239104$ configurations. Combined with the order-independent arguments for the other arcs, all $2349088560$ geometric classes above miss the literal upper neighbor by one flip. Run the two companion Python files in the same directory; the 90-vertex script imports the 93-vertex script's small-block and cover helpers.

## 4. Independent one-apex bound

A strong one-apex presentation of $N$ (the vertex deletion is transitive) would require exactly 122 spine vertices. Its normal contact word starts with 1 and ends with 0. The first-to-last inversion alone contributes $2^{w-2}$ at spine width $w\ge3$, and $N<2^{121}$ forces $w\le122$. Conversely swapping adjacent $01$ to $10$ increases the weighted inversion read: the two bits contribute a new inversion and all affected earlier or later old inversions move to nondecreasing distance weights. Thus at a given width the maximum is achieved on a two-run word, where the exact value is $2^{a+b}-2^a-2^b+3$. Its maximum at $w=121$ is $2^{121}-3\cdot2^{60}+3<N$. The maximum is nondecreasing with width, so no shorter spine works.

The latter result does not prove an unrestricted one-apex exclusion. The finite scanner checks all compositions of width 122 into at most eight alternating runs, modulo reverse complement. Its accelerated final-two-run inversion is exact: if current read and carrier are $H,P$ after an even number of runs, and the remaining two widths are $a,b$ with sum $w$, the final count is
\[
 H+(P+1)2^w+2-(P+1)2^a-2^b.
\]
Here $P$ is even. To hit $N$, set $K=H+(P+1)2^w+2-N$. If $a\ne b$, then $v_2(K)=\min(a,b)$, giving at most two candidate assignments; if $a=b$ it is separately checked. The exact replay found no hit up to eight runs. More runs remain open.
