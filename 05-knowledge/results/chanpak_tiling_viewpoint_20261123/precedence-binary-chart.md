# A binary nonmultiplicative chart inside transitive precedence

Let M be a nonempty tournament module, with D=H(M), Q its number of ordered two-slot Hamilton-path covers allowing empty slots, and transitive precedence tournament. Put M after a transitive spine prefix v_0,…,v_(L−1), so all prefix vertices beat M, and add apex a with every M vertex beating a. An O contact is a→v_i; an I contact is v_i→a. Read B=Σ_(i:O) 2^(L−1−i), binary O tags in reverse positional weight, and put
\[
E(B)=\sum_{p>q:\,p,q\text{ both set bits of }B}k_{p-q},
\quad k_1=2,\quad k_d=2^{d-1}\ (d\ge2).
\]
Then the exact ordinary height is `H=(1−E(B)) D+B Q`.

Proof: if M is transitive of order t, a Hamilton path splits into ascending words before and after the apex. If its first vertex after the apex lies before its last vertex before the apex in the spine, their tags are O,I at distance d, giving \(2^{d-1}\) free assignments of intermediate vertices. Otherwise the two endpoints are consecutive with tags I,O. Paths starting or ending at the apex, together with those consecutive I,O choices, number one plus the number of consecutive O,I choices (by the change of tag from first to last). This gives the weights \(k_d\) and a constant one. For a fixed prefix O, the Q-leading contribution from I's inside M is 2^(L−1−i)(2^t−1), plus one when it is adjacent to the block. Full possible contribution to later prefix positions is 2^(L−1−i), except the final position has none and its adjacency bonus compensates. The only lost contributions are the O,O pairs requested by E. Thus H=1−E+B·2^t. For arbitrary M, every Hamilton path assigns the module to at most two runs around the apex; the fixed external contacts are uniform and allow ordinary paths inside each occupied run. One occupied run contributes D, two Q−2D. Therefore the height is linear in D,Q, fixed by t=1,2 and equal to the claimed chart.

The coefficient of D equals 1 exactly when B is zero or has one set bit, since all E summands are positive. Therefore this natural single-apex chart cannot independently tune an arbitrary Q coefficient while keeping the D read additive and unchanged.

## Sufficient spacing for transitive precedence

Assume |M|≥3, all O positions are separated by at least 4, and at least three prefix I's follow the final O. Cross cycles through M all include a and a nonempty directed module path. If the observation pair includes at least one module vertex, the selected module vertices are mandatory and at most two chosen; the number of possible supersets is an even power of 2, and each has odd Hamilton-path count. Thus those cycles do not change original precedence across blocks or within the M fiber. For a pair of external vertices, the parity of a nonempty module choice instead equals one (2^|M|−1 choices); the external parity is exactly what obtains on replacing M by a single terminal I singleton.

Now in this contracted one-apex word compare each O in isolation with all other sites I. Its active component is one reverse outer arc in a transitive order, with at least five vertices and natural precedence order by the [outer-arc formula](PRECEDENCE-FIBER.md#exponentially-many-classes-with-the-same-parity-order). Summing these individual cycle changes and those raw reversed O arcs recovers the multiple-O tournament, except for spurious cycles from an earlier O to a later O as their maximum site. For any observation pair each such spurious type still has ≥1 unforced intermediate vertex: the O distance is ≥4 and at most two other sites are forced. Its cycle parity is zero. Thus the external precedence order stays natural. Combining the arguments gives exact precedence order `prefix ⇒ D(M) ⇒ a`, so the output is rigid. Taking the first prefix position to be O makes the whole output strong: a reaches first then all later prefix/module vertices, every module vertex returns to a, and every prefix vertex reaches M.

## Consequences

For the H=25 seed M=J_4⇒J_4 one has Q=32²=1024. The chart contains `25+1024·2^q` for q≥0, as well as the two-digit rail `25+1024(2^p+2^q)−25·2^(p−q−1)` for q≥0,p−q≥3 (use the exact criterion below), and all sparse-B outputs. These exact strong nonmultiplicative extensions are not claimed to cover an interval.

More generally there is a direct strong CRT construction at all three primes without extra seeds. For a modulus 2^a3^b11^c take q≥max(a,3), g=5·3^max(b−1,0)·11^max(c−1,0), h=g+4, and
\[
 B=2^q(1+2^g)(1+2^h).
\]
Its four bits satisfy the spacing restriction; B is divisible by the modulus by LTE (or elementary lifting of 2^g+1). With x=2^g,y=2^h the six gap contributions yield 2E=2x+2y+xy+y/x. Since x=−1 mod33, E=−1 mod33 and 1−E=2 mod33; E is also even by spacing. Thus A=1−E is a unit for the full modulus and the chart gives H≡A D. Choose the established ordinal bank word M to realize D≡A^(−1)r (padding it if needed); output is strong, rigid and in transitive precedence, now at a common order of these words and every larger order by padding *inside M*. This proof is independent of the alternative capped-join mechanism.

## Exact wrap operator

A simpler nonmultiplicative transformation wraps M between two singleton guards a→M→b→a. It has H'=D+Q. Modulo two, directed cross-cycles using both guards correspond to one nonempty module subset, each with odd path count. For pairs inside M or between M and a guard, those extra cycle counts are even provided |M|≥3; for the guard pair they are odd, reversing it in precedence. Thus the full precedence order is a⇒D(M)⇒b. For |M|=2 its sole internal pair also reverses, still transitive. At size1 this is a directed triangle. For the ordered k-slot reads the exact closure is
\[
 Q'_k=k^2Q_k+kQ_{k+1}−kQ_{k-1}.
\]
Assign guards to distinct slots (k(k−1) choices) and each contributes no extra path restriction in its own slot; if they share a slot the induced formula is H(S)+Q_2(S) with a subtraction of one for S empty, giving k(Q_k+Q_{k+1}−Q_{k-1}). In particular Q'_2=4Q+2R−2D.

## Exact precedence criterion: gaps three suffice and are necessary
For |M|≥3 the precedence order is transitive if and only if D(M) is transitive and no two O tags are within distance two. There is no restriction on how many I's follow the final O. The sufficiency proof above still works with gap3: a spurious cycle between two O endpoints at distance3 contributes odd parity precisely when both chosen observation vertices are the two intermediate I's. This swaps their adjacent positions in the otherwise natural prefix; distinct such swaps cannot overlap. In a single isolated-O contracted component of active order ≥5 the J_s result applies; at orders3 and4 its restriction to external sites still has natural order (order3 is a directed triangle with its terminal contracted-I singleton removed; order4 only swaps that singleton with the preceding interior). Thus small trailing gaps cause no problem.

For necessity, write o_i=[O at i], extending by zero on the left and using the contracted terminal I singleton on the right. Directed cycles through a and prefix vertex i make an odd contribution only when i is a min/max of an adjacent O,I pair, or i is the mandatory middle of a length-two O,I gap. Consequently the orientation a→i in D(T) is the majority bit of (o_(i−1),o_i,o_(i+1)): over F₂ it is the sum of the three pair products. If two O's lie within distance two, some such bit is 1. Together with the unchanged cross orientations i→M→a this gives a directed triangle in D(T). The entire iff statement was independently replayed by direct Hamilton words for every prefix length≤7, with transitive M of order3.

For a fixed module this complete binary chart still cannot cover a dense family: E≤B by summing geometric weights for each O, so H≥D+(Q−D)B. Binary numbers with no two 1's within distance two among n positions satisfy the recurrence f_n=f_(n−1)+f_(n−3) and grow as O(ρ^n), with ρ^3=ρ²+1. Hence only O(X^(log₂ρ)) heights through X, exponent about 0.551. By contrast exact sparse rails include B=1 at D25,Q1024, yielding the prime 1049, and B=2^10 or 2^12 yields the further primes 1048601 and 4194329 (ordinary finite trial division), all rigid with transitive precedence and strong.

## Literal H=25 seed
The exact seed is the ordinal join of two copies of the canonical four-site J4 rows (6,8,10,1), i.e. (246,248,250,241,96,128,160,16). Each component has (D,Q)=(5,32), so their products are (25,1024). The replay checks full masks and both reads. In the exact gap-three criterion, the canonical precedence order is itself a Hamilton path iff no gap equals three and the module's canonical order is itself Hamiltonian: each gap-three swap puts adjacent original I sites backwards relative to their physical forward spine edge. All prefix→module→apex junctions are already forward.

Replay: `node precedence-binary-chart.js`.


### Strongification preserving any prescribed modulus
For any integer modulus m=2^a r, r odd, take r set-bit positions q+jg for j=0,…,r−1, with q≥a and g≥max(4,a+1) a multiple of the order of 2 modulo r (take order 1 at r=1). Their sum B is divisible by m. At the odd modulus each gap summand in E is 1/2, so E=binom(r,2)/2=0 modulo r; at 2^a each gap weight vanishes since g≥a+1. Thus the strong chart satisfies H≡D mod m exactly and preserves transitive (or physically Hamiltonian) precedence whenever the module has it. This transports an existing modular class without needing an inverse unit correction.

## A strongification that preserves the residue at any modulus

For any modulus \(m=2^a r\) with \(r\) odd, choose \(r\) set bits at \(q+jg\), \(0\le j<r\), where \(q\ge a\), \(g\ge\max(4,a+1)\), and \(2^g\equiv1\bmod r\). For \(r=1\) the last condition is empty. Then \(B\equiv0\bmod m\). The gaps give \(E\equiv0\bmod 2^a\); modulo \(r\), every pair contributes the inverse of two, so \(E\equiv\binom r2/2\equiv0\). Thus appending this collar to any module with transitive precedence preserves its height modulo \(m\) and makes it strong, rigid and still transitive in precedence. Pad an original module of order less than three with singleton factors before applying the collar. This assertion requires no surjectivity claim for the input module at primes outside the original three.
