# Exact perfect midpoints with three transitive ports

**PROVED.** For every integer exponent \(p\ge5\), let \(E_p=2^{p-1}(2^p-1)\) and \(v=2^{p-3}\). There are two strong tournaments on \(2p\) vertices differing by one arc with Hamilton path counts
\[
 E_p\pm(6v^2+2v+1).
\]
The two shorter port widths can be *any* positive split of \(p-3\). At exponent 61 the halfgap is 498460498419343452915388399908552705. A finite classification below proves that this is the only possible halfgap at this exponent for three cap-two ports on a strong source of order at most nine when the flipped arc has two ordinary endpoints. At order ten a second, wider endpoint pair appears. The original halfgap stays minimal throughout the same architecture through order ten, including arbitrary target-dependent widths.

## A literal one-arc pair

Number the nine quotient vertices from zero. Use the outgoing bitmasks
\[
 (462,508,504,368,481,449,128,264,64)
\]
and reverse \(5\to8\). At the ordered ports \((1,2,6)\), substitute transitive tournaments of widths \((u,w,z)\), all ordinary vertices being singletons. The source and flipped quotient are strong. For example both contain the directed Hamilton cycle \(0,1,2,4,8,6,7,3,5,0\). The path-count formulas in \(X=2^u,Y=2^w,Z=2^z\) are
\[
 H_-=-1-5XY-Z+26XYZ,\qquad H_+=1-XY-Z+38XYZ .
\]
The three marked vertices each visit at most twice in compressed paths. A single visit to a transitive port has one filling, and two visits have \(2^u-2\) fillings for width \(u\). The coefficients follow by enumerating compressed words oriented by the displayed masks with every ordinary vertex used once and each port used once or twice. Ports 1 and 2 are consecutive elements of a single transitive module in the explicit source, which explains why only their product \(XY\) occurs. [Independent compressed-word method](../chanpak_support_birth_grammar_20261122/check_perfect_port.py)

Set \(u+w=z=p-3\). Taking the average gives \(32v^2-4v=E_p\), and taking half the difference gives \(6v^2+2v+1\). The order is \(6+2(p-3)=2p\). Neither the construction nor the classification needs \(2^p-1\) to be prime.

In fact this width rule is **necessary for this pair at every exponent**. Write \(A=XY=2^a\), \(B=Z=2^b\). If \(a\ne b\) the 2-adic valuation of \(32AB-3A-B\) is \(\min(a,b)\), since the term of least valuation is unique. This would force \(\min(a,b)=p-1\). Then \(AB\ge2^{2p-1}\) and, because \(A,B\ge1\), \(3A+B\le4AB\); the midpoint is already too large. At equality of the widths the midpoint is \(2^{a+2}(2^{a+3}-1)\), forcing \(a=p-3\). The 2-adic target thus determines the two effective widths without a search.

## Exhausting arbitrary widths at exponent 61

**FINITE-EXACT.** Consider every strong unlabeled tournament of orders five through ten. A configuration consists of three marked vertices that have at most two visits in every compressed spanning word, and one arc between two ordinary vertices whose reversal preserves all three caps. The flipped quotient may be reducible. Ports may have width one. Orders below five have too few ordinary vertices for this move.

The two-run condition has an elementary exact test. At a marked vertex \(p\), each interval between two visits contains an arc from an outneighbor of \(p\) to an inneighbor of \(p\), since the interval begins with the former and ends with the latter. Record the set of ordinary endpoints of every such arc. A return with no ordinary endpoints can be repeated. Two returns with disjoint sets of ordinary endpoints can be concatenated at \(p\), giving three visits with all ordinary symbols still distinct. Insert any missing quotient symbols by the ordinary insertion argument for paths in a tournament. Conversely a compressed path visiting \(p\) three times exhibits two returns with disjoint ordinary supports. Thus for each mark all the supports must be nonempty and pairwise intersect.

| Order | Strong sources | Marked triples passing caps | Ordinary reversals preserving caps | Distinct ordered tensor pairs |
|---|---:|---:|---:|---:|
| 5 | 6 | 15 | 15 | 15 |
| 6 | 35 | 84 | 220 | 169 |
| 7 | 353 | 458 | 2,232 | 1,476 |
| 8 | 6,008 | 2,674 | 21,130 | 13,603 |
| 9 | 178,133 | 18,649 | 225,567 | 154,226 |
| 10 | 9,355,949 | 180,852 | 3,252,636 | 141,212 retained by the exact filters below |

For each pair of cubical path-count polynomials solve their sum \(P(X,Y,Z)=2E_{61}\) exactly. Widths require no heuristic cutoff. In a strong source a triangle goes through each port, giving a compressed path with that port used twice and every other port once by inserting the missing symbols. The source height is at least \(2^w-2\) for each port width \(w\). Since \(2^{122}-2>2E_{61}\), every width is at most 121. For each possible sorted assignment of the port widths \(x\le y\le z\), \(2^x\mid2E_{61}-P(0,0,0)\). After substituting \(X=2^x\), write \(P=C+AY+BZ+DYZ\). Then \(2^y\mid2E_{61}-C\). For each remaining candidate,
\[
  2^z=\frac{2E_{61}-C-A2^y}{B+D2^y}.
\]
The denominator is positive by the same two-visit word for the last port. Test for an exact positive integral power of two and the required order, breaking ties by port index. This calculation searches every width triple, including target-dependent, unbalanced and constant widths.

Orders five through eight give no midpoint. Order nine gives 456 pairs of an ordered tensor pair and a width tuple; the table deduplicates tensor pairs separately for each order. After forgetting port order and the direction of the flip, just two polynomial pairs survive. Writing \(A\) for the product of two port powers and \(B\) for the remaining power, their forms are
\[
 (-1-5A-B+26AB,\ 1-A-B+38AB),\qquad
 (-1-A-5B+26AB,\ 1-A-B+38AB).
\]
All 456 have \(A=B=2^{58}\), order 122 and the two heights
\[
 2159995493150488291739304216045289471,\qquad
 3156916489989175197570081015862394881.
\]
The width multisets are exactly \(\{r,58-r,58\}\), \(1\le r\le29\). This rules out an adjacent-unit straddle throughout the stated universe and shows exactly where a nominal three-port rule collapses to two variables. The enumeration measures distinct tensors with widths, not the number of isomorphism classes.

### Order ten gives exactly one further pair

Of the 3,252,636 physical reversals at order ten, 94,560 have opposite endpoint constant coefficients, giving 58,158 distinct ordered tensor pairs. Their unrestricted width search yields 904 certificates, all giving the same numerical endpoints as order nine. Among nonzero-constant pairs the exact filters in the next section retain 83,054 distinct pairs. These yield six certificates with widths \(\{1,57,58\}\), all giving the two heights
\[
 1246151246048358630126743178633543679,\qquad
 4070760737091304859182642053274140673.
\]
This is the complete numerical list together with the order-nine pair. In particular the second pair is wider and the minimum halfgap is unchanged after examining all strong order-ten quotients.

The second pair has a second uniform construction. Take the strong quotient with outgoing masks \((922,1020,1017,1008,608,961,897,784,16,256)\), mark vertices \((1,2,3)\) with widths \((p-3,1,p-4)\), and flip \(8\to4\). The flip makes vertex 8 a sink. With \(X=2^{p-3}\), \(Z=2^{p-4}\), the two compressed path polynomials after inserting the width-one middle port are
\[
 1+19X-49Z+98XZ,\qquad -1+5X-15Z+30XZ.
\]
At \(Z=X/2\) they give \(E_p\pm(17X^2-3X/2+1)\) on \(2p+1\) vertices for every \(p\ge5\). Both known midpoint families at 61 are therefore reductions to two varying widths, even with three available. The surviving order-ten tensors with opposite constants come in four types up to port permutation and reversal, and the six nonzero-constant successes come in three types. Orders eleven and up were not included.

## Why large exponents require a vanishing constant

**PROVED for any fixed multiaffine coefficient tensor.** Suppose a fixed integer multiaffine polynomial in finitely many port powers takes the value \(2^p(2^p-1)\) for unbounded exponents. Along a subsequence, order the widths increasingly and pass to a subsequence on which each width is constant or tends to infinity. Freeze the constant widths. The remaining polynomial must have zero constant term. Indeed its other monomials are divisible by the power of two of the least diverging width, while the target is divisible by \(2^p\). Both valuations tend to infinity. At least two widths must diverge: with only one the target would be a fixed coefficient times a power of two, whose odd part is bounded. A bounded schedule on a single short port can therefore help only through exact cancellation in the remaining constant. This argument applies even to schedules chosen adaptively from each integer target.

A finite version accelerates the order-ten search. Write \(P\) for the sum of the two endpoint polynomials and assume its coefficients have absolute value less than \(2^{14}\). For a nonzero constant term, the smallest width \(x\) is at most 13. If the next constant \(C=P(0,0,0)+P_X(0,0,0)2^x\) is nonzero, its absolute value is below \(2^{28}\), forcing \(y\le27\). After substituting both widths the remaining polynomial has form \(K+D2^z\), with \(|K|,|D|<2^{56}\). Matching \(2^{61}(2^{61}-1)\) forces \(z\ge62\), which would give \(K\equiv2^{61}\pmod{2^{62}}\), an impossibility. Thus a nonzero initial constant must cancel exactly when the shortest port is frozen. The enumerator's optional `small` pass keeps all such cancellations, and retains every nonzero-constant pair exceeding the coefficient bound for a direct check.

## The full width space is a union of affine lattice slices

**PROVED.** For any fixed finite set of paired cap-two tensors, all choices of positive widths and exponents realizing a specified perfect midpoint form an effectively computable finite union of affine lattice slices in the positive orthant. The same statement holds for the joint midpoint and halfgap-one problem. More generally the prescribed heights may be any finite integer sums of powers of two whose exponents are integer affine functions of the parameters and nonnegative in the domain. This describes the actual width vectors as well as the projected exponents.

Expand either required equality as \(\sum_{i=1}^m c_i 2^{L_i}=0\), deleting zero coefficients. Partition a solution into minimal zero subsums. Write \(K=\lceil\log_2(1+\sum_i|c_i|)\rceil\). In a minimal subsum, two consecutive ordered exponents cannot differ by more than \(K\). If they did, the absolute value of the sum of all terms before the gap would be strictly less than the smallest power of two after it. Divisibility of that partial sum by the latter power would then make it zero, violating minimality. Thus all exponent differences within a subsum have absolute value at most \((m-1)K\). Enumerate the set partitions and these bounded integer offsets, retaining just the offset lists whose normalized exponential sums vanish. Each choice imposes only linear equalities \(L_i-L_j=d_{ij}\), and each solution of those equalities gives a genuine zero sum. Intersecting the two finite families handles simultaneous endpoint constraints.

Projecting onto the exponent gives an effectively ultimately periodic set by integer linear elimination in Presburger arithmetic. This complements the repo’s [context-free exponent-language theorem](../chanpak_decimal_kernel_20261122/EXPONENT-LANGUAGES.md) and its [unit-equation argument for different recurrence roots](../chanpak_tiling_viewpoint_20261123/TARGET.md). Here the common base permits an elementary, explicit search for the relative exponents, and the target-selected schedules need satisfy no regularity condition. It supplies a decision procedure for each fixed list of sources. A proof of Chan–Pak would still require a family covering arbitrary admissible integer targets.

To reproduce, run in this directory with nauty's `gentourng` and a C++17 compiler:

```sh
python3 -B center_three.py --witness
c++ -O3 -std=c++17 center_three.cpp -o /tmp/chanpak-center-three
for n in 5 6 7 8 9; do
  gentourng -cq "$n" | /tmp/chanpak-center-three "$n" all > "/tmp/chanpak-center-three-$n.txt"
done
python3 -B center_three.py --arbitrary /tmp/chanpak-center-three-[5-9].txt
gentourng -cq 10 | /tmp/chanpak-center-three 10 > /tmp/chanpak-center-three-ten-zero.txt
gentourng -cq 10 | /tmp/chanpak-center-three 10 small > /tmp/chanpak-center-three-ten-small.txt
python3 -B center_three.py --arbitrary /tmp/chanpak-center-three-ten-*.txt
```

[Enumerator](center_three.cpp) · [Exact width solver](center_three.py)
