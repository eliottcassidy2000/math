# Pollock's five: Petersen edges, golden cancellation, and what the center forgets

2026-11-26. **Proved:** simultaneous local lifting for an entire bounded anchor stencil; a geometric description of all linear cubic-cancellation charts; nonlinear integer families outside the rational charts; a density-zero barrier, including a growing range of cubic defects. **Finite exact:** the largest prime admitting a cubic-balanced representation is 421. None of these restricts unrestricted Pollock representations by the same amount; [Pollock's five-tetrahedral conjecture](https://conjectures.io/problems/pollocksconjecture-pollock-tetrahedral) stays open.

Use $F(x)=(x^3-x)/6$ for integer $x$. The positive-centered convention is $x\ge1$, so $F(1)=0$ and $T_k=F(k+1)$. Padding with $F(1)$ handles fewer than five summands.

## 1. The old defect stencil cannot be excluded by congruences alone

The [June defect-pair descent](../../../07-reflections/pollock-tetrahedral-defect-pair-descent.md) starts with $N=T_k+r$ and looks for a four-sum at $r$ or $r+\binom{k+1}{2}$. The [Sierpinski/carry scout](../../../07-reflections/pollock-sierpinski-carry-and-defect-pairs.md) proposed upgrading marginal residue tests to joint ones. Here is the exact limit: **every residue modulo every modulus is a sum of two tetrahedra. Moreover an entire prescribed finite run of remainders can use the same three contacts and change only the fourth by unit index steps.** This is a statement about residue lifts, not lifts within the size of the target.

For a prime $p\ge5$ put $f(x)=x^3-x$ over $\mathbb F_p$. Its range has cardinality
\[
|f(\mathbb F_p)|=(2p+\chi_p(-3))/3,
\]
and every value has a regular preimage. Indeed a critical root of $f(x)-v$ is double (it is not zero), and the remaining root is distinct and simple. With $t=\#\{x:3x^2=1\}$ the number of two-point fibers is exactly $t$. An ordered collision of two distinct roots satisfies $x^2+xy+y^2=1$. This conic has $p-\chi_p(-3)$ points: write $(2x+y)^2+3y^2=4$ and count either in the split product or via the norm of $\mathbb F_{p^2}^\times$. Its diagonal contributes $t$. With $n_3$ the number of three-point fibers, $2t+6n_3=p-\chi_p(-3)-t$ and the image is $p-t-2n_3$, giving the formula. The image has more than half the residues, so any $6r$ is a sum of two image values. Lift one regular root while holding the second fixed. For $p=2,3$, a single summand suffices: Hensel-lift $f(x)=6r$ on the branch $x=0\bmod p$, working modulo $p^{e+1}$ to get $F(x)=r\bmod p^e$. Take the second summand zero. CRT on the two indices uses modulus $p^e$ for $p\ge5$ and $p^{e+1}$ for $p=2,3$. Positive lifts exist by adding multiples of $6M$.

For any $M,k,r$, fixed contact $g\ge1$, and stencil $j=0,\ldots,J\le k$, choose two common indices $x,y$ representing $r+T_k-F(g)\bmod M$. Use $z_j=z_0+j$ with $z_0\equiv-(k+1)\pmod{6M}$. Since $F$ is odd,
\[
F(x)+F(y)+F(g)+F(z_j)\equiv r+T_k-T_{k-j}\pmod M
\]
holds on the entire stencil. The three common contacts really are the same, and their index phases lift simultaneously even across different prime powers. In fact they can all be made positive and under a shared index cap when the cap leaves one complete period (plus $J$ steps for $z$). What is not controlled is the exact quotient after division by $M$, or the available positive index window when the modulus must grow. Thus the observed absence of long *exact* defect pairs cannot follow from a local obstruction to these four-sums, even with this strong common-contact prescription. Finite congruence obstructions that also encode height, or global arguments about true defects, are not ruled out.

This is the precise reverse transfer to the [five-core Chan–Pak shell notebook](../five_core_web_20261111/fivefold/DYADIC_SHELLS.md): even the correct joint phase of a whole stencil can leave the absolute target unknown. A better invariant is a signed quotient or a height condition on a common contact, not another marginal coverage table.

## 2. An actual pentagon, an actual Petersen graph

For five centered inputs $x_i=c+y_i$ with $\sum y_i=0$, exact expansion gives
\[
\sum_i F(x_i)=5F(c)+\frac c2\sum_i y_i^2+\frac16\sum_i y_i^3.
\]
So complete cubic cancellation is precisely the **Clebsch diagonal cubic** $\sum y_i=\sum y_i^3=0$ in the sum-zero projective three-space of the five coordinates. This is the quadratic-carrier strategy of the [tournament cycle-spectrum notebook](../../../01-canon/theorems/THM-498-cycle-spectrum-pollock-completeness-hierarchy.md) applied one degree higher: it works only on a cubic constraint, not on the whole score-like hyperplane.

The line classification is classical (see [Huybrechts, Example 2.15](https://www.math.uni-bonn.de/people/huybrech/Notes.pdf)); the following gives its Pollock and graph interpretations. There are fifteen rational lines, all permutations of $(a,-a,b,-b,0)$. To mark them, put a vertex at each of the ten projective points $[e_i-e_j]$ indexed by unordered pairs of coordinate positions. The rational line through two such marks exists exactly when their supports are disjoint: disjoint supports give the displayed form, while for a one-coordinate overlap the mixed cubic term does not vanish identically. Therefore the rational lines are **the fifteen edges of the Petersen graph** $KG(5,2)$ on those ten marks. This identifies the marked incidence, not every intersection among all 27 lines.

For the other twelve, order the five labeled positions around an unoriented pentagon $C$. Take the span of the sine and cosine modes of frequency one around $C$. The sum of the coordinates vanishes, and the cubic sum does too: a sum of three frequencies from $\{+1,-1\}$ is never zero modulo five. These are twelve different real planes. Their orthogonal projectors have off-diagonal entries $(2/5)\cos(2\pi(i-j)/5)$ in cyclic positions, whose signs recover adjacency in $C$. These planes have no nonzero rational vector: their eigenvalue for the integer adjacency matrix of $C$ is $(\sqrt5-1)/2$. Galois conjugation gives the other eigenplane, the first mode for the complementary pentagon; there are six such conjugate pairs. The rational planes and pentagon planes thus supply 27 distinct lines. The cubic is smooth, since a singular point on the hyperplane would have all $y_i^2$ equal and an odd number of signed equal coordinates cannot sum to zero unless all vanish. The classical 27-line theorem gives exhaustion.

This is the [fivefold boundary and conjugate projection](../fivefold_boundary_20261018g.md) with an exact cubic meaning. On the rational side, keeping one coordinate fixed at zero restricts the marked graph to three possible pairings on the other four positions, exactly the [tetrahedral pair charts](../tetrahedral_pair_refinement_20261004.md). On the real side the golden chart cancels perfectly but has no nontrivial integral contact. In particular every two-dimensional **rational linear offset space** that cancels the cubic gives only the paired carrier
\[
5F(c)+c(a^2+b^2),\qquad c\ge 1+\max(|a|,|b|).
\]
A rational one-dimensional space on the cubic but off those lines is also allowed; linear chart classification does not exhaust integer points on the surface. Nor is either 27-line symmetry or this Petersen incidence an identification of graph edges with tournament arc responses.

### A nonlinear chart that escapes all fifteen paired carriers

Take $u=(1,-1,t,-t,0)$ and $v=(1,1,-1,0,-1)$. On the chord $\lambda u+\mu v$ the mixed terms are $3\lambda\mu(\lambda(2-t^2)+\mu t)$. Its third intersection, choosing $\lambda=t$, $\mu=t^2-2$, is
\[
y(t)=((t-1)(t+2),(t-2)(t+1),2,-t^2,2-t^2).
\]
Directly its first and third moments vanish and its half square norm is $Q(t)=2t^4-5t^2+8$. For every integer $c\ge1-\min_i y_i(t)$,
\[
5F(c)+cQ(t)=\sum_i F(c+y_i(t))
\]
is an actual nonnegative-index construction. If $t\equiv4\bmod121$, then $Q(t)\equiv77\bmod121$, hence $v_{11}(Q(t))=1$. Since 11 is inert in the Gaussian two-square problem ($11\equiv3\bmod4$), no rational paired chart with integer offsets gives that same correction at this center. This uses the **actual prime factor**, not the parity of the total number of factors. For a small example halve the $t=4$ vector, obtaining $(9,5,1,-8,-7)$ of half norm 110. Center $c=9$ gives $1590=F(18)+F(14)+F(10)+F(1)+F(2)$, where 110 cannot be the paired correction. The claim is at the same center, not uniqueness of the target's representations.

## 3. A global obstruction to forcing exact or near cubic balance

For unrestricted positive integer coordinates $x_1,\dots,x_5$ let $S=\sum x_i$, $Q_2=\sum x_i^2$, and define the integer cubic defect
\[
D=\frac15\sum_i(5x_i-S)^3=25\sum_i x_i^3-15S Q_2+2S^3 .
\]
No integrality of the mean has been assumed. For $N=\sum F(x_i)$,
\[
150N-D=S(15Q_2-2S^2-25),\qquad D\equiv2S^3\pmod5.
\]
This is the missing quotient. At $D=0$, $5\mid S$, so $c=S/5$ is an integer and $c\mid6N$ by the centered expansion. If $N$ is a prime, $c<N$ by convexity ($N\ge5F(c)$); hence $c\mid6$. The complete search of nondecreasing positive tuples of sum $5c$ for $c=1,2,3,6$ gives maximum prime **421**, attained at $(1,2,6,10,11)$. Thus large primes force a genuinely unbalanced choice even if all nonlinear cancellation charts are permitted.

The obstruction is much broader than primes. From positivity and Jensen,
\[
(6N)^{1/3}\le S,\qquad S^3/25-S\le6N,
\]
so $S\le(300N)^{1/3}$ when $S\ge8$. Consequently for $N\asymp X$, fixed $D$ puts a divisor $S\asymp X^{1/3}$ in $150N-D$. Apply Ford's [divisor-in-a-dyadic-interval theorem, (1.1)](https://www.ford126.web.illinois.edu/wwwpapers/hxy2y_2012.pdf), with
\[
L(X)=(\log X)^\delta(\log\log X)^{3/2},\qquad
\delta=1-(1+\log\log2)/\log2 .
\]
For $0\le B\le\sqrt X$, the number of $N\le X$ that have any five-term representation with $|D|\le B$ is
\[
O\big(\sqrt X+(2\lfloor B\rfloor+1)X/L(X)\big).
\]
To see the uniformity, discard $N<\sqrt X$, split the rest into descending dyadic intervals $Z/2<N\le Z$, and fix $D$. The positive integers $150N-D\le151Z$ are distinct, and their eligible divisor interval is covered by a fixed number of Ford intervals $(y,2y]$ with $y\asymp Z^{1/3}$. This gives $O(Z/L(Z))$ for that $D$, then sum the geometric series and the possible $D$. This upper bound includes all unmarked permutations and centers. In particular exact balance has density zero, and even a fixed finite defect library, or a slowly growing band $B=o(L(X))$, cannot hit almost every target. This is a no-go theorem for a particular proof strategy, **not** a claim that these targets lack five-term representations.

A constructive next target is an exact nonzero-defect selector: choose $S$ and prescribe a feasible $D\equiv150N\bmod S$, then solve the resulting joint moment system for positive indices. The residue lift in §1 does not select that quotient; the paired/taxicab charts in §2 provide honest starts but cannot be the whole bank. In the language of the old weighted-cut sessions, both a desired response and the source configuration have to survive the same move.

Replay: `python replay.py`. It checks regular image counts, prime-power lifting and CRT, marked Petersen/pentagon counts, nonlinear identities, the finite prime maximum, and the defect factorization.
