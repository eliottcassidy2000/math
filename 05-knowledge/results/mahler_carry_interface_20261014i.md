# A finite rational-height obstruction for matching Mahler and Collatz codes

2026-10-14. **PROVED:** the finite separation bound, irrationality consequence,
and source-dependent delay family below. **FINITE-EXACT:** the declared
checker universe. The two infinite problems remain **OPEN**: this packet
proves neither existence/nonexistence of Mahler Z-numbers nor Collatz
convergence.

## 1. Inheritance and the actual map

The closest proved mechanism is [THM-2228, carry tails and integral
stabilization](../../01-canon/theorems/THM-2228-mahler-three-halves-carry-tail-and-integral-stabilization.md).
It gives separate native source coordinates and real-tail conditions. Its
denominator-19 rational pseudo-source is the hostile: arbitrarily long finite
safety is possible without an ordinary positive source at infinite depth.
The corrected near miss is **MISTAKE-588** in the maintained
[mistakes ledger](../../01-canon/MISTAKES.md): a generic giant-source runway
comparison was not an identified Mahler map. The least-used sidecar here is
the **ordinary numerator/denominator of the transported Collatz source**.
Our concept board is binary addresses / two different affine carriers /
real safety / ordinary-source stabilization / rational height / finite
periodicity.

There is already a genuine block bridge in
[THM-4469, adjacent supercritical block pairs](../../01-canon/theorems/THM-4469-mahler-bridge-adjacent-block-pairs.md).
Its smallest pair has length 10, seven odd letters, multiplier 2187/1024,
and adjacent carries 4726/4727. That theorem concerns a specified generalized
Mahler interval. It does not identify ordinary Mahler multiplication by 3/2
with a positive Collatz orbit. The present map uses the **same binary code**
in two different native encodings and proves a restriction on the resulting
source.

Write, on the 2-adic integers or rational numbers of odd denominator,

\[
M(x)=\frac{3x+(x\bmod2)}2,\qquad
T(x)=\begin{cases}x/2,&x\equiv0\pmod2,\\(3x+1)/2,&x\equiv1\pmod2.\end{cases}
\tag{1}
\]

For positive integers, $M(x)=\lceil3x/2\rceil$. For a binary word
$w=(e_0,\ldots,e_{L-1})$, the two affine formulas are

\[
M^L(x)=\frac{3^Lx+C_M(w)}{2^L},\qquad
T^L(x)=\frac{3^{|w|_1}x+C_T(w)}{2^L},                 \tag{2}
\]

where adjoining a bit $e$ updates the two carries by
$C_M(we)=3C_M(w)+e2^{|w|}$ and
$C_T(we)=3^eC_T(w)+e2^{|w|}$. In each case, the numerator's divisibility
by $2^L$ is equivalent to the whole prescribed parity word: peeling the
first parity gives an induction, since the multiplicative coefficient is
odd. Thus each code has one native residue modulo $2^L$ for each map.
These compatible finite bijections give an isometry

\[
\Psi:\mathbb Z_2\longrightarrow\mathbb Z_2,\qquad
\operatorname{code}_T(\Psi(A))=\operatorname{code}_M(A),\qquad
\Psi M=T\Psi.                                         \tag{3}
\]

The source and target are two marked 2-adic dynamical systems. The map
preserves the full ordered code and its precision; it need not preserve
positivity, ordinary integrality, rationality, size, or a stopping target.
The API `transport(A,L)` retains $A$, $L$, the common word and **both**
native residues. Its auditor requires the supplied original $A$.

The classical primary sources are Mahler's
[1968 paper, reprinted by EMS](https://ems.press/books/dms/252/4994), for
the Z-number definition, and Rozier's
[Parity sequences of the 3x+1 map on the 2-adic integers and Euclidean
embedding, INTEGERS 19 (2019), A8](https://www.ipgp.fr/~rozier/pub/Parity.pdf),
§§1–2, for the Collatz parity coordinate and rational periodic examples.
The finite bound below is derived here; no analytic distribution result is
used. The current [Mahler frontier](../../00-navigation/MAHLER-THREE-HALVES-FRONTIER-2026-08-23.md)
retains the stronger previously audited finite-safety and moving-source
boundaries.

## 2. Rational Collatz sources with no consecutive odd steps are periodic

For any rational $x=u/v$ in lowest terms, with $v>0$ odd, the orbit of $T$
stays in the fixed lattice $v^{-1}\mathbb Z$. A reduction of the displayed
fraction may shrink the denominator, but never invalidates this lattice.
On the three two-bit words that avoid 11,

\[
T^2(x)=x/4,\quad(3x+2)/4,\quad(3x+1)/4,
\qquad |T^2(x)|\le\tfrac34|x|+\tfrac12.             \tag{4}
\]

Define exact integers

\[
B=\max\{|u|,2v\},\qquad K=\left\lfloor\frac{3B+v}{2}\right\rfloor,
\qquad N=2K+1.                                       \tag{5}
\]

As long as a **finite** prefix avoids 11, its even-clock numerators have
absolute value at most $B$, by (4); the intervening numerators are at most
$K$, by one application of (1). Therefore states at times $0,\ldots,N$
belong to the $N$-element set $\{-K,\ldots,K\}/v$, provided the first $N$
bits avoid 11. Some actual states $x_i=x_j$ repeat, with
$0\le i<j\le N$. Determinism now makes the actual rational orbit periodic
from $i$ onward. This last conclusion uses a repeated actual state, not
an assumption that a second orbit continues to match it.

Consequently an all-depth no-11 rational Collatz code is eventually
periodic. This small constrained result needs neither Collatz convergence
nor the general rational Periodicity Conjecture. The exact finite observer
stops at the second bit of 11 or at a repeated state; it checks the 11
boundary before declaring a cycle, so the wraparound is retained.

## 3. A total finite separation certificate

**Theorem.** Let $A\ge1$ be an ordinary integer and $x=u/v$ a reduced
rational of odd positive denominator. Put $B,K,N$ as in (5), and
$b=\operatorname{bitlength}(A+1)$. Within the first

\[
H=2N+b                                                   \tag{6}
\]

bits, either the Mahler code of $A$ differs from the Collatz code of $x$,
or their common code contains 11. Indices begin at zero, so the witnessed
mismatch or second 1 has index strictly less than $H$.

*Proof.* Suppose the first $H$ bits agree and contain no 11. Apply the
finite-prefix argument of §2 to obtain $x_i=x_j$, $j\le N$. Write
$d=j-i$, $w$ for its period, $C=C_M(w)$ and $D=3^d-2^d$. The unique
Mahler source having the purely periodic code $w^\infty$ is

\[
z=-C/D,\qquad 0\le C\le D.                           \tag{7}
\]

The carry bound follows by comparing every bit to 1; for the all-one word
the carry is exactly $3^d-2^d$. Let $a=M^i(A)>0$. By the native-cell
bijections, the exact number of matching bits between $a$ and $z$ is

\[
h=v_2(a-z)=v_2(Da+C),                                 \tag{8}
\]

because $D$ is odd and $Da+C$ is a nonzero positive integer. The rational
Collatz orbit really has the periodic continuation from $i$, so the two
codes disagree at index $i+h$, irrespective of any eventual Mahler behavior.
Moreover $M(a)+1\le\tfrac32(a+1)$ for positive integral $a$. Hence

\[
0<Da+C\le D(a+1)<3^d(3/2)^i(A+1)<2^{2d+i+b}.
\]

It follows that $i+h<2(i+d)+b=2j+b\le H$, contradicting the assumed
finite agreement. This proves the theorem. No infinite matching hypothesis
has entered the argument. $\square$

The implementation `compare` is a total, source-bound finite consumer.
`audit_separation` authenticates the supplied two literal words, source,
horizon, first-event index and event type. For observers that already
exhibit a rational cycle, `cycle_mismatch_bound` independently computes
the exact index in (8) and compares it with literal iteration.

If the rational height satisfies $\max\{|u|,v\}\le Q$, then
$B\le2Q$, $K\le\lfloor7Q/2\rfloor$, and $N\le7Q+1$. Thus the
uniform height-budget version is

\[
H_Q=14Q+2+\operatorname{bitlength}(A+1).               \tag{9}
\]

A common no-11 prefix of length $H_Q$ excludes **every** rational source
of height at most $Q$. This is a finite exclusion of a specified rational
height, not a bound on the complexity or magnitude of an arbitrary 2-adic
source.

In particular, if $z_A=\Psi(A)$ has an all-depth no-11 code, then every
reduced $u/v\in\mathbb Z_2$ of height at most $Q$ satisfies

\[
v_2(z_A-u/v)<14Q+2+b,\qquad
|z_A-u/v|_2>2^{-(14Q+2+b)}.                            \tag{9a}
\]

This is an explicit exponential-in-height separation. It is much weaker
than a polynomial-in-height lower bound or a finite irrationality exponent;
it supplies neither of those stronger approximation conclusions.

## 4. The all-depth consequence and its missing coordinate

For a hypothetical Z-number, THM-2228 supplies a positive integer $A$
whose Mahler carries satisfy

\[
Y_n=\sum_{j\ge0}e_{n+j}(2/3)^{j+1}<1\quad\hbox{for all }n.
\tag{10}
\]

Two consecutive 1s contribute $2/3+4/9=10/9$, so strict safety implies
no 11. Therefore:

> If a positive ordinary Mahler source has an all-depth no-11 itinerary,
> its same-code Collatz source $\Psi(A)$ is an irrational 2-adic integer.
> In particular every hypothetical Mahler Z-number's integer part has
> this property.

One can also prove this consequence qualitatively: §2 would make the
rational Collatz code eventually periodic; the unique Mahler source of
the periodic tail is nonpositive by (7), impossible for a positive Mahler
iterate. Equations (6)–(9) strengthen that argument to finite receipts.

No-11 is only necessary for real safety. The periodic code 01 has no 11,
but its largest Mahler tail is $6/5$. Its Collatz source is 2, whereas its
Mahler source is $-2/5$. Thus even the familiar positive Collatz ROOT cycle
does not supply a positive Mahler candidate through this code map.
Conversely, irrational transported source alone proves neither safety
nor positive stabilization.

The companion [triangular quotient packet](recursion_carry_bridges_20261014i.md)
forgets a different coordinate: it identifies a carry word with its
complement, and restores it with one marked sheet. Here the code itself
is fully retained; ordinary rational-source membership is what changes.
Both interfaces therefore require their stated source coordinates.

## 5. Arbitrarily long strictly safe matches, with a changing source

For every $t\ge0$ put

\[
L_t=18t+8,\qquad A_t=\frac{2^{L_t}-9}{19}.             \tag{11}
\]

Since $2^{18}\equiv1\pmod{19}$ and $2^8\equiv9\pmod{19}$,
$A_t$ is a positive integer. The period 100 has Mahler source $-9/19$
and Collatz source $1/5$; the latter actual orbit is
$1/5\to4/5\to2/5\to1/5$. Now

\[
v_2(A_t+9/19)=v_2(19A_t+9)=L_t.                       \tag{12}
\]

The two maps' codes at $(A_t,1/5)$ agree for exactly $L_t$ positions,
then differ. The three infinite Mahler tails of 100 are
$18/19,8/19,12/19$, all strictly below 1. Every finite suffix of the
matching prefix is therefore strictly below $18/19$.

The first two ordinary sources are 13 and 3,532,045. This is a family
whose source increases with the requested matching horizon. It proves
that no finite comparison bound depending only on the rational observer
can replace (6). It does not give one fixed positive source with all-depth
safety. Allowing the nonpositive rational Mahler source $-9/19$ itself
does give infinite agreement with $1/5$, showing why the ordinary
positive-integer premise is essential.

## 6. Exact reproduction and next obligation

The checker covers both native cell maps for all words of lengths 0–10,
with three positive integer lifts each; all 510 nonempty periods of length
at most 8; all 193 reduced rationals with numerator -15…15 and odd
denominator 1…15, each compared with all Mahler sources 1…128; and six
members of (11), through 98 matching bits. It separately checks 1,946
exact cycle mismatch indices from (8). Malformed type/source/precision
packets are rejected before equality can identify Boolean or float aliases.

```text
python -B -X utf8 04-computation/experiments/mahler_carry_interface_20261014i.py
python -B -O -X utf8 04-computation/experiments/mahler_carry_interface_20261014i.py
```

The useful remaining obligation is to constrain the actual **irrational**
Collatz 2-adic sources selected by positive Mahler integers, while retaining
real strict safety and source stabilization. A finite-state graph or a
rational-cycle classifier alone cannot certify those predicates. The
bounded observer is still useful: it provides exact finite rational-height
exclusions with a terminating verifier and no convergence oracle.

Normal and optimized runs match the saved JSON at **65,814 exact checks**.
LF-normalized SHA256:
`ff62a36185f5d7a2c79021461c14e83548330ebdd7026e137c27f451cceaf63e`.
