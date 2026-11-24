# 19 as a Keith hit, an Ulam miss, and a Kaprekar fork

**PROVED** and **FINITE-EXACT**, scoped as below. This merges the
[Collatz-27 filler result](FILLER-COLLATZ-27.md) with the repo's
[Kaprekar decimal packet](../chanpak_kaprekar_decimal_20261123/README.md),
which records the physical \(6173\leftrightarrow6175\) flip and
separates the raw Kaprekar map from its normalized gap clock.
There is no transfer of an asymptotic property of Ulam integers to
ordinary Hamilton counts.

## The three integer rules disagree at 19

The supplied sequence is \(U(1,2)\): after \(1,2\), admit an integer iff
it has exactly one unordered representation as the sum of two distinct
already admitted numbers. The prefix through 48 is reproduced by the
replay. In particular \(19=1+18=3+16=6+13=8+11\) has **four** witnesses
and is rejected. See [Steinerberger's Ulam study](https://arxiv.org/abs/1507.00267)
for the distinct-summand rule and the experimental angular signal.

By contrast, under [Keith's repfigit rule](https://cadaeic.net/keithnum.htm),
keep the ordered decimal digits as seeds and sum the previous
digit-width terms: \(19\) returns from \((1,9)\) as \(1,9,10,19\).
All two-digit decimal Keith numbers are \(14,19,28,47,61,75\).
With Fibonacci numbers \(F_1=F_2=1\), a two-digit number \(10a+b\)
matches term \(k\ge3\) (one-based) precisely if
\((10-F_{k-2})a=(F_{k-1}-1)b\). For \(k=3\) there is no solution,
for \(k=4,\ldots,8\) this gives respectively \(19;\ 14,28;\ 47;\ 75;\ 61\),
and for \(k\ge9\) the left side is negative. In the Ulam sequence
28 and 47 pass the uniqueness test; 61, like 19, fails it (three
representations for 61). Thus neither rule determines the other.

The actual marked tournament witnesses in [the two-edit deck](two-port-witnesses.js)
send the three two-digit **prime** Keith numbers along
\(19\to77\to133,\ 47\to119\to329,\ 61\to93\to427\).
The inherited witness heights are used here, not a claim that the
digit rules are physical edits.

| Keith start | Last two digits of its upper endpoint | First width-two Kaprekar images of those two blocks |
|---|---|---|
| 19 | 33 | 72 and 00 |
| 47 | 29 | 27 and 63 |
| 61 | 27 | 45 and 45 |

In particular \(19\to72\to45\) under the padded decimal subtraction
enters the nonzero 5-cycle \(09\to81\to63\to27\to45\to09\).
But the paired two-digit block 33 collapses to zero immediately.
For 61 the two clock states instead merge after one step.
These are comparisons of **suffix blocks**, not assertions about
Kaprekar on full-width large counts.

## Why the 19 fork has a repdigit on just one end

In base \(B\), the only two-digit Keith number that hits at the
second generated term is \((1,B-1)_B=2B-1\): the equation
\(Ba+b=a+2b\) forces \(b=(B-1)a\), hence \(a=1\).
For an integer multiplier \(m\ge2\) and \(B\ge m\), the last two digits
of \(m(2B-1)\) are \((2m-1\bmod B,\ B-m)\). They coincide iff
\(B\mid3m-1\). With \(m=7\) the only such bases are \(10,20\).
Of these the starting Keith values are \(19,39\), with 19 the only
prime. This explains the 19/33 pair in the sevenfold decimal chart
rather than merely observing it.

There is also a complementary **carry-free** family. For every
integer \(b\ge1\), put \(B=7b+5\) and \(h=8b+5=(1,b)_B\).
Then \(7h=(7,7b)_B\) is still two-digit, and the term
\(5a+8c\) from seeds \((a,c)\) hits both \(h\) and \(7h\) in the same
base and at the same index. The repo's
[A-frame transitive family](TWO-APEX-SCALABLE.md) has \(h=6\cdot2^t+5\),
so \(b=3\cdot2^{t-2}\) for \(t\ge2\) supplies infinitely many literal
two-edit pairs of this common-base form. The base varies and the
starts are not asserted prime. In fact every sufficiently large
\(h\equiv5\bmod8\) qualifies for this variable-base interpretation,
so by the dyadic result one quarter of the dimension-19 bank, and
one quarter of the sharper sieve below, admit it.

A raw width-two Kaprekar move is much stricter. For an eligible prime
\(p\ge5,\ p\ne7\), the only \(p\to7p\) move is the padded one-digit
case \(p=5,B=8\). Indeed for \(p=aB+b\) with \(a\ge1\), the equation
\((B-1)|b-a|=7p\) forces \(b>a\), and writing \(\Delta=b-8a\)
gives \(B=8+63a/\Delta\), \(1\le\Delta\le7\) (use \(b<B\)), and
\(p=(7a+\Delta)(9a+\Delta)/\Delta\). Since \(\Delta\) divides the
numerator and both factors exceed \(\Delta\), distribute its factors
between the two numerators to get a nontrivial integer factorization.
In reverse, \(7p\to p\) requires \(B=p+1\) and gap 1 by primality.
For \(p\ge11\), the digits of \(7p\) are \((6,B-7)\), so
\(|12-p|=1\) selects \(p=11,13\). At \(p=5,B=6\) the digits coincide.
Thus \(77\to11\) in base 12 and \(91\to13\) in base 14 are the
only reverse moves, both physically liftable via the inherited
two-edit witnesses in the opposite direction. The former endpoints
11 and 77 are both Ulam; the latter endpoints 13 and 91 are both
Keith in a **different** common base, 12, with seeds \((1,1),(7,7)\).
No common-base raw move here has both endpoints two-digit Keith:
its smaller endpoint has just one digit in the Kaprekar radix.

## Ulam ambiguity of filler splits can be read from physical heights

Multiplicative independence of the third readings in FILLER-JOINT
makes \(\log R(e)=\sum_j e_j\log R_j\) a faithful encoding of the
nonnegative lattice. For a nonzero multiplicity vector \(e\),
the number of unordered unequal splits into two nonempty fillers is
\[
 {\cal N}(e)=\frac{\prod_j(e_j+1)-\mathbf1_{e\ {\rm all\ even}}}{2}-1.
\]
Here the object counted is a pair of multiplicity profiles, not a
pair of labelled tournaments. The indicator deletes the equal pair
and the minus one deletes the empty split. Exactly one split requires
two different generators once each, or a single generator to power
3 or 4. Exactly four splits, the analogue of 19's four Ulam witnesses,
requires a single generator to power 9 or 10, or exponents \(4,1\)
on two different generators.

The genuine recursive Ulam rule tests only splits into previously
accepted profiles. Starting from all elementary vectors and weighting by
\(\log R\), it is the basis-vector Ulam set transported from the
lattice, as in [Kravitz–Steinerberger](https://arxiv.org/abs/1705.01883).
On the squarefree cube, only supports of size one or two pass:
triples have three singleton/pair splits, quadruples three pair/pair
splits, and higher supports have none by induction. Off the cube, the
static split count above does **not** decide recursive Ulam membership.

There is an exact physical proxy for balanced splits. Fix \(0<\theta<1/2\)
and let \(X=h(e)\) be large. For positive profiles \(v,w\) with
\(h(v),h(w)\ge X^\theta\),
\[
 v+w=e\quad\Longleftrightarrow\quad
 \left| \frac{h(v)h(w)}{20X}-1\right|<X^{-\theta/8}.
\]
Indeed \(h=20R(1+O(R^{-1/7}))\), so a true product split makes the
relative error \(O(X^{-\theta/7})\). Conversely a near product puts
\(R(v)R(w)\asymp R(e)\asymp X\), with all exponent differences
\(O(\log X)\). The two scalar log-products cannot be that close
if distinct, by the same fixed-log lower bound as in FILLER-JOINT:
their separation is \(\gg(\log X)^{-C}\), whereas the displayed
errors are power-decaying. Then independence gives the exact split.
Thus keeping real scale and enough precision can restore uniqueness
that would be lost by retaining only a decimal suffix. This is not
an equality \(h(v)h(w)=20h(e)\), nor a test on arbitrary tournaments.


## 47 closes a mixed Kaprekar–Collatz loop

The number 47 is both Ulam (uniquely \(11+36\)) and Keith
(\(4,7,11,18,29,47\)). It also closes a hybrid loop:
\[
 47 \xrightarrow{\text{Kaprekar}} 27
 \xrightarrow{\text{seven ordinary Collatz steps}}47,
\]
via \(27,82,41,124,62,31,94,47\). In particular this connects 19's
rejected Ulam profile to 47's accepted profile and to the earlier
Collatz-27 experiment, without identifying any of the three rules.

Extend the hybrid map by replacing **only the last two decimal digits**
in the Kaprekar step. For \(m\ge1\) and \(z\ge0\), a start
\[
 h=47+25\cdot16^m z
 \quad\longmapsto\quad
 47+25\cdot27^m z
\]
follows \(m\) consecutive traversals of the mixed loop. Indeed after
\(j\) traversals the displacement from 47 is
\(25\cdot27^j16^{m-j}z\). For \(j<m\) it is divisible by 100,
preserving the Kaprekar suffix, and by 16, preserving the Collatz
parities. The seven-step affine slope is \(27/16\) (three odd steps
and four halvings). Conversely these two conditions force the
initial displacement to be a multiple of \(25\cdot16^m\) if
all \(m\) traversals are to use exactly these suffix and parity reads.
The exact repeat is at 47; all nonzero positive lifts expand.

The replay gives an explicit **865-vertex** start in the first
13-alternative support \(S_1\) below, with multiplicities
\((9,12,7,2,8,3,1,10,12,8,8,9,6)\) in its displayed order.
It has \(h\equiv47\bmod(25\cdot16^4)\) and \(h\equiv12197\bmod15015\).
The physical sevenfold square exists by the common frame. The
Kaprekar and Collatz steps are arithmetic observations on its count,
not additional physical flips. In the sharp five-prime-sieved
bank below, for any fixed \(m\) the fraction that shadows at least
\(m\) traversals tends to \(1/(20\cdot2^{4m-1})\): the decimal uniformity
projects onto the units modulo \(25\cdot16^m\).

## Additional literal overlaps

Both 69 and 483 are Ulam and are Keith in common base 61,
with digits \((1,8)\) and \((7,56)\). Moreover a common ten-vertex
square has corners \(69,341,277,483\), with the route through 341
entirely Ulam. The unique sums for its three checkpoints are
\(16+53,\ 2+339,\ 69+414\), respectively. The alternative
intermediate 277 is not Ulam. This tests branch continuation on
actual heights, rather than assuming matching endpoint labels
force matching intermediate labels.

A second square of arithmetic observations uses base 26:
\(55=(2,3)_{26}\) and \(385=(14,21)_{26}\) are Keith at the same
term \(8a+13b\). Raw Kaprekar sends them to 25 and 175, respectively:
it commutes with multiplication by seven here, because the input
digits scale without carries. The replay checks independent literal
two-flip squares for \(55\to385\) and \(25\to175\); no common host or
literal horizontal Kaprekar moves are claimed.

## Sieve before counting decimal phases: 19 drops to 12

This fixed-bank theorem concerns **possible prime starts**, not prime
density. The exponent-19 sieve by 3 from FILLER-COLLATZ-27 is only
the first pass. Let \(P=3\cdot5\cdot7\cdot11\cdot13\). Among starts
of the full 28-alternative bank in \([X,2X]\), the number coprime to
\(P\) is asymptotic to
\[
 (c_{S_1}+c_{S_2})(\log X)^{12},\qquad
 c_S=\frac{\log2}{12!\prod_{j\in S}\log R_j}.
\]
In particular prime starts in this fixed bank number at most
\(O((\log X)^{12})\) per shell. Label the first/second alternatives
of zero-based type \(i\) by \(ia,ib\). There are exactly two
maximal-size survivor supports:
\[
\begin{split}
 S_1&=\{0a,0b,1a,2a,2b,3a,3b,6a,6b,7a,7b,11a,11b\},\\
 S_2&=\{0a,1a,2a,2b,3a,3b,6a,6b,7a,7b,11a,13a,13b\}.
\end{split}
\]
Here a support is the set of alternatives used at least once.
A term \(cZ\) of \(h=11D+42Q+20R\) vanishes modulo a prime if that
prime divides \(c\) or any included factor's corresponding reading.
Thus supports for which all three terms vanish modulo any one of the
five primes must be thrown away. This is a downward-closed test on
supports, not a primality test. The short exhaustive certificate
branches on the three coordinate choices forced by each prime:
to avoid that prime's structural zero, retain one nonzero summand,
deleting all factors that kill it. Every surviving support is contained
in a terminal branch. The replay returns exactly \(S_1,S_2\) at size
13 and no larger support. On either full support exactly one summand
actually survives for each prime, in the order \(D,Q,R,Q,(Q\text{ or }R)\).
There can be no accidental cancellation there. The other supports
have size at most 12 and yield \(O((\log X)^{11})\) starts. Each
13-support contributes the sharp thin-strip count from FILLER-JOINT
with its own 13 log weights, and the eventual injectivity there
prevents collisions even across the two supports.

For every fixed \(k\ge1\), these sieve survivors are **asymptotically
uniform over the units modulo \(10^k\)**. This same decimal law
already holds after discarding only nonunits of \(3\cdot5\), whose
shell exponent is 17. To see the latter exponent: after deleting
the eight generators killed by the 3-sieve, the two alternatives
of type 5 automatically make \(5\mid h\) by their \(D,Q\) readings.
The 18 remaining entries all have \(Q,R\) units modulo 5 and
\(4\mid Q\). On their full support \(D\equiv0\bmod5^k\) for all but
a lower-order family: if not, the multiplicity of \(0a\), with
\(5\mid D_{0a}\), is bounded in terms of \(k\); fixing it reduces
the thin-strip dimension by one. The same argument holds for both
13-supports, since both contain \(0a\).
Also \(Q\equiv0\bmod2^k\) throughout any sufficiently large shell.
Consequently the height is generically \(11D+20R\) on the 2-side
and \(42Q+20R\) on the 5-side; in particular all generic long
18-filler starts are units modulo 30.

For clarity, here is the lift proof of equal decimal fibers.
The seed coordinates are \((D,R\bmod16;\ Q,R\bmod5)\).
In each of the three supports, their finite generated group has
128 states, already generated by first alternatives of types
\(0,1,2,11\). The replay gives the four seed states per each pair of
odd height residue mod 16 and unit height residue mod 5.
Above level 16 on the binary side, powers \(4\cdot5^{k-1}\) of the
first alternatives of types 0 and 2 freeze the entire mod-\(5^k\)
side and supply independent binary lifts by the square calculation in
FILLER-COLLATZ-27. They generate \((1+16\mathbb Z)^2\) modulo \(2^k\),
and \(11D+20R\) has equal fibers there by solving for \(D\).
For the 5-side, put \(e=2^{\max(k,4)-2}\) and take type \(0a\)
to exponent \(e\), freezing the binary side while preserving all
seed phases. Write \(A=Q_{0a}^e,\ B=R_{0a}^e\).
Since \(34084^4\equiv11\bmod25\), we have \(v_5(A-1)=1\);
also \(B\equiv1\bmod5\). For different indices \(t,s\) the difference
of \(42QA^t+20RB^t\) and \(42QA^s+20RB^s\) has exact valuation
\(1+v_5(t-s)\), by lifting powers and the extra factor 5 in \(20\).
Thus the index cycle supplies each available 5-adic height lift once.
Uniformity on the finite product group follows by these equal-fiber
actions, then by counting exponent congruence classes using the
thin strip. Small \(k<4\) also follows by projection.

Consequences on the actual \(P\)-sieved bank: any one of the suffixes
19, 47, 61 holds asymptotically on \(1/40\) of its starts.
For width-two Kaprekar, the image multiplicities over these **40 unit
suffix classes** (not only the three Keith suffixes) are
\[
\begin{array}{c|rrrrrrrrrr}
\text{image}&00&09&18&27&36&45&54&63&72&81\\\hline
\text{classes}&4&7&6&5&4&4&4&3&2&1
\end{array}
\]
as counted by \(9|a-b|\). In particular the 19-versus-33 basin split
appears on an infinite physical family retaining the full exponent 12
of this sharper sieve. Neither a Keith suffix nor a Kaprekar state
says the entire physical count is a Keith number or is prime.

Replay: python3 filler-keith-ulam-check.py.
