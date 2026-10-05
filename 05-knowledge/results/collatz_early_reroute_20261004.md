# Reroute earlier: a smaller branch clock, new critical cancellations, and a ternary-height boundary

2026-10-04. **PROVED:** the exact least clock and complete finite applicability
criterion for the specified early-reroute grammar; two positive-density
extensions of the previous family bank; Mersenne exponent repairs; and an
infinite ternary-height obstruction. **FINITE-EXACT:** the declared independent
replays and supplied-child demonstrations. **OPEN:** a cover of every positive
source, and home certificates for the free-parameter children. No priority claim
beyond the explicitly compared repository banks.

## 1. Inheritance and what changes

The closest mechanisms are the [incoming branch-toll and graph rank](collatz_branch_toll_rank_20261004.md),
the [composed complement construction](collatz_complement_routing_20261004.md),
the [binary16 and ternary sibling bank](collatz_binary_ternary_guard_fusion_20261004.md),
and the [reset switch](checked_switch_phase19_20261004.md). We reuse their
literal odd-map identities, not an assumption of Collatz convergence.
The [sibling grammar](creative_sibling_20260925.md) already warns that pure
inverse ports do not enlarge a family that is closed under forward routes.
The gain below is consequently stated relative to named, guarded dispatch
banks. It is not a claim of new coverage beyond every earlier proof grammar.

Canonical hostiles are 7,27,703; all still escape this new selector. The
corrected near miss is to demand payment from an intermediate sibling rather
than from the final child against the original source. The least-used sidecar
is the position at which rerouting starts. Our live board is **original source /
reroute position / real branch toll / ternary precision / critical rank /
supplied child / residual cofactor**. Anchor: enlarge exact critical cancellation.
Niche: turn the choice of reroute position into a complete bounded selector.
Wildcard: repair some Mersenne classes and then test the opposite prime tower.

The previous construction followed the whole initial run and its reset before
rerouting. Starting before the reset, or at the supplied source itself, needs
fewer ternary digits. The result is an exact reroute with a smaller child; its
child remains a proof obligation unless a certificate is supplied.

## 2. A uniform construction at any initial-ones checkpoint

Let \(U(n)=\operatorname{oddpart}(3n+1)\) on positive odd integers, and

\[
S(x)=4x+1,\qquad S^k(x)=4^kx+(4^k-1)/3.
\]

Fix \(n>1\), and an actual initial string of \(r\ge0\) valuation-one steps.
Thus \(r\le v_2(n+1)-1\), and

\[
x=U^r(n)=(3/2)^r(n+1)-1.
\]

For \(k\ge1\), define the new clock

\[
D(k)=\min\{d\ge1:3^d>2^{d+2k-1}\}.
\tag{1}
\]

It begins \(2,6,9,12,16,19,23,26,30,33,\ldots\). Unlike the inherited
clock \(E(k)=\min\{d:3^d>2^{d+2k}\}\), it pays only \(4^k/2\): a final
half-child switch pays the extra factor two. Always \(E(k)-D(k)\in\{1,2\}\).

Suppose \(3^r\mid k\). Put

\[
\begin{aligned}
d&=D(k),& \ell&=r+d,& M&=3^d,\\
C_r(k)&=\frac{2^{r+1}(4^k-1)}{3^{r+1}},&
G_r(k)&=4^k-C_r(k),&
P&=2^{d+2k-1}.
\end{aligned}
\]

These are integers, and \(0<G_r(k)<4^k\). The exact source guard is

\[
4^k(n+1)\equiv C_r(k)\pmod{3^d},
\quad\text{equivalently}\quad 4^kn+G_r(k)\equiv0\pmod{3^d}.
\tag{2}
\]

On this guard define

\[
h+1=\frac{2^{d-1}\bigl(4^k(n+1)-C_r(k)\bigr)}{3^d}.
\tag{3}
\]

**PROVED.** The child is a positive odd integer, \(0<h<n\), with

\[
v_2(h+1)=r+d.
\tag{4}
\]

Writing \(a=v_2(3x+1)\), the actual words give

\[
n\xrightarrow{(1^r,a)}J
\xleftarrow{(1^{r+d-1},2,a+2k-2)}h.
\tag{5}
\]

In particular no unspecified valuation is rounded or replaced by its minimum.
The final exponent is positive because \(k\ge1\).

**Proof.** Start with the checked intermediate

\[
m=2^\ell(S^k(x)+1)/3^\ell-1,\qquad h=(m-1)/2.
\]

Formula (2) is exactly \(3^\ell\mid S^k(x)+1\); algebra gives (3).
The quotient \((S^k(x)+1)/3^\ell\) is positive and has valuation exactly one,
since \(S^k(x)\equiv1\pmod4\). Hence \(h+1\) has valuation exactly
\(\ell\), proving positivity, oddness and (4). After the child's first
\(\ell-1\) ones its value is \((S^k(x)-2)/3\). Its next valuation is exactly
two and it reaches \(S^{k-1}(x)\); the final step has valuation \(a+2k-2\)
and reaches \(J=U(x)\). This also proves every intermediate is positive.

For the size test write

\[
h+1=\rho(n+c),\quad
\rho=\frac{2^{d+2k-1}}{3^d}<1,\quad
c=1-(2/3)^{r+1}(1-4^{-k})\in(0,1).
\]

Then \(h<\rho(n+1)-1<n\). None of the child's preceding nodes can be1:
the initial ones grow and the penultimate node is \(S^{k-1}(x)>1\).
The common endpoint may be1, in which case it is the first hit. QED.

### Exact necessity and the finite selector

Consider this same grammar with any positive inverse depth \(\ell\), without
presupposing payment or divisibility of \(k\). A positive integer child smaller
than \(n\) requires \(d=\ell-r\ge D(k)>0\).
Indeed, the preceding affine formula remains valid over the rationals. If
\(\rho\ge1\), then

\[
h-n=(\rho-1)(n+1)-\rho(1-c)>-2/3,
\]

because \(1-c<2/3\) and the expression is increasing in \(\rho\ge1\).
The integer difference cannot be negative. Thus payment is equivalent to
the strict coefficient condition within legal rows of this grammar.

With \(\ell>r\), inverse integrality forces \(3^r\mid k\): the two terms
in \(2^r(S^k(x)+1)\) have ternary valuations at least \(r\) and exactly
\(v_3(k)\), using \(v_3(4^k-1)=1+v_3(k)\). A smaller valuation in the
second term cannot cancel. Conversely that divisibility makes \(C_r(k)\)
integral, and (2) is sufficient.

For fixed \(r,k\), guards at larger \(d\) are nested inside the guard at
\(D(k)\). Therefore existence of any paid depth is equivalent to the one
minimum-depth test (2). Further, positivity and (4) give

\[
2^{r+D(k)}<n+1,\qquad D(k)>2k-1.
\tag{6}
\]

Consequently the complete applicability test is finite: take
\(0\le r\le v_2(n+1)-1\), \(k\ge1\) divisible by \(3^r\), and

\[
r+D(k)\le\lfloor\log_2n\rfloor,
\qquad k\le\frac{\lfloor\log_2n\rfloor-r}{2}.
\]

Check (2) for exactly these candidates. Failure means no member of this entire
early-reroute grammar applies, including all larger inverse depths. It does
not mean no other construction or later checkpoint can work.

The isometric address is

\[
R_r(k)=-1+(2/3)^{r+1}(1-4^{-k}).
\]

On indices \(k=3^ru\), its difference valuation is
\(v_3(R_r(k)-R_r(k'))=v_3(u-u')\). Thus it retains the inherited ternary
odometer while keeping the real clock and position cost separate.

## 3. Two new arithmetic cells at the original source

Choose \(r=0,k=2,d=6\), and the binary cylinder \(n\equiv123\pmod{128}\).
The source guard is \(n\equiv273\pmod{729}\). CRT gives

\[
n_t=89211+93312t,\qquad h_t=62655+65536t,\qquad t\ge0.
\tag{7}
\]

Their exact words are \(n_t\xrightarrow{(1)}J_t\xleftarrow{(1^5,2,3)}h_t\).
Every source has first four actual valuations \(1,2,1,2\), and all four
iterates exceed the source. This excludes immediate descent, the first-reset
at-least-three switch, and all sixteen inherited binary-debt rows: the only
matching first boundary in that bank requires the fourth exponent six.
It also excludes the elementary smaller-predecessor guards \(n=2\pmod3\)
and \(n=4\pmod9\), because every source is \(3\pmod9\).

The old sibling rows \(j=0,1\) are incompatible. Its \(j=2\) guard has
depth seven and selects exactly \(t\equiv2\pmod3\). Retain the other two
parameter branches:

\[
\begin{array}{c|c}
\text{source}&\text{child}\\\hline
89211+279936t&62655+196608t\\
182523+279936t&128191+196608t
\end{array}
\tag{8}
\]

For each branch, the remaining old rows start at \(E(3)=11\), whereas its
own ternary precision is seven. Since \(E(j+1)\ge E(j)+3\), their entire
relative occupancy is at most

\[
\sum_{j\ge3}3^{7-E(j)}\le\frac{3^{-4}}{1-3^{-3}}=\frac1{78}.
\tag{9}
\]

This is a natural-density statement, not countable additivity of arbitrary
sets: old successful rows at sources at most \(X\) have only \(O(\log X)\)
possible indices by their positivity cutoff. Counting each residue class
adds only \(O(\log X)\) endpoint errors. Finite truncations and the geometric
tail prove existence and the bound after division by \(X\).

The comparison also includes the **whole previous composed schema** applicable
to these sources. Their maximal initial one-run is one, so its only relevant
indices are \(k=1,10,19,\ldots\). The first two are incompatible (mod3 and
mod9 respectively). The rest start at its precision62 for \(k=19\), and each
increment by9 raises precision by at least30. Their extra relative occupancy
is at most

\[
\varepsilon=\frac{3^{7-62}}{1-3^{-30}}.
\]

The previous schema has the same finite-height cutoff, so the same density
argument applies. Subtracting all these earlier domains produces two actual
disjoint new dispatch classes, each occupying at least
\(77/78-\varepsilon>0.98717\) of its progression. Their combined density
among positive odd integers is at least

\[
\frac{2(77/78-\varepsilon)}{64\cdot3^7}
=\frac{172212682662848571118682075}{12208653583262207315470753502976}.
\tag{10}
\]

This is new smaller-dependency coverage against those named banks. It is not
a density of completed Collatz proofs or a comparison against every rule ever
studied in the repository.

Every source in (8) is \(27\pmod{96}\), hence a critical state of the incoming
rank \(R(n)=(3^{v_2(n-1)}((n-1)/2^{v_2(n-1)})^2,v_2(n-1))\).
Both source and child have \(v_2(n-1)=1\); \(h<n\) therefore also proves
\(R(h)<R(n)\). These are all-height cancellations inside that exact critical
class. Each source is divisible by3 and has no odd predecessor of its own;
the checked common future is essential.

## 4. Moving one step still adds coverage

Choose \(r=1,k=3,d=9\). Then \(R_1(3)=3690\pmod{19683}\), and the same
binary cylinder gives

\[
n_t=2424699+2519424t,\qquad h_t=2018303+2097152t.
\]

The words are \(n_t\xrightarrow{(1,2)}J_t\xleftarrow{(1^9,2,6)}h_t\).
This progression is incompatible with the first four rows \(k=0,1,2,3\)
of the **strengthened original-source bank** (including its inherited \(k=0\)
inverse-one rule). Its remaining rows begin at \(D(4)=12\). Thus at least

\[
1-\frac{3^{9-12}}{1-3^{-3}}=\frac{25}{26}
\]

lies outside that whole strengthened bank, again in relative natural density.
It is also distinct from the previous displayed \(k=10\) complement family:
the latter's residue at depth nine is9000, not3690. No assertion that these
children return to the same source class is made.

At \(r=0\), the new bank retains precisely the old sibling addresses, with
the shorter clock \(D(k)\), together with the old \(k=0\) inverse-one class.
The new \(k=1,d=2\) class is exactly the inherited \(n=4\pmod9\) predecessor
\((8n-5)/9\); it receives no novelty credit. Exact prefix-union arithmetic
through \(k=31\), with the proved tail after that index, gives

\[
\begin{aligned}
\delta_{\rm new}&\in[0.4458180914749391,0.4458180914749392],\\
\delta_{\rm old+12}&\in[0.4449019034751894,0.4449019034751895].
\end{aligned}
\]

These are densities of ternary dispatch domains among odd sources; they are
not added to the overlapping binary coverage. The general tail is
\(27/(26\,3^{D(K+1)})\), and (6) supplies the natural-density justification.

The incoming branch toll is a lower bound on any path with a specified first
branch excess. It does not prove this new bank is globally shortest. For
example (7) has first reverse excess two and delay six; the unrestricted
branch lower bound is four. What is sharp here is the least clock **within
the displayed inverse-run plus one reset-switch grammar**. The extra inverse
valuation two changes both the real cost and the ternary address.

Both displayed source families lie in the critical residue27 mod96, and both
children have root-distance valuation one. Thus the same strict root-rank
cancellation holds for the position-one family too. Modulo19, retain the row
instead of requiring the previous phase-identity normalization. If its child
formula is \(h=(Pn+I)/M\), then

\[
h\equiv(PM^{-1})n+IM^{-1}\pmod{19}.
\]

This affine map is a permutation because \(P\) and \(M\) are powers of2 and3.
The source periods are also units modulo19, so every phase occurs. The program
checks all19 source phases and inverse transport for both rows. The identity
phase may be lost, but the retained row makes the phase transformation lossless;
none of these congruences supplies the size inequality or a child certificate.

## 5. Repairing sparse exponent families inside the old Mersenne hostile

The earlier after-reset schema misses every \(2^{2j+1}-1\), \(j\ge1\).
For the original-source row \(k\not\equiv0\pmod3\), \(R_0(k)+1\) is a
ternary unit. Since

\[
\operatorname{ord}_{3^d}(2)=2\cdot3^{d-1},
\]

there is a unique exponent residue \(p_0\) such that
\(2^{p_0}\equiv R_0(k)+1\pmod{3^d}\). This order follows from
\(4=1+3\) and the elementary valuation identity
\(v_3(4^{3^s}-1)=s+1\). One of three exponent lifts works at each new
ternary digit; the script verifies that uniqueness directly.

The following entire Mersenne exponent progressions therefore have the new
smaller dependency:

| \(k,d\) | Exponents \(p\) | Residue mod6 |
|---|---|---|
| \(1,2\) | \(p=5+6t\) | 5; inherited small inverse word |
| \(4,12\) | \(p=75379+354294t\) | 1 |
| \(7,23\) | \(p=36969036897+62762119218t\) | 3 |

Here \(t\ge0\). These are sparse progressions within the indicated mod6
classes, not a cover of those classes. The 75,379-bit least source is expanded
and independently replayed through the thirteen child edges. The third row
is kept symbolic: modular powering verifies its source identity, and the
all-height theorem proves its dependency without allocating \(2^p\).

## 6. An infinite hostile remains: ternary towers

**PROVED.** No positive power \(n=3^a\), \(a\ge1\), belongs to the entire
adaptive initial-ones reroute grammar, at any position, index or inverse depth.

Use the complete minimum-depth criterion. Its guard is
\(4^kn+G_r(k)\equiv0\pmod{3^d}\), with \(0<G_r(k)<4^k\).
If \(d\le a\), it forces \(3^d\mid G_r(k)\), impossible since the paid
clock gives \(3^d>4^k>G_r(k)\). If \(d>a\), it forces
\(n=3^a\mid G_r(k)\), impossible because (6), oddness, and \(D(k)\ge2k\) give
\(n>2^{r+d}\ge4^k>G_r(k)\). Larger depths cannot repair the failed minimum
guard. This is a complete all-height failure mechanism for this grammar.

More generally every successful source obeys the necessary real/ternary bound

\[
3^{v_3(n)}<2(3/2)^{\lfloor\log_2n\rfloor}.
\tag{11}
\]

If \(d\le v_3(n)\), the previous contradiction already applies. Otherwise
the guard forces \(3^{v_3(n)}\mid G_r(k)\), and
\(G_r(k)<4^k<2(3/2)^d\le2(3/2)^{\lfloor\log_2n\rfloor}\).
The implementation checks the failure of (11) by exact integer arithmetic.

Writing \(n=3^au\) with odd ternary unit \(u\), an explicit excluded wedge is

\[
9^a\ge32u^3.
\tag{12}
\]

Indeed \((3/2)^5<2^3\), so the right side of (11) is at most
\(2n^{3/5}\). Condition (12) implies \(3^a\ge2n^{3/5}\), contradicting
(11). Thus for **every fixed odd cofactor \(u\)**, all sufficiently deep
ternary lifts are outside this schema. It is not just an isolated least-source
exception or a fixed finite inverse-depth obstruction.

### Divide this hostile among other constructions

The failure is attached to the schema, not to convergence. Other inherited
operations handle the following powers \(3^a\):

* Even \(a\): \(n\equiv1\pmod4\), so the immediate odd step is smaller.
* \(a\equiv1\pmod4\): the first valuations are \(1,b\) with \(b\ge3\),
  so the inherited half-child reset switch applies. At \(a=1\), the child
  is ROOT1 and first-hit trimming is required.
* \(a\equiv7\pmod8\): the first valuations are \(1,2,b\), \(b\ge2\).
  Powers of3 modulo32 give \(27\cdot3^a+23\equiv0\pmod{32}\), and
  \(U^3(n)\le(27n+23)/32<n\).

The survivor of these three specific dispatch tests is exactly
\(a\equiv3\pmod8\), whose sources satisfy \(n\equiv27\pmod{96}\).
This is not an assertion that none of the many other existing rules apply.
It identifies the next narrow target: later checkpoints that can discharge the
original source despite large ternary valuation. The current selector's
residual criterion is exact; a universal replacement for it remains open.

## 7. Evidence, API and proof-memory boundary

The [program](../../04-computation/experiments/collatz_early_reroute_20261004.py)
provides `row(r,k)`, `apply(n,row)`, and `candidates(n)`; the last returns all
minimum-depth possibilities permitted by the proved source-height cap.
Rows are checked against their canonical integer coefficients before use.
Malformed types and forged rows are rejected; a forged-intercept regression
preserves the demonstrated boundary. Typed caches keep Boolean arguments from
aliasing an already-cached integer key.

`attach_supplied_child` verifies the supplied AST's exact child value, consumes
the actual child prefix to the join, and prepends the checked source prefix.
It performs no encoder or orbit search. Its three demonstration premises are
explicitly constructed by bounded child iteration, cap5000; source iteration
is used only afterward as an independent audit. No demonstration provides the
free-parameter children for the global theorem.

Finite universes are fixed:36 family controls (\(r=0..3\), \(k=3^ru\),
\(u=1,2,4\), three parameters);2,059 legal shallow-clock comparisons among
\(n=3..399\), \(r=0..2\), \(k=1..12\), \(\ell\le12\);64 parameters from
each new refined progression;32 ternary address indices for each density
prefix;32 powers of3; six cofactor thresholds with four heights each;64
power-of-three exponent classification controls; the three symbolic Mersenne
rows and one literal large source; and three supplied-child AST demonstrations.
All checks use exceptions and remain active under `-O`.

Run `python -B 04-computation/experiments/collatz_early_reroute_20261004.py`
and repeat with `python -O -B`. The [saved stdout](collatz_early_reroute_20261004.out)
records the exact constants, directed rational density rounding, universes,
and limits. These controls support the separately proved quantifiers; they
do not replace the proofs or imply universal convergence.
