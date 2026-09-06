---
id: THM-4447
title: "LRC14 composite-clock capacity and small-clock reduction"
status: >
  PROVED ELEMENTARY + PROVED RELATIVE TO CITED LRCUpTo13 +
  VERIFIED-EXACT + INDEPENDENTLY AUDITED. A tail on a q-sheet pack fibre
  has an exact floor/ceiling bad-label count and sharp capacity
  gcd(q,t) ceil(q/[7 gcd(q,t)]). Optimizing after absorbing tails along
  every divisor closes every primitive ten-pack clock at least five and
  leaves only explicit clocks two, three, and four. That small-clock
  corollary independently recovers a prior unnumbered effective-clock
  result; the pointwise chamber formula and general divisor packaging are
  the added content. LRC(14) remains open.
source: pattern111_next + root + clean-room referee, 2026-09-06
depends_on:
  - LRCUpTo13
related:
  - THM-737-pack-clock-sampling-measure-dispatch
  - THM-761-multi-exception-sheet-covering-bound
  - THM-765-safe-component-tooth-deck-and-hereditary-primitivity
  - THM-2060-crt-tail-coset-saturation
  - THM-2064-multitail-sheet-capacity-and-dyadic-seam
  - THM-4004-lrc14-three-detuned-divisor-comb-profile
  - THM-4446-lrc14-primitive-ten-pack-descent-and-dilation-rays
prior_result: 05-knowledge/results/lrc14_effective_clock_empty_core_sep06.md
primary_script: 04-computation/lrc14_composite_clock_capacity_thm4447.py
primary_output: 05-knowledge/results/lrc14_composite_clock_capacity_thm4447.out
primary_script_sha256: 06dfd8d0da921e680d690a295c5baaa6d2ca1fc85c212dd387393cf6daf9da72
primary_output_sha256: 3006a5bd34b4cef5dd4d3c64a170bc51664f335d86dcd4b8cf90e9376bd83061
independent_script: 04-computation/lrc14_composite_clock_capacity_thm4447_independent.py
independent_output: 05-knowledge/results/lrc14_composite_clock_capacity_thm4447_independent.out
independent_script_sha256: a749488996d5870a2582fbb2a749b6ff440ef407bc6d9150d1ec87ba4a851f02
independent_output_sha256: c7daf0a767e17a6427eb2ee9987412ec3ffb1a2c0498b9fa3b8b68073559a94c
report: 05-knowledge/results/lrc14_composite_clock_capacity_thm4447.md
audit: 05-knowledge/results/lrc14_composite_clock_capacity_thm4447_independent_audit.md
hash_basis: raw LF repository bytes
---

# THM-4447 -- LRC14 composite-clock capacity and small-clock reduction

**PROVED ELEMENTARY + PROVED RELATIVE TO CITED LRCUpTo13 +
VERIFIED-EXACT + INDEPENDENTLY AUDITED.** The pointwise count and divisor
operation below are new canonical packaging. The final three-clock list is
an independent proof of a conclusion already recorded in the unnumbered
effective-clock result linked in the frontmatter. No priority is claimed for
that corollary, and LRC(14) remains **OPEN**.

For a finite positive integer set \(R\), put
\[
 G_R=\{y\in\mathbb R/\mathbb Z:\|ry\|\ge1/14
                         \text{ for every }r\in R\}.
\]

## 1. Exact pointwise clock count

Let \(q\ge2\), \(y\in G_R\), and consider the complete fibre
\[
                    x_j={y+j\over q},\qquad j\in\mathbb Z/q\mathbb Z.
                                                               \tag{1}
\]
Every point in (1) remains safe for \(qR\). For a tail \(t>0\), define
\[
 K_t(y)=\{j:\|t(y+j)/q\|<1/14\}.                              \tag{2}
\]
Put
\[
 g=\gcd(t,q),\qquad m=q/g,\qquad
 \delta={ty\over q},\qquad A=-m\delta-{m\over14}.             \tag{3}
\]
Multiplication by \(t/g\) permutes the order-\(m\) grid, and every grid
point has exactly \(g\) labelled preimages. Since the relevant open interval
has length \(m/7<m\), reduction modulo \(m\) loses no multiplicity. Hence
\[
\boxed{\quad
 |K_t(y)|
   =g\,\#\bigl(\mathbb Z\cap(A,A+m/7)\bigr)
   =g\left(\left\lceil A+{m\over7}\right\rceil
                 -\lfloor A\rfloor-1\right).
\quad}                                                        \tag{4}
\]
In particular the sharp uniform capacity is
\[
             B_q(t)=g\left\lceil{m\over7}\right\rceil
                   =g\left\lceil{q\over7g}\right\rceil.        \tag{5}
\]
Every value in (5) is attained by some phase translate.

For tails \(T=\{t_1,\ldots,t_r\}\), (4) gives the pointwise bound
\[
 \#\{\hbox{common safe labels above }y\}
 \ge q-\left|\bigcup_iK_{t_i}(y)\right|
 \ge q-\sum_i|K_{t_i}(y)|.                                   \tag{6}
\]
Consequently
\[
                \sum_i B_q(t_i)<q                             \tag{7}
\]
leaves a common safe label above every \(y\in G_R\). Fibre integration also
gives the measure estimate
\[
 \mu(G_{qR\cup T})\ge
 \mu(G_R)\max\left(0,1-{1\over q}\sum_iB_q(t_i)\right).        \tag{8}
\]

## 2. Low-count chambers, walls, and equality

Write \(m/7=k+s\), where \(k=\lfloor m/7\rfloor\), \(0\le s<1\), and
write \(f=\{A\}\). Formula (4) has exactly two possible orbit counts:

- if \(s=0\), the count is \(k-1\) at \(f=0\), and \(k\) otherwise;
- if \(0<s<1\), the count is \(k+1\) when \(f>1-s\), and \(k\) when
  \(f\le1-s\).

Thus the labelled count is either \(B_q(t)\) or \(B_q(t)-g\). A resonance
wall can initiate a drop of exactly \(g\) labels because strict danger
excludes a grid point at equality. When \(s\ne0\), however, the lower count
persists through the whole chamber \(f\in[0,1-s]\); it is not exclusively
an endpoint phenomenon.

If \(\sum_iB_q(t_i)=q\), the cap has no conclusion. A full cover occurs
exactly when every tail is cap-tight and the sets \(K_{t_i}(y)\) form a
disjoint partition of \(\mathbb Z/q\mathbb Z\). Any overlap or any
low-count tail leaves a safe label. Exact partition controls are
\[
\begin{array}{c|c|c}
q&T&y\\ \hline
3&(1,4,5)&23/112\\
4&(1,11,14)&86/539.
\end{array}                                                   \tag{9}
\]
These are hostile pack phases, not non-lonely full rows.

## 3. Divisor absorption

Let \(d\mid q\), \(d\ge2\). Repackage the same physical row by absorbing
every \(d\)-divisible tail:
\[
\begin{aligned}
 qC\cup T&=dR_d\cup D_d,\\
 R_d&=(q/d)C\cup\{t/d:t\in T,\ d\mid t\},\\
 D_d&=\{t\in T:d\nmid t\}.                                  \tag{10}
\end{aligned}
\]
No speed or label is discarded. Whenever \(G_{R_d}\ne\varnothing\),
the sufficient capacity is
\[
 \operatorname{Cap}_d(T):=
 \sum_{t\in D_d}\gcd(t,d)
       \left\lceil{d\over7\gcd(t,d)}\right\rceil<d.            \tag{11}
\]
The strongest certificate minimizes
\(\operatorname{Cap}_d(T)-d\) over the divisors for which the absorbed pack
is known safe. The operation preserves the physical row and common branch
label; it changes only the chosen pack/tail decomposition.

## 4. Primitive ten-pack corollary

**Recovered corollary / independent proof.** Let \(P\) contain ten distinct
positive speeds, let \(T=\{t_1,t_2,t_3\}\) contain three further distinct
positive speeds, and suppose \(S=P\cup T\) is a primitive thirteen-speed
row. Let
\[
                         c=\gcd(P)>1.                           \tag{12}
\]
Then \(S\) is \(1/14\)-lonely unless one of the following necessary
small-clock signatures holds:
\[
\begin{array}{c|l}
c&\hbox{signature not closed by the theorem}\\ \hline
2&\hbox{zero or one even tail},\\
3&3\nmid t_1t_2t_3,\\
4&\hbox{exactly one tail has }v_2(t)=1
       \hbox{ and the other two are odd}.
\end{array}                                                   \tag{13}
\]
In particular every exact ten-pack clock \(c\ge5\) closes.

**Proof.** Write \(P=cC\), so \(\gcd(C)=1\). For every divisor \(d>1\)
of \(c\), primitivity prevents all three tails from being \(d\)-divisible.
Thus \(R_d\) in (10) has at most twelve distinct speeds. Cited LRC through
thirteen total runners gives a phase of clearance at least \(1/13>1/14\),
so (11) applies.

If a prime \(p\ge5\) divides \(c\), take \(d=p\). At most two tails are
absorbed; every residual tail has capacity \(\lceil p/7\rceil\). Therefore
\[
                  (3-r)\lceil p/7\rceil<p                     \tag{14}
\]
for \(0\le r\le2\): check \(p=5,7\), and use
\(3(p+6)/7<p\) for \(p\ge11\).

It remains that \(c=2^a3^b\). If \(a,b>0\), use \(d=2\) when at least two
tails are even. Otherwise, if any tail is divisible by three, use \(d=3\):
one or two such tails are absorbed and at most two unit-capacity tails
remain. If none is divisible by three, use \(d=6\); with at most one even
tail the total capacity is at most \(2+1+1=4<6\).

For \(c=3^b\), a three-divisible tail is handled by \(d=3\). If there is no
such tail, \(c=3\) is equality, while \(d=9\) gives total capacity
\(3\lceil9/7\rceil=6<9\) for \(b\ge2\).

For \(c=2^a\), two even tails close at \(d=2\). With no even tail, direct
use of \(d=c\) closes every \(c\ge4\). With exactly one even tail, valuation
at least two closes at \(d=4\); valuation one closes at \(d=8\) for
\(c\ge8\). The remaining cases are exactly (13). \(\square\)

The same proof gives a nonintrinsic sufficient form: it is enough that a
chosen \(q>1\) divide all ten pack speeds and that the full row be primitive.
Using the exact gcd in (12) makes the residual signature intrinsic. For a
coherent row \(qC\cup T\), every \(q\ge5\) therefore closes; at smaller
displayed \(q\), any additional body gcd should first be folded into \(c\).

## 5. Prior result and the live frontier

The unnumbered result
\[
 \texttt{05-knowledge/results/lrc14\_effective\_clock\_empty\_core\_sep06.md}
\]
already obtained the same four effective-order profiles using THM-761/765,
hereditary deletion gcds, and pair lcm constraints. Section 4 is retained
because divisor absorption gives an independent direct proof and places that
conclusion under a theorem identifier. The novelty claim is limited to the
pointwise formula (4), its complete two-chamber endpoint description, and
the general optimization (10)--(11).

At clock two, the one-even case is exactly the old eleven-pack/two-odd-tail
dyadic seam after absorption; the zero-even case is a ten-pack/three-odd-tail
cover. Clock three is the component-address problem treated on bounded
bodies by THM-4442 and on nonprimitive packs by THM-4446. Clock four reduces
at \(d=2\) to the same two-odd-tail equality when its unique even tail has
valuation one. Exact pack phases show all three scalar boundaries can cover,
so further progress must use body component placement or mask intersections,
not a stricter sum of individual capacities.

## 6. Verification and scope

The primary verifier checks every translated order grid through \(140\),
literal clocks through \(60\), and all \(43{,}991\) primitive gcd-signature
triples through clock \(210\), including arbitrary three-adic tails. It
replays strict endpoint, equality-partition, and divisor-two controls.

The clean-room verifier imports no primary or repository geometry code. It
rederives the count from independent cyclic distances and checks every
primitive gcd signature through clock \(512\), together with symbolic case
guards and the displayed hostiles. Normal and optimized runs agree.

The finite sweeps audit the proof; they do not supply its all-clock
quantifier. Nothing here empties any signature in (13), proves arbitrary
entry into a coherent ten-pack chart, or proves LRC(14).
