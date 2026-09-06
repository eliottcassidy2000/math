---
id: THM-4448
title: "LRC14 general shore attachment and decoder-pair cones"
status: >
  PROVED ELEMENTARY + PROVED RELATIVE TO CITED LRCUpTo13 +
  FINITE-EXACT + INDEPENDENTLY AUDITED. A lower-dimensional witness
  produces a protected closed quotient arc carrying one persistent inverse
  sheet; a scaled r-speed shore attaches whenever that arc is at least as
  long as the largest open shore-danger component. For two-speed THM-3818
  decoder shores the sharp finite-atlas constant is 29/196, while the
  filter-free bounded-ratio constant is 15/98. Cofinal hostile families
  prove that cross height cannot select an arbitrary body component.
  LRC(14) remains open.
source: component-address continuation + root + clean-room referee, 2026-09-06
depends_on:
  - LRCUpTo13
related:
  - THM-1007-weak-target-single-killer-closure
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
  - THM-4032-lrc14-d3-affine-defect-lattice-boundary
  - THM-4052-lrc14-affine-component-width-escape-cones
  - THM-4066-lrc14-diagonal-intercept-pullback-and-exact-affine-ray-closure
  - THM-4148-first-window-width-universal-odd-tail-lrc14-transfer
  - THM-4335-lrc14-owner-permutation-component-address-and-minority-renewal
  - THM-4442-lrc14-bounded-ten-body-parity-free-scale-three-completion
  - THM-4444-lrc14-signed-112-sharp-one-ray-classification
  - THM-4446-lrc14-primitive-ten-pack-descent-and-dilation-rays
primary_script: 04-computation/lrc14_general_shore_attachment_thm4448.py
primary_output: 05-knowledge/results/lrc14_general_shore_attachment_thm4448.out
primary_script_sha256: f2070a2b77d17b2635c73efe53cbddf450e94fc71161a71c2fbc12a55d41488b
primary_output_sha256: 9dfabe053f72d88a3115d2dff171ab085527c9ca21d001cad19c9325737b9667
independent_script: 04-computation/lrc14_general_shore_attachment_thm4448_independent.py
independent_output: 05-knowledge/results/lrc14_general_shore_attachment_thm4448_independent.out
independent_script_sha256: 989ed36ca551e0e31d29071579e13e4a5ba696d2ea3776f8a61d4683f374fb26
independent_output_sha256: b02e4827a56d453df05836bc8d6ed13eacfb70253d0ffa742d5d96bc0aa5f916
report: 05-knowledge/results/lrc14_general_shore_attachment_thm4448.md
audit: 05-knowledge/results/lrc14_general_shore_attachment_thm4448_independent_audit.md
hash_basis: raw LF repository bytes
---

# THM-4448 -- LRC14 general shore attachment and decoder-pair cones

**PROVED ELEMENTARY + PROVED RELATIVE TO CITED `LRCUpTo13` +
FINITE-EXACT + INDEPENDENTLY AUDITED.** The attachment theorem and hostile
families below are analytic. The tail cells and pair maxima have two
independent exact implementations; the displayed bounded head is a primary
finite-exact control.
The theorem abstracts a protected-arc move already present in
THM-4052/4066/4148; the general \(r\)-shore formulation, exact decoder-pair
constants, and the explicit obstruction to prescribed-component entry are
the added packaging. No priority is claimed for the protected-arc idea, and
LRC(14) remains **OPEN**.

For a finite positive-integer set \(R\), write
\[
 G_R=\{y\in\mathbb R/\mathbb Z:\|ry\|\ge1/14
                         \text{ for every }r\in R\}.          \tag{1}
\]
All safe sets are closed and all danger sets use the strict inequality.

## 1. A protected arc from a lower-dimensional witness

Let \(1\le r\le9\). Let \(B\) contain \(10-r\) distinct positive speeds,
let \(U\) contain \(r\) distinct positive speeds, and let
\(h\in\mathbb Z_{\ge1}\), with \(B\cap hU=\varnothing\). Let \(T\) be any
three distinct positive speeds; tails divisible by three are allowed. Define
\[
 D_U=\{z:\|uz\|<1/14\text{ for some }u\in U\}.                \tag{2}
\]
Assume \(G_U\ne\varnothing\), so \(D_U\) is a proper open subset of the
circle, and let \(\delta(U)\) be the maximum length of one of its open
components. The length is at most one. In the present range \(r\le9\), the
cited lower-dimensional Lonely Runner theorem also supplies a positive safe
arc, but only properness and this weak bound are needed below.

The physical base \(3B\cup T\) has at most \(13-r\) distinct speeds, hence
at most \(14-r\) runners after adjoining the stationary runner. Cited LRC
through thirteen total runners gives \(x_0\in\mathbb R/\mathbb Z\) such
that
\[
 \|3b x_0\|,\ \|t x_0\|\ge {1\over14-r}
 \quad(b\in B,\ t\in T).                                    \tag{3}
\]
Put \(y_0=3x_0\) and
\[
 \rho_*(x_0)=\min\left\{
  \min_{b\in B}{\|3b x_0\|-1/14\over b},\
  \min_{t\in T}{3(\|t x_0\|-1/14)\over t}\right\}.           \tag{4}
\]
This radius is positive. Moreover,
\[
 \rho_*(x_0)\ge\rho_r:=
 \min\left\{{r\over14(14-r)\max B},\
            {3r\over14(14-r)\max T}\right\}.                 \tag{5}
\]

Let
\[
 J_*=\{y:\operatorname{dist}(y,y_0)\le\rho_*(x_0)\}.         \tag{6}
\]
Since a body term in (4) is at most \(3/7<1/2\), this is a
proper closed circle arc. There is therefore a unique continuous local
inverse branch
\[
 \sigma:J_*\longrightarrow\mathbb R/\mathbb Z,\qquad
 3\sigma(y)=y,\quad \sigma(y_0)=x_0.                         \tag{7}
\]
For all \(y\in J_*\),
\[
 \|by\|\ge1/14\quad(b\in B),\qquad
 \|t\sigma(y)\|\ge1/14\quad(t\in T).                         \tag{8}
\]
Indeed \(\|\,\cdot\,\|\) is one-Lipschitz: displacement \(s\) in quotient
phase costs at most \(bs\) for \(b\), while (7) costs at most \(ts/3\)
for a tail. Equality at an endpoint remains safe. Thus \(J_*\) is not just
a scalar mass estimate: it retains one physical inverse sheet throughout.

## 2. General \(r\)-shore attachment

If
\[
                      2h\rho_*(x_0)\ge\delta(U),             \tag{9}
\]
then
\[
                       G_{\,3(B\cup hU)\cup T}\ne\varnothing.\tag{10}
\]
In particular, the witness-independent condition
\[
                           2h\rho_r\ge\delta(U)              \tag{11}
\]
suffices.

To prove this, observe that the packet-danger set in quotient coordinates
is
\[
 D_{hU}=\{y:hy\in D_U\}.                                    \tag{12}
\]
If a component of \(D_U\) has length \(\ell\le1\), each component of its
degree-\(h\) pullback has length \(\ell/h\). Thus all components of (12)
have length at most \(\delta(U)/h\). Were the connected closed arc \(J_*\)
contained in \(D_{hU}\), it would lie in one open component. That is
impossible under (9). The equality case is included: a closed arc cannot
fit in an open arc of the same length, since its two endpoint gaps cannot
both be positive. Hence some \(y\in J_*\) is safe for \(hU\), and
\(\sigma(y)\) is simultaneously safe for \(3B\), \(3hU\), and \(T\).

For the normalized \(r=1\) packet \(U=\{1\}\), (11) gives the familiar
sanity bounds \(h\ge13\max B\) and \(3h\ge13\max T\). Single-speed attachment was already
available through THM-1007 and the component machinery cited above.

## 3. Exact two-speed cones

For coprime \(1\le p<q\), let
\[
 \delta(p,q)=\max\{|I|:I\text{ is an open component of }
 D_{\{p,q\}}\}.                                             \tag{13}
\]
Specialize the theorem to an eight-speed set \(A\), the packet
\(h\{p,q\}\), and
\[
 M=\max A,\qquad E=\max T,\qquad
 \rho=\min\{1/(84M),1/(28E)\}.                              \tag{14}
\]
Then
\[
                         2h\rho\ge\delta(p,q)                \tag{15}
\]
implies
\[
                     G_{\,3(A\cup h\{p,q\})\cup T}
                         \ne\varnothing.                     \tag{16}
\]
No primitivity or ternary-unit hypothesis is used in this implication.

The finite THM-3818 decoder atlas consists of the reduced pairs
\[
 \gcd(p,q)=1,\quad p<q,\quad p+q\le356,                      \tag{17}
\]
for which every prime dividing \(p+q\) is \(2\bmod3\) and occurs to
exponent at most two. There are exactly \(5{,}855\) such pairs. Exact
strict-open interval union on denominator \(14pq\), independently replayed
with rational endpoints, gives
\[
 \boxed{\max_{\text{THM-3818 pairs}}\delta(p,q)=29/196,}
 \qquad\text{uniquely at }(p,q)=(1,28).                     \tag{18}
\]
At the leader, the speed-one tooth has length \(1/7\), and the two
speed-28 teeth at its endpoints add \(1/392\) on each side.
Consequently the uniform decoder-pair cone is
\[
                         14h\ge87M,\qquad14h\ge29E.          \tag{19}
\]
For the two THM-4444 critical packets, its tail inequalities are
\[
\begin{array}{c|c}
T&\text{tail-scale condition}\\ \hline
\{1,5,11\}&h\ge23,\\
\{2,11,20\}&h\ge42.
\end{array}                                                  \tag{20}
\]

The filter in (17) is load-bearing. Over all \(19{,}314\) coprime
\(p<q\) with \(p+q\le356\), the finite-exact maximum is
\[
 \boxed{\max\delta(p,q)=15/98,}
 \qquad\text{uniquely at }(p,q)=(1,14),                     \tag{21}
\]
and the corresponding filter-free cone is
\[
                           7h\ge45M,\qquad7h\ge15E.          \tag{22}
\]
Open components must not be merged through a safe equality wall:
\((1,13)\) is the regression control that catches this error.

## 4. Exact component addresses for the critical tails

This section assumes \(3\nmid t\) for every \(t\in T\). For
\(y\in\mathbb R/\mathbb Z\), let
\[
 \operatorname{Fib}(y)=\{x:3x=y\},\qquad
 K_t(y)=\{x\in\operatorname{Fib}(y):\|tx\|<1/14\},           \tag{23}
\]
and
\[
 F_T=\{y:\bigcup_{t\in T}K_t(y)=\operatorname{Fib}(y)\}.     \tag{24}
\]
Each \(K_t(y)\) contains at most one sheet. Represent \(y\) in \([0,1)\)
and label its fibre by \(x_j=(y+j)/3\), \(j\in\mathbb Z/3\mathbb Z\).
If \(\|ty\|<3/14\) and \(n_t(y)\) is the unique nearest integer to \(ty\),
the killed sheet has label
\[
                         j=-n_t(y)t^{-1}\pmod3.              \tag{25}
\]
At equality the tail is safe on that sheet.

Exact subdivision at \(y=n/t\pm3/(14t)\) gives
\[
\begin{array}{c|c|c}
T&\text{open component of }F_T&\text{owner labels}\\ \hline
(1,5,11)&(25/154,31/154)&(0,1,2)\\
        &(123/154,129/154)&(2,1,0)\\ \hline
(2,11,20)&(5/56,3/28)&(0,1,2)\\
         &(123/280,129/280)&(1,2,0)\\
         &(151/280,157/280)&(1,0,2)\\
         &(25/28,51/56)&(2,1,0).
\end{array}                                                  \tag{26}
\]
Their total measures are \(6/77\) and \(11/140\), respectively, recovering
the scalar values of THM-4444 while retaining the lost cyclic addresses.

Let \(C\) be any nonempty finite positive speed set. Since \(0\notin G_C\),
every connected component has an ordinary interval representative
\(J=[L,R]\subset(0,1)\), including possible singleton components. Let \(J\)
run over all of them. Then the following are equivalent:

1. \(G_{3C\cup T}=\varnothing\);
2. \(G_C\subset F_T\);
3. every \(J\) is strictly contained in an open component of \(F_T\);
4. for every \(J\), there are integers \(n_t\) such that
   \[
   {n_t\over t}-{3\over14t}<L\le R<
   {n_t\over t}+{3\over14t}\quad(t\in T),                    \tag{27}
   \]
   and
   \[
   \{-n_t t^{-1}\bmod3:t\in T\}=\mathbb Z/3\mathbb Z.        \tag{28}
   \]

The map \(x\mapsto3x\) proves the first equivalence. A connected closed
component contained in the finite open union \(F_T\) lies in one of its
open components. On that component an active tail cannot change its
nearest tooth without crossing an equality-safe wall, so (27)--(28) are
exactly the remaining conditions. This also proves the converse. The strict
endpoint inequalities are essential.

As a finite control, all \(43{,}758\) ten-subsets of \([18]\) escape
\(F_T\) for each packet in (26); through height thirteen the quotient test
also agrees with literal construction of the physical safe set. This is
**FINITE-EXACT** signal only, not an extrapolation to arbitrary height.

## 5. Cofinal obstruction to prescribed-component entry

Let
\[
 A=\{1,2,\ldots,8\},\qquad C_N=A\cup\{N,4N\}.                \tag{29}
\]
For \(T=(1,5,11)\), take
\[
            N=53+2310k,\qquad y_*=2/11,\qquad
            k\in\mathbb Z_{\ge0}.                           \tag{30}
\]
For \(T=(2,11,20)\), take
\[
            N=121+210k,\qquad y_*=1/10,\qquad
            k\in\mathbb Z_{\ge0}.                           \tag{31}
\]
In each family, \(C_N\) is primitive, every body speed is strictly safe at
\(y_*\), and the distinguished pair \(\{N,4N\}\) has cross height \(N\)
against \(A\). In (30), \(y_*\) is the centre of
\((25/154,31/154)\), at wall distance \(3/154\); in (31), it lies in
\((5/56,3/28)\), at nearest-wall distance \(1/140\).

The component of \(G_{C_N}\) through \(y_*\) is contained in a single
speed-\(N\) safe interval and therefore has diameter at most \(6/(7N)\).
The strict comparisons
\[
 {6\over7\cdot53}<{3\over154},\qquad
 {6\over7\cdot121}<{1\over140}                              \tag{32}
\]
put the entire marked component inside the corresponding failure cell.
Thus every one of its three physical lifts is spoiled, while the
distinguished cross height tends to infinity.

These are **not** non-lonely rows: cone (19) completes them. They refute
only the tempting stronger claim that large cross height certifies an
arbitrary or otherwise preselected body component. A valid entry
argument must existentially select a component or retain the sheet/address
data carried by the protected arc.

## 6. Consequence and remaining wedge

Within a live THM-3818 \(d=3\) chart whose ten-body is
\(A\cup h\{p,q\}\), the theorem removes the scale-dominant region (19).
For the two critical packets the honest residual is contained in
\[
\begin{array}{c|c}
T&\text{necessary residual disjunction}\\ \hline
(1,5,11)&14h<87\max A\ \text{or}\ h\le22,\\
(2,11,20)&14h<87\max A\ \text{or}\ h\le41.
\end{array}                                                  \tag{33}
\]
This must be intersected with THM-4446's inherited gcd-one and
cross-height-above-\(91^6\) conditions; none implies another. The next
sharp problem is opposite-scale entry: use the actual pair-specific
\(\delta(p,q)\) and the cyclic endpoint word to manufacture one protected
component in the wedge (33), without selecting a component in advance.

Tournament Analysis is intentionally not used. The two-speed packet has no
intrinsic orientation; its lawful datum is the unoriented reduced pair
together with scale and component address. Reciprocal reflection preserves
\(\delta(p,q)\), while an oriented rational edge would discard precisely
the scale and address needed by (15).

## 7. Verification and scope

The primary exact implementation checks (18), (21), (26), the bounded head,
strict walls, and both cofinal families. The independent implementation
reconstructs the event intervals and pair components without importing the
primary code, rechecks arbitrary-tail attachment arithmetic and topology,
and audits the analytic constants and hostile progressions. Normal and
optimized runs agree.

The result proves a cofinal shore-attachment cone and a precise
prescribed-component obstruction. It does not prove arbitrary component
entry, close the opposite-scale wedge, or settle LRC(14).
