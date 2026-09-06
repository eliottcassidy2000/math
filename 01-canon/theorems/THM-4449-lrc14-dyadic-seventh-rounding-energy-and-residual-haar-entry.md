---
id: THM-4449
title: "LRC14 dyadic seventh-rounding energy and residual Haar entry"
status: >
  PROVED ELEMENTARY + FINITE-EXACT + INDEPENDENTLY AUDITED. The dyadic
  two-tail cross-comb has an exact seventh-rounding coordinate. Three odd
  tails have sharp physical failure mass at most 214/1449, with projective
  equality shape (1,9,23); for odd 3-unit tails the sharp value is 72/539,
  with equality shape (1,7,11). The associated clock-two Haar entries are
  sufficient gates, not universal body floors. LRC(14) remains open.
source: composite-clock residual continuation + two independent exact audits, 2026-09-06
depends_on:
  - THM-2060-crt-tail-coset-saturation
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-4153-third-tier-haar-finite-exception-pool43-odd-tail-transfer
related:
  - THM-4041-lrc14-d2-affine-defect-edge-boundary
  - THM-4136-fixed-body-universal-odd-tail-lrc14-completion
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4152-second-tier-haar-finite-exception-pool40-odd-tail-transfer
  - THM-4414-lrc14-six-separated-contact-capacity-collapse
  - THM-4447-lrc14-composite-clock-capacity-and-small-clock-reduction
primary_script: 04-computation/lrc14_dyadic_seventh_rounding_thm4449.py
primary_output: 05-knowledge/results/lrc14_dyadic_seventh_rounding_thm4449.out
primary_script_sha256: 88954edd6a39ee544d2548514044f49cb42d8e923df47e4a606bc776f321add3
primary_output_sha256: fffa16096b12ec941bd196e1ace42b099f8acfa7fc7537fcb8e42a5fe102acb6
sharp_script: 04-computation/lrc14_dyadic_physical_union_sharp_thm4449.py
sharp_output: 05-knowledge/results/lrc14_dyadic_physical_union_sharp_thm4449.out
sharp_script_sha256: eaceaeb497a3327ecfcfb9df5446c6e2f22934546497363c0387e8c003ff77fb
sharp_output_sha256: 780724feeb66641f6a165ca2c47c424895c38ad31ef4866d61b3ccef503dc15f
independent_script: 04-computation/lrc14_dyadic_physical_union_sharp_thm4449_independent.py
independent_output: 05-knowledge/results/lrc14_dyadic_physical_union_sharp_thm4449_independent.out
independent_script_sha256: aaad9c327e2246f96da5e2b7eb317ad524e26f1e2c2c5975376256aed9a053e5
independent_output_sha256: b5c495fc57ca2523d6bf774a0798db1964cd2a130767cc3497cf04d011901d05
report: 05-knowledge/results/lrc14_dyadic_seventh_rounding_thm4449.md
audit: 05-knowledge/results/lrc14_dyadic_seventh_rounding_thm4449_independent_audit.md
hash_basis: raw LF repository bytes
---

# THM-4449 -- LRC14 dyadic seventh-rounding energy and residual Haar entry

**PROVED ELEMENTARY + FINITE-EXACT + INDEPENDENTLY AUDITED.** The infinite
reductions are analytic. Two independent exact interval implementations audit
every finite residue/exposure box, equality case, and endpoint coordinate.
The pair caps and five-ray localization recover THM-4150/4152/4153
specializations; the seventh-rounding packaging, pair-energy theorem, sharp
three-tail physical caps, and exposure compiler are new relative to the
searched repository truth surface. No literature-priority claim is made. No
universal ten-body Haar floor follows, and LRC(14) remains **OPEN**.

Put
\[
 D_s=\{x\in\mathbb R/\mathbb Z:\|sx\|<1/14\},\qquad
 G_A=(\mathbb R/\mathbb Z)\setminus\bigcup_{a\in A}D_a.       \tag{1}
\]
Danger sets are open and safe sets are closed. An **odd 3-unit** is a
positive odd integer not divisible by three.

## 1. Seventh rounding

For distinct positive odd $a,b$, let $\Sigma_{a,b}$ be the quotient
phases on which the two tails collectively kill both lifts:
\[
 \Sigma_{a,b}=\{y:\forall j\in\{0,1\},
 \min(\|a(y+j)/2\|,\|b(y+j)/2\|)<1/14\}.                    \tag{2}
\]
Write $a=tp,b=tq$, where $t=\gcd(a,b)$, $0<p<q$, and
$\gcd(p,q)=1$. Define
\[
 \alpha=(p+q)/2,\quad\beta=(q-p)/2,
\]
\[
 d=\alpha-7\lfloor\alpha/7+1/2\rfloor,\qquad
 e=\beta-7\lfloor\beta/7+1/2\rfloor.                       \tag{3}
\]
Then $d,e\in\{-3,\ldots,3\}$ and
\[
 \boxed{\mu(\Sigma_{a,b})={2\over49}
       \left(1+{e^2-d^2\over pq}\right).}                  \tag{4}
\]

For proof, let $A(x)=1_{\{\|x\|<1/14\}}$. The two disjoint physical cross
terms form the doubling pullback of $\Sigma_{p,q}$, while
\[
 \widehat A(0)=1/7,\qquad \widehat A(k)={\sin(\pi k/7)\over\pi k}.
\]
Orthogonality leaves $(qk,-pk)$; the half-translate contributes
$(-1)^k$. Product-to-sum and the Bernoulli cosine series give
\[
 \mu(\Sigma_{p,q})={2\over49}+{2\over pq}
 \left[B_2(\{1/2+\beta/7\})-B_2(\{1/2+\alpha/7\})\right].   \tag{5}
\]
Nearest-seventh reduction and $B_2(1/2+z)=z^2-1/12$ prove (4). Finally,
$\Sigma_{tp,tq}=m_t^{-1}(\Sigma_{p,q})$, and $m_t$ preserves Haar
measure. QED.

The coordinate splits additive half-sum/difference residues from the
multiplicative decay $pq$. Reciprocal reflection reverses a projective
rational ray but preserves the squared-residue formula. Common scale is
invisible to mass and must remain as an address sidecar.

Since $e^2-d^2\le9$, finite substitution below the resulting product
cutoffs gives the sharp recovered caps
\[
 \mu(\Sigma_{a,b})\le4/63                                      \tag{6}
\]
for all odd tails, uniquely at primitive ratio $(1,9)$, and
\[
 \mu(\Sigma_{a,b})\le4/77                                      \tag{7}
\]
for odd 3-units, uniquely at $(1,11)$. For (7), $pq\ge33$ is immediate;
the possible $pq=33$ boundary exits the 3-unit domain, and the exact
$pq<33$ list finishes the proof.

## 2. Pair energy is not physical union mass

For an odd-3-unit triple $T$, define
\[
 E(T)=\sum_{\{a,b\}\subset T}\mu(\Sigma_{a,b}).             \tag{8}
\]
Then
\[
 \boxed{E(T)\le124/847,}                                      \tag{9}
\]
with equality exactly at permutations and common odd-3-unit dilates of
$(1,11,121)$. Indeed, (4) bounds every pair with $pq\ge117$ by $4/91$.
The complete high-edge alphabet is
\[
 (1,11),(1,23),(5,11),(1,37),(1,25).                       \tag{10}
\]
Two $(1,11)$ edges force the equality shape. Otherwise normalize one high
edge, test the reciprocal attachments of the other four rays, and, when no
$(1,11)$ edge occurs, test the $4\times4$ pairs of high-ray types. The
exact compatibility table reconstructed by the audit has every entry
strictly below $124/847$. This is a finite residue calculation, not a
height cutoff.

For any odd triple put
\[
 U_T=\bigcup_{s\in T}D_s,\qquad P_T=U_T\cap(U_T+1/2).        \tag{11}
\]
If $F_T$ is its quotient failure set, then
$m_2^{-1}(F_T)=P_T$, hence $\mu(F_T)=\mu(P_T)$. One odd
tail cannot kill both lifts, so
\[
 F_T=\bigcup_{\{a,b\}\subset T}\Sigma_{a,b}.               \tag{12}
\]
Let $\Omega_T$ be the mixed locus where one lift has one owner and the
other has the remaining two. Pointwise owner counting yields
\[
 E(T)=\mu(F_T)+\mu(\Omega_T).                              \tag{13}
\]
Thus this is an undirected two-colour owner cut with ties, not a tournament;
an arbitrary orientation destroys exactly the multiplicity correction. At
the energy equality shape the two terms are $108/847$ and $16/847$.

## 3. Sharp three-tail physical caps

For every three distinct positive odd tails,
\[
 \boxed{\mu(F_T)\le214/1449,}                               \tag{14}
\]
with equality exactly at permutations and common positive odd dilates of
$(1,9,23)$. For odd 3-units the sharp refinement is
\[
 \boxed{\mu(F_T)\le72/539,}                                 \tag{15}
\]
with equality exactly at permutations and common odd-3-unit dilates of
$(1,7,11)$.

Write the defect in (4) as \(\delta=(e^2-d^2)/(pq)\). Then
\[
 E(T)=6/49+(2/49)\sum\delta.                                \tag{16}
\]
For (14), the neutral edge level is $128/621$, because
\[
 214/1449=6/49+(2/49)3(128/621).                            \tag{17}
\]
If every edge is neutral, (13) proves the bound. Equality there is impossible:
it would force $128\mid(e^2-d^2)$, whose positive value is at most nine.
A strict exceptional edge has $pq<44$, giving exactly
\[
 (1,9),(1,11),(1,23),(3,11).                               \tag{18}
\]
For (15), the neutral level is $1/11$. Simultaneous equality of all three
edges is excluded by exact substitution, and a strict edge has $pq<99$,
giving precisely (10).

Fix an exceptional ray $(p,q)$, present its tails as $(tp,tq)$, and call
the third tail $c$. Dividing the triple gcd preserves mass and leaves
$\gcd(t,c)=1$. With
\[
 A=D_p\cup D_q,\qquad R=(A+1/2)\setminus A,                 \tag{19}
\]
the exact disjoint owner expansion is
\[
 \mu(F_{tp,tq,c})=\mu(\Sigma_{p,q})
             +2\mu(m_t^{-1}R\cap D_c).                     \tag{20}
\]
If an a.e.-minimal interval representative of $R$ has $N$ essential circle
components, bounded variation and the coprime
Fourier diagonals $(ck,-tk)$ give
\[
 \left|\mu(m_t^{-1}R\cap D_c)-\mu(R)/7\right|\le {N\over3tc}. \tag{21}
\]
Specifically, $\operatorname{Var}(1_R)=2N$, $\operatorname{Var}(1_{D_1})=2$, and
$\sum_{k\ne0}k^{-2}=\pi^2/3$.

Here `N` is a BV/a.e. quantity: isolated deleted endpoints are filled before
counting. It is not the number of actual connected components of the strict
open representative. THM-4451 records why the distinction matters.

For target $M$, let
\[
 \eta=(M-\mu(\Sigma_{p,q}))/2-\mu(R)/7.                    \tag{22}
\]
Exact interval geometry gives positive margins. Equations (20)--(21) are
strict beyond $tc>N/(3\eta)$. The nonintegral cutoffs are
\[
\begin{array}{c|rrrr}
\text{odd ray}&(1,9)&(1,11)&(1,23)&(3,11)\\ \hline
N/(3\eta)&30429/128&409101/1822&10143/26&409101/1822
\end{array}                                                \tag{23}
\]
and
\[
\begin{array}{c|rrrrr}
\text{3-unit ray}&(1,11)&(1,23)&(1,25)&(1,37)&(5,11)\\ \hline
N/(3\eta)&5929/15&86779/135&309925/489&219373/243&35035/96.
\end{array}                                                \tag{24}
\]
Literal rational interval enumeration below them checks 1,704 full-odd and
2,801 3-unit primitive presentations. There are no violations and only the
equality shapes in (14)--(15). An independent implementation lifts all
endpoints to denominator $\operatorname{lcm}(14pqt,14c)$ and agrees exactly.

At full-odd equality the physical pullback has 16 components of maximum
width $10/483$, while the quotient has 8 of maximum width $20/483$. At
3-unit equality the corresponding data are 12, $1/77$ and 6, $2/77$.
Doubling preserves mass but changes component addresses, counts, and widths.

## 4. Residual-clock consequences

For a ten-body $C$, failure of $2C\cup T$ implies $G_C\subset F_T$.
The former is compact and the latter proper open, so containment cannot hold
at equal Haar mass. Thus the sufficient clock-two gates, including equality,
are
\[
 \mu(G_C)\ge214/1449\quad\text{for odd tails},\qquad
 \mu(G_C)\ge72/539\quad\text{for odd 3-unit tails}.          \tag{25}
\]

The remaining dyadic signatures retype as
\[
 2C\cup\{2r,a,b\}=2(C\cup\{r\})\cup\{a,b\},\qquad
 4C\cup\{2r,a,b\}=2(2C\cup\{r\})\cup\{a,b\}.             \tag{26}
\]
The elementary transfers
\[
 \mu(G_{C\cup\{r\}})\ge\mu(G_C)-1/7,qquad
 \mu(G_{2C\cup\{r\}})\ge\mu(G_C)/2                        \tag{27}
\]
combine with (6)--(7). For odd-3-unit pairs, sufficient original-body levels
are $15/77$ and $8/77$, respectively.

At absorbed mass at least $4/91$, inherited THM-4153 localization leaves
only the five rays (10), with largest primitive component widths
\[
 2/77,\ 2/161,\ 9/385,\ 2/259,\ 2/175                       \tag{28}
\]
in the same order. If $L_H$ is the longest positive-length component of
$G_H$, failure at common scale $t$ requires $tL_H$ to be smaller than
the corresponding width. The crude original-body localization levels from
(27) are $17/91$ and $8/91$.

These are conditional entry gates, not asserted body floors. No residual
clock is declared closed. The exact reproduction commands and raw-LF hashes
are recorded in the result note and independent audit.
