---
id: THM-4446
title: "LRC14 primitive ten-pack descent and dilation rays"
status: >
  PROVED RELATIVE TO CITED LRCUpTo13, THM-3818, THM-4052, AND THM-4442
  + VERIFIED-EXACT + INDEPENDENTLY AUDITED. A primitive clock-three row with
  a ten-pack of nontrivial gcd is lonely. Consequently all 286 bounded
  THM-4442 bodies remain closed on every integer dilation ray. In the live
  THM-3818 d=3 chart a survivor must instead have primitive ten-pack,
  distinguished cross height above 91^6, and tail maximum below eleven
  times the pack maximum. This is an entry reduction, not a proof of LRC(14).
source: entry-bridge + root continuation, 2026-09-06
depends_on:
  - LRCUpTo13
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
  - THM-4052-lrc14-affine-component-width-escape-cones
  - THM-4442-lrc14-bounded-ten-body-parity-free-scale-three-completion
related:
  - THM-4004-lrc14-three-detuned-divisor-comb-profile
  - THM-4032-lrc14-d3-affine-defect-lattice-boundary
  - THM-4335-lrc14-owner-permutation-component-address-and-minority-renewal
primary_script: 04-computation/lrc14_primitive_ten_pack_descent_thm4446.py
primary_output: 05-knowledge/results/lrc14_primitive_ten_pack_descent_thm4446.out
primary_script_sha256: 41a862b59835c4aa3f9e99192bf6a832f63ede8c66c39299719979c13fada545
primary_output_sha256: f3f295cdd2eabbd506f8a707fad4b9658924d8393f9eae4e62a593225fbeebe5
independent_script: 04-computation/lrc14_primitive_ten_pack_descent_thm4446_independent.py
independent_output: 05-knowledge/results/lrc14_primitive_ten_pack_descent_thm4446_independent.out
independent_script_sha256: 740382cda20a151ed474cbc191f75cb342e40bcde0e5855f3e4b4937b0c868ce
independent_output_sha256: bdfcac0262a25615b1d83d6a307a081607aa9f71991ed593e50030fcce2a8e2f
report: 05-knowledge/results/lrc14_primitive_ten_pack_descent_thm4446.md
audit: 05-knowledge/results/lrc14_primitive_ten_pack_descent_thm4446_independent_audit.md
hash_basis: raw LF repository bytes
---

# THM-4446 -- LRC14 primitive ten-pack descent and dilation rays

**PROVED RELATIVE TO CITED LRCUpTo13, THM-3818, THM-4052, AND
THM-4442 + VERIFIED-EXACT + INDEPENDENTLY AUDITED.** The result removes
common body dilation as a possible source of a clock-three obstruction. It
does not supply arbitrary entry into the chart, and LRC(14) remains
**OPEN**.

For a finite positive integer set \(A\), write
\[
 G_A=\{x\in\mathbb R/\mathbb Z:\|ax\|\ge1/14
                         \text{ for every }a\in A\}.
\]

## 1. Open-arc grid cap

Let \(D=\{z:\|z\|<1/14\}\), the strict danger arc of length \(1/7\).
For every \(m\ge1\) and every translate of the order-\(m\) grid,
\[
 \#\{j\in\mathbb Z/m\mathbb Z:\alpha+j/m\in D\}
       \le \left\lceil {m\over7}\right\rceil .                  \tag{1}
\]
Indeed, \(r\) cyclically consecutive grid points in an open arc of length
\(1/7\) have span \((r-1)/m<1/7\), which gives (1). When \(7\mid m\),
openness excludes a spurious extra endpoint: the bound is \(m/7\), not
\(m/7+1\).

## 2. Primitive ten-pack descent

**Proposition.** Let \(C\) be a set of ten distinct positive integers and
let \(T=\{a,b,c\}\) be three distinct positive integers prime to three.
Suppose the thirteen-speed row
\[
                              S=3C\cup T                         \tag{2}
\]
is primitive. If \(\gcd(C)>1\), then \(G_S\ne\varnothing\).
More strongly, for any prime \(p\mid\gcd(C)\) and every
\(y\in G_{C/p}\), at least one point in the complete fibre of
\(x\mapsto3px\) above \(y\), labelled by
\[
                       x_j={y+j\over3p},\qquad 0\le j<3p,        \tag{3}
\]
is safe for the full row.

**Proof.** Write \(C=pC'\). Every lift (3) preserves the body:
\[
 \|(3pc')x_j\|=\|c'(y+j)\|=\|c'y\|\ge1/14.                    \tag{4}
\]
Fix \(t\in T\). On the labels in (3), its phase orbit has order
\[
 m={3p\over\gcd(t,3p)}
\]
and each orbit point occurs \(\gcd(t,3p)\) times. Because \(3\nmid t\),
there are two cases:
\[
\begin{array}{c|c|c}
 p\mid t & m=3 & \text{at most }p\text{ labels killed},\\
 p\nmid t & m=3p & \text{at most }B_p=\lceil3p/7\rceil
                                      \text{ labels killed}.
\end{array}                                                     \tag{5}
\]
The first case is impossible when \(p=3\), which only improves the count.
Primitivity of (2) forbids all three tails from being divisible by \(p\),
so at most two tails use the first row of (5). Moreover
\[
 B_p\le {3p+6\over7}<p\qquad(p\ge2).                            \tag{6}
\]
Hence the union of all bad label sets has size at most
\[
                         2p+B_p<3p.                              \tag{7}
\]
One common label survives. Finally, the cited LRC theorem for the ten-speed
set \(C'\) supplies a \(y\) with clearance at least \(1/11\), hence a
\(y\in G_{C'}\). This proves the proposition. \(\square\)

The union-bound margin is attained. For
\[
 C'=(1,2,3,4,6,7,8,9,10,11),\quad T=(1,4,10),\quad y=23/56,
\]
the primitive row \(6C'\cup T\) leaves exactly label \(j=3\) among the six
lifts. At \(p=3\), no tail is \(p\)-divisible and the actual uniform margin
is three. The proof uses cardinality ten only to obtain the cited pack phase.
The same argument is therefore valid, at threshold \(1/14\), for a primitive
row made from any pack of at most twelve distinct speeds covered by
LRCUpTo13 and at most three ternary-unit tails.

## 3. Unbounded completion of all THM-4442 dilation rays

**Corollary.** Let \(C_0\in\binom{\{1,\ldots,13\}}{10}\), let \(h\ge1\),
and let \(T\) be any three distinct positive ternary units such that the
displayed union has thirteen speeds. Then
\[
                          G_{\,3hC_0\cup T}\ne\varnothing .      \tag{8}
\]

**Proof.** Every ten-subset of \([13]\) has gcd one, since a prime has at
most six multiples there. Let \(g=\gcd(3hC_0\cup T)\). The tails imply
\(3\nmid g\); since the body gcd is \(3h\), it follows that \(g\mid h\).
Divide the row by \(g\). This preserves lonely-time existence and gives the
primitive row
\[
                    3(h/g)C_0\cup(T/g).                           \tag{9}
\]
The divided tails remain distinct ternary units. If \(h/g=1\), THM-4442
applies. If \(h/g>1\), the pack in (9) has gcd \(h/g>1\), so the proposition
applies. These cases exhaust the ray. \(\square\)

Thus every tail-universal ten-pack base chart at clock three propagates
along all integer body-dilation rays. The base case and larger-scale descent
are different mechanisms; normalization in (9) is essential.

## 4. Exact restriction in the THM-3818 d=3 chart

Work conditionally in THM-3818/4004's exact rank-eleven 11+2 branch
\[
 n=s(u_1,\ldots,u_{11})\mathbin\sqcup t(p,q),\qquad
 s\in\{1,2\},\quad\gcd(s,t)=\gcd(p,q)=1,\quad p+q\le356,        \tag{10}
\]
with \(W=V_{\rm dec}\). Suppose \(3\mid t\) and exactly eight \(u_i\) are
divisible by three. Define
\[
\begin{aligned}
 C&=\{su_i/3:3\mid u_i\}\cup\{(t/3)p,(t/3)q\},\\
 T&=\{su_j:3\nmid u_j\},\\
 P&=\{(t/3)p,(t/3)q\}\subset C.
\end{aligned}                                                     \tag{11}
\]
Then \(n=3C\cup T\). If this row is non-lonely, it necessarily satisfies
\[
 \gcd(C)=1,\qquad
 \kappa(C;P)>Q=91^6=567869252041,\qquad
 \max T<11\max C,                                                \tag{12}
\]
where
\[
 \kappa(C;P)=
 \min_{x\in P,\ z\in C\setminus P}
       \max\!\left({x\over\gcd(x,z)},{z\over\gcd(x,z)}\right).    \tag{13}
\]

The first condition is the contrapositive of the proposition. For every
cross pair \(x,z\), the primitive support-two relation between actual speeds
\(3x,3z\) has height equal to the corresponding term in (13) and meets both
decoder components. THM-3818 forbids such a crossing row through height
\(Q\), proving the second condition. Finally, THM-4052 closes the \(d=3\)
width cone \(\max T\ge11\max C\), proving the strict third condition. In
particular, \(\max C>Q\).

This separates the bounded chart from actual entry: a survivor is not a
large dilation of a bounded body. It must instead be primitive and have
enormous reduced ratios across its distinguished \(2+8\) bipartition.

## 5. The missing component-address coordinate

For \(y\in G_C\), tail \(t\), and sheet \(j\in\mathbb Z/3\mathbb Z\), put
\[
 K_t(y)=\{j:\|t(y+j)/3\|<1/14\},\qquad
 F_T=\{y:K_a(y)\cup K_b(y)\cup K_c(y)=\mathbb Z/3\mathbb Z\}.   \tag{14}
\]
The three-to-one circle map gives the exact equivalence
\[
              G_{3C\cup T}=\varnothing
                    \quad\Longleftrightarrow\quad G_C\subseteq F_T. \tag{15}
\]
THM-4032 detects when the ambient fully-spoiled set \(F_T\) is nonempty,
but does not retain the placement of \(G_C\) needed for (15). The cheapest
faithful sidecar is the cyclic word of tail-sheet enter/exit events on one
component of \(G_C\), including its endpoint owners. Equivalently, it is the
piecewise three-bit safe-sheet mask used in THM-4442.

The next open implication is therefore not a scalar mass estimate:
\[
 \gcd(C)=1,\quad\kappa(C;P)>Q
   \quad\Longrightarrow\quad
 \text{some component-address cell of }G_C\text{ has a safe sheet}. \tag{16}
\]
Statement (16) is **OPEN**. The signed \((1,1,1)\) family in THM-4445 has
cofinal fully-spoiled mass above \(6/77\), so a finite-head Haar improvement
cannot prove it.

## 6. Typed structural comparison

The proved mechanism is a genuine additive--multiplicative duality:

- **source:** a prime divisor of the body gcd and three tail divisibility
  bits;
- **map:** multiplication by \(p\) creates the additive orbit
  \(\mathbb Z/(3p)\mathbb Z\) of common lifts;
- **preserved:** every body clearance and one common label;
- **destroyed:** the component address of the chosen pack phase, which the
  uniform count does not need;
- **sidecar:** primitivity excludes the tail word 111.

For (13), the intrinsic rational object is instead a bipartite reduced-ratio
graph. Reciprocal reflection \(x/z\leftrightarrow z/x\) preserves its height
and swaps shores. Orienting it by numerical order would create a cosmetic
tournament and discard the component address required by (15).

## 7. Exact audit and scope

The primary verifier exhausts every translated order grid through order 91,
including endpoint controls \(m=7,21\). It also tests all 286 ten-subsets of
\([13]\), eight primes, and every possible count \(0,1,2\) of
\(p\)-divisible tails: 6,292 literal primitive rows.

The clean-room verifier imports no primary or repository geometry code. It
checks grid orders through 210, all primes through 1000, every admissible
ternary-unit tail event cell through 20 for \(p=2,3,5,7\), all bounded-body
gcds, and 292,500 normalizations: 2,189,781 explicit checks. Normal and
optimized runs agree on both paths.

The computations audit the proof's finite-orbit mechanism; they are not the
source of its all-height quantifier. Nothing here proves arbitrary entry
into (11), the open implication (16), or LRC(14).
