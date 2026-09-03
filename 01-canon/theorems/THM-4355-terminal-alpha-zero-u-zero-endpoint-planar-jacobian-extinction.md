---
id: THM-4355
title: "Terminal-alpha-zero U-zero endpoint planar Jacobian extinction"
status: >
  PROVED RELATIVE TO THM-4230/4327/4341/4344/4352 + VERIFIED-EXACT +
  INDEPENDENTLY HOSTILE-AUDITED. The inherited exact-weight-twelve gate
  Z=beta_11=zeta_3=W=xi_10=U=alpha_11=0 with K!=0 is extinct. Its 64 exact
  supports have twelve fans and exactly four collision walls: odd A4/A2,
  even A3, and two-branch A15. Every carrier and collision tail of positive
  genus has positive form order. With THM-4354 this closes the full U=0,K!=0
  endpoint. Seam entry, JC(2), and DC(2) are not claimed.
source: root + jc-deep-u0 + independent A4/A2/A15 referees / next-sharp session, 2026-09-02
depends_on:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4327-generic-exact-weight-twelve-endpoint-wall-extinction
  - THM-4341-odd-self-similar-cusp-reciprocal-tail-duality
  - THM-4344-clean-cubic-infinity-exit-planar-jacobian-extinction
  - THM-4352-even-self-similar-cusp-reciprocal-parity-and-attachment-law
related:
  - THM-4353-simultaneous-zero-endpoint-planar-jacobian-extinction
  - THM-4354-first-normal-owner-u-zero-endpoint-planar-jacobian-extinction
mistake_firewall:
  - MISTAKE-522
  - MISTAKE-531
  - MISTAKE-540
primary_script: 04-computation/jc2_m12_terminal_alpha_zero_u_zero_endpoint_extinction_thm4355.py
primary_output: 05-knowledge/results/jc2_m12_terminal_alpha_zero_u_zero_endpoint_extinction_thm4355.out
primary_script_sha256: 9e27fb396184e86d7ee397accd38d34a9a1b6737ce49c6cef33336d603d12c98
primary_output_sha256: f81477c19bdcb590a5b02ea453714ea8efc80bfe0d75ee686364c1536952b9a7
referee_script: 04-computation/jc2_m12_terminal_alpha_zero_u_zero_endpoint_extinction_independent_referee_thm4355.py
referee_output: 05-knowledge/results/jc2_m12_terminal_alpha_zero_u_zero_endpoint_extinction_independent_referee_thm4355.out
referee_script_sha256: b9a5b3f00e66949972e56bcb5f12382ee9f2a6bc34dade97b4b5b02018bee4ea
referee_output_sha256: 78e1fedef14ce56278aabeec6eb7196e6d89acab78ce262b780457c8f3a0730b
hash_basis: raw LF bytes
audit: >
  PASS. The 769-check primary and independently written 774-check referee
  reconstruct the literal source, all 64 exact supports and twelve fans,
  separate hostile atlases, eighteen faces, clean graph/genus ledgers, four
  collision charts, critical depths, primitive target-compatible form
  orders, attachment addresses, and the reversible natural index. Normal and
  optimized executions byte-match both frozen outputs.
---

# THM-4355 -- Terminal-alpha-zero U-zero endpoint planar Jacobian extinction

**PROVED RELATIVE TO THM-4230/4327/4341/4344/4352 + VERIFIED-EXACT +
INDEPENDENTLY HOSTILE-AUDITED. JC(2), DC(2), AND SEAM ENTRY REMAIN OPEN.**

## 1. Statement and inheritance

Work over $\mathbb C$ in the inherited $b=d=0$ reduced $(2,3)$,
exact-weight-twelve seam. Put $y=sp$, $e=-1376/135$, $u=\upsilon_5$,
$\alpha=\alpha_{11}$, and

~~~text
H=-3p+(8/3)p^2+e p^3+K y^2+Phi p^2y+Delta p^4
  +Theta py^2+eta p^3y+zeta_3 y^3+u p^5
  +xi_10 p^2y^2+alpha p^4y+beta_11 py^3
  +U p^6+W p^3y^2+Z y^4,

K=2848/45-(7/6)Delta,
F_Q=(s^2-p)(1-QH)-Q s^2/2.                              (1)
~~~

Impose

~~~text
Z=beta_11=zeta_3=W=xi_10=U=alpha=0,       K!=0.          (2)
~~~

No nonautomorphic planar Keller pair lies on (1)--(2), with the remaining
coefficients arbitrary subject to the seam relation. The conclusion is
relative to the inherited toroidal, good-target, and proper-flat interfaces.
Together with THM-4354, it closes the full endpoint

~~~text
Z=beta_11=zeta_3=W=xi_10=U=0,             K!=0.          (3)
~~~

Neither statement supplies seam entry or proves JC(2) or DC(2).

The inheritance pass is:

- closest proved mechanism: THM-4354, but its eleven nodes
  $\alpha S^{11}=1$ leave the finite torus when $\alpha$ vanishes;
- canonical hostile: the two ends of an even cusp need not land in one
  connected complement;
- corrected near miss: the A15 tail has two marked ends but adds no graph
  cycle because those ends join different complementary components;
- least-used sidecar: the endpoint-to-component address map.

The live concept board was

~~~text
endpoint owner | odd versus even cusp | persistent delta | attachment address
source base versus target base | support rank | labelled multigraph.         (4)
~~~

## 2. Exact support quotient and natural address

On (2), the three multiply visible support coefficients include

~~~text
(2,3,1)=K-e,       (2,4,1)=Theta-Delta,       (2,5,1)=-u. (5)
~~~

The first vanishes exactly at

~~~text
Delta_*=3968/63,                 K=e=-1376/135.          (6)
~~~

For $\Theta\ne0$, five coupled $(\Delta,\Theta)$ support classes may be
represented by

~~~text
(0,1), (Delta_*,Delta_*), (Delta_*,2Delta_*), (1,1), (1,2).
~~~

Crossing these with the presence bits of $u,\Phi,\eta$ gives $5\cdot8=40$
supports. For $\Theta=0$, the three $\Delta$ classes
$0,\Delta_*,1$ give $3\cdot8=24$. Thus (2) has exactly 64
coefficient-realizable supports.

The two strata have respectively five and seven lower fans:

| stratum | supports | fan | Pick $(2A,B,I)$ | $(V,E,b_1)$ | carrier genus | total genus |
|---|---:|---|---:|---:|---:|---:|
| $u\Theta\ne0$ | 20 | BU,FUT,T | $(32,12,11)$ | $(8,16,9)$ | 2 | 11 |
| $u=0,\Delta\eta\Theta\ne0$ | 8 | C4e,BETA,FET,T | $(29,11,10)$ | $(5,12,8)$ | 2 | 10 |
| $u=\eta=0,\Delta\Theta\ne0$ | 8 | BDT,T | $(26,10,9)$ | $(3,9,7)$ | 2 | 9 |
| $u=\Delta=0,\eta\Theta\ne0$ | 2 | C3e,BETA,FET,T | $(28,10,10)$ | $(5,12,8)$ | 2 | 10 |
| $u=\eta=\Delta=0,\Theta\ne0$ | 2 | CTH,BTH,T | $(24,10,8)$ | $(4,10,7)$ | 1 | 8 |
| $\Theta=0,u\eta\ne0$ | 6 | BU,FUE,EETA | $(31,11,11)$ | $(8,16,9)$ | 2 | 11 |
| $\Theta=\eta=0,u\ne0$ | 6 | BU,EU | $(30,10,11)$ | $(7,15,9)$ | 2 | 11 |
| $\Theta=u=0,\Delta\eta\ne0$ | 4 | C4e,BETA,EETA | $(28,10,10)$ | $(4,11,8)$ | 2 | 10 |
| $\Theta=u=\eta=0,\Delta\ne0$ | 4 | BDELTA,EDELTA | $(24,10,8)$ | $(6,12,7)$ | 1 | 8 |
| $\Theta=u=\Delta=0,\eta\ne0$ | 2 | C3e,BETA,EETA | $(27,9,10)$ | $(4,11,8)$ | 2 | 10 |
| $\Theta=u=\eta=\Delta=0,\Phi\ne0$ | 1 | CPHI,BPHI,EDELTA | $(21,9,7)$ | $(4,9,6)$ | 1 | 7 |
| $\Theta=u=\eta=\Delta=\Phi=0$ | 1 | BDEEP | $(18,8,6)$ | $(2,6,5)$ | 1 | 6 |

The primary hostile atlases independently toggle support deletions inside the
two honest $\Theta$ strata. They have $128$ keys, $96$ supports, nine fans
and $128$ keys, $72$ supports, eleven fans. Their tagged union has 256 keys,
168 supports, and twenty fans. The referee instead builds one
support-decoupled 128-support, seventeen-fan atlas. Both contain every exact
support; neither asserts realizability of the extra states.

There is a reversible natural address. Number the five $\Theta\ne0$
coefficient classes by $c=0,\ldots,4$ and set

~~~text
n=1+8c+4[ u!=0 ]+2[ Phi!=0 ]+[ eta!=0 ],       1<=n<=40. (7)
~~~

For $\Theta=0$, number the three $\Delta$ classes by $d=0,1,2$ and set

~~~text
n=41+8d+4[ u!=0 ]+2[ Phi!=0 ]+[ eta!=0 ],      41<=n<=64. (8)
~~~

The certificates check both inverses. The odd-square display retains the
ordinal because

~~~text
((2n-1)^2-1)/8=n(n-1)/2=T(n-1).                         (9)
~~~

## 3. Exhaustive face and order ledger

For a supporting plane $(a,b,c)$, the target-compatible order is
$L(5/6-a-b-c)$, where $L$ clears both the source plane and the target
sixth root.

| face | plane | initial equation, up to a torus monomial | geometry | base/order |
|---|---|---|---|---:|
| BU | $(1/10,1/5,-1/5)$ | $(P-S^2)(uP^5-1)$ | six rational, 10 nodes | $30/22$ |
| FUT | $(1/5,1/5,-2/5)$ | $S^2(1-uP^5-\eta SP^4-\Theta S^2P^3)$ | genus 2 off L5 | $30/25$ |
| T | $(1/2,0,-1)$ | $S^2(1-S^2P^2(K+\Theta P))$ | rational | $6/8$ |
| C4e | $(0,1/4,-1/4)$ | $P(\Delta P^4+\eta SP^4-1)$ | rational | $12/10$ |
| BETA | $(1/9,2/9,-2/9)$ | $(P-S^2)(\eta SP^4-1)$ | two rational, 9 nodes | $18/13$ |
| FET | $(1/5,1/5,-2/5)$ | $S^2(1-\eta SP^4-\Theta S^2P^3)$ | genus 2 | $30/25$ |
| BDT | $(1/8,1/4,-1/4)$ | $(P-S^2)(\Delta P^4+\Theta S^2P^3-1)$ | rational plus genus 2 off L15 | $24/17$ |
| C3e | $(-1/3,1/3,-1/3)$ | $P(eP^3+\eta SP^4-1)$ | rational | $6/7$ |
| CTH | $(0,1/3,-1/3)$ | $P(eP^3+\Phi SP^3+\Theta S^2P^3-1)$ | genus 1 off L3 | $6/5$ |
| BTH | $(1/8,1/4,-1/4)$ | $(P-S^2)(\Theta S^2P^3-1)$ | two rational, 8 nodes | $24/17$ |
| FUE | $(1/5,1/5,-2/5)$ | $S^2(1-uP^5-\eta SP^4)$ | rational | $30/25$ |
| EETA | $(1/3,1/6,-2/3)$ | $S^2(1-\eta SP^4-KS^2P^2)$ | genus 2 | $6/6$ |
| EU | $(3/10,1/5,-3/5)$ | $S^2(1-uP^5-KS^2P^2)$ | genus 2 | $30/28$ |
| BDELTA | $(1/8,1/4,-1/4)$ | $(P-S^2)(\Delta P^4-1)$ | five rational, 8 nodes | $24/17$ |
| EDELTA | $(1/4,1/4,-1/2)$ | $S^2(1-\Delta P^4-\Phi SP^3-KS^2P^2)$ | genus 1 off L4 | $12/10$ |
| CPHI | $(0,1/3,-1/3)$ | $P(eP^3+\Phi SP^3-1)$ | rational | $6/5$ |
| BPHI | $(1/7,2/7,-2/7)$ | $(P-S^2)(\Phi SP^3-1)$ | two rational, 7 nodes | $42/29$ |
| BDEEP | $(1/6,1/3,-1/3)$ | $(P-S^2)(eP^3+KS^2P^2-1)$ | rational plus genus 1, 6 nodes | $6/4$ |

The positive-genus normalization identities are

~~~text
FUT:     Y^2=P(4Theta+(eta^2-4uTheta)P^5),
FET:     Y^2=P(4Theta+eta^2P^5),
BDT-C:   Y^2=Theta P(1-Delta P^4),
CTH:     Y^2=(Phi^2-4eTheta)P^4+4Theta P,

EETA:    Y^2=eta^2P^6+4K,
EU:      Y^2=K(1-uP^5),
EDELTA:  Y^2=(Phi^2-4KDelta)P^4+4K,
BDEEP-C: Y^2=K(1-eP^3).                                (10)
~~~

On the displayed clean strata these branch polynomials are squarefree and
give exactly the genera in the fan table.

## 4. Collision completeness

Every toric edge is reduced on its owner gate except one variable factor in
each of four faces:

| face | variable edge word | wall | type |
|---|---|---|---|
| FUT | $\Theta+\eta X+uX^2$ | L5: $\eta^2=4u\Theta$ | A4 |
| CTH | $\Theta+\Phi X+eX^2$ | L3: $\Phi^2=4e\Theta$ | A2 |
| BDT | $(X-1)(\Delta X+\Theta)$ | L15: $\Delta+\Theta=0$ | A15 |
| EDELTA | $\Delta+\Phi X+KX^2$ | L4: $\Phi^2=4K\Delta$ | A3 |

The exact edge audit checks all face-edge words, not only this table.
Coefficientwise collection occurs before specialization: the deletions in
(5) are support walls, while L5, L3, L15, and L4 are the complete
singularity list.

## 5. Odd A4 and A2 walls

On L5 with $u\Theta\ne0$, use the source-minimal parameter $\delta$:

~~~text
Q=delta^5, s=delta^-1 S, p=delta^-1 P, G=delta^2 F_Q,
A(v)=u+eta v+Theta v^2,       B(v)=Delta+Phi v+Kv^2,
x=P^-1, v=S/P, rho=delta*x.
~~~

The literal primitive FUT chart gives

~~~text
x^7G=(v^2-rho)[x^5-A(v)-rho B(v)-e rho^2-(8/3)rho^3+3rho^4]
     -rho^5v^2/2.                                         (11)
~~~

Writing $A(v)=\Theta(v-a)^2$ makes $a\ne0$, so the prefactor is a unit.
Formal Morse preparation yields $y^2=x^5-\psi(\rho)$. If
$q=2Ka+\Phi$, its first selectors are

~~~text
c1=Delta+Phi a+Ka^2,
c2=e-q^2/(4Theta),
c3=8/3+Kq^2/(4Theta^2),
c4=-3-K^2q^2/(4Theta^3).                                (12)
~~~

If $c_1=c_2=c_3=0$, then $K/\Theta=45/172$ and
$c_4=-99/43$. Hence the depth is exactly one of $1,2,3,4$; the primary
also gives seam-compatible witnesses for every row.

For the common target base put $\delta=t^6$, so $Q=t^{30}$. The primitive
A4 packets are

| depth | tail genus | persistent delta | $(\operatorname{ord}t,\operatorname{ord}x,\operatorname{ord}y)$ | form order |
|---:|---:|---:|---:|---:|
| 1 | 2 | 0 | $(4,6,15)$ | 115 |
| 2 | 1 | 1 | $(1,4,10)$ | 35 |
| 3 | 1 | 1 | $(2,18,45)$ | 95 |
| 4 | 0 | 2 | $(1,24,60)$ | 85 |

Each tail is one-ended, so its one attachment adds no graph cycle.

At L3, put $A(S)=e+\Phi S+\Theta S^2$ and $B(S)=8/3+KS^2$ and use
$Q=\delta^3$, $s=S$, $p=\delta^{-1}P$, $G=\delta F_Q$, and
$\rho=\delta x$. The exact CTH chart is

~~~text
x^4G=(rho S^2-1)[x^3-A(S)-rho B(S)+3rho^2]-rho^4S^2/2. (13)
~~~

Again the prefactor is a unit. The two critical selectors are

~~~text
c1=8/3+Ka^2,             c2=-3-K^2a^2/Theta.            (14)
~~~

If $c_1=0$, then $c_2=-99/43$, so only depths one and two occur. On the
target base $\delta=t^2$, $Q=t^6$, their primitive $(t,x,y)$ rows are
$(2,2,3)$ and $(1,4,6)$, with form orders 13 and 11.
The first tail is elliptic; the second is rational with persistent delta
one. Both are one-ended.

## 6. The even A3 wall

On $\Theta=u=\eta=\alpha=0$, $\Delta\ne0$, L4 doubles the EDELTA
quadratic. Use $Q=\delta^4$, $s=\delta^{-1}S$, $p=\delta^{-1}P$,
$G=\delta^2F_Q$, and put $x=P^{-1}$, $v=S/P$, $\rho=\delta x$. The exact
chart is

~~~text
x^6G=(v^2-rho)[x^4-(Delta+Phi v+Kv^2)-e rho
               -(8/3)rho^2+3rho^3]-rho^4v^2/2.          (15)
~~~

The double root is interior and the prefactor is a unit. The critical value
begins

~~~text
psi(rho)=e rho+(8/3)rho^2-3rho^3+(1/2)rho^4+O(rho^5),  (16)
~~~

so depth one is forced. On the target base $\delta=t^3$, so $Q=t^{12}$,
the primitive row is

~~~text
x=tX, y=t^2Y,
Y^2=X(X^3-e).                                           (17)
~~~

This elliptic tail has form order 12 and two marked ends. Both ends attach
to the same connected BDELTA complement. The collision complement has
$(V,E,b_1)=(7,12,6)$; adding the tail and two attachments gives
$(8,14,7)$. Thus $b_1+g_{\rm tail}=7+1=8$, exactly the fan genus.

## 7. The two-branch A15 wall

On $\alpha=u=\eta=0$, $\Delta\Theta\ne0$, impose
$\Theta=-\Delta$. In the target-compatible BDT chart

~~~text
Q=sigma^24, s=sigma^-3S, p=sigma^-6P,
G=sigma^6 F_Q, x=S^-1, w=x^2P-1, rho=sigma^3x,
Frec=-x^10G,
~~~

the literal source becomes

~~~text
Frec=w[x^8-Delta w(1+w)^3-Phi rho(1+w)^3
       -rho^2{K(1+w)^2+e(1+w)^3}
       -(8/3)rho^4(1+w)^2+3rho^6(1+w)]
     +rho^8/2.                                           (18)
~~~

At $\rho=0$ the two smooth branches meet with quotient
$\mathbb C[[x,w]]/(w,x^8)$, hence contact length eight and type A15.
After normalizing the contact, one branch belongs to the rational component
$R:P=S^2$ and the other to the genus-two component $C$, already joined to
the rational T face. The complement therefore has two components, $R$ and
$C$--T.

If $x_\pm$ denote the two discriminant branches, put
$\chi_\pm=x_\pm^8$. Then

~~~text
chi_plus =Phi rho+(K+e)rho^2+(8/3+lambda)rho^4+O(rho^5),
chi_minus=Phi rho+(K+e)rho^2+(8/3-lambda)rho^4+O(rho^5),
lambda^2=-2Delta.                                       (19)
~~~

The exact depths are consequently $1,2,4$ according as $\Phi\ne0$;
$\Phi=0,K+e\ne0$; or $\Phi=0,K+e=0$. In the last row the seam forces

~~~text
Delta=2048/45, K=1376/135, Theta=-2048/45,
(8/3+lambda)(8/3-lambda)=1472/15!=0,                    (20)
~~~

so there is no deeper case.

| depth | normalized packet | persistent delta | graph $b_1$ | form order |
|---:|---|---:|---:|---:|
| 1 | two rational signs, 7 mutual nodes | 1 | 6 | 116 |
| 2 | two rational signs, 6 mutual nodes | 2 | 5 | 16 |
| 4 | one smooth genus-three tail | 4 | 0 | 14 |

At depth four the tail is

~~~text
Y^2=(X^4-8/3)^2+4096/45,                                (21)
~~~

whose degree-eight branch polynomial is squarefree. One marked end lands on
$R$ and the other on $C$. The two attachments join the complementary
components and create no graph cycle. In every row,

~~~text
2+g_tail+b1+persistent_delta=9.                          (22)
~~~

The exact residue is a unit times $\sigma^{17}x^6dx/y$; the primitive
orders in the table are $17\operatorname{ord}(\sigma)-\operatorname{ord}(x)$.
This direct calculation is necessary because the raw $x^6$ numerator does
not reach the A15 conductor $x^8$.

## 8. Attachment rank and representation loss

If finite graphs $P$ and $C$ receive $k$ attachment edges without vertex
identification, Euler's identity gives

~~~text
b1(final)=b1(P)+b1(C)+k-c(P)-c(C)+c(final).              (23)
~~~

For A3, both packet and complement are connected and $k=2$, so (23) adds
one. For A15, the complement has two components and the final graph is
connected, so it adds zero. The endpoint address, not parity alone, decides
the graph contribution.

This also explains why a simple tournament is the wrong carrier. Collapsing
parallel nodes or forgetting which complementary component receives an end
can erase graph cycles. The precise cross-problem research analogy with LRC
owner/arrival states is recorded in the associated reflection; it supplies
no LRC implication.

## 9. Proper-flat extinction and audit

After a finite common base change and toric regularization, every special
component on (2) is rational or one of the positive-genus carriers and tails
listed above. Rational curves map constantly to the good elliptic target;
positive form order makes every higher-genus component map constantly.
Keeping actual positive component multiplicities on one common dominating
base, the inherited proper-flat identity gives

~~~text
deg(phi_generic^*L)=sum_Gamma m_Gamma deg(phi_Gamma^*L)=0,             (24)
~~~

contradicting the positive generic response degree of a nonautomorphic Keller
pair. This proves (2), and THM-4354 supplies the $\alpha\ne0$ half of (3).

The 769-check primary and 774-check referee independently reconstruct the
support/fan census, face polygons and normalizations, clean graphs, all four
local charts, depths, primitive target orders, attachment identities, and
genus conservation. The referee explicitly tests the false A15
same-complement rule and the coefficientwise-specialization and
source-base/target-base hostiles. Both ordinary and optimized runs byte-match
the frozen outputs in the frontmatter.

The remaining endpoint is $U=K=0$. Seam entry, JC(2), and DC(2) remain open.
