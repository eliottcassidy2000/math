---
id: THM-4356
title: "Double-zero endpoint planar Jacobian extinction"
status: >
  PROVED RELATIVE TO THM-4230/4327/4341/4344/4352 + VERIFIED-EXACT +
  INDEPENDENTLY HOSTILE-AUDITED. The inherited exact-weight-twelve gate
  Z=beta_11=zeta_3=W=xi_10=U=K=0, hence Delta=5696/105, is extinct. Its 48
  exact supports have fifteen fans and twelve face planes. Exactly three
  within-owner collision divisors occur: an odd A4, a two-branch A15, and a
  smooth boundary tangency. Every positive-genus carrier or collision tail
  has positive primitive form order. Together with THM-4350/4351/4353/4354/
  4355 this closes the full endpoint Z=beta_11=zeta_3=W=xi_10=0. Seam entry,
  JC(2), and DC(2) are not claimed.
source: root + jc-u0-k0-scout + independent clean-room referee / next-sharp session, 2026-09-02
depends_on:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4327-generic-exact-weight-twelve-endpoint-wall-extinction
  - THM-4341-odd-self-similar-cusp-reciprocal-tail-duality
  - THM-4344-clean-cubic-infinity-exit-planar-jacobian-extinction
  - THM-4352-even-self-similar-cusp-reciprocal-parity-and-attachment-law
related:
  - THM-4350-first-normal-owner-cubic-infinity-exit-planar-jacobian-extinction
  - THM-4351-double-normal-owner-zero-cubic-infinity-exit-planar-jacobian-extinction
  - THM-4353-simultaneous-zero-endpoint-planar-jacobian-extinction
  - THM-4354-first-normal-owner-u-zero-endpoint-planar-jacobian-extinction
  - THM-4355-terminal-alpha-zero-u-zero-endpoint-planar-jacobian-extinction
mistake_firewall:
  - MISTAKE-522
  - MISTAKE-531
  - MISTAKE-540
primary_script: 04-computation/jc2_m12_double_zero_endpoint_extinction_thm4356.py
primary_output: 05-knowledge/results/jc2_m12_double_zero_endpoint_extinction_thm4356.out
primary_script_sha256: c96e4a4de6d58c29c87061fc072e9f0f5bf476ae215cd6320690ccb308337d9c
primary_output_sha256: 5fb0a8dc2e3900ac655380cf4315b650fd01647231b8abd239cdc8fcbcb88b6b
referee_script: 04-computation/jc2_m12_double_zero_endpoint_extinction_independent_referee_thm4356.py
referee_output: 05-knowledge/results/jc2_m12_double_zero_endpoint_extinction_independent_referee_thm4356.out
referee_script_sha256: 558f4331e9e8860e2bea199073ca18e8867f30439a8dbac5ae91590a430f5118
referee_output_sha256: f05c0d23b39a608a6fc24ab3472c68f2d3f445ab1bd0ac1bb8fc18965d0a5e1a
hash_basis: raw LF bytes
audit: >
  PASS. The 625-check primary and independently written 956-check referee
  reconstruct the literal source, all 48 exact supports and fifteen fans,
  separate hostile atlases, every face and edge, nineteen global polygon/
  packet ledgers, carrier normalizations, graphs, primitive orders, the A4
  and A15 critical-depth packets, and the smooth N-wall blowups. The referee
  explicitly audits coefficientwise specialization, source-minimal versus
  target-compatible bases, and the A15 endpoint-to-component addresses.
  Normal and optimized executions byte-match both frozen outputs.
---

# THM-4356 -- Double-zero endpoint planar Jacobian extinction

**PROVED RELATIVE TO THM-4230/4327/4341/4344/4352 + VERIFIED-EXACT +
INDEPENDENTLY HOSTILE-AUDITED. THE DISPLAYED DOUBLE-ZERO GATE IS EXTINCT.
SEAM ENTRY, `JC(2)`, AND `DC(2)` REMAIN OPEN.**

## 1. Statement and inheritance

Work over `C` in the inherited `b=d=0` reduced `(2,3)`, exact-weight-twelve
seam. Put `y=sp`, `e=-1376/135`, `u=upsilon_5`,
`alpha=alpha_11`, and

```text
H=-3p+(8/3)p^2+e p^3+K y^2+Phi p^2y+Delta p^4
  +Theta py^2+eta p^3y+zeta_3 y^3+u p^5
  +xi_10 p^2y^2+alpha p^4y+beta_11 py^3
  +U p^6+W p^3y^2+Z y^4,

K=2848/45-(7/6)Delta,
F_Q=(s^2-p)(1-QH)-Q s^2/2.                              (1)
```

Impose

```text
Z=beta_11=zeta_3=W=xi_10=U=K=0,
Delta=D=5696/105.                                        (2)
```

Then no nonautomorphic planar Keller pair lies on `(1)--(2)`, with
`alpha,Theta,Phi,eta,u` arbitrary. The conclusion is relative to the
inherited toroidal, good-target, and proper-flat interfaces. With
THM-4350/4351/4353/4354/4355, this closes every coefficient choice on

```text
Z=beta_11=zeta_3=W=xi_10=0.                              (3)
```

It does not supply entry into this endpoint, into the exact-weight-twelve
seam, or prove `JC(2)` or `DC(2)`.

The inheritance pass is:

- closest proved mechanisms: THM-4353's `K=0,U!=0` owner tree and smooth
  `N` tangency, and THM-4354/4355's `U=0,K!=0` carrier and cusp packets;
- canonical hostile: on `alpha=u=eta=0,Theta=-D`, eight intersections run
  to one boundary point, but the two stable ends land on different
  normalized complementary components;
- corrected near miss: setting `U=0` or `K=0` in an old face equation does
  not preserve face ownership, because roots can escape and components can
  merge;
- least-used sidecar: the normalized-component address of every tail end.

The live board was

```text
double owner deletion | support quotient | face merger | boundary-root motion
repeated edge | carrier genus/order | attachment address | natural rank.     (4)
```

## 2. Literal support quotient and reversible rank

Use lifted support coordinates `(deg_s,deg_p,deg_Q)`. Collecting equal
monomials before specialization gives exactly twenty possible points:

```text
(0,1,0):-1; (2,0,0):1; (2,0,1):-1/2;
(0,2,1):-3; (2,1,1):3;
(0,3,1):8/3; (2,2,1):-8/3;
(0,4,1):e; (2,3,1):-e;
(0,5,1):D; (2,4,1):Theta-D; (4,3,1):-Theta;
(0,6,1):u; (2,5,1):-u;
(1,4,1):Phi; (3,3,1):-Phi;
(1,5,1):eta; (3,4,1):-eta;
(1,6,1):alpha; (3,5,1):-alpha.                         (5)
```

Thus `Theta` has three support states: `0`, `D`, and nonzero unequal to
`D`. The other four parameters contribute independent presence bits, so

```text
3*2^4=48                                                     (6)
```

distinct supports occur. They split into `24` with `alpha!=0`, `16` with
`alpha=0,Theta!=0`, and `8` with `alpha=Theta=0`.

Let `c=0,1,2` label the three `Theta` support states. The address

```text
n=1+16c+8[alpha!=0]+4[eta!=0]+2[Phi!=0]+[u!=0]          (7)
```

bijects the supports with `{1,...,48}`; both certificates check its inverse.
An odd-square display retains exactly the same ordinal because

```text
((2n-1)^2-1)/8=n(n-1)/2=T(n-1).                         (8)
```

The rank does not encode a continuous collision divisor; the collision word
remains a separate sidecar.

## 3. The fifteen exact fans

The 48 supports have exactly fifteen fan words:

| exact face word | supports |
|---|---:|
| `BU` | 1 |
| `BDT/BDelta` | 5 |
| `C4a,Balpha` | 1 |
| `C5,Balpha` | 1 |
| `C4e,Beta` | 1 |
| `BU,F5` | 9 |
| `BU,EUPhi` | 1 |
| `BDelta,EDPhi` | 1 |
| `C4a,Balpha,E11` | 8 |
| `C4a,Balpha,N` | 3 |
| `C5,Balpha,E11` | 8 |
| `C5,Balpha,N` | 3 |
| `C4e,Beta,F5` | 4 |
| `C4e,Beta,N` | 1 |
| `BU,F5,N` | 1 |

Here `BDT` and `BDelta` share one plane, as do `FUT,FET,FUE`, denoted
collectively by `F5`. Consequently the fifteen words use twelve distinct
planes. The primary hostile atlas independently toggles five source rows and
two aggregate cancellations. It has `128` keys, `96` distinct supports,
`23` fans, and `17` planes. The referee instead toggles lifted points and
obtains `1024` supports, `52` fans, and `24` planes. Every exact support
embeds in both over-atlases; neither atlas asserts realizability of its extra
states.

## 4. Complete face, component, and order ledger

For a supporting plane `(a,b,c)`, let `L` clear both its source denominator
and the target sixth root. The target-compatible form order is
`L(5/6-a-b-c)`.

| face | plane | initial equation, up to a torus monomial | geometry | base/order |
|---|---|---|---|---:|
| `C4a` | `(-1/4,1/4,-1/4)` | `P(DP^4+alpha SP^5-1)` | rational | `12/13` |
| `C5` | `(0,1/5,-1/5)` | `P(uP^5+alpha SP^5-1)` | rational | `30/25` |
| `C4e` | `(0,1/4,-1/4)` | `P(DP^4+eta SP^4-1)` | rational | `12/10` |
| `Balpha` | `(1/11,2/11,-2/11)` | `(P-S^2)(alpha SP^5-1)` | two rational, 11 nodes | `66/49` |
| `BU` | `(1/10,1/5,-1/5)` | `(P-S^2)(uP^5-1)` | six rational, 10 nodes | `30/22` |
| `Beta` | `(1/9,2/9,-2/9)` | `(P-S^2)(eta SP^4-1)` | two rational, 9 nodes | `18/13` |
| `BDT` | `(1/8,1/4,-1/4)` | `(P-S^2)(DP^4+Theta S^2P^3-1)` | rational plus genus 2, 8 nodes off `L15` | `24/17` |
| `BDelta` | same | `(P-S^2)(DP^4-1)` | five rational, 8 nodes | `24/17` |
| `FUT` | `(1/5,1/5,-2/5)` | `S^2(1-uP^5-eta SP^4-Theta S^2P^3)` | genus 2 off `L5` | `30/25` |
| `FET` | same | `S^2(1-eta SP^4-Theta S^2P^3)` | genus 2 | `30/25` |
| `FUE` | same | `S^2(1-uP^5-eta SP^4)` | rational | `30/25` |
| `EDPhi` | `(1/4,1/4,-1/2)` | `S^2(1-DP^4-Phi SP^3)` | rational | `12/10` |
| `E11` | `(2/7,1/7,-4/7)` | `S^2(1-alpha SP^5-Theta S^2P^3)` | genus 3 | `42/41` |
| `EUPhi` | `(2/5,1/5,-4/5)` | `S^2(1-uP^5-Phi SP^3)` | rational | `30/31` |
| `N` | `(1,0,-2)` | `S^2[1-SP^3(Phi+eta P+alpha P^2)]` | rational | `6/11` |

The only positive-genus normalization identities are

```text
E11:     Y^2=P(alpha^2 P^7+4Theta),
FUT:     Y^2=P(4Theta+(eta^2-4uTheta)P^5),
BDT-C:   Y^2=Theta P(1-DP^4).                           (9)
```

They are squarefree on their clean owner cells and give genera `3,2,2`.
The `u=0,eta!=0` specialization `FET` remains automatically squarefree.
The other displayed faces have an explicit rational parametrization or
factor into the stated rational components. In particular, their polygon
interior counts are not substituted for normalization genus.

Away from the collision divisors, the distinct graph ledgers are:

| stratum | `(V,E,b1)` | carrier genus | total genus |
|---|---:|---:|---:|
| `alpha*Theta!=0` | `(4,13,10)` | 3 | 13 |
| `alpha!=0,Theta=0`, without/with `N` | `(3,12,10)` / `(4,13,10)` | 0 | 10 |
| `alpha=0,Theta*u!=0` | `(7,15,9)` | 2 | 11 |
| `alpha=u=0,Theta*eta!=0` | `(4,11,8)` | 2 | 10 |
| `alpha=u=eta=0,Theta!=0` | `(2,8,7)` | 2 | 9 |
| terminal `u*eta!=0`, without/with `N` | `(7,15,9)` / `(8,16,9)` | 0 | 9 |
| terminal `u!=0,eta=0,Phi!=0` | `(7,15,9)` | 0 | 9 |
| terminal `u!=0,eta=Phi=0` | `(6,10,5)` | 0 | 5 |
| terminal `u=0,eta!=0`, without/with `N` | `(3,10,8)` / `(4,11,8)` | 0 | 8 |
| terminal `u=eta=0,Phi!=0` | `(6,12,7)` | 0 | 7 |
| terminal `u=eta=Phi=0` | `(5,8,4)` | 0 | 4 |

For every support, `E-V+1=b1` and `b1+sum(g_i)` equals the global Newton
polygon's interior count. There are nineteen coefficient-sensitive global
polygons/packets, and each also satisfies

```text
sum_j(packet_j-1)=2I-2.                                  (10)
```

The large graph genera are literal: `Balpha`, `BU`, `Beta`, and `BDelta`
retain respectively `11,10,9,8` separately addressed transverse nodes.
Collapsing parallel nodes to one simple edge would destroy the proof ledger.

## 5. Collision completeness

The exact edge audit leaves precisely three non-owner repeated words:

| face | variable edge word | divisor | local type |
|---|---|---|---|
| `FUT` | `Theta+eta X+uX^2` | `L5=eta^2-4uTheta` | odd `A4` |
| `BDT` | `(X-1)(DX+Theta)` | `L15=D+Theta` | two-branch `A15` |
| `N` | `Phi+eta X+alpha X^2` | `LN=eta^2-4alpha Phi` | smooth boundary tangency |

All other nonlinear edge discriminants are nonzero away from deletion of
their named owner. The three rows occupy disjoint owner cells, so there is no
unlisted simultaneous collision.

## 6. The `A4` wall at `K=0`

Assume `alpha=0`, `uTheta!=0`, and `L5=0`. Put

```text
a=-eta/(2Theta),             u=Theta a^2,             a!=0. (11)
```

In the source-minimal `FUT` chart

```text
Q=sigma^5, s=sigma^-1 S, p=sigma^-1 P, G=sigma^2 F_Q,
x=P^-1, v=S/P, rho=sigma x,
```

the literal equation is

```text
x^7G=(v^2-rho)[x^5-A(v)-rho B(v)-e rho^2
                -(8/3)rho^3+3rho^4]-rho^5v^2/2,
A(v)=Theta(v-a)^2,             B(v)=D+Phi v.             (12)
```

The prefactor is a unit at the collision, and Morse preparation gives
`y^2=x^5-psi(rho)`. Its successive critical selectors are

```text
c1=D+Phi a,          c2=e-Phi^2/(4Theta),          c3=8/3. (13)
```

Thus only depths `1,2,3` occur; the depth-four row from the neighboring
`K!=0` problem cannot survive. Every depth is nonempty. On the common target
base `Q=t^30`, the stable packets are

| depth | tail genus | persistent delta | primitive `(ord t,ord x,ord y)` | form order |
|---:|---:|---:|---:|---:|
| 1 | 2 | 0 | `(4,6,15)` | 115 |
| 2 | 1 | 1 | `(1,4,10)` | 35 |
| 3 | 1 | 1 | `(2,18,45)` | 95 |

The reduced branch polynomials are respectively
`X(X^4-c)`, `X^3-c`, and `X(X^2-c)`, with `c!=0`. Every tail is one-ended,
so it creates no graph cycle, and

```text
9+g_tail+delta_persistent=11.                            (14)
```

## 7. The two-branch `A15` wall

Assume

```text
alpha=u=eta=0,                 Theta=-D.                  (15)
```

The two `BDT` components are

```text
R: P=S^2,               C: 1-DP^4+DS^2P^3=0.             (16)
```

The second is smooth of genus two. Off `(15)`, their eight intersections
satisfy `(D+Theta)S^8=1`, so all eight run to `S=infinity` at the wall. In
the target-compatible chart

```text
Q=tau^24, s=tau^-3 S, p=tau^-6 P, G=tau^6 F_Q,
x=S^-1, z=x^2P-1, rho=tau^3x,
```

the exact reciprocal equation is

```text
Frec=z[x^8-Dz(1+z)^3-Phi rho(1+z)^3-e rho^2(1+z)^3
        -(8/3)rho^4(1+z)^2+3rho^6(1+z)]+rho^8/2.         (17)
```

At `rho=0`, this is `z[x^8-Dz(1+z)^3]`; the branch intersection length is
eight, hence type `A15`. If `lambda^2=-2D`, the two critical `x^8` branches
begin

```text
chi_+=Phi rho+e rho^2+(8/3-lambda)rho^4+O(rho^5),
chi_-=Phi rho+e rho^2+(8/3+lambda)rho^4+O(rho^5).        (18)
```

Since `e!=0`, the exact depth is one when `Phi!=0` and two when `Phi=0`.
The exceptional equations and graph ledgers are

| depth | exceptional equation | persistent delta | graph `(V,E,b1)` | form order |
|---:|---|---:|---:|---:|
| 1 | `D Z^2-X(X^7-Phi)Z=0` | 1 | `(4,9,6)` | 116 |
| 2 | `D Z^2-X^2(X^6-e)Z=0` | 2 | `(4,8,5)` | 16 |

Each packet has two rational sign components. They meet at respectively
seven and six simple mutual nodes. Direct reciprocal-chart substitution
identifies the zero sign with `R` and the nonzero sign with `C`; the two ends
therefore attach to different complementary components. The attachments
connect those components and add no graph cycle. In both rows,

```text
2+delta_persistent+b1=9.                                 (19)
```

The raw residue has only an `x^6` numerator, below the `A15` conductor
`x^8`; the positive orders `116,16` come from the direct target-compatible
valuation, not a conductor shortcut.

## 8. The `N` wall is smooth

This wall occurs only when `Theta=0,alpha*Phi!=0`. Write

```text
J(P)=Phi+eta P+alpha P^2=alpha(P-a)^2.                   (20)
```

The exact source-minimal reciprocal chart is

```text
H=R B-P^3J(P)C,
C=1-sigma^2 P R^2,
B=C(1-sigma D0(P))-sigma/2,
D0(P)=-3P+(8/3)P^2+eP^3+DP^4+uP^5.                     (21)
```

At `(sigma,R,P)=(0,0,a)`, `H_R=1`. Hence the curve has one smooth branch
`R=unit*(P-a)^2`, not an `A_k` collision tail. For the first ordinary
blowup, the complementary exceptional chart is empty; after the second, the
complementary chart contains one reduced transverse point. Both exceptional
divisors are rational. For both adjacent caps the puncture packet changes

```text
(10,5,4,2,2,1) -> (10,5,4,3,1),                        (22)
```

while its defect remains `18`. This is the one neighboring local model that
specializes without changing its mechanism.

## 9. Ownership under specialization

The literal reconstruction, rather than coefficientwise substitution in an
old fan, is essential:

- the six roots of `1-UP^6` run to `P=infinity` as `U` vanishes;
- the old `D6` and `E10` equations merge with the non-`R` component of
  `Balpha`;
- `Eeta` merges with `Beta`, and `EU` with the five vertical components of
  `BU`;
- the old rational `T` leaf is pruned when `K` vanishes;
- `EDelta` becomes rational `EDPhi`, or merges with `BDelta` at `Phi=0`;
- the old `A3` equation `Phi^2=4KDelta` becomes the owner wall `Phi=0`, and
  the old `A2` row requiring `Delta=0` is impossible because `D!=0`.

Thus the preserved object is the normalized labelled multigraph with root
and endpoint addresses, not the list of specialized face equations and not a
simple tournament.

## 10. Proper-flat extinction and audit

After one finite common base change and toric regularization, every special
component on `(2)` is rational or one of the carriers and tails listed above.
Rational curves map constantly to the good elliptic target. Every
positive-genus component has strictly positive pullback-form order and also
maps constantly. Retaining the actual positive multiplicities on one common
dominating base, the inherited proper-flat identity gives

```text
deg(phi_generic^*L)=sum_Gamma m_Gamma deg(phi_Gamma^*L)=0,             (23)
```

contradicting the positive generic response degree of a nonautomorphic
Keller pair. This proves `(2)`.

The 625-check primary and 956-check clean-room referee independently
reconstruct the source, supports, fans, faces, edge completeness, graph and
packet ledgers, all local charts and depths, primitive target orders, and
the `R/C` attachment addresses. Their deliberately different hostile atlases
also test the false coefficientwise-specialization and same-complement rules.
Normal and optimized executions byte-match the frozen outputs named in the
frontmatter.

Together with the four complementary `(U,K)` cases supplied by
THM-4350/4351, THM-4353, and THM-4354/4355, this closes `(3)`. Exact-seam
entry, `JC(2)`, and `DC(2)` remain open.
