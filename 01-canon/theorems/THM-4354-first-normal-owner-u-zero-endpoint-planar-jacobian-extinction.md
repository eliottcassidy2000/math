---
id: THM-4354
title: "First-normal-owner U-zero endpoint planar Jacobian extinction"
status: >
  PROVED RELATIVE TO THM-4230/4327/4344 + VERIFIED-EXACT + INDEPENDENTLY
  HOSTILE-AUDITED. The inherited exact-weight-twelve gate
  Z=beta_11=zeta_3=W=xi_10=U=0 with K*alpha_11!=0 is extinct. Its 64 exact
  supports have six fans; the reducible middle face has eleven transverse
  nodes and graph genus ten, while the unique smooth genus-three carrier has
  positive form order. Seam entry, JC(2), and DC(2) are not claimed.
source: root + jc-deep-u0 / planar-Jacobian next-sharp session, 2026-09-02
depends_on:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4327-generic-exact-weight-twelve-endpoint-wall-extinction
  - THM-4344-clean-cubic-infinity-exit-planar-jacobian-extinction
related:
  - THM-4340-u-zero-cubic-edge-planar-jacobian-extinction
  - THM-4350-first-normal-owner-cubic-infinity-exit-planar-jacobian-extinction
  - THM-4353-simultaneous-zero-endpoint-planar-jacobian-extinction
mistake_firewall:
  - MISTAKE-522
  - MISTAKE-531
  - MISTAKE-540
primary_script: 04-computation/jc2_m12_first_normal_owner_u_zero_endpoint_extinction_thm4354.py
primary_output: 05-knowledge/results/jc2_m12_first_normal_owner_u_zero_endpoint_extinction_thm4354.out
primary_script_sha256: 630909a2819ef809a049ce84effd23390d7e1b1709961becb275632ee173cdd0
primary_output_sha256: f8c5c28b52e97a0d2e8e782e31fbbbea10bf70a5b22c00337919bb942aa018f2
referee_script: 04-computation/jc2_m12_first_normal_owner_u_zero_endpoint_extinction_independent_referee_thm4354.py
referee_output: 05-knowledge/results/jc2_m12_first_normal_owner_u_zero_endpoint_extinction_independent_referee_thm4354.out
referee_script_sha256: b324a2b6562d0b4c060aabf7111f80e7ca62274c8c04cb4c722a29ca9218b5c7
referee_output_sha256: 32657943b9cf7055c739ae0b2a8d067e9f6037c64e89d36c3733363e00323af0
hash_basis: raw LF bytes
audit: >
  PASS. The 880-check primary and independently written 870-check referee
  reconstruct the literal source, exact and hostile support atlases, all
  primitive faces and form orders, the eleven-node middle multigraph, both
  carrier normalizations, edge reductions, global Pick/packet/graph ledgers,
  and the closed natural-address bijection. Normal and optimized runs
  byte-match both frozen outputs.
---

# THM-4354 -- First-normal-owner U-zero endpoint planar Jacobian extinction

**PROVED RELATIVE TO THM-4230/4327/4344 + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED. `JC(2)`, `DC(2)`, AND SEAM ENTRY REMAIN OPEN.**

## 1. Statement and inheritance

Work over `C` in the inherited `b=d=0` reduced `(2,3)`, exact-weight-twelve
seam, with the complete sixteen-row source

```text
H=-3p+(8/3)p^2-(1376/135)p^3+K y^2+Phi p^2y
  +Delta p^4+Theta py^2+eta p^3y+zeta_3 y^3
  +u p^5+xi_10 p^2y^2+alpha p^4y+beta_11 py^3
  +U p^6+W p^3y^2+Z y^4,

K=2848/45-(7/6)Delta,          u=upsilon_5,
alpha=alpha_11.                                           (1)
```

Impose

```text
Z=beta_11=zeta_3=W=xi_10=U=0,       K*alpha!=0.           (2)
```

No nonautomorphic planar Keller pair lies
on `(1)--(2)`, with `Delta,Theta,Phi,eta,u` arbitrary subject to the seam
relation. The conclusion is relative to the inherited toroidal, good-target,
and proper-flat interfaces. It proves neither seam entry nor `JC(2)` or
`DC(2)`.

The inheritance pass is:

- closest proved mechanism: THM-4350's first-normal-owner replacement faces;
- canonical hostile: the six roots of `1-UP^6` do not remain finite when
  `U` vanishes, so the old `M` face and its six attachments cannot be
  specialized coefficientwise;
- corrected near miss: a face polygon with ten interior points is not a
  genus-ten component when its initial polynomial factors;
- least-used sidecar: the eleven separately addressed intersections of the
  two rational middle components.

The Anchor / Niche / Wildcard board is

```text
U-zero cap | first normal owner | reducible middle | mu_11 nodes
carrier replacement | graph genus | fan grid | natural support address.     (3)
```

## 2. Exact support quotient and the `3 x 2` fan grid

Put

```text
F_Q=(s^2-p)(1-QH)-Q*s^2/2.                              (4)
```

On `(2)`, the three multiply visible support coefficients are

```text
(2,3,1)=K+1376/135,
(2,4,1)=Theta-Delta,
(2,5,1)=-u.                                             (5)
```

The first one vanishes exactly at

```text
Delta_*=3968/63,             K=-1376/135.               (6)
```

The pair `(Delta,Theta)` therefore has eight support classes, represented by

```text
(0,0),(0,1),
(Delta_*,0),(Delta_*,Delta_*),(Delta_*,2Delta_*),
(1,0),(1,1),(1,2).                                     (7)
```

Here the last three pairs represent the generic nonzero-`Delta` classes;
the numerical representative `1` is not an extra coefficient restriction.
Combining `(7)` with the independent presence bits of `u,Phi,eta`, while
`alpha` is forced present, gives exactly

```text
8*2^3=64                                                (8)
```

distinct realizable supports.

The lower fan depends only on the first live left cap owner

```text
L=5 if u!=0;       L=4 if u=0,Delta!=0;
L=3 if u=Delta=0,                                      (9)
```

and on whether `Theta` is nonzero. The exact six-fan census is:

| left cap | `Theta!=0` | `Theta=0` |
|---|---:|---:|
| `C5`, `u!=0` | 20: `C5,B,E11,T` | 12: `C5,B,E10` |
| `C4`, `u=0,Delta!=0` | 16: `C4,B,E11,T` | 8: `C4,B,E10` |
| `C3`, `u=Delta=0` | 4: `C3,B,E11,T` | 4: `C3,B,E10` |

A hostile atlas forces the `alpha` row present but independently toggles the
other five optional rows and three multiply visible deletions. Its 256 keyed
configurations collapse to 168 distinct supports. All 64 exact supports embed
literally, and all hostile supports have the same six fans, in keyed counts
`64,64,32,32,32,32`. This is a fan-stability over-atlas, not a coefficient
realizability census.

## 3. Exact faces, primitive charts, and form orders

For a supporting plane `(a,b,c)`, the complete face list, up to invertible
torus monomials, is:

| face | plane | polygon | initial equation | geometry | target base/order |
|---|---|---|---|---|---:|
| `C5` | `(0,1/5,-1/5)` | `(0,1),(1,6),(0,6)` | `P(uP^5+alpha SP^5-1)` | rational | `30/25` |
| `C4` | `(-1/4,1/4,-1/4)` | `(0,1),(1,6),(0,5)` | `P(Delta P^4+alpha SP^5-1)` | rational | `12/13` |
| `C3` | `(-2/3,1/3,-1/3)` | `(0,1),(1,6),(0,4)` | `P(eP^3+alpha SP^5-1)` | rational | `6/9` |
| `B` | `(1/11,2/11,-2/11)` | `(0,1),(2,0),(3,5),(1,6)` | `(P-S^2)(alpha SP^5-1)` | two rational components, 11 nodes | `66/49` |
| `E11` | `(2/7,1/7,-4/7)` | `(2,0),(4,3),(3,5)` | `S^2(1-alpha SP^5-Theta S^2P^3)` | smooth genus 3 | `42/41` |
| `T` | `(1/2,0,-1)` | `(2,0),(4,2),(4,3)` | `S^2(1-S^2P^2(K+Theta P))` | rational | `6/8` |
| `E10` | `(3/8,1/8,-3/4)` | `(2,0),(4,2),(3,5)` | `S^2(1-alpha SP^5-KS^2P^2)` | smooth genus 3 | `24/26` |

Here `e=-1376/135`. The density formula

```text
5/6-(a+b+c)                                             (10)
```

gives the positive orders `25,13,9,49,41,8,26` on the
displayed common target bases.

The literal primitive source charts are

```text
C5:  Q=sigma^5,  s=S,          p=sigma^-1 P, G=sigma F_Q;
C4:  Q=sigma^4,  s=sigma S,    p=sigma^-1 P, G=sigma F_Q;
C3:  Q=sigma^3,  s=sigma^2 S,  p=sigma^-1 P, G=sigma F_Q;
B:   Q=sigma^11, s=sigma^-1 S, p=sigma^-2 P, G=sigma^2 F_Q;
E11: Q=sigma^7,  s=sigma^-2 S, p=sigma^-1 P, G=sigma^4 F_Q;
T:   Q=sigma^2,  s=sigma^-1 S, p=P,          G=sigma^2 F_Q;
E10: Q=sigma^8,  s=sigma^-3 S, p=sigma^-1 P, G=sigma^6 F_Q. (11)
```

Direct substitution of the complete source `(1)` gives exactly the seven
initials above; the cap factor `P` is a torus monomial and not a second cap
component.

The sole positive-genus carrier has one of the exact normalizations

```text
E11: Y=2Theta*S*P^2+alpha*P^4,
     Y^2=P(alpha^2P^7+4Theta);                           (12)

E10: Y=2K*S*P+alpha*P^4,
     Y^2=alpha^2P^8+4K.                                (13)
```

The branch polynomials in `(12)--(13)` are squarefree under their displayed
owner gates. Both smooth projective normalizations have genus three. No
carrier collision divisor remains inside `(2)`.

## 4. Eleven middle nodes, edges, and genus conservation

Write the two `B` components as

```text
R: P=S^2,                  C: alpha*S*P^5=1.             (14)
```

Their intersection is

```text
alpha*S^11=1.                                            (15)
```

Because `alpha!=0` and the characteristic is zero, `(15)` has eleven
distinct roots. Differentiating along `R` gives `11alpha*S^10!=0` at every
root, so all eleven intersections are transverse nodes. Hence the two
rational components of `B` alone have graph genus

```text
11-2+1=10.                                              (16)
```

A complete reduced edge word is

```text
C5:  alpha*X-1; alpha+uX;       u-X^5;
C4:  alpha*X-1; alpha+Delta X;  Delta-X^4;
C3:  alpha*X-1; alpha+eX;       e-X^3;
B:   X-1; 1-alpha X; alpha(X-1); alpha-X;
E11: 1-Theta X; -Theta-alpha X; X-alpha;
T:   1-KX^2; -K-Theta X; X-Theta;
E10: 1-KX^2; -K-alpha X; X-alpha.                       (17)
```

Every polynomial in `(17)` is reduced on its exact owner gate. The cap meets
`B` in one node, `B` meets the carrier in one node, and, when `Theta!=0`,
`E11` meets the rational `T` leaf in one node. Therefore

```text
Theta!=0: (V,E,b1)=(5,14,10),
Theta=0:  (V,E,b1)=(4,13,10).                           (18)
```

Adding the unique genus-three carrier gives total genus `13` in every fan.
The six independent global ledgers are:

| cap / right owner | global polygon | `(2A,B,I)` | packet |
|---|---|---:|---:|
| `C5 / Theta` | `(0,1),(2,0),(4,2),(4,3),(3,5),(1,6),(0,6)` | `(36,12,13)` | `(10,8,5,3,2,2,1)` |
| `C5 / 0` | `(0,1),(2,0),(4,2),(3,5),(1,6),(0,6)` | `(35,11,13)` | `(10,10,5,2,2,1)` |
| `C4 / Theta` | `(0,1),(2,0),(4,2),(4,3),(3,5),(1,6),(0,5)` | `(35,11,13)` | `(10,8,5,3,2,2,1)` |
| `C4 / 0` | `(0,1),(2,0),(4,2),(3,5),(1,6),(0,5)` | `(34,10,13)` | `(10,10,5,2,2,1)` |
| `C3 / Theta` | `(0,1),(2,0),(4,2),(4,3),(3,5),(1,6),(0,4)` | `(34,10,13)` | `(10,8,5,3,2,2,1)` |
| `C3 / 0` | `(0,1),(2,0),(4,2),(3,5),(1,6),(0,4)` | `(33,9,13)` | `(10,10,5,2,2,1)` |

Every packet has defect `24=2*13-2`. Thus Pick, the dual graph, carrier
normalization, and Riemann--Hurwitz give four independent views of the same
genus ledger.

## 5. Proper-flat extinction

After a finite common base change and toric regularization, every special
component on `(2)` is either rational or the smooth genus-three carrier
`E11/E10`. A rational curve maps constantly to the good elliptic target, and
the positive form order `41/26` makes the carrier map constant. Retaining
actual positive component multiplicities on one common dominating base, the
inherited proper-flat identity gives

```text
deg(phi_generic^*L)=sum_Gamma m_Gamma deg(phi_Gamma^*L)=0,           (19)
```

contradicting the positive generic response degree of a nonautomorphic Keller
pair. This proves the gate extinct relative to the inherited interfaces.

## 6. Unexpected connections and natural indexing

The middle face is a sharp warning against a simple tournament model. If its
two rational components are treated as vertices, a tournament can retain only
one directed edge, while `(15)` supplies eleven separately addressed parallel
nodes and graph genus ten. The faithful carrier is a labelled multigraph (or
an edge torsor over `mu_11`), not a tournament. The map

```text
node -> its root in {S:alpha*S^11=1}                   (20)
```

preserves incidence and transversality; quotienting by the cyclic action
destroys ten graph cycles. The required sidecar is the root address.

The 64 exact support patterns also have a closed natural address. Order the
eight coupled classes in `(7)` by `c=0,...,7`, and put

```text
(w,p,eeta)=(1_[u!=0],1_[Phi!=0],1_[eta!=0]),
n=1+8c+4w+2p+eeta.                                    (21)
```

Then `(21)` bijects exact supports with `{1,...,64}`. One may equivalently
call the support the `n`th odd-square address `(2n-1)^2` while storing `n`.
The six-fan coordinate `(L,1_[Theta!=0])` is a coarser `3 x 2` quotient: it
forgets `Phi,eta`, two support cancellations, and the individual `mu_11`
attachments. A lawful compressed record is

```text
(n; cap owner, right owner, plane word,
    edge word, mu_11 attachment labels, primitive bases).           (22)
```

This also explains what happened to the six old `M` attachments of
THM-4350. As `U` tends to zero, their roots have valuation
`v(P)=-v(U)/6` and leave through `P=infinity`; they are replaced by the
eleven nodes `(15)`, not specialized into them. What survives is total genus,
not an identity of labelled edges.

## 7. Audit and scope

The primary certificate deterministically reconstructs `(4)--(19)` in 880
checks. The independently written 870-check referee rebuilds support coupling,
the factorized middle face, all eleven nodes, rational parametrizations, edge
reductions, carrier squarefreeness, all six packet ledgers, and the
primitive/common-base distinction. Both ordinary and optimized executions
byte-match the frozen outputs recorded in the frontmatter.

The complementary `alpha=0` endpoint strata, seam entry, `JC(2)`, and `DC(2)`
remain outside this theorem.
