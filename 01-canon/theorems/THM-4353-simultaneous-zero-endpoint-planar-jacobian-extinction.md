---
id: THM-4353
title: "Simultaneous zero-endpoint planar Jacobian extinction"
status: >
  PROVED RELATIVE TO THM-4230/4327/4344 + VERIFIED-EXACT + INDEPENDENTLY
  HOSTILE-AUDITED. The entire inherited exact-weight-twelve gate
  Z=beta_11=zeta_3=W=xi_10=K=0, U!=0 is extinct. Its 48 exact supports have
  eight fans. The only positive-genus faces are automatically smooth
  genus-three carriers of positive form order; all Theta=0 faces are
  rational. The sole repeated quadratic edge is one smooth boundary-tangent
  branch, whose two index-two punctures merge to index three and whose two
  ordinary blowups add only rational components. Seam entry, JC(2), and
  DC(2) are not claimed.
source: root + jc-deep-k0 / planar-Jacobian next-sharp session, 2026-09-02
depends_on:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4327-generic-exact-weight-twelve-endpoint-wall-extinction
  - THM-4344-clean-cubic-infinity-exit-planar-jacobian-extinction
related:
  - THM-4161-y-only-double-top-root-planar-jacobian-exclusion
  - THM-4342-clean-cubic-zero-exit-planar-jacobian-extinction
  - THM-4350-first-normal-owner-cubic-infinity-exit-planar-jacobian-extinction
  - THM-4351-double-normal-owner-zero-cubic-infinity-exit-planar-jacobian-extinction
mistake_firewall:
  - MISTAKE-522
  - MISTAKE-531
  - MISTAKE-540
primary_script: 04-computation/jc2_m12_simultaneous_zero_endpoint_extinction_thm4353.py
primary_output: 05-knowledge/results/jc2_m12_simultaneous_zero_endpoint_extinction_thm4353.out
primary_script_sha256: 35b61f1765f7c6b4f143941e77faa205f728f0ae9fb79c94996a5dde9cf5350f
primary_output_sha256: 1837e48346606cd179fe51b73b37314fe26e1be8366597857b69c188eb7c16d4
referee_script: 04-computation/jc2_m12_simultaneous_zero_endpoint_extinction_independent_referee_thm4353.py
referee_output: 05-knowledge/results/jc2_m12_simultaneous_zero_endpoint_extinction_independent_referee_thm4353.out
referee_script_sha256: 6ea7db7d4902edb90100ec189bfa7fda6dbb7655c4ff5c054619cd58f723fd23
referee_output_sha256: cfd0e2ab8106896dac836bd43661410ecbf467d9f4b3fa7bf2fe85216c18ddd0
hash_basis: raw LF bytes
audit: >
  PASS AFTER PACKET REPAIR. The 431-check primary rebuilds the literal source,
  48 exact supports, 128-key/96-support hostile atlas, all seven primitive
  faces, carrier discriminants, orders, edges, Pick/graph ledgers, the exact
  reciprocal differential, and both charts of both boundary blowups. The
  independently written 637-check referee reconstructs the same data and
  catches the required specialization from packet (11,6,4,2,2,1) to
  (11,6,4,3,1). Normal and optimized streams byte-match both frozen outputs.
---

# THM-4353 -- Simultaneous zero-endpoint planar Jacobian extinction

**PROVED RELATIVE TO THM-4230/4327/4344 + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED. THE ENTIRE DISPLAYED GATE IS EXTINCT. `JC(2)`, `DC(2)`, AND
SEAM ENTRY REMAIN OPEN.**

## 1. Statement and inheritance

Work over `C` in the inherited `b=d=0` reduced `(2,3)`, exact-weight-twelve
seam. Retain the complete sixteen-row source

```text
H=-3p+(8/3)p^2-(1376/135)p^3+K y^2+Phi p^2y
  +Delta p^4+Theta py^2+eta p^3y+zeta_3 y^3
  +u p^5+xi_10 p^2y^2+alpha p^4y+beta_11 py^3
  +U p^6+W p^3y^2+Z y^4,

K=2848/45-(7/6)Delta,                 u=upsilon_5,
alpha=alpha_11.                                           (1)
```

Impose

```text
Z=beta_11=zeta_3=W=xi_10=K=0,          U!=0.              (2)
```

Thus `Delta=5696/105`. Then no nonautomorphic
planar Keller pair lies on `(1)--(2)`, with `Theta,alpha,eta,Phi,u`
arbitrary. The conclusion is relative to the inherited toroidal, good-target,
and proper-flat interfaces. It proves neither entry into this seam nor
`JC(2)` or `DC(2)`.

The inheritance pass is:

- closest proved mechanism: THM-4350's owner merger and positive-order
  replacement carriers;
- canonical hostile: THM-4351's repeated edge is a two-branch `A5` contact
  and creates a tail, whereas the repeated edge below is one smooth branch
  tangent to the toric boundary and creates no tail;
- corrected near miss: a repeated outer root changes the puncture packet even
  when it changes neither geometric genus nor graph genus;
- least-used sidecar: THM-4161's rule that tangency multiplicity adds
  ramification defects rather than puncture indices.

The Anchor / Niche / Wildcard board is

```text
K-zero fan | first normal owner | support quotient | rational pruning
boundary tangency | puncture defect | owner tree | natural address.          (3)
```

## 2. Exact supports and the eight-fan owner tree

For

```text
F_Q=(s^2-p)(1-QH)-Q*s^2/2,                              (4)
```

the gate `(2)` has 22 active-or-optional lifted support points. The only
coefficient coupling relevant to the support classification is

```text
(2,4,1)=Theta-5696/105.                                 (5)
```

Consequently `Theta` has three support states:

```text
Theta=0;        Theta=Delta;        Theta notin {0,Delta}. (6)
```

Combining `(6)` with the four independent presence bits of
`alpha,eta,Phi,u` gives exactly

```text
3*2^4=48                                                  (7)
```

distinct realizable supports. Literal lower-hull enumeration gives exactly
eight fans:

| conditions | exact supports | fan |
|---|---:|---|
| `Theta!=0, alpha!=0` | 16 | `M,D6,E11` |
| `Theta!=0, alpha=0` | 16 | `M,E01` |
| `Theta=0, alpha!=0, (Phi,eta)!=(0,0)` | 6 | `M,D6,N` |
| `Theta=0, alpha!=0, Phi=eta=0` | 2 | `M,D6` |
| `Theta=alpha=0, Phi*eta!=0` | 2 | `M,Eeta,N` |
| `Theta=alpha=Phi=0, eta!=0` | 2 | `M,Eeta` |
| `Theta=alpha=eta=0, Phi!=0` | 2 | `M,EPhi` |
| `Theta=alpha=eta=Phi=0` | 2 | `M` |

The factor two in each `Theta=0` terminal is the support bit of `u`; it does
not create another face.

A hostile Boolean atlas independently toggles the five optional source rows
and the two multiply visible lifted points. Its `2^7=128` keyed configurations
collapse to 96 distinct supports. Every exact support embeds literally, and
the hostile atlas has the same eight fans. This is a stability probe, not a
claim that the extra 48 distinct supports are coefficient-realizable.

## 3. Every face, normalization, and form order

For a plane `(a,b,c)`, lifted height is `a*r+b*l+c`. Up to an invertible
torus monomial, the exhaustive face table is:

| face | plane | polygon | initial equation | genus | target base/order |
|---|---|---|---|---:|---:|
| `M` | `(1/12,1/6,-1/6)` | `(0,1),(2,0),(2,6),(0,7)` | `(P-S^2)(UP^6-1)` | seven rational components | `12/9` |
| `D6` | `(1/6,1/6,-1/3)` | `(2,0),(3,5),(2,6)` | `S^2(1-UP^6-alpha SP^5)` | 0 | `6/5` |
| `E11` | `(2/7,1/7,-4/7)` | `(2,0),(4,3),(3,5)` | `S^2(1-alpha SP^5-Theta S^2P^3)` | 3 | `42/41` |
| `E01` | `(1/4,1/6,-1/2)` | `(2,0),(4,3),(2,6)` | `S^2(1-UP^6-Theta S^2P^3)` | 3 | `12/11` |
| `Eeta` | `(1/3,1/6,-2/3)` | `(2,0),(3,4),(2,6)` | `S^2(1-UP^6-eta SP^4)` | 0 | `6/6` |
| `EPhi` | `(1/2,1/6,-1)` | `(2,0),(3,3),(2,6)` | `S^2(1-UP^6-Phi SP^3)` | 0 | `6/7` |
| `N` | `(1,0,-2)` | `(2,0),(3,3),(3,5)`; endpoint deletions give `(2,0),(3,4),(3,5)` or `(2,0),(3,3),(3,4)` | `S^2(1-SP^3(Phi+eta P+alpha P^2))` | 0 | `6/11` |

Here `target base/order` uses the least common base needed to compare the
source face with the target differential. The source chart of `N` itself has
base one. The inherited density formula is

```text
5/6-(a+b+c),                                             (8)
```

and gives the positive orders `9,5,41,11,6,7,11` in table order.

The faces `D6,Eeta,EPhi,N` are rational because their displayed equations
solve linearly for `S`. The only positive-genus normalizations are

```text
E11: y=2Theta*S*P^2+alpha*P^4,
     y^2=P(4Theta+alpha^2P^7),
     disc=-53971714048*Theta^8*alpha^12;                 (9)

E01: y=Theta*S*P^2,
     y^2=Theta*P(1-UP^6),
     disc=46656*Theta^12*U^5.                           (10)
```

On their exact owner gates these degree-eight or degree-seven branch
polynomials are squarefree, so both carriers are smooth of genus three. There
is no carrier discriminant left to resolve anywhere on `Theta!=0`.

The literal primitive source charts are

```text
M:    Q=sigma^12, s=sigma^-1 S, p=sigma^-2 P, G=sigma^2 F_Q;
D6:   Q=sigma^6,  s=sigma^-1 S, p=sigma^-1 P, G=sigma^2 F_Q;
E11:  Q=sigma^7,  s=sigma^-2 S, p=sigma^-1 P, G=sigma^4 F_Q;
E01:  Q=sigma^12, s=sigma^-3 S, p=sigma^-2 P, G=sigma^6 F_Q;
Eeta: Q=sigma^6,  s=sigma^-2 S, p=sigma^-1 P, G=sigma^4 F_Q;
EPhi: Q=sigma^6,  s=sigma^-3 S, p=sigma^-1 P, G=sigma^6 F_Q;
N:    Q=sigma,    s=sigma^-1 S, p=P,          G=sigma^2 F_Q. (11)
```

Direct substitution of the complete source `(1)` gives exactly the seven
initial equations above; no genus or owner assertion is inferred from the
polygon alone.

## 4. Edge schemes, graph genus, and puncture packets

All fixed outer and internal edge schemes are reduced in characteristic zero
under their displayed owner conditions and `U!=0`. The moving right-boundary
words and internal words are:

| owner case | moving right-boundary word | internal word |
|---|---|---|
| `Theta*alpha!=0` | `1-Theta X; Theta+alpha X; alpha+UX` | `1-UX^6; 1-alpha X` |
| `Theta!=0, alpha=0` | `1-Theta X; Theta+UX` | `1-UX^6` |
| `Theta=0, alpha*Phi!=0` | `1-Phi X; Phi+eta X+alpha X^2; alpha+UX` | `1-UX^6; 1-alpha X` |
| `Theta=Phi=0, alpha*eta!=0` | `1-eta X; eta+alpha X; alpha+UX` | `1-UX^6; 1-alpha X` |
| `Theta=Phi=eta=0, alpha!=0` | `1-alpha X; alpha+UX` | `1-UX^6` |
| `Theta=alpha=0, Phi*eta!=0` | `1-Phi X; Phi+eta X; eta+UX` | `1-UX^6; 1-eta X` |
| only `eta!=0` | `1-eta X; eta+UX` | `1-UX^6` |
| only `Phi!=0` | `1-Phi X; Phi+UX` | `1-UX^6` |
| no normal owner | `1-UX^6` | none |

The only possibly nonreduced scheme is

```text
J_edge(X)=Phi+eta X+alpha X^2,
disc(J_edge)=eta^2-4alpha Phi.                           (12)
```

Section 5 treats its double-root wall. Away from that wall, the global
polygon, source-infinity packet, and contracted normalized graph are:

| conditions | global `(2A,B,I)` | packet | `(V,E,b1)` | carrier genus |
|---|---:|---|---:|---:|
| `Theta*alpha!=0` | `(37,11,14)` | `(11,8,6,5,1)` | `(9,19,11)` | 3 |
| `Theta!=0, alpha=0` | `(36,10,14)` | `(13,11,5,1)` | `(8,18,11)` | 3 |
| `Theta=0, alpha*Phi!=0` | `(32,12,11)` | `(11,6,4,2,2,1)` | `(9,19,11)` | 0 |
| `Theta=Phi=0, alpha*eta!=0` | `(31,11,11)` | `(11,6,5,2,1)` | `(9,19,11)` | 0 |
| `Theta=Phi=eta=0, alpha!=0` | `(30,10,11)` | `(11,6,6,1)` | `(8,18,11)` | 0 |
| `Theta=alpha=0, Phi*eta!=0` | `(31,11,11)` | `(11,7,4,2,1)` | `(9,19,11)` | 0 |
| only `eta!=0` | `(30,10,11)` | `(11,7,5,1)` | `(8,18,11)` | 0 |
| only `Phi!=0` | `(30,10,11)` | `(11,8,4,1)` | `(8,18,11)` | 0 |
| no normal owner | `(24,14,6)` | `(11,1,1,1,1,1,1,1)` | `(7,12,6)` | 0 |

The `M` face consists of seven rational components with twelve internal
nodes, hence graph genus six. Every first right face meets `M` in the six
simple roots of `1-UX^6`; adding it changes `(V,E)` by `(1,6)` and raises
`b1` to eleven. When `N` is present, it is a rational leaf joined by one
reduced internal edge, changing `(V,E)` by `(1,1)`. This proves the graph
ledger without using Pick. Each displayed packet independently satisfies
`sum(e_i-1)=2I-2`.

In particular, setting `K=0` in the nonzero-owner fans of THM-4350 prunes a
rational terminal leaf:

```text
M-D6-E11-T -> M-D6-E11,
M-E01-T    -> M-E01.                                   (13)
```

Both `V` and `E` fall by one, so `b1=11` is conserved. This is the exact
zero/infinity endpoint relation; it is not an import of THM-4342's
`W!=0` elliptic zero-exit face.

## 5. Repeated outer edge: one smooth tangent branch

Assume

```text
Theta=0,       alpha*Phi!=0,       eta^2=4alpha Phi.
```

Write

```text
J(P)=Phi+eta P+alpha P^2=alpha(P-a)^2,
a=-eta/(2alpha)!=0.                                    (14)
```

The exact primitive `N` chart from `(11)` is

```text
G_N=(S^2-sigma^2P)(1-SP^3J(P)-sigma D(P))-sigma S^2/2,

D(P)=-3P+(8/3)P^2-(1376/135)P^3
     +(5696/105)P^4+uP^5+UP^6.                         (15)
```

At the boundary `S=infinity`, put `R=1/S` and multiply by `R^3`. Then

```text
H=(R-sigma^2PR^3)(1-sigma D(P))
  -(1-sigma^2PR^2)P^3J(P)-sigma R/2.                   (16)
```

The exact identities

```text
H(0,P,sigma)=-alpha P^3(P-a)^2,
H_R(0,a,sigma)=1-sigma D(a)-sigma/2                    (17)
```

show that `H_R` is a DVR unit. Formal implicit preparation therefore gives

```text
R=(P-a)^2 V(P-a,sigma),              V(0,0)!=0.         (18)
```

Thus the curve is smooth. There are not two branches, no vertical `A_k`
singularity, and no ramified collision tail.

For completeness, absorb `V` into the boundary coordinate and write
`q=P-a`, `r=R/V`, so the reduced pair is `r=q^2` against `r=0`. The first
ordinary blowup has chart

```text
r=q r_1:       C': r_1=q,       B':r_1=0,       E_1:q=0. (19)
```

The three reduced components meet at one point. Blowing up that point gives

```text
r_1=q r_2:     C'':r_2=1,       B'':r_2=0,      E_2:q=0,
q=r_1 q_2:     C'':q_2=1,       E_1':q_2=0,     E_2:r_1=0. (20)
```

The reduced total transform is now simple normal crossings. Both exceptional
components are rational; no base change was used. The first blowup alone
makes the two strict transforms transverse, but leaves a triple point, which
is why `(20)` records the second blowup.

The repeated root does change the puncture packet. On `H=0`, the source form
in this chart is, up to the fixed toric unit,

```text
dS/G_P=-R*dR/H_P=R*dP/H_R.                             (21)
```

At a simple root of `J`, `ord(R)=1` and the puncture index is two. On `(14)`,
`ord_q(R)=2`, so the two former index-two punctures merge to one index-three
puncture:

```text
(11,6,4,2,2,1)  ->  (11,6,4,3,1).                     (22)
```

Both sides have defect 20, equal to `2*11-2`. The genus and graph ledgers
therefore remain unchanged. Freezing the generic packet on this wall would
violate the support/multiplicity firewall.

## 6. Proper-flat extinction

After a finite common base change and toric regularization, every special
component on `(2)` is either

1. a rational component from `M,D6,Eeta,EPhi,N`, a toric chain, or the
   boundary resolution `(19)--(20)`; or
2. the smooth genus-three carrier `E11` or `E01`, on which the pulled-back
   Keller form has positive order `41` or `11`.

Every map from a rational projective curve to the good elliptic target is
constant. Positive form order makes each genus-three carrier map constant as
well. Retaining the actual positive component multiplicities on one common
dominating base, the inherited proper-flat identity gives

```text
deg(phi_generic^*L)=sum_Gamma m_Gamma deg(phi_Gamma^*L)=0            (23)
```

for a positive-degree target line bundle `L`. This contradicts the positive
generic response degree of a nonautomorphic Keller pair. Therefore the gate
`(1)--(2)` is extinct. **QED.**

## 7. Unexpected connections and lawful natural indexing

The support patterns themselves admit a closed natural address. Put

```text
c=0  if Theta=0,
c=1  if Theta=Delta,
c=2  if Theta notin {0,Delta},

(a,e,p,w)=(1_[alpha!=0],1_[eta!=0],1_[Phi!=0],1_[u!=0]),
n=1+16c+8a+4e+2p+w.                                   (24)
```

Then `(24)` is a bijection from the 48 exact support patterns to
`{1,...,48}`. Equivalently, the pattern may be called the `n`th odd-square
address `(2n-1)^2`, while storing only `n`. This indexes support, not the
continuum of coefficient values and not the collision divisor `(14)`.

The eight fans form a rooted first-owner decision tree, not a tournament:
some pairs of states have no intrinsic directed comparison. The lawful
compressed record is

```text
(n; owner-tree leaf, plane word, edge word,
    disc(J_edge) status, puncture packet, primitive bases).          (25)
```

The comparison with THM-4351 and THM-4352 isolates the decisive extra bit.
A two-branch collision can carry delta and create a two-ended tail/graph
cycle; the one-branch boundary tangency `(18)` carries neither, but merges
two puncture defects. Forgetting branch count makes these configurations look
combinatorially alike and gives the wrong geometry. The cheapest decisive
test is the transverse derivative `H_R`: a unit selects the tangent branch;
its vanishing would reopen the collision-tail hierarchy.

## 8. Audit and scope

The primary exact certificate checks `(4)--(23)` in 431 assertions. An
independently written reconstruction repeats the source, supports, fans,
primitive charts, edges, carriers, orders, graphs, reciprocal residue, packet
change, and both blowups in 637 assertions. Both normal/optimized pairs
byte-match the frozen outputs named in the frontmatter.

What is proved is only the gate `(1)--(2)`, relative to the inherited
toroidal, good-target, and proper-flat interfaces. Seam entry, other cells,
`JC(2)`, and `DC(2)` remain open.
