---
id: THM-4414
title: "LRC14 six-separated contact capacity collapse"
status: >
  PROVED ANALYTICALLY RELATIVE TO THM-4409 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. A sharp sparse-interval lemma makes every degree-zero THM-4409
  max-flow equal its edgewise-minimum envelope at arbitrary height. The
  envelope has an exact raw-carrier projection formula and a crossing-hinge
  loss. THM-4434 now proves its universal 6/77 ceiling. Entry,
  synchronization, and LRC(14) remain open.
source: root + cross_frontier_bridge + network_universal / LRC14 continuation session, 2026-09-05
depends_on:
  - THM-4409-lrc14-third-sheet-component-network-certificate
related:
  - THM-2216-residual-capacity-hinge-gram-law
  - THM-2658-balanced-lift-helly-circular-arc-gain-nerve-and-wrap-boundary
  - THM-2789-interval-gram-tomography-and-graceful-gap-tail-quadratic-detector
  - THM-4392-lrc14-raw-carrier-box-spline-fourier-poisson-duality
  - THM-4434-lrc14-universal-scale-three-network-projection-bound
primary_script: 04-computation/lrc14_six_separated_contact_capacity_collapse_thm4414.py
primary_output: 05-knowledge/results/lrc14_six_separated_contact_capacity_collapse_thm4414.out
primary_script_sha256: 591c2a2ff540d5eb95e4baf7fb3c02cd5828be6ba8e617c7f92182bb9f82cedd
primary_output_sha256: cc70d23dac4ddea4e05f31426b84819a817552bf1e172d16f4501e15c53fb65e
independent_audit_script: 04-computation/lrc14_six_separated_contact_capacity_collapse_thm4414_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_six_separated_contact_capacity_collapse_thm4414_independent_audit.out
independent_audit_script_sha256: f29e52f223c9c422563e2c7dec1083c5ac2b78c7dcf3332968a66151a238dc64
independent_audit_output_sha256: 0f5b3fb45e0cbf4e0e48773e1b368dfcdd85c69d3331476d1d4c4dabd8a9847a
hash_basis: raw LF bytes
audit: >
  PASS. The primary verifies six-separation, edgewise feasibility, hinge and
  quantitative-Helly identities on every height-79 sheet graph. A separate
  raw-y implementation retains the two owner-orientation classes, derives the
  carrier envelope over one common integer denominator, and independently
  reproduces the height-79 counts while extending the finite signal to height
  149. Normal and optimized outputs byte-match; no assertion is optimized out.
---

# THM-4414 -- LRC14 six-separated contact capacity collapse

**PROVED ANALYTICALLY RELATIVE TO THM-4409 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED.** This removes a max-flow optimization from the degree-zero local
three-speed quotient. THM-4434 subsequently proves the arithmetic ceiling.
This theorem does not prove arbitrary chart entry, synchronization, or
`LRC(14)`, which remains **OPEN**.

## 1. Sharp sparse-interval lemma

Let `A` and `B` be finite families of pairwise disjoint bounded open intervals
on a line. Suppose that there are `s_A,s_B>0` such that

```text
|I|<=s_A for I in A,     consecutive A-gaps >=2s_A,
|J|<=s_B for J in B,     consecutive B-gaps >=2s_B.       (1)
```

Join `I` to `J` when their intersection has positive length, and give every
vertex its interval length as capacity. Then

```text
maxflow(A,B)=sum_(I~J) min(|I|,|J|).                      (2)
```

Indeed, the right side is an edgewise upper bound. Assign each edge the load

```text
x_(I,J)=min(|I|,|J|).                                    (3)
```

If `I` meets `d>=2` intervals of `B`, it spans the `d-1` intervening gaps, so

```text
|I|>2(d-1)s_B>=d s_B
   >=sum_(J~I)|J|>=sum_(J~I)x_(I,J).                     (4)
```

For `d<=1` the capacity inequality is immediate. Interchanging the families
checks every `B` vertex, so `(3)` is feasible and proves `(2)`.

The factor two is sharp. For any `rho<2`, two unit intervals separated by a
gap `rho` admit one interval meeting both but having length below two, so its
two edge-min loads exceed its capacity. The exact audit uses `rho=3/2`, edge
sum `2`, and vertex capacity `17/10`. At the threshold `rho=2`, a positive-
contact control has edge sum `2` and spanning capacity `11/5`.

The contact graph in `(2)` is a forest: the union of the two internally
disjoint interval families is an interval graph and bipartite, hence has no
cycle. Forest structure alone does not prove `(2)`; the factor-two spacing is
what eliminates star competition.

## 2. Intersection preserves sparsity

Call a disjoint family `rho`-sparse at scale `s` when its components have
length at most `s` and consecutive gaps at least `rho*s`. If `A` and `B` are
`rho`-sparse at scales `s_A,s_B`, the nonempty components of `A intersect B`
are `rho`-sparse at scale `min(s_A,s_B)`.

Each new component lies in one parent from each family, giving the length
bound. Between consecutive new components at least one parent changes, and
the new gap contains that parent's gap. Thus it is at least
`rho*min(s_A,s_B)`.

For speed `w`, a fixed sheet of the THM-4409 danger set has pieces of length

```text
s_w=1/(7w)                                                (5)
```

and gaps `6s_w`. Cutting a wrapping interval at zero only shortens the two
end pieces and creates no new internal adjacency. Therefore, on every fixed
sheet permutation `pi`, both

```text
P_pi=D_(w_i,pi(i)) intersect D_(w_j,pi(j)),
K_pi=D_(w_k,pi(k))                                       (6)
```

are six-sparse at their respective scales. Applying `(2)` gives, at every
height,

```text
kappa_pi=sum_(I~J) min(|I|,|J|),
U_net(i,j;0,0)=sum_pi kappa_pi.                          (7)
```

The sheet qualifier is essential. In raw `y=3x` coordinates, the two cyclic
orientations of the three owner labels must remain separate endpoint classes.
Collapsing them can manufacture false vertex competition. Once retained,
they introduce no multiplicity in the final raw-carrier sum below.

## 3. Exact raw-carrier projection formula

Let `w1<w2<w3` and use the THM-4413 margins

```text
p_i(C)=3(w_j+w_k)-14|C_i|.                               (8)
```

Write

```text
Lambda_w^live={C in Z^3:
  C dot w=0,
  every C_i is nonzero modulo 3,
  every p_i(C)>0}.                                      (9)
```

For the coordinate pair `{j,k}` with exact third coordinate `i`, define

```text
E_i(w)=sum_(C in Lambda_w^live)
  min(3/(7w3), p_i(C)/(14w_jw_k)).                       (10)
```

Then the complete six-sheet degree-zero capacity is exactly

```text
U_net(j,k;0,0)=E_i(w).                                   (11)
```

To see this, a pair component has length

```text
min(3/(7w_j),3/(7w_k),p_i/(14w_jw_k)),                  (12)
```

while the opposite component has length `3/(7w_i)`. Their edgewise minimum
is `(10)` because `3/(7w3)` is the shortest of all three caps. Each raw
`y`-component has three inverse `x` branches, each contributing one third of
this minimum; their sum is one raw term. The two owner orientations partition
the raw carriers and cause no factor of two or three.

The physical mass is the full three-facet minimum

```text
mu(F_w)=sum_(C in Lambda_w^live) min_(i=1,2,3)
  min(3/(7w3),p_i(C)/(14w_jw_k)).                        (13)
```

Thus THM-4409's remaining universal degree-zero target is no longer a flow
problem. It is the explicit projection inequality

```text
min_(i=1,2,3) E_i(w)<=6/77.                              (14)
```

Equation `(14)` is now **PROVED** by [THM-4434](THM-4434-lrc14-universal-scale-three-network-projection-bound.md) in the odd ternary-unit domain, with sole primitive equality `(1,5,11)`.

## 4. Exact crossing-hinge loss

For one contact edge put

```text
p=|I|, q=|J|, z=|I intersect J|.
```

Then

```text
min(p,q)-z=min(|I minus J|,|J minus I|).                 (15)
```

Consequently

```text
network envelope - physical mass
 =sum_over_contacts min(directed exclusive left tail,
                        directed exclusive right tail). (16)
```

The graph quotient has no Hall-competition loss at degree zero. Its entire
slack is the crossing hinge `(16)`, which vanishes exactly when each contact
is nested. For `(1,5,11)`, every pair choice has six nested matching edges and
both sides equal `6/77`. For the first strict THM-4409 graph hostile,

```text
w=(1,19,79), pair=(1,79),
envelope=8/553=152/10507,
physical=108/10507,
hinge=44/10507.                                          (17)
```

There is also a cheapest exact repair. If the pair component is
`I=A_j intersect A_k` and it contacts `J=A_i`, quantitative interval Helly
gives

```text
|A_1 intersect A_2 intersect A_3|
 =min(|A_1 intersect A_2|,
      |A_1 intersect A_3|,
      |A_2 intersect A_3|).                              (18)
```

Hence the two cross-pair overlap lengths, or just their minimum, restore the
exact edge length without full endpoint geometry. In raw-carrier language,
these are precisely the other two facets omitted from `(10)`.

## 5. Exact audits and finite frontier

The primary exhausts all `2,910` primitive distinct positive odd ternary-unit
triples through height 79, all 8,730 pair rows and all six sheet assignments.
It checks 137,682 contacts and 8,529,915 live predicates. It reproduces

```text
pair successes:       2818, 2855, 2859
selected pair counts:  400,  533, 1977
selected exact rows:  1747
physical below/on/above 6/77: 2909, 1, 0.                (19)
```

The independent raw-carrier implementation extends the finite signal to
height 149:

```text
eligible speeds: 50; primitive triples: 19,429;
best-envelope failures: 0;
unique equality: (1,5,11);
strongest row beyond height 79: (1,67,133), 60/931.      (20)
```

It also rejects the tempting unweighted average claim
`E_1+E_2+E_3<=18/77`; the first exact hostile is `(1,5,7)`, with sum
`66/245>18/77`. Selection among the three projections is essential.

Reproduce with

```powershell
python -B 04-computation/lrc14_six_separated_contact_capacity_collapse_thm4414.py 79
python -B -O 04-computation/lrc14_six_separated_contact_capacity_collapse_thm4414.py 79
python -B 04-computation/lrc14_six_separated_contact_capacity_collapse_thm4414_independent_audit.py --height 149 --audit-height 79
python -B -O 04-computation/lrc14_six_separated_contact_capacity_collapse_thm4414_independent_audit.py --height 149 --audit-height 79
```

The sparse lemma is all-height analytic; `(20)` is only **FINITE-EXACT**.
For nonconstant weights, length separation need not control mass competition,
so THM-4409's finite-degree weighted max-flows do not collapse without an
additional oscillation bound. THM-4434 now closes `(14)` by a zonotope slice
integral and exact head. The sharp remaining problems are the quantitative
body-safe-set floor, arbitrary entry, synchronization, and `LRC(14)`.
