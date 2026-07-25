---
id: THM-2305
title: "Canonical blocker-word handoff hypergraph at prescribed expiration"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. On the
  prescribed-return arm of THM-2296, the terminal blocker word is
  canonically one of {a}, {b}, and {a,b}. Either a pure exclusive-owner
  transfer or a genuine two-target fork receives at least one third of the
  return mass. On that same exact word stratum, every nonzero C_13
  character has energy strictly larger than one ninth of the squared
  return floor, and covariant endpoint-Prony supplies a signed
  current-service coefficient with gauge index at most 4S-1. More
  generally, fork mass d and the larger pure transfer obey the sharp trade
  max(p_a,p_b)>=(rho-d)/2. If all three owner labels independently enter
  pure return branches, their transfer graph contains a directed 2-cycle or
  3-cycle; this is not yet an orbit cycle because edge composition loses
  the incoming/outgoing subset intersection. No scalar row is excluded and
  LRC(14) remains open.
source: codex-2026-07-25-canonical-handoff-hypergraph
depends_on:
  - THM-2255-valuation-separated-pair-cap-and-exclusive-owner-mass
  - THM-2296-prescribed-expiration-return-or-bounded-ancestry-resonance
  - THM-2299-rooted-current-service-energy-and-base-phase-no-go
related:
  - THM-2269-marked-expiration-root-spectrum-and-branch-state-no-go
  - THM-2294-anchored-plucker-tournament-and-kakeya-address-bank
script: 04-computation/lrc14_canonical_handoff_hypergraph_thm2305.py
output: 05-knowledge/results/lrc14_canonical_handoff_hypergraph_thm2305.out
script_sha256: 43e45d9048bced40e12ef99033e5420bb007bdd98234b1e92b9813070710e536
output_sha256: de4924c8889c04b4d50ab501f1b3bd8118e5ff7495e903dfde7c7176214fc771
hash_basis: working-tree bytes (LF)
---

# THM-2305 -- the prescribed handoff is a blocker-word hyperedge

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2299 assigns overlap points to one of two target blockers in a fixed
order. That is sufficient for a named service label, but the order is a
gauge: a point dangerous for both remaining blockers has not made an
exclusive owner transfer.

The intrinsic terminal observable is instead the complete nonempty blocker
word among the other two labels. This produces an exact alternative:

```text
one owner expires
  -> a pure transfer to one other exclusive owner,
     or a fork simultaneously dangerous for both other blockers.       (1)
```

The distinction is load-bearing. Pure transfers are directed graph edges.
Forks are directed hyperedges and must not be cosmetically oriented.

## 1. Exact blocker-word partition

Use the scalar-cover notation of THM-2296:

```text
A_0=C_H minus union_(i=1)^5 D_(q_i),

E_j
 =A_0 intersection D_(c_j)
       minus union_(ell!=j)D_(c_ell),

R_j=A_0 minus D_(c_j),

c_j=13^(lambda_j)u_j,
k_j=lambda_j+1.                                      (2)
```

Fix the selected owner `j`, and call the other two labels `a,b`. Define
the three exact terminal blocker-word strata

```text
Q_(j,{a})
 =A_0 intersection D_(c_a)
       minus (D_(c_j) union D_(c_b))
 =E_a,

Q_(j,{b})
 =A_0 intersection D_(c_b)
       minus (D_(c_j) union D_(c_a))
 =E_b,

Q_(j,{a,b})
 =A_0 intersection D_(c_a) intersection D_(c_b)
       minus D_(c_j).                                (3)
```

The scalar cover says that, away from its null exceptional set,

```text
R_j
 =Q_(j,{a}) disjoint_union Q_(j,{b})
                     disjoint_union Q_(j,{a,b}).     (4)
```

This partition has no ordering choice. The first two strata are literally
the already defined exclusive-owner sets of the other labels; the third is
the simultaneous two-blocker stratum.

Put

```text
G_j=P^(lambda_j)1_(E_j),

rho_(j,sigma)
 =integral_(Q_(j,sigma)) P G_j,                     (5)

sigma in {{a},{b},{a,b}}.
```

The transfer identity and (4) give

```text
rho_j
 =measure(E_j intersection T^(-k_j)R_j)
 =sum_sigma rho_(j,sigma).                           (6)
```

All three terms are nonnegative.

## 2. The sharp pure-edge/fork trade

Write

```text
p_a=rho_(j,{a}),
p_b=rho_(j,{b}),
d=rho_(j,{a,b}).                                    (7)
```

Then

```text
max(p_a,p_b)>=(rho_j-d)/2.                           (8)
```

This elementary inequality is sharp for every fixed pair `(rho_j,d)`: take
the two pure masses equal. Consequently,

```text
d>=rho_j/3,

or

max(p_a,p_b)>rho_j/3.                                (9)
```

Thus at least one canonical word stratum has mass at least `rho_j/3`.
More precisely, equality in the fork branch of (9) is allowed, while
failure of that branch makes the pure edge strict.

On THM-2296's quantitative return arm this gives the exact mass floors

```text
strict:
 max-word mass
 >=39002430583/160481782761300;

repeated-first:
 max-word mass
 >=13560199813/160481782761300.                     (10)
```

For a pure stratum, the transfer identity has the literal form

```text
measure(E_j intersection T^(-k_j)E_t)
 =rho_(j,{t})>0.                                    (11)
```

It is therefore an intrinsic time-directed exclusive-owner edge

```text
j -> t.                                             (12)
```

For the double stratum it is instead the hyperedge

```text
j -> {a,b}.                                         (13)
```

No choice of one head in (13) preserves the exact target word.

## 3. Every root character survives on the same word

Retain THM-2299's thirteen predecessor values

```text
v_r(y)=G_j((y+r)/13),

M_k(y)=sum_(r=0)^12 v_r(y)zeta^(-kr),
                                      1<=k<=12.     (14)
```

Because `support(G_j) subset D_(u_j)`, each root fibre has at most two
positive entries. THM-2299's sharp odd-root chord bound gives

```text
|M_k(y)|>P G_j(y)                                   (15)
```

whenever `P G_j(y)>0`, simultaneously for all twelve nonzero characters.

Let `sigma` be any terminal word. Cauchy--Schwarz and
`measure(Q_(j,sigma))<=1` give

```text
integral_(Q_(j,sigma)) |M_k|^2
 >rho_(j,sigma)^2.                                  (16)
```

Select the word supplied by (9). On the quantitative return arm,

```text
strict:
 integral_(Q_(j,sigma)) |M_k|^2
 >1521189591381733719889
   /25754402598245085852777690000;

repeated-first:
 integral_(Q_(j,sigma)) |M_k|^2
 >183879018968485234969
   /25754402598245085852777690000.                  (17)
```

These inequalities hold for every `k=1,...,12` on one and the same exact
blocker word.

## 4. A bounded signed coefficient on the exact word

Gauge-correct as in THM-2299:

```text
N_k(y)=exp(-2*pi*i*k*y/13)M_k(y),

W_(k,sigma)=1_(Q_(j,sigma))N_k.                     (18)
```

The function `N_k` is periodic. Equation (16) makes `W_(k,sigma)` nonzero.
Both `G_j` and every indicator in (3) use only the nine original scalar
boundary banks. Hence

```text
#Jump(G_j)<=2S,
#Jump(1_(Q_(j,sigma)))<=2S.                         (19)
```

The step amplitude in (18) has at most `4S` jumps. Covariant
endpoint-Prony therefore gives

```text
0<=h<=4S-1
```

with

```text
(W_(k,sigma))_hat(h)
 =integral_(Q_(j,sigma)) M_k(y)
    exp(-2*pi*i(h+k/13)y)dy
 !=0.                                               (20)
```

Thus the pure-edge/fork alternative preserves not only energy but a bounded
signed coefficient with the exact terminal word still present.

## 5. The honest graph consequence

There is a useful conditional three-owner statement. Suppose each of the
three exclusive owners independently has a positive prescribed-return arm
and, for each owner, the selected large word in (9) is pure. Selecting one
large pure target defines

```text
f:{1,2,3}->{1,2,3},
f(j)!=j.                                            (21)
```

Every finite functional digraph contains a directed cycle. Since (21) has
three vertices and no loops, the cycle has length two or three. Therefore
the pure branch contains

```text
j <-> t,

or

1 -> 2 -> 3 -> 1                                   (22)
```

after relabelling.

This is a legitimate directed relation because its pairwise observable is
positive-measure exclusive-owner transfer at the source's actual expiration
clock. It is not generally a tournament: both directions may occur, a
vertex may have two positive targets, and a fork must remain a hyperedge.

Most importantly, (22) is not yet an orbit cycle. An incoming transfer
lands on a positive subset of `E_t`, while an outgoing transfer starts from
another positive subset of `E_t`; the two subsets can be disjoint. The
missing sidecar is their actual intersection, or a mixing/continuation
mechanism that forces one. Bare blocker labels discard it.

## 6. Connection and stopping ledger

```text
source:
  THM-2296's quantitative prescribed-time return and THM-2299's rooted
  current-service word;

map:
  stratify the current residual by its complete terminal blocker word
  before selecting a large target;

preserved:
  source owner, actual clock, exact target word, pure-versus-double
  multiplicity, transfer mass, every nonzero root character, and a bounded
  signed gauge coefficient;

new object:
  a directed graph edge j->t for a pure exclusive-owner handoff, or the
  irreducible directed hyperedge j->{a,b} for simultaneous service;

destroyed by the bare hypergraph:
  terminal component, complex base phase, root sheet, and the intersection
  needed to compose consecutive positive-measure transfers;

next decisive split:
  charge the double-target fork by a two-blocker intersection mechanism,
  or force incoming/outgoing overlap on a pure owner cycle.             (23)
```

The theorem is conditional only on the already stated scalar-cover and
return-arm hypotheses. The three-owner cycle paragraph is explicitly
conditional on all three positive pure returns; the current LRC reduction
guarantees only one quantitatively selected owner. No profile is eliminated.

## 7. Exact verification

The companion uses exact integer and `Fraction` arithmetic. It checks the
three blocker-word truth tables for every selected source label, the sharp
one-third trade, both branch floors and energy floors, the `4S-1` jump
ledger, and all loop-free functional graphs on three labels. Reproduce with

```bash
python3 04-computation/lrc14_canonical_handoff_hypergraph_thm2305.py
python3 -O 04-computation/lrc14_canonical_handoff_hypergraph_thm2305.py
```

Every load-bearing check raises explicitly, so optimized mode runs the same
audit. QED.
