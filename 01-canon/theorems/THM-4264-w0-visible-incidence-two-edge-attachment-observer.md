---
id: THM-4264
title: "W=0 visible-incidence two-edge attachment observer"
status: >
  PROVED COROLLARY TO THM-4249/4253/4258 + VERIFIED-EXACT. On every one of
  the current 1,512 visible map-ratio incidences, two consecutive attachment
  differences vanish iff all twelve vanish. Thus the concrete observer table
  has 3,024 group-value rows, improving THM-4258's ungated three-edge table of
  4,536. After the proof reaches the two-torsion sidecar, the visible-projector
  relation cuts the ambient 64 recurrence words to exactly 16, with period
  profile {1:1,3:3,6:12}; one edge is module-theoretically insufficient.
  No incidence is evaluated or excluded. W=0, M=12, normal lift, seam entry,
  JC(2), and DC(2) remain OPEN.
source: root/cross-frontier-overnight/2026-08-27
depends_on:
  - THM-4249-w0-cyclic-projector-missing-eigenline-attachment-squeeze
  - THM-4253-w0-degree-forty-two-norm-three-profile-exclusion
  - THM-4258-w0-three-sample-attachment-recurrence-and-two-torsion-sidecar
  - THM-4259-w0-explicit-hlambda-normalization-and-glue-dictionary
mistake_firewall:
  - MISTAKE-521
script: 04-computation/jc23_w0_visible_two_edge_observer_thm4264.py
output: 05-knowledge/results/jc23_w0_visible_two_edge_observer_thm4264.out
script_sha256: 7aa38c69f60fd854d2c51b7b32036e13dd164544dcec68d772aa470cc9bc97d4
output_sha256: f941b7d18567a0e62660ec8b7cc3cee0649364d491422cfc75beb3c898b614a4
independent_script: 04-computation/jc23_w0_visible_two_edge_rank_audit_thm4264.py
independent_output: 05-knowledge/results/jc23_w0_visible_two_edge_rank_audit_thm4264.out
independent_script_sha256: 57bdd4c44bdf1cea8a4836d3ec735a867eb75a3fd2773a09f496e45f47271de5
independent_output_sha256: 4ae1bc8d76d17c68d26c9e14425d51dba9c15a694a47c6805e1fdd3be0c2d124
hash_basis: raw bytes
audit: >
  PASS. One implementation enumerates all 64 F4 cubic-recurrence words and
  the 16 projected words, including both sharp hostiles. A second path builds
  and row-reduces the 24-coordinate binary constraint matrices directly. It
  obtains ambient/projected nullities 6/4 and projected two-prefix nullity 0.
  Normal and optimized streams are byte-identical for both scripts.
---

# THM-4264 -- `W=0` visible-incidence two-edge attachment observer

**PROVED COROLLARY TO THM-4249/4253/4258 + VERIFIED-EXACT. NO CURRENT
INCIDENCE IS EXCLUDED BY THIS THEOREM.**

## 1. Exact incidence gate

Use THM-4249's `W=0` curve, elliptic target, and normalized Hom basis:

```text
C_0:x^6+y^4=1,              E_0:Y^2=X^3+1,
O=Z[omega],                 omega^2+omega+1=0,
Q_j=tau^j Q_0,              D_j=[Q_(j+1)-Q_j],
m=b f+c g+d h in M_0={a_u=0},
delta_j=m(D_j) in E_0.                                  (1)
```

Let `T` be precomposition by `tau`, put

```text
pi=omega^2-1,       P=v(Q_0),
V=P_v(T)m=dv,       H=P_L(T)m=2ell.                    (2)
```

Every one of the current `864+648=1,512` map-ratio incidences retains its
actual visible ideal `(d)` and satisfies

```text
[d*pi]P=O.                                             (3)
```

This is stronger than membership in a maximal torsion envelope after the
actual ideal has been forgotten. Replacing `d` or `P` by associates preserves
`(3)`. Source shifting rotates the difference word and target-unit scaling
multiplies every entry by a unit, so observer zero is invariant under the
source-target orbit quotient.

## 2. Two-edge theorem

For every current incidence in `(1)--(3)`,

```text
delta_j=O for all j mod 12
iff delta_0=delta_1=O.                                 (4)
```

Equivalently, the three point images `m(Q_0),m(Q_1),m(Q_2)` are equal iff all
twelve attachment images are equal. This is a **two-difference observer on
the visible incidence relation**, not a two-difference observer on all of
`M_0`.

### Visible projector collapse

THM-4249 gives

```text
v(Q_j)=[omega^(2j)]P.                                  (5)
```

Hence `(3)` implies, for every `j`,

```text
V(D_j)=[d*omega^(2j)*(omega^2-1)]P=O.                 (6)
```

Thus `V(Q_j)` is constant. If `A(T)=sum_k a_kT^k` and `S` shifts a difference
word, covariance gives the typed identity

```text
(A(T)m)(D_j)=sum_k[a_k]m(D_(j+k))=A(S)delta_j.         (7)
```

Since `P_v(T)m=V`, equation `(6)` supplies the extra relation

```text
P_v(S)delta=0.                                         (8)
```

### Two zero differences force two-torsion

On `M_0`, THM-4249 gives

```text
2m=V+H,                   T^2H=-omega H.               (9)
```

Write `p_j=m(Q_j)` and `H_j=H(Q_j)`. If
`delta_0=delta_1=O`, then `p_0=p_1=p_2`. By `(6)` and `(9)`,

```text
H_0=H_1=H_2.                                           (10)
```

But `H_2=[-omega]H_0`, so `[1+omega]H_0=O`. The element
`1+omega=-omega^2` is a unit, hence `H_0=H_1=O`; the order-two recurrence
kills every `H_j`. Together with visible collapse, `(9)` now yields

```text
[2]delta_j=O for every j.                              (11)
```

### Confluent recurrence shortens under the visible gate

THM-4258 proves

```text
P_0(S)delta=0,
P_0(X)=(X-omega^2)(X^2+omega).                        (12)
```

On `E_0[2]~=O/(2)=F_4`, set `N=S+omega^2`. Then

```text
P_0(S)=N^3,
P_v(S)=-omega^2(S+omega)(S^2+omega)
      =unit*(1+N)N^2.                                 (13)
```

On `N^3=0`, the factor `1+N` is invertible with inverse `1+N+N^2`.
Equations `(8)` and `(13)` imply

```text
N^2 delta=0,
delta_(j+2)=[omega]delta_j.                            (14)
```

The two initial zeros propagate through all twelve positions, proving `(4)`.

## 3. Exact sidecar and sharp gates

Once the proof has reached `(11)` and both projected maps collapse, the
conditional sidecar universe is

```text
{delta in E_0[2]^12:P_0(S)delta=0,
                       P_v(S)delta=0,
                       sum_j delta_j=O}.               (15)
```

Equation `(14)` gives

```text
delta_(2k)=omega^k delta_0,
delta_(2k+1)=omega^k delta_1.                          (16)
```

Thus `(15)` has exactly `4^2=16` words, and their period profile is

```text
{period 1:1, period 3:3, period 6:12}.                 (17)
```

All telescope. One edge is insufficient even in this projected module: the
nonzero word

```text
(0,1,0,omega,0,omega^2,0,1,0,omega,0,omega^2)         (18)
```

has `delta_0=0` and satisfies both polynomial relations. Geometric
realization of `(18)` is not asserted.

Conversely, the ambient `64=4^3` module in THM-4258 cannot be shortened
without `(8)`. For nonzero `R in E_0[2]`,

```text
delta_j=(omega^2)^j binom(j,2)R                        (19)
```

satisfies `P_0(S)delta=0` and has `delta_0=delta_1=0`, `delta_2!=0`, but its
visible-projector word is nonzero. For `R=1` these two words are

```text
delta=(0,0,omega,1,0,0,1,omega^2,0,0,omega^2,omega),
P_v delta=(omega,1,omega^2,omega,1,omega^2,
           omega,1,omega^2,omega,1,omega^2).           (20)
```

This identifies the exact load-bearing gate rather than weakening the
ambient theorem.

## 4. Exact workload

THM-4258's intrinsic difference map is

```text
F_m(Q)=m([tau Q-Q]),          delta_j=(T^jF_m)(Q_0).   (21)
```

The visible-gated table therefore needs

```text
degree 34: 864*2=1,728,
degree 42: 648*2=1,296,
total:             3,024.                              (22)
```

This saves 1,512 rows over THM-4258's three-difference table and 15,120 over
the naive twelve-difference table. If curve-map values are evaluated instead,
the same check uses the three point images `p_0,p_1,p_2`; those are not three
difference rows. THM-4258's observer-degree ceilings `120/150` remain valid.

THM-4259 supplies the explicit `H_lambda` normalization and glue dictionary,
so `(22)` is concretely evaluable without a two-torsion halving ambiguity.
No row has been evaluated here; consequently no incidence is removed.

## 5. Independent audits and scope

The recurrence companion enumerates all 64 ambient and 16 projected words.
The independent path constructs the `24`-coordinate binary constraint
matrices directly and finds

```text
ambient rank/nullity             18/6,
projected rank/nullity           20/4,
ambient + two zero prefixes      22/2,
projected + one zero prefix      22/2,
projected + two zero prefixes    24/0.                 (23)
```

Reproduction:

```bash
python3 -B 04-computation/jc23_w0_visible_two_edge_observer_thm4264.py
python3 -B -O 04-computation/jc23_w0_visible_two_edge_observer_thm4264.py
python3 -B 04-computation/jc23_w0_visible_two_edge_rank_audit_thm4264.py
python3 -B -O 04-computation/jc23_w0_visible_two_edge_rank_audit_thm4264.py
```

The observer is tangential after specialization to `W=0`; it carries no
`W`-normal information. A relative Hom section and enough transverse Hasse
jets are still needed for a normal lift. No current predicate-preserving map
transfers `(21)` to seam entry. Thus the theorem proves neither entry nor
`W=0` emptiness, exact `M=12`, `JC(2)`, or `DC(2)`. **QED.**
