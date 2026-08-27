---
id: THM-4275
title: "Opposite-parity attachment observers and confluent sample matroids"
status: >
  PROVED RELATIVE TO THM-4258/4264/4272 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. On every one of
  THM-4264's inherited 1,512 visible map-ratio incidences, any two attachment
  samples of opposite parity have exactly the same fibres as the full
  twelve-edge word; all 36 unordered pairs work, including gaps three and
  nine. The ambient/visible abstract recurrence sidecars have exact sample
  matroids governed respectively by residues modulo four and parity. A
  one-edge visible quotient has perfectly uniform fibres but cannot recognize
  collapse. At THM-4272's rank-twelve fat contact, generic Fourier point
  values lose rank 12-to-1 on collision while the full normal-Hasse observer
  retains rank twelve. Simple special-fibre roots impose no uniform
  transverse jet depth. No new incidence is excluded and JC(2)/DC(2) remain
  open.
source: codex-creative-frontiers-20260827
depends_on:
  - THM-4258-w0-three-sample-attachment-recurrence-and-two-torsion-sidecar
  - THM-4264-w0-visible-incidence-two-edge-attachment-observer
  - THM-4272-lambda-zero-a23-contact-and-e0-infinity-jet-obstruction
related:
  - THM-4263-moving-multigraph-filtered-jet-and-finite-factor-density-transport
  - THM-4274-confined-confluent-observer-fibre-and-density-transport
script: 04-computation/jc23_confluent_sample_matroids_thm4275.py
output: 05-knowledge/results/jc23_confluent_sample_matroids_thm4275.out
script_sha256: 28df3d1549e3112b44eda7fcb4f377db1bfedebfaf6e3d4a47cba7dc060a01df
output_sha256: cc76f44caaffaece2ea85f59f0d077f9789f9e39835e500f6cdaae3e9a98759f
independent_script: 04-computation/jc23_confluent_sample_matroids_rank_audit_thm4275.py
independent_output: 05-knowledge/results/jc23_confluent_sample_matroids_rank_audit_thm4275.out
independent_script_sha256: d24525bae78053809cafefeb3014bff2f267c152e8020f73eb8784b19c686ef6
independent_output_sha256: cc9d6fe1acbbe01ebc5879f7e5c973c49cd2778afc8cf401edb90e9aabe02cfa
hash_basis: raw bytes
audit: >
  PASS. A dependency-free standard-library path enumerates all 64 ambient
  and 16 visible F4 recurrence words, all 4,096 coordinate subsets in each
  module, every pair/triple rank and regular fibre, all 36 faithful
  opposite-parity pairs and 30 sharp same-parity failures. It also checks
  1,728 Fourier inversion cells across the twelve nonzero F13 Kummer
  parameters, the coalesced-value/Hasse hostile, and 512 transverse
  truncation fibre pairs. An independent path constructs only the coordinate
  functionals and row-reduces them over F4/F13, without enumerating words or
  fibres, and reproduces every rank census. Normal, optimized, and
  fixed-hash-seed streams are byte-identical for both paths.
---

# THM-4275 -- opposite-parity attachment observers and confluent sample matroids

**PROVED RELATIVE TO THM-4258/4264/4272 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED. NO NEW INCIDENCE IS EXCLUDED; `JC(2)` AND `DC(2)` REMAIN OPEN.**

## 1. Opposite parity preserves the full attachment-word fibre

Fix any one of THM-4264's `1,512` inherited visible map-ratio incidences. Keep
its notation

```text
E_0:Y^2=X^3+1,                  O=Z[omega],
Q_j=tau^jQ_0,                  D_j=Q_(j+1)-Q_j,
delta_j(m)=m(D_j),             j in Z/12.                         (1)
```

Let `C_I` be any candidate subset in that fixed incidence universe: its
elements lie in the same linear module `M_0` and satisfy the inherited visible
relation with the fixed incidence data. Define

```text
Phi(m)=(delta_0(m),...,delta_11(m)),
Obs_(r,s)(m)=(delta_r(m),delta_s(m)).                              (2)
```

Then, for every `r,s in Z/12` of opposite parity,

```text
Phi(m)=Phi(m') iff Obs_(r,s)(m)=Obs_(r,s)(m')
for all m,m' in C_I.                                               (3)
```

Thus every opposite-parity pair is a complete observer of the **full fibre
partition**, not merely the zero/collapse predicate. There are

```text
6*6=36                                                               (4)
```

unordered such pairs. By minimum cyclic gap they split as

```text
gap 1:12,                  gap 3:12,                  gap 5:12.     (5)
```

In particular, gaps `3` and `9` work even though they do not generate
`Z/12`; parity orbits, not additive generation by the gap, control the proof.

### Proof

Let `x=m-m'` and write `epsilon_j=x(D_j)`. Every input identity used below is
linear, so it survives passage to this difference. THM-4264 supplies the
visible collapse and hidden recurrence

```text
V_x(D_j)=0,
2x=V_x+H_x,                  T^2H_x=-omega H_x.                    (6)
```

Put `eta_j=H_x(D_j)`. Equations `(6)` give

```text
eta_j=2epsilon_j,            eta_(j+2)=-omega eta_j.                (7)
```

If the two samples in `(3)` agree, then `epsilon_r=epsilon_s=0` for one even
and one odd index. Since `-omega` is a unit, `(7)` propagates the first zero
through the even parity orbit and the second through the odd orbit. Hence

```text
eta_j=0 and 2epsilon_j=0 for every j.                              (8)
```

The difference word now lies in `E_0[2]^12`. The two inherited operator
relations of THM-4264 are

```text
P_0(S)epsilon=0,              P_v(S)epsilon=0.                     (9)
```

With `N=S+omega^2`, THM-4264 factors them as

```text
P_0(S)=N^3,
P_v(S)=unit*(1+N)N^2.                                             (10)
```

On `ker N^3`, the operator `1+N` is invertible, so `(9)--(10)` imply

```text
N^2epsilon=0,
epsilon_(j+2)=omega epsilon_j.                                    (11)
```

The zero at `r` kills its parity orbit under `(11)`, and the zero at `s`
kills the other. Therefore `epsilon=0`, proving the forward implication in
`(3)`; the reverse implication is immediate. `square`

This proof does not compare candidates belonging to different incidence
data, and it does not assert that the abstract hostile words below are
geometrically realized.

## 2. Exact coordinate-sample matroids

The following classifications concern the conditional abstract recurrence
modules already proved in THM-4258/4264.

### 2.1 Ambient cubic module

THM-4258 parametrizes its 64 words over `F_4` as

```text
delta_j=alpha^j(c_0+binom(j,1)c_1+binom(j,2)c_2),
(c_0,c_1,c_2) in F_4^3,              alpha=omega^2.                 (12)
```

The unit factor `alpha^j` does not change coordinate-observer rank. Modulo
two, the four distinct evaluation rows, indexed by `j mod 4`, are

```text
0:(1,0,0),       1:(1,1,0),       2:(1,0,1),       3:(1,1,1).     (13)
```

Any three are independent and all four sum to zero. Therefore, for a sample
set `J subset Z/12`,

```text
rank_F4 Obs_J=min(|{j mod 4:j in J}|,3).                            (14)
```

This is `U_(3,4)` with each element replaced by three parallel positions.
The exact subset census is

```text
all 4,096 subsets: rank 0:1, rank 1:28, rank 2:294, rank 3:3,773;
66 pairs:          rank 1:12, rank 2:54;
220 triples:       rank 1:4,  rank 2:108, rank 3:108.               (15)
```

Thus a three-coordinate ambient observer is complete exactly when it meets at
least three residue classes modulo four; adjacency is irrelevant.

### 2.2 Visible-gated quadratic module

THM-4264 parametrizes its 16 words as

```text
delta_(2k)=omega^k d_0,
delta_(2k+1)=omega^k d_1,             (d_0,d_1) in F_4^2.           (16)
```

Therefore

```text
rank_F4 Obs_J=|{j mod 2:j in J}|.                                  (17)
```

The sample matroid has two rank-one elements, each with six parallel copies.
Its exact census is

```text
all 4,096 subsets: rank 0:1, rank 1:126, rank 2:3,969;
66 pairs:          rank 1:30, rank 2:36;
220 triples:       rank 1:40, rank 2:180.                           (18)
```

The 36 rank-two pairs in `(18)` are exactly those proved geometrically
complete on the inherited incidences in Section 1. Same-parity insufficiency
is sharp on the full abstract sidecar. Representative nonzero words missed by
all even, respectively all odd, samples are

```text
(0,1,0,omega,0,omega^2,0,1,0,omega,0,omega^2),
(1,0,omega,0,omega^2,0,1,0,omega,0,omega^2,0).                     (19)
```

### 2.3 Perfect density can still destroy collapse

Give the 16 visible words uniform measure and retain one coordinate. The
linear map

```text
e_j:K_vis -> F_4,                 delta |-> delta_j                 (20)
```

is surjective with four equal fibres. Its normalized fibre weight in
THM-4274's notation is identically one, so it transports every uniform target
event with exact measure. Nevertheless, `e_j^(-1)(0)` contains the zero word
and three nonzero words. Full collapse is not saturated under this quotient.

In the exact observer ledger,

```text
full certificate:       the twelve-edge word in K_vis;
complete observer:      one even and one odd coordinate in F_4^2;
quotient:               forget either coordinate;
preserved:              uniform pushforward and one parity seed;
destroyed:              the other parity seed and collapse-fibre saturation;
minimal restoring data: one sample from the other parity.                         (21)
```

This is a concrete witness that measure preservation and predicate
preservation have different boundaries.

## 3. A rank-twelve collision turns Fourier values into Hasse channels

Let `K` be a field whose characteristic does not divide `12` and which
contains a primitive twelfth root of unity `zeta`. Put

```text
A_Lambda=K[b]/(b^12-Lambda).                                      (22)
```

After the Kummer base change `Lambda=a^12`, take

```text
f(b)=sum_(k=0)^11 c_k b^k,
y_j=f(zeta^j a),                       0<=j<12.                     (23)
```

Fourier orthogonality gives, for every `k`,

```text
(1/12)sum_(j=0)^11 zeta^(-jk)y_j=c_k a^k.                          (24)
```

Hence the twelve labelled point values have rank twelve when `a!=0`. At
`a=0`, all sections coalesce and the value map has rank one: it retains only
`c_0`. Equivalently, the discriminant of `b^12-Lambda` is

```text
unit * 12^12 * Lambda^11.                                         (25)
```

On the fat fibre

```text
A_0=K[b]/(b^12),                                                   (26)
```

the normal-Hasse observer is

```text
f |-> (D_b^[k]f(0))_(0<=k<12)=(c_0,...,c_11),                     (27)
```

and remains rank twelve. The minimal hostile is `f=b`: all twelve coalesced
point values vanish, while `D_b^[1]f(0)=1`.

THM-4272 identifies `(26)` with the `Lambda=0` rank-twelve contact on the
`W=0` face. Thus the lawful information transition there is

```text
Lambda!=0: twelve labelled reduced evaluations after Kummer splitting;
Lambda=0:  one fat point with twelve coefficient/Hasse channels.           (28)
```

Selecting a few coalesced reduced values is not a wall observer. One must
retain the full finite-flat restriction, or equivalently all normal-Hasse
channels after choosing the target local coordinate. This statement supplies
no missing descent from a resolved Keller map to the raw fat contact.

## 4. Finite transverse depth and its sharp moving-shell failure

Let `C` be a nonempty finite set and let `h:C->K[[W]]`. If not all values are
equal, define

```text
R_min=1+max_(h(c)!=h(c')) ord_W(h(c)-h(c'));                        (29)
```

if all are equal, put `R_min=0`. Then, for every integer `R>=0`,

```text
c |-> h(c) mod W^R has exactly the same fibres as h
iff R>=R_min.                                                       (30)
```

Indeed, unequal series collide modulo `W^R` exactly when their difference
has order at least `R`. Formula `(30)` is THM-4274's pair-difference test in
the `W`-adic filtration and is effective once the finite shell is explicit.

No bound follows from simple roots on the special fibre. For every `N>=2`
and any `t_0`, put

```text
D_N(W,t)=t-t_0,
E_N(W,t)=t-t_0-W^N.                                                (31)
```

Both special-fibre roots are simple, but their root sections differ by
exact order `N`; every transverse jet through order `N-1` collides. Thus
THM-4265's simple reduced roots do not bound an off-fibre observer depth.
A coherent off-fibre dictionary and an effective pair-contact bound remain
open prerequisites.

## 5. Reproduction and scope

Run

```text
python -B 04-computation/jc23_confluent_sample_matroids_thm4275.py
python -B -O 04-computation/jc23_confluent_sample_matroids_thm4275.py
python -B 04-computation/jc23_confluent_sample_matroids_rank_audit_thm4275.py
python -B -O 04-computation/jc23_confluent_sample_matroids_rank_audit_thm4275.py
```

The program enumerates the abstract sidecars and finite algebra controls. The
proof in Section 1, not enumeration of the 16 sidecar words, transports the
opposite-parity result to the inherited incidence hypotheses.

This theorem does **not** assert realization of every abstract word, evaluate
or exclude a new map-ratio incidence, build an off-`W=0` Hom dictionary,
descend a resolved Keller map to the raw `A_23` contact, cross `U=0` or `Z=0`,
prove exact-`M=12` entry, or establish `JC(2)`/`DC(2)`.
