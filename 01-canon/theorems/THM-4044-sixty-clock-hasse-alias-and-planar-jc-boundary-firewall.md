---
id: THM-4044
title: "Sixty-clock Hasse alias and planar-JC boundary firewall"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Over every field containing
  60 distinct sixtieth roots of unity, the depth-k Hasse observer on those
  roots has kernel ((P^60-1)^k). On the pure-P planar-JC residual ideal P^2K[P]
  its kernel is P^2(P^60-1)^kK[P]; hence it is injective through degree 60k+1
  and first fails, in one dimension, at degree 60k+2. The sharp alias
  Delta_k=P^2(P^60-1)^k preserves the full torus observer, the endpoint
  polynomial, and its first P-Hasse jet, but changes [P^2] by (-1)^k. Thus no
  fixed finite 60-clock depth can certify THM-3997's mandatory [p^2]R!=0 for
  an unbounded residual; the exact missing sidecar is the second boundary
  Hasse jet or a degree cap. The same theorem makes THM-4035's broad 4D
  Kakeya spine a cubic interpolation frame, not an unrestricted residual
  evaluator. Conversely, one clock determines Sun's degree-at-most-eight
  centered atoms, and three clock jets determine THM-4034's degree-178
  conductor after scalar extension. No residual alias is asserted to satisfy
  the Keller equations, and JC(2) remains open.
source: jc2-merged-frontiers-0824 / exact clock-kernel audit, 2026-08-24
audit: >
  PASS. The primary certificate realizes the theorem in F_61 with phi=44,
  checks Delta_k through k=6, exhausts all 487635 broad four-phase
  Vandermonde minors, verifies the half-clock node-address involution, and
  checks the bounded Sun and conductor controls. An independent confluent
  evaluation-matrix path obtains ranks (60,60,61), (120,120,121), and
  (180,180,181) before failure, at first failure, and after adjoining the
  second boundary jet. Normal and optimized streams match both frozen
  outputs.
depends_on:
  - THM-3997-reduced-two-three-hasse-repair-and-zero-residual-no-go
  - THM-4035-sixty-clock-separation-and-finite-kakeya-spine
related:
  - THM-4010-confluent-consecutive-hasse-observer-kernel-index-and-smith-firewall
  - THM-4034-exceptional-quartic-global-conductor-degree-178
  - THM-4037-centered-binomial-parity-and-singular-fibres
  - THM-4038-ap-deficit-holonomic-sixty-phase-law
script: 04-computation/jc2_sixty_clock_hasse_alias_thm4044.py
output: 05-knowledge/results/jc2_sixty_clock_hasse_alias_thm4044.out
script_sha256: cc49cd7024fdeaab6c0d668c9ca497ee113ea6c45b3c9735d836322781629898
output_sha256: 4e73f5d3fed3bec3966d914f4981777f1a3b78892321bd443ed83f2799185fc0
independent_audit_script: 04-computation/jc2_sixty_clock_hasse_alias_thm4044_independent_audit.py
independent_audit_output: 05-knowledge/results/jc2_sixty_clock_hasse_alias_thm4044_independent_audit.out
independent_audit_script_sha256: e5c97fdcff0a2f0824dd4c70582228bc7acfd5e55dccf8850a2dc21bddf169bf
independent_audit_output_sha256: fbb7ab4429fcf6ade9d912f482f3b99aee71de01539713fb204c39952625f056
independent_audit_semantic_sha256: 00f526270644924d75acf32de7f638f9b9dffa04ddda83d5320c76b7a93b4e94
hash_basis: raw LF bytes
---

# THM-4044 -- the exact `C_60` confluent kernel

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.** A sixty-phase clock is a
lossless address and a sharp bounded-degree interpolator. It is not a
complete observer of an unbounded planar-Jacobian residual. The first lost
coordinate is exactly the coefficient which THM-3997 forces to be nonzero.

## 1. Complete confluent kernel

Let `K` be a field whose characteristic does not divide `60` and which
contains the group

```text
mu_60={t in K:t^60=1},                    |mu_60|=60.    (1)
```

For `f in K[P]`, define its Hasse derivatives by

```text
f(P+T)=sum_(j>=0) D_P^[j]f(P) T^j.                    (2)
```

For `k>=1`, the depth-`k` clock observer is

```text
O_k(f)=(D_P^[j]f(t))_(t in mu_60, 0<=j<k).             (3)
```

At a fixed `t`, all entries in `(3)` vanish precisely when `(P-t)^k` divides
`f`. The factors are pairwise coprime, so

```text
ker O_k = product_(t in mu_60)(P-t)^k K[P]
        = (P^60-1)^k K[P].                              (4)
```

This proves the all-`k` statement; it is not a finite matrix extrapolation.
Because `P` is coprime to `P^60-1`, intersecting `(4)` with the pure residual
ideal gives

```text
ker(O_k|P^2K[P])=P^2(P^60-1)^k K[P].                  (5)
```

Consequently:

```text
O_k is injective on K[P]_(degree<60k),
O_k is injective on P^2K[P]_(degree<=60k+1),            (6)
```

and its first kernel on the latter space is one-dimensional:

```text
Delta_k=P^2(P^60-1)^k,                deg Delta_k=60k+2.
                                                               (7)
```

The boundary expansion is sharp:

```text
D_P^[0]Delta_k(0)=D_P^[1]Delta_k(0)=0,
D_P^[2]Delta_k(0)=[P^2]Delta_k=(-1)^k !=0.             (8)
```

Thus the second boundary Hasse jet, rather than another nonzero phase, is
the first sidecar which separates the minimal alias.

## 2. The planar-JC residual firewall

Use `P,Y` here to distinguish the abstract observer variables from the
source variables in THM-3997. Extend the characteristic-zero base field if
necessary so that it contains `mu_60`; scalar extension does not change a
coefficient's vanishing. Apply `(3)` coefficientwise in `Y` and strengthen it
with both available lower boundary layers:

```text
E_k(R)=
  ((D_P^[j]R)(t,Y))_(t in mu_60,j<k),
  R(0,Y),
  (D_P^[1]R)(0,Y).                                     (9)
```

For every residual `R in (P^2,Y)` and every `c in K`, equations `(7)--(8)`
give

```text
R_c=R+c Delta_k in (P^2,Y),
E_k(R_c)=E_k(R),
[P^2]R_c=[P^2]R+c(-1)^k.                              (10)
```

The final scalar in `(10)` can be prescribed arbitrarily, in particular
made zero. Therefore the observer `(9)` cannot certify

```text
[P^2]R !=0.                                             (11)
```

This is exactly THM-3997's first mandatory pure residual coordinate, which
equals `(D_P^[2]R)(0,0)`. A lawful complete observer must retain that boundary
jet, impose `deg_P R<=60k+1`, or use an independently target-preserving
sidecar.

Equation `(10)` is an **observer counterexample**, not a Keller pair. Adding
`Delta_k` need not preserve the source-normal equations, the all-row Hasse
image condition, weighted face, or generic-fibre map. Hence the theorem
refutes a proposed finite-clock inference, not the forced coefficient or
`JC(2)`.

## 3. What the four-dimensional Kakeya spine really preserves

In the exact `F_61` model of THM-4035, `phi=44` has order `60`, so

```text
t_r=phi^r,                         {t_r}=F_61^*.         (12)
```

The broad one-clock rows are

```text
B(r)=[1:t_r:t_r^2:t_r^3].                              (13)
```

They are precisely evaluation rows on `F_61[P]_(<=3)`. Four distinct rows
have the Vandermonde determinant

```text
product_(i<j)(t_j-t_i) !=0,                            (14)
```

so all `binom(60,4)=487635` phase quartets reconstruct a cubic. More
generally, `O_k` reconstructs the residue class modulo `(P^60-1)^k`, and
becomes an actual polynomial decoder only under `(6)`.

For four selected phases, the hostile

```text
P^2 product_(i=1)^4(P-t_i)^k                           (15)
```

has all selected Hasse jets below `k` equal to zero while its `P^2`
coefficient is nonzero. The full-clock version of `(15)` is `(7)`. Thus
four-wise transversality is a bounded interpolation statement; it neither
reconstructs an unrestricted JC residual nor supplies a Keller determinant.

There is a precise boundary-chart match. The torus `{t_r}` omits `P=0`, and
the lost coefficient in `(8)` is a second normal coordinate at that omitted
boundary. THM-4035 likewise omits projective directions with a zero
coordinate. This is a common quotient loss, not a transfer from Euclidean
Kakeya geometry to `JC(2)`.

## 4. The half-clock is a node address, not an owner

Since `phi^30=-1`, the half-period acts by

```text
t_(r+30)=-t_r.                                         (16)
```

For each `r`, put in `F_61`

```text
a_r=-(2/3)t_r^2.                                       (17)
```

The THM-3992 deleted-line normalization

```text
nu_a(X)=(X^2+a, X^3+(3a/2)X)                          (18)
```

sends `X=t_r` and `X=t_(r+30)` to the same node

```text
(-a_r/2,0).                                             (19)
```

Thus the half-clock labels the two normalization addresses of all `30`
split nonzero node parameters over `F_61`. It does not force the equation
`X^2=-3a/2`, label nonsplit nodes, choose an address owner, or prove the
global balance/nonproperness alternatives of THM-3996. Antipodality by
itself is not a collision theorem.

## 5. Why Sun atoms and the degree-178 conductor are positive controls

For even `m`, the centered binomial polynomial

```text
B_m(T)=binom((T+m-1)/2,m)                              (20)
```

is even: `B_m(-T)=B_m(T)`. The clock shift `(16)` therefore realizes the
reflection `n -> m-1-n` for `m=2,4,6,8`. Each polynomial has degree at most
eight, strictly below `60`, so `O_1` determines it exactly. The missing phase
`T=0` is the unique fixed point of the reflection; torus data do not by
themselves restore the parity/fixed-point information isolated in THM-4037.

The Fibonacci/triangular scalar quotient gives a prior address loss. By
THM-4035,

```text
(F_r mod 10,T_r mod 30)
```

has twelve doubleton fibres, while adjoining `r mod 2` recovers the full
phase. Even after that phase repair, `(7)` gives the independent evaluator
loss. Hence there are two separate gates:

```text
scalar state + parity  -> phase address,
phase address + degree/boundary sidecar -> polynomial evaluator.           (21)
```

There is also a higher-degree positive control. THM-4034's conductor
polynomial has degree `178`, and

```text
178 < 3*60.                                              (22)
```

After extending its coefficient field to contain `mu_60`, `O_3` determines
that conductor uniquely by confluent interpolation. This does not make
internal multiplication recover the three cotangent lines lost in
THM-3737: interpolation and the multiplier algebra are different observers.

## 6. Connection ledger and strict stopping boundary

The exact connection contract is

```text
source:       the pointed clock C_60;
target:       mu_60 subset A^1 and K[P]/(P^60-1)^k;
map:          r -> zeta^r followed by confluent evaluation;
preserved:    residue class, and the polynomial itself under degree <60k;
destroyed:    multiples of (P^60-1)^k and the P=0 normal coefficient;
sidecar:      a degree cap or the second boundary Hasse jet;
hostile:      Delta_k=P^2(P^60-1)^k.                   (23)
```

The consecutive-node factor in THM-4010/4038 and the stable coefficient
label `J_6` in the Russell lane must also remain separated. One is an
interpolation polynomial; the other is an order index in a pulled-back
Jacobian series. A common numeral `6` supplies no map between them.

This theorem proves no upper degree bound on `R`, no new all-row Hasse
identity, no weighted-face exclusion, no actual Russell lift, no Kakeya
dimension bound, and no representation or minimality theorem for Sun's
2-4-6-8 problem. `JC(2)` remains **OPEN**.

## 7. Reproduction

```bash
python3 -B 04-computation/jc2_sixty_clock_hasse_alias_thm4044.py
python3 -B -O 04-computation/jc2_sixty_clock_hasse_alias_thm4044.py
python3 -B 04-computation/jc2_sixty_clock_hasse_alias_thm4044_independent_audit.py
python3 -B -O 04-computation/jc2_sixty_clock_hasse_alias_thm4044_independent_audit.py
```

The normal and optimized streams match the frozen transcripts. The primary
semantic hash is
`0c35001112cfda56ec6a6495f7907d84bf84077a1948128a412f8dfbfadd2855`;
the independent hash is recorded in the frontmatter. **QED.**
