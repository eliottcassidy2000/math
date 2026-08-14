---
id: THM-3360
title: "Uniform reflected high-pair floor by admissible affine tails"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT.  Every reflected high pair in the
  4,044 upper-median physical contexts has overlap at least 1/294.  The last
  22,890 disconnected-low affine rays are discharged by an exact rational
  limit atlas, a uniform O(1/n) averaging lemma, a p>=273 many-turn exit, and
  239,233,858 literal residual evaluations with zero failures.  This sharpens
  the weighted horn-tree proof of THM-3355: its four centered-grid horns are
  defects of that lower bound, not exceptions to the physical floor.  As a
  redundant corollary it also closes disconnected-low reflected assignments.
  Arbitrary-residue k<=1, projected k=2,3, lower-body physical entry, the
  six-body/seven-tail rung, and LRC(14) remain open.
source: codex-major-frontiers-2026-08-14
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-3350-connected-low-full-tree-atlas-dense-closure-and-uniform-tail
  - THM-3352-connected-low-all-head-universal-physical-forest-closure
  - THM-3355-disconnected-low-affine-tail-and-reflected-branch-closure
related:
  - THM-3200-fixed-lrc-channel-cleared-overlap-quasipolynomial-and-mass-recurrence-boundary
  - THM-3211-uniform-lrc-channel-limit-bernoulli-cubic-and-sharp-floor
tail_wrapper: 04-computation/lrc14_disconnected_low_affine_tail_kps_s171.py
tail_compiler: 04-computation/lrc14_disconnected_low_affine_tail_kps_s171.cpp
tail_output: 05-knowledge/results/lrc14_disconnected_low_affine_tail_kps_s171.out
tail_wrapper_sha256: 6c3d0d704b79813b25dd4da832bd047dd0c62304797533ffca5a919d5625285d
tail_compiler_sha256: 0804e5067a0bbaa3c96f6f3d69915b7baf55a4a8a6b4a35792ec0b4b9ed75e64
tail_output_sha256: f711b71da9f410847fb6d8fa65893f1e4502fe17152d52e80ab7432483db2f5d
hash_basis: LF-normalized bytes
---

# THM-3360 -- uniform reflected high-pair floor by admissible affine tails

**PROVED + FINITE-EXACT + VERIFIED-EXACT.**

## 1. Statement

Use the reflected physical-overlap semantics of THM-3352.  Fix any of its
`4,044` feasible upper-median contexts `(L,j,e,f)`.  For positive unequal raw
levels `p<q`, write

```text
(p,q)=g(P,Q),                 gcd(P,Q)=1.                (1)
```

If the primitive pair is high,

```text
P+Q>=8,                                                     (2)
```

then its exact physical overlap satisfies

```text
I_(L,j,e,f)(p,q) >= 1/294.                                 (3)
```

Consequently, let six distinct positive reflected levels be assigned to the
six labels of any upper-median body.  Join two levels when their primitive
pair has `P+Q<=7`.  If this low graph is disconnected, the assignment has a
five-edge high spanning tree and closes by Hunter's inequality.  If it is
connected, THM-3352 closes it.  Hence every reflected assignment with six
distinct levels closes.

This gives an independent closure of the full reflected branch.  THM-2941 already closes all
assignments on `2,442` bodies.  Its same-level good-pair graph is `K_6` on
`3,001` of the `3,003` bodies; the only two exceptions belong to the earlier
robust-complete block and are also closed there by their located low-channel
sidecar.  Therefore the remaining `561` bodies have a complete same-level
good-pair graph.  Any repeated level supplies a one-edge Bonferroni closure,
while the all-distinct case is the preceding paragraph.  Thus the entire
matched-residue reflected six-drift branch closes, redundantly with the
weighted horn-tree proof of THM-3355.

## 2. Reduction to the affine bank

The already verified reductions partition all pairs in `(2)` as follows.

1. A first-clock primitive-oscillation estimate gives `(3)` whenever
   `q>=8p`.
2. For `q<8p`, all `1,514` contexts with `L>=4592` satisfy `(3)` at every
   dilation.
3. On the remaining `2,530` small-ruler contexts, every non-`3:5` channel
   with `g>=4` satisfies `(3)`, and the `3:5` channel satisfies it for every
   `g>=1`.
4. The exact head scanner covers every remaining `g<=3` pair with `p<264`.
   Its exhaustive global minimum is

   ```text
   158/46397 = 1/294 + 55/13640718.                       (4)
   ```

For `p>=264`, the generalized Dirichlet reduction supplies

```text
1<=d<=8,  0<=a<=7d,  1<=|c|<=46,
p=p0+dn,  q=q0+(d+a)n,  9|c|<=p,
dq-(d+a)p=c,                                               (5)
```

where `1<=p0<=d`.  The last inequality is the inherited Dirichlet
admissibility condition; a formal point before it is not part of the selected
witness cover.  The zero-resonance case cannot occur here when `g<=3`.
There are exactly `22,890` canonical rows in `(5)`, containing `17,206`
distinct `(d,a,c)` triples.  This is a cover rather than a partition, which
is harmless because every row is proved safe.  An independent exact quotient
audit maps it onto `14,168` primitive carrier rays; the proof below does not
need to identify duplicate rows.

## 3. Exact homogenized limit

Put `D=d+a`, `k=gcd(d,D)`, `P=d/k`, `Q=D/k`, and

```text
alpha=p0-e/L,          beta=q0-f/L,
u=(ej mod L)/L,        v=(fj mod L)/L,
A=Lc+De-df.                                                 (6)
```

Body safety makes the physical overlap exactly

```text
I_n = integral_0^1 chi((knP+alpha)x-u)
                    chi((knQ+beta)x-v) dx,                 (7)
```

where `chi` is the circular radius-`1/14` indicator.  Define the periodized
tent and primitive-pair fibre by

```text
T_R(theta)=sum_m [R-|theta+m|]_+,
Phi_(P,Q)(theta)=
 [T_((P+Q)/14)(theta)-T_(|P-Q|/14)(theta)]/(PQ).            (8)
```

Then the exact ray limit is

```text
J = integral_0^1 Phi_(P,Q)(theta0-A t/(kL)) dt,
theta0=(-D(ej mod L)+d(fj mod L))/(kL).                     (9)
```

Every breakpoint in `(9)` is rational, so integrating the piecewise-quadratic
tent primitive is exact rational arithmetic.  Over the `79` context lanes not
already removed by the midpoint envelope, the compiler evaluates

```text
17,206*79 = 1,359,274                                      (10)
```

limits.  Their minimum is

```text
J_min=709/48048,
J_min-1/294=1273/112112>0,                                 (11)
```

at `(d,a,c;L,j,e,f)=(3,8,-1;168,90,12,1)`.

## 4. Effective averaging lemma

Set

```text
N=kn,   epsilon=alpha/P,   B=P beta-Q alpha=A/(kL),
delta=B/[P(N+epsilon)].                                     (12)
```

The substitution `y=(N+epsilon)x` rewrites `(7)` as

```text
I_n = 1/(N+epsilon) integral_0^(N+epsilon)
      chi(Py-u) chi(Qy+delta y-v) dy.                       (13)
```

For

```text
h(s)=integral_0^1 chi(Py-u)chi(Qy+s-v)dy                   (14)
```

translation of one circular interval gives

```text
|h(s)-h(t)|<=2|s-t|.                                       (15)
```

There is one additional variable-frequency issue inside a unit block.
Interpolate the second clock by

```text
chi((Q+lambda delta)y+s-v),       0<=lambda<=1.             (16)
```

Both endpoint frequencies are at least one.  For either boundary family at
frequency `R>=1`, write its roots in `[0,1]` as `(m+theta)/R`.  If
`R=M+r`, direct summation in the cases `theta<=r` and `theta>r` gives

```text
sum_(0<=m+theta<=R)(m+theta)<=R^2.                          (17)
```

Distributional differentiation of `(16)` therefore costs at most
`|delta|` per boundary family, hence `2|delta|` per block.  Endpoint roots
follow by one-sided approximation.

Split `(13)` after its first `N` unit blocks; `epsilon` need not be below one.
The first clock has mean `1/7` and mean-zero primitive oscillation
`6/(49P)`.  Combining normalization and the final segment, the within-block
bound `(17)`, the `2`-Lipschitz left Riemann sum for `(14)`, and the endpoint
change from `delta N` to `B/P` proves

```text
|I_n-J| <=
 [2epsilon/7+6/(49P)+|B|(3+epsilon)/P]/(N+epsilon).         (18)
```

With `z0=Lp0-e` and `z=Lp-e`, this is equivalently

```text
|I_n-J| <=
 [2z0/7+6L/49+3|A|/k+|A|z0/(dL)]/z.                       (19)
```

Thus `(11)` and `(19)` give an explicit finite convergence head on every
ray/context row.

## 5. Exact residual and the repaired turn boundary

The compiler retains the actual moving gcd by residue classes modulo `|c|`;
it never treats gcd as a ray coordinate.  It also applies the exact midpoint
bound before literal evaluation.  A second analytic exit uses the generalized
many-turn theorem when

```text
floor(p/d)|A|/z>=5.                                        (20)
```

For the edgewise target `(3)`, this exit is used only at `p>=273` and
`9|c|<=p`.  The
weaker theorem at `264<=p<273` beats `D_max/5` and would suffice only after
summing five edges.  The entire admissible strip is instead sent through the
midpoint or literal exact certificate.

The frozen census is

```text
ray/context rows             1,808,310
convergence exits            1,089,868
many-turn exits                718,442
maximum finite ray head           1,149
literal mass evaluations   239,233,858.                  (21)
```

Its exact literal minimum is

```text
29565457/1887366115
 at (d,a,c,p0,q0;L,j,e,f;n,p,q)
 = (8,16,-8,8,23;168,90,12,1;111,896,2687),               (22)
```

strictly above `1/294`.  The C++ engine uses checked signed `__int128`
floor moments.  Ninety-nine deterministic probes agree with the canonical
Python physical-mass engine and an independent rational implementation of
`(9)`.  Ordinary and optimized wrappers reproduce the same pinned summary;
the hashes are in the front matter.

## 6. Hunter corollary and exact scope

For six distinct levels, the universal singleton debt from THM-3350 is

```text
D_max=186636088362/11773143757375.                         (23)
```

If the low graph is disconnected, its high complement is connected and has
a five-edge spanning tree.  Equations `(3)` and `(23)` give the strict credit

```text
5/294-D_max
 =570672686921/494472037809750>0.                          (24)
```

so Hunter closes the assignment.  THM-3352 handles the complementary
connected-low case.

This independently removes disconnected component scales and, together with
THM-2941's same-level graph, makes the `561`-body reflected ledger empty.  The
primary closure remains THM-3355; the new content here is the stronger uniform
physical edge floor, including its four centered-grid horn lanes.  The conclusion
still uses the reflection relation: the six drifts are `Lq_e-e`, with residue
matched to the body label.  It does not cover arbitrary drifts `Lq_i+b_i`,
where `b_i` is independent of the body label.  That arbitrary-residue `k=1`
sector, projected `k=2,3`, lower-body physical entry, and the six-body/
seven-tail rung remain separate obligations.  In particular, THM-3360 is not
a proof of LRC(14).
