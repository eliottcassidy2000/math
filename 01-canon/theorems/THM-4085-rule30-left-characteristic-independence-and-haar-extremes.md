---
id: THM-4085
title: "Rule-30 left-characteristic independence and Haar extremes"
status: >
  PROVED universal left-characteristic independence theorem + PROVED
  rational-slope entropy floor + PROVED Rule-30 marked-block temporal
  independence and Haar extreme law + VERIFIED-EXACT + INDEPENDENT EXACT
  AUDIT. For every radius-one left-permutive cellular automaton, outputs at
  finitely many spacetime cells with distinct left-characteristic addresses
  are independent uniform symbols under Haar initial data. On a Rule-30
  rational ray this gives full iid traces at nonpositive slopes and an exact
  min-entropy floor at positive subluminal slopes. THM-4050's radius-r marked
  blocks are mutually independent whenever their center times differ by at
  least 2r-1; the threshold is uniformly sharp, since radius two at gap two
  has joint mass 3/128 instead of 1/64. Retaining the finite-line cemetery
  class, the capped maximum stopping address satisfies M_N/log_4(N)->1
  almost surely, and limsup Z_k/log_4(k)=1 almost surely after the finitely
  many possible cemetery values. These are Haar-random-row theorems, not
  distinguished-single-seed statements; every Rule-30 prize remains OPEN.
source: codex-frontier-synthesis-creative-20260825f / moving-observer niche
audit: >
  PASS. The primary direct-orbit path exhausts all sixteen binary
  left-permutive elementary rules on 131 distinct-address point sets, six
  rational rays, three separated marked-block packets, the exact gap-two
  hostile, and cemetery classes through depth six. It performs 374164 exact
  assignments and 20876 consequence gates. The no-import audit instead
  builds Boolean algebraic normal forms through depth six, verifies that the
  extreme-left variable occurs as its lone pivot monomial, replays 11226
  triangular inverses, symbolically reduces the hostile second block, and
  solves both cemetery systems. It performs 26340 monomial gates and 217418
  evaluation gates. Normal and optimized streams match the frozen outputs;
  both scripts contain zero assert nodes and zero float literals.
depends_on:
  - THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary
  - THM-4050-rule30-half-arc-marked-cylinder-and-radius-nine-hostile
related:
  - THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
  - THM-4048-rule30-periodicity-balance-and-model-firewalls
  - THM-4064-rule30-cyclotomic-kernel-character-and-c60-alias-obstruction
  - THM-4065-rule30-temporal-cylinder-transfer-infinite-zero-columns-and-next-zero
script: 04-computation/rule30_left_characteristic_independence_thm4085.py
output: 05-knowledge/results/rule30_left_characteristic_independence_thm4085.out
script_sha256: aee4cae0ab8f4aaad58bbe410b01fc24a6bcf48d4c5fdc2e5f658e5c7ff6fb7f
output_sha256: 3e9dd6cc40ff3cdaa688a486c5b6bf180a0a3fd172bcd898b043607c8ad372c2
independent_audit_script: 04-computation/rule30_left_characteristic_independence_thm4085_independent_audit.py
independent_audit_output: 05-knowledge/results/rule30_left_characteristic_independence_thm4085_independent_audit.out
independent_audit_script_sha256: 816865aa2a750eb0640b4a734ac9fe9d5a65e93a24df9371db3759db93fa2682
independent_audit_output_sha256: 31bef89c56340b4cb4b8f5625a5d7b2c78ab667182b7fc39071bbc9c41735954
hash_basis: raw LF bytes
---

# THM-4085 -- left-characteristic independence and Haar extremes

**PROVED in the universal and Haar scopes stated below.** THM-3456 makes the
vertical center trace of a left-permutive cellular automaton a free Haar
coordinate. The actual fresh coordinate is not vertical time: it is the
leftmost characteristic address `j-t` of a spacetime cell. Distinct addresses
remain independent even when their full dependency cones overlap. This
observation both handles rationally moving observers and turns THM-4050's
one-time marked-cylinder law into an exact sparse temporal independence law.

The consequence is an almost-sure logarithmic extreme law for the Haar
stopping address. It is sharply separated from the distinguished Rule-30
single seed. No random-row typicality is transferred to that one named orbit.

## 1. The universal left-characteristic theorem

Let `A` be a finite alphabet of size `Q>=2`, and let

```text
F:A^Z -> A^Z,
(F(x))_j=f(x_(j-1),x_j,x_(j+1))                     (1)
```

be radius one and **left permutive**: for every fixed `(c,r)`, the map

```text
l |-> f(l,c,r)                                       (2)
```

is a permutation of `A`. For a spacetime cell `v=(t,j)` put

```text
lambda(v)=j-t,               rho(v)=j+t.             (3)
```

Thus `(F^t(x))_j` depends only on the initial interval
`[lambda(v),rho(v)]`.

### Theorem 1.1 (finite characteristic-coordinate bijection)

Let `v_1,...,v_m` be any finite family of spacetime cells whose addresses

```text
lambda_i=lambda(v_i)                                 (4)
```

are pairwise distinct. Fix every initial coordinate except

```text
x_(lambda_1),...,x_(lambda_m).                       (5)
```

Then the map

```text
(x_(lambda_i))_(i=1)^m
  |-> ((F^(t_i)(x))_(j_i))_(i=1)^m                  (6)
```

is a bijection `A^m -> A^m`.

Consequently, under the uniform Bernoulli product measure on `A^Z`, the
outputs in `(6)` are independent uniform `A`-valued symbols. More strongly,
their conditional law remains uniform after every nonpivot initial
coordinate is fixed. In the binary case every nonempty Walsh character of
the output vector therefore has expectation zero.

### Proof

For one cell `v=(t,j)`, the output `(F^t(x))_j` is a permutation of the
extreme input `x_(j-t)` when all other inputs are fixed. At `t=0` this is the
identity. For the induction step,

```text
(F^t(x))_j
 =f((F^(t-1)(x))_(j-1),
    (F^(t-1)(x))_j,
    (F^(t-1)(x))_(j+1)).                             (7)
```

Only the first argument in `(7)` depends on `x_(j-t)`. By induction it is a
permutation of that variable, and `(2)` composes it with another permutation.
This is THM-3456's extreme-input mechanism translated to an arbitrary cell.

Now order the cells so that

```text
lambda_1>lambda_2>...>lambda_m.                      (8)
```

The first output depends on its own pivot and on coordinates strictly to its
right, but on no later pivot in `(8)`. It therefore chooses a unique
`x_(lambda_1)` for any prescribed first output. After that choice, the second
output uniquely chooses `x_(lambda_2)`, and so on. This triangular inverse
constructs exactly one pivot vector for every target vector, proving `(6)`.
Uniform independent pivots and a conditional bijection give the Haar claim.

Dependency intervals may overlap arbitrarily in this proof. The load-bearing
hypothesis is distinctness of the **extreme pivots**, not disjointness of the
full cones. `QED`

## 2. Rationally moving observers

Let `b,p` be integers, let `q>=1`, and sample the moving path

```text
j_t=b+floor(pt/q),
Y_t=(F^t(X))_(j_t),                 t>=0,            (9)
```

from Haar initial data `X`. Restrict first to physical slopes
`-q<=p<=q`.

If `-q<=p<=0`, then `j_(t+1)-j_t` is `-1` or `0`, so

```text
lambda_(t+1)-lambda_t
  =j_(t+1)-j_t-1 in {-2,-1}.                        (10)
```

Every address is distinct. The complete moving trace `(Y_t)` is therefore an
iid uniform process. This includes the vertical center at `p=0`, but also
every rational left-moving observer in the light cone.

For `0<=p<=q`, the addresses are nonincreasing and change by `0` or `-1`.
Among `t=0,...,N-1` they take exactly

```text
d_N=N-floor(p(N-1)/q)                               (11)
```

distinct values. Choose one time over each value. The corresponding `d_N`
outputs are iid, so for every word `w in A^N`,

```text
P((Y_0,...,Y_(N-1))=w) <= Q^(-d_N).                 (12)
```

Thus the base-`Q` min-entropy and Shannon entropy are at least `d_N`, and the
directional entropy-rate floor is

```text
liminf_(N->infinity) H_infinity(Y_[0,N))/N
 >= 1-p/q.                                          (13)
```

This is a floor, not an equality. Repeated-address outputs can still acquire
additional entropy from right-hand coordinates.

### The sharp repeated-address hostile

Distinctness cannot be omitted. For Rule 30, use the cells `(0,0)` and
`(1,1)`, both with left address zero. If the initial bits are `x_0,x_1,x_2`,
then

```text
Y_0=x_0,
Y_1=x_0+x_1+x_2+x_1 x_2=x_0+(x_1 or x_2).          (14)
```

Their exact joint counts over the eight input triples are

```text
(00,01,10,11) : (1,3,3,1).                          (15)
```

Both marginals are fair, but

```text
E((-1)^(Y_0+Y_1))=-1/2.                             (16)
```

The shared pivot cancels in the character, exposing the biased right-tail
contrast `x_1 or x_2`. This is the smallest address-collision sidecar. It also
marks the right-characteristic slope-one boundary: one free Haar symbol per
distinct address is guaranteed, not one per observed cell.

## 3. Marked blocks become temporally independent on an exact lattice

Now specialize to Rule 30 from a Haar-random initial row `X`. Write

```text
a_t(j)=(F^t(X))_j.                                   (17)
```

Retain THM-4050's terminal line and finite stopping address:

```text
T_k(h)=a_(k+h)(h),                  -k<=h<=0,
Z_k=min{a in {1,...,k}:T_k(-a)=0},
Z_k=infinity if this finite line contains no zero.  (18)
```

The cemetery value in `(18)` is mandatory. MISTAKE-492 records the omission
of that class from an earlier mass partition.

For `1<=r<=k`, define the length-`2r-1` marked block

```text
B_(k,r)
 =(a_(k-r)(-r),a_(k-r)(-r+1),...,a_(k-r)(r-2)).     (19)
```

THM-4050's right-characteristic lemma gives the equivalence

```text
Z_k>r                     iff
B_(k,r)=1 0^(2r-2).                                  (20)
```

The cells in `(19)` have the consecutive left-address interval

```text
I_(k,r)
 ={-k,-k+1,...,2r-2-k}.                              (21)
```

### Theorem 3.1 (exact separated-block independence)

Fix `r>=1`. Let `K` be any finite set of indices `k>=r` satisfying

```text
|k-k'|>=2r-1                 for distinct k,k' in K. (22)
```

Then the vectors

```text
(B_(k,r))_(k in K)                                  (23)
```

are mutually independent and each is uniform on `{0,1}^(2r-1)`. In
particular, the events `{Z_k>r}`, `k in K`, are iid Bernoulli with parameter

```text
p_r=2^(-(2r-1)).                                    (24)
```

Indeed, `(22)` makes the address intervals `(21)` pairwise disjoint. Every
cell in the concatenated blocks therefore has a distinct left address, and
Theorem 1.1 applies to all `(2r-1)|K|` cells at once. Notice that their full
initial dependency cones usually overlap; no spatial mixing estimate is
being inserted.

### The uniform threshold is sharp

Take `r=2`. Then each marked word is `100` and the separation in `(22)` is
three. At `k=2` the event fixes

```text
(x_-2,x_-1,x_0)=(1,0,0).                            (25)
```

At `k=4`, put

```text
(a,b,c,d)=(x_-4,x_-3,x_1,x_2).                      (26)
```

Two exact Rule-30 steps reduce the second marked block to

```text
(a,b,1+c+d+cd)=(a,b,not(c or d)).                   (27)
```

For it also to equal `100`, one needs `a=1`, `b=0`, and `(c,d)!=(0,0)`.
Exactly three of the sixteen free assignments work. Hence

```text
P(Z_2>2 and Z_4>2)=3/128,
P(Z_2>2)P(Z_4>2)=1/64=2/128.                        (28)
```

The gap-two intervals overlap at the pivot `-2`, and the events are not
independent. At gap three they are independent by Theorem 3.1, with exact
joint mass `1/64=8/512`. Therefore the uniform separation `2r-1` cannot be
reduced to `2r-2`, already at radius two.

## 4. The Haar stopping radius has logarithmic extremes

THM-4050 gives, for every `1<=r<=k`,

```text
P(Z_k>r)=2^(-(2r-1)),
P(Z_k=infinity)=2^(-(2k-1)).                         (29)
```

The second series in `(29)` is summable. The first Borel--Cantelli lemma
therefore shows that almost surely only finitely many `Z_k` equal infinity.
To state a deterministic maximum before that random last cemetery time, put

```text
bar Z_k = Z_k       if Z_k<infinity,
bar Z_k = k+1       if Z_k=infinity,
M_N=max_(1<=k<=N) bar Z_k.                           (30)
```

### Theorem 4.1 (extreme-value exponent)

Almost surely,

```text
M_N/log_4(N) -> 1.                                  (31)
```

On the same probability-one event of eventual finiteness,

```text
limsup_(k->infinity) Z_k/log_4(k)=1.                 (32)
```

### Upper bound for `(31)`

Fix `epsilon>0`, let `N_m=2^m`, and put

```text
u_m=ceil((1+epsilon)log_4(N_m)).                     (33)
```

If `k<u_m`, then `bar Z_k<=k+1<=u_m`. For `k>=u_m`, `(29)` and the union
bound give

```text
P(M_(N_m)>u_m)
 <= 2 N_m 4^(-u_m)
 <= 2 N_m^(-epsilon).                               (34)
```

The right side is summable in `m`. Borel--Cantelli proves the dyadic upper
bound almost surely, and monotonicity of `M_N` interpolates between successive
dyadic `N_m`. Letting `epsilon` decrease through positive rationals gives

```text
limsup M_N/log_4(N)<=1.                              (35)
```

### Lower bound for `(31)`

Fix `0<epsilon<1` and put

```text
r_m=floor((1-epsilon)log_4(N_m)),
L_m=2r_m-1.                                         (36)
```

Choose an arithmetic progression `K_m` inside `[N_m/2,N_m]` with step
`L_m`. For all large `m` it has

```text
|K_m|>=N_m/(2L_m)-1.                                (37)
```

By Theorem 3.1, the events `{Z_k>r_m}`, `k in K_m`, are independent. Their
common probability is

```text
p_m=2^(-L_m)=2*4^(-r_m)
    >=2 N_m^(-(1-epsilon)).                          (38)
```

Thus, for one absolute `c>0` and all large `m`,

```text
P(no k in K_m has Z_k>r_m)
 =(1-p_m)^|K_m|
 <=exp(-|K_m|p_m)
 <=exp(-c N_m^epsilon/log(N_m)).                    (39)
```

The final series is summable in `m`. Borel--Cantelli shows that almost every
row has such an exceedance in every sufficiently late dyadic block. Hence

```text
liminf M_N/log_4(N)>=1.                              (40)
```

again using monotonic interpolation and then positive rational `epsilon`.
Equations `(35)` and `(40)` prove `(31)`.

### Proof of `(32)`

For the upper bound, set

```text
v_k=floor((1+epsilon)log_4(k)).                      (41)
```

For all large `k`, `(29)` gives

```text
P(Z_k>v_k)
 =2*4^(-v_k)
 <=8 k^(-(1+epsilon)).                              (42)
```

This series is summable, so almost surely only finitely many such exceedances
occur. For the lower bound, `(39)` supplies infinitely many
`k_m in [N_m/2,N_m]` with

```text
Z_(k_m)>r_m=(1-epsilon+o(1))log_4(N_m)
           >=(1-epsilon+o(1))log_4(k_m).             (43)
```

The cemetery values are already eventually absent. Letting `epsilon` tend to
zero proves `(32)`.

The proof gives only the leading extreme exponent. The overlapping events at
gaps below `2r-1`, including `(28)`, prohibit silently replacing the process
by a fully independent one; no Gumbel law or complete temporal mixing theorem
is asserted.

## 5. Preservation, loss, and prize boundary

| source -> target | preserves | destroys / fails to control | needed sidecar / boundary |
|---|---|---|---|
| spacetime cell -> left address | its extreme permutive input | the rest of its cone | full nonpivot boundary for inversion |
| distinct cells -> pivot vector | exact Haar joint law | contrasts inside a repeated address | collision character such as `(16)` |
| rational ray -> distinct-address selector | `d_N` iid outputs and entropy floor | entropy from repeated addresses | within-address response |
| marked block -> address interval | exact word and exceedance event | overlaps at shorter time gaps | separation `2r-1` |
| fixed-radius marginal -> extremes | tail exponent and sparse iid trials | full dependence law | all-time union bound plus separated trials |
| Haar initial row -> named single seed | local Rule-30 equations only | random fresh pivots and almost-sure conclusions | THM-3456's inverse seed boundary |

This theorem is not another fixed-column periodicity statement. Its vertices
are arbitrary spacetime cells and rationally moving observers; its operation
is conditioning on characteristic pivots. Conversely, it does not touch the
single-seed center except as a firewall. The one-seed initial row has no Haar
freshness, and none of `(12)`, `(24)`, `(31)`, or `(32)` transfers to it.

Therefore the following remain **OPEN**:

1. boundedness or unboundedness of THM-4050's addresses for the single seed;
2. center non-eventual-periodicity and limiting balance for that seed;
3. every fixed-seed query-complexity lower bound and every Rule-30 prize;
4. a full dependence or limiting extreme distribution for overlapping Haar
   stopping events.

THM-4064's cyclotomic character and THM-4065's temporal-cylinder transfer
live on fixed eventually periodic left-front words. The Walsh collision in
`(16)` is a different, Haar spacetime character. Their common lesson is typed:
a quotient can retain every marginal while losing the character selected by
a moving observer.

## 6. Exact verification

Run from the repository root:

```bash
python3 04-computation/rule30_left_characteristic_independence_thm4085.py
python3 -O 04-computation/rule30_left_characteristic_independence_thm4085.py
python3 04-computation/rule30_left_characteristic_independence_thm4085_independent_audit.py
python3 -O 04-computation/rule30_left_characteristic_independence_thm4085_independent_audit.py
```

The primary script evolves finite rows directly. It exhausts all sixteen
binary left-permutive elementary rules on the bounded point atlas, verifies
the rational-ray consequence rather than only address counts, makes the
separated marked blocks jointly uniform, derives `(27)` by direct Rule-30
evolution, and retains both cemetery masses from MISTAKE-492.

The independent audit never imports those distributions. It computes Rule
30 in the Boolean algebraic-normal-form quotient, checks the unique pivot
monomial through depth six, constructs bounded triangular inverses, obtains
`(27)` by symbolic substitution, and solves the finite cemetery systems. The
universal claims are the proofs above; the finite replays are hostile controls
and implementation audits, not extrapolations.

**QED.**
