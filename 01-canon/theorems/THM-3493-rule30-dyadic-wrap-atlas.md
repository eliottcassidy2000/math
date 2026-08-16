---
id: THM-3493
title: "Rule 30 dyadic period floor, wrap-prefix atlas, and Mersenne hard face"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  The persistent left boundary pair sharpens the seed-period lower bound to
  the dyadic floor and excludes every two-cycle center endpoint.  Each dyadic
  depth block consequently has either no wrapped depths or one initial wrap
  prefix whose center word is zeros followed by one.  At the complementary
  hard Mersenne endpoint, the pointed profile and current functionals require
  respectively (n+1)/2 and n independently load-bearing individual Hasse
  coordinates on their declared ambient modules.  A 2^27-step exact scout
  certifies that depths 5 through 2^28-1 are hard.  No Rule 30 prize is
  claimed.
source: root-rule30-next-targets-20260816
audit: >
  An independent hostile audit rederived the period floor, innovation
  valuation, block atlas and density reductions, both ambient adaptive-query
  adversaries, and the compiler scope.  A separate two-limb 68-bit scout
  reproduced all 28 committed valuations; ordinary and optimized runs equal
  the stored transcript byte-for-byte: ACCEPT.
depends_on:
  - THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
  - THM-3480-rule30-staircase-transducer-entropy-and-nonrectangular-macroblock-compiler
  - THM-3489-rule30-packed-restart-and-pointed-pascal-face
related:
  - THM-3463-rule30-mealy-section-suffix-parity-current-and-complexity-boundary
  - THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum
script: 04-computation/rule30_dyadic_wrap_atlas_thm3493.py
output: 05-knowledge/results/rule30_dyadic_wrap_atlas_thm3493.out
script_sha256: ef80ad966d70631a4a28937f177655c9e4280dafb9ac16bd32232cc30f31f78c
output_sha256: 0dac3f198615ece230ba2e9b5dd9df74119f19e64a5833d71c6e507f6fc5754b
hash_basis: raw bytes
---

# THM-3493 -- Rule 30 dyadic period floor, wrap-prefix atlas, and Mersenne hard face

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-3489 closes the center whenever its depth has reached the lower-width
seed period, but it deliberately leaves the frequency of such wraps open.
The present theorem extracts the exact dyadic geometry of that question.
The new input is elementary but load-bearing: the two cells on the physical
left boundary are always `11`.  That one extra boundary bit rules out the
only equality case left by THM-3458's period lower bound.

The result has two faces.  On the fixed Rule 30 seed it gives an exact
wrap-prefix atlas and strong conditional density reductions.  On ambient
phase modules it gives a sharp individual-coordinate query obstruction at
the hard Mersenne endpoint.  These scopes must not be identified: the latter
is not a time lower bound for the distinguished seed.

## 1. Conventions, inheritance, and live objects

Use THM-3458's inward right-edge packing

```text
R_0=1,
R_(t+1)=Phi(R_t),
Phi(x)=x xor ((2x) or (4x)).                           (1)
```

If `a_t(j)` is the line evolution from the single seed, then

```text
R_t=sum_(k=0)^(2t) a_t(t-k)2^k,
c_t=a_t(0)=bit_t(R_t).                                (2)
```

Let `P_k` be the least positive seed return time modulo `2^k`, and put

```text
epsilon_k=bit_k(R_(P_k)),
P_(k+1)=2^(epsilon_k)P_k.                              (3)
```

Thus every `P_k` is a power of two.  Enumerate the innovation depths by

```text
kappa_1<kappa_2<...,
epsilon_(kappa_r)=1.                                  (4)
```

They are infinite because the period tower is unbounded, and counting the
earlier innovations gives

```text
P_(kappa_r)=2^(r-1),
P_(kappa_r+1)=2^r.                                    (5)
```

Call a depth

```text
wrapped: P_k<=k,
hard:    P_k>k.                                       (6)
```

The inheritance pass is:

1. closest proved mechanism: THM-3489's packed restart and pointed Pascal
   face;
2. canonical hostile: scalar period towers can have wrap density zero or
   one while obeying the same monotone dyadic inequalities;
3. corrected near miss: the possible endpoint `k=2P_k` must be tested
   against the second physical left-boundary bit, not only the top bit; and
4. least-used sidecar: the valuation of `R_(2^m)-1`, which is the next
   innovation depth and hence the exact cut inside a dyadic block.

The live board is the boundary pair, period tower, innovation valuation,
dyadic wrap cut, hard Mersenne endpoint, pointed Hasse face, and staircase
fallback compiler.

## 2. The boundary pair and the dyadic period floor

### Lemma 2.1 (persistent left boundary pair)

For every `t>=1`,

```text
a_t(-t)=a_t(-t+1)=1.                                  (7)
```

Equivalently,

```text
bit_(2t)(R_t)=bit_(2t-1)(R_t)=1.                      (8)
```

### Proof

At time one the two cells are one.  If the old left pair is `11`, then the
new leftmost cell and its neighbor are respectively

```text
f(0,0,1)=1,
f(0,1,1)=1.                                          (9)
```

Induction proves (7).  Under the packing (2), sites `-t` and `-t+1` have
indices `2t` and `2t-1`, proving (8).

### Theorem 2.2 (dyadic period floor)

For every `k>=1`,

```text
boxed:
P_k >= 2^floor(log_2 k).                              (10)
```

In particular,

```text
boxed:
k != 2P_k                                             (11)
```

for every depth `k`.

### Proof

Put `q=2^floor(log_2 k)`, so `q<=k<2q`.  If `k>q`, THM-3458 gives

```text
P_k>=ceil(k/2)>q/2.                                   (12)
```

Because `P_k` is a power of two, (12) implies `P_k>=q`.

It remains to handle the left endpoint `k=q`.  The old lower bound permits
only one additional case, `P_q=q/2`.  Were that equality true, the seed
return would give

```text
R_(q/2)=1 mod 2^q.                                    (13)
```

But Lemma 2.1 at time `q/2` says that bit `q-1` of this row is one, whereas
bit `q-1` of the seed is zero.  This contradicts (13).  The case `q=1` is
`P_1=1` directly.  Equation (10) follows.

If `k=2P_k`, then `k` is a power of two and `P_k=k/2`, precisely the excluded
left-endpoint case.  This proves (11).

Thus THM-3489's abstract wrap trichotomy sharpens for the actual Rule 30
seed to

```text
k<P_k:       hard pointed terminal arc;
P_k<=k<2P_k: c_k=epsilon_k.                           (14)
```

The possible `k=2P_k` endpoint is not merely unobserved; it cannot occur.

## 3. Innovation valuations and the exact dyadic wrap atlas

For `m>=0`, put

```text
q=2^m,
B_m={q,q+1,...,2q-1},
v_m=kappa_(m+1).                                      (15)
```

### Lemma 3.1 (dyadic row valuation)

For every `m>=0`,

```text
boxed:
v_m=nu_2(R_q-1),                                     (16)
1<=v_m<=2q-1.                                        (17)
```

### Proof

At `k=kappa_(m+1)`, equation (5) says `P_k=q`.  Hence

```text
R_q=1 mod 2^k.                                       (18)
```

The lift bit at this depth is one, so bit `k` of `R_q` is one.  Therefore
`R_q-1` is divisible by `2^k` but not `2^(k+1)`, proving (16).

The packed row is odd.  Subtracting one clears its bit zero without changing
any higher bit.  Lemma 2.1 gives a one at bit `2q-1`, so the least positive
set bit occurs no later than `2q-1`.  This proves (17).  Notice that the
generic top-support argument would give only `v_m<=2q`; the second boundary
bit is exactly what makes the inequality strict.

Applying (17) at the preceding scale also gives, for `m>=1`,

```text
kappa_m<=q-1<q.                                      (19)
```

### Theorem 3.2 (dyadic wrap-prefix atlas)

For every `m>=0`, exactly one of the following holds.

#### Hard block

If `v_m<q`, then

```text
B_m subset H.                                        (20)
```

#### Wrap-prefix block

If `q<=v_m<=2q-1`, then

```text
W intersect B_m={q,q+1,...,v_m},
H intersect B_m={v_m+1,...,2q-1},                    (21)
```

and the complete wrapped center word is

```text
boxed:
(c_q,c_(q+1),...,c_(v_m))=0^(v_m-q)1.                (22)
```

An empty interval in (21) is omitted.  Thus wrapped depths form an initial
prefix of each dyadic block, and every nonempty prefix contains exactly one
center one, at its final depth.

### Proof

The case `m=0` is direct: `R_1=7`, so `v_0=kappa_1=1`, while
`P_1=1` and `epsilon_1=1`.  Thus `B_0={1}` is a one-letter wrap prefix
with center word `1`.

Assume now that `m>=1`.  Then (19) applies.  The `m`th innovation has
already occurred before `q`, and the period is exactly `q` from depth
`kappa_m+1` through depth `kappa_(m+1)=v_m`.

If `v_m<q`, the next innovation has also occurred before the block begins.
Hence `P_k>=2q>k` throughout `B_m`, proving (20).

If `v_m>=q`, then `P_k=q` for `q<=k<=v_m`.  These depths are wrapped, and
the innovation word is zero there except for `epsilon_(v_m)=1`.  Equation
(14) therefore gives (22).  Immediately after `v_m` the period is at least
`2q`, so every remaining `k<2q` is hard.  This proves (21).

### Corollary 3.3 (exact block ledgers)

Define the positive part `x_+=max(x,0)`.  The wrap length and wrapped-one
count in `B_m` are exactly

```text
r_m=(v_m-2^m+1)_+,
u_m=1_[v_m>=2^m].                                    (23)
```

Consequently, for every `M>=1`,

```text
|W intersect [1,2^M-1]|
  =sum_(m=0)^(M-1)r_m,                               (24)

sum_(k in W, 1<=k<=2^M-1)c_k
  =sum_(m=0)^(M-1)u_m
  <=M.                                               (25)
```

For an arbitrary cutoff `N>=1`, a final partial block can still contain at
most one wrapped one, so

```text
boxed:
sum_(k in W, k<=N)c_k <= floor(log_2 N)+1.            (26)
```

Equation (24) is an exact law for wrap frequency.  It does not by itself
bound that frequency away from either zero or one, because the unknown
Rule 30 data are precisely the valuations `v_m`.

## 4. Hard-core density reductions

Put

```text
A(N)=sum_(1<=k<=N)c_k,
A_H(N)=sum_(k in H, k<=N)c_k,
h(N)=|H intersect [1,N]|.                            (27)
```

Equation (26) gives the exact reduction

```text
0<=A(N)-A_H(N)<=floor(log_2 N)+1,                    (28)
A(N)<=h(N)+floor(log_2 N)+1.                         (29)
```

Therefore:

1. if `A(N)/N` has limit `delta`, then

   ```text
   liminf_(N->infinity) h(N)/N >= delta;             (30)
   ```

2. in particular, center balance would force lower hard-depth density at
   least `1/2`; and
3. if `h(N)/N->0`, then the center has one-density zero and is not eventually
   periodic.

For the last statement, (29) first gives `A(N)/N->0`.  There must be
infinitely many nonempty wrap blocks: otherwise all sufficiently large
blocks would be entirely hard, contradicting `h(N)/N->0`.  Theorem 3.2 puts
one center one in every nonempty wrap block, so the center has infinitely
many ones.  An eventually periodic binary sequence of density zero is
eventually zero, a contradiction.

These are conditional prize reductions.  No unconditional hard-depth
density, center density, or nonperiodicity statement has been proved.

## 5. Mersenne easy-block / hard-face dichotomy

Fix `m>=0`, put

```text
q=2^m,
n=2q-1.                                              (31)
```

At the final depth of `B_m`, Theorem 3.2 becomes a sharp dichotomy:

```text
v_m=2q-1:
  (c_q,...,c_n)=0^(q-1)1;                            (32)

v_m<=2q-2:
  n is hard.                                         (33)
```

In the hard branch let `p=P_n`.  Then `p>n` and `p` is a power of two, so

```text
p>=2q.                                               (34)
```

Retain THM-3489's complete terminal profile and current

```text
F_n=A_n Q_n,
M_j(f)=xor_(0<=h<p)binom(h,j)f(h).                   (35)
```

### Theorem 5.1 (profile-coordinate query bound)

In the hard branch,

```text
boxed:
c_n=xor_(a subseteq_bits n-1) M_(p-n+a)(F_n).        (36)
```

This face contains exactly

```text
2^popcount(n-1)=2^m=q=(n+1)/2                       (37)
```

coordinates.  Every displayed order is odd and every coordinate is
independently load-bearing on `im(A_n)=V_p`.

The parity follows because `p-n` is odd while every submask of the even
number `n-1` is even.  Since `n` is odd, THM-3489's image formula has
repeated-root loss `2^nu_2(n)-1=0`; hence `A_n` is invertible and its image is
all of `V_p`.

Declare the following precise query model.  An algorithm receives oracle
access to the individual coordinates `M_j(f)` of an arbitrary
`f in im(A_n)`, may choose orders adaptively, and must output `f(p-n)` for
every such profile.  Its worst-case query count is at least `q`.

Indeed, run it on the zero profile.  If fewer than `q` coordinates are read,
some displayed face order `j_0` is omitted.  The profile `Y^(j_0)`, with
`Y=X+1`, lies in `im(A_n)`, changes only Hasse coordinate `j_0`, and toggles
the marked point because

```text
binom(j_0,p-n)=1 mod 2.                              (38)
```

The adaptive transcript is unchanged, so the same output cannot be correct
on both profiles.

### Theorem 5.2 (fixed-holonomy current-coordinate query bound)

The equivalent current formula is

```text
boxed:
c_n=xor_(a proper-submask-of n)
 M_(1+((p-n-1) bit_or a))(Q_n).                      (39)
```

It contains exactly

```text
2^popcount(n)-1=2^(m+1)-1=n                         (40)
```

individual current coordinates.  All their indices are at least one, so
they exclude the holonomy coordinate

```text
M_0(Q_n)=epsilon_n.                                  (41)
```

Fix either value of `M_0`.  In the query model where an algorithm may read
individual current Hasse coordinates and must compute the terminal
functional for every current in that affine holonomy slice, its worst-case
query count is at least `n`.  The same zero-transcript argument applies after
choosing any base current with the prescribed `M_0`: adding `Y^(j_0)` at an
unqueried displayed order fixes `M_0`, changes no other Hasse coordinate,
and toggles (39).

Theorems 5.1 and 5.2 are ambient-input coordinate lower bounds.  They do not
apply to arbitrary linear-combination queries, do not assert that actual
Rule 30 currents fill the ambient module, and are not word-RAM, circuit, or
time lower bounds for the single seed.

## 6. A conditional staircase batch compiler

THM-3480's staircase algorithm advances an explicitly stored finite row and
then extracts its center.  Stop the same computation at time `q`, retain the
emitted final row `R_q`, and scan its words for the least positive set bit.
This computes

```text
v_m=nu_2(R_q-1)                                      (42)
```

in

```text
O(q^2/log^2 q) word operations,
O(q/log q) words.                                    (43)
```

The final scan costs only `O(q/log q)` word operations.  No new transition
table or uncharged preprocessing is used.

The atlas then emits the entire wrap prefix in `B_m` for free beyond output
cost.  In the easy Mersenne branch `v_m=2q-1`, one height-`q` computation
therefore compiles all `q` center bits in (32).  In every other branch it
certifies the remaining suffix as hard; the original THM-3480 compiler is
still the fixed-seed fallback for a requested hard center.

This is an adaptive branch compiler, not a better worst-case exponent.
Determining `v_m` costs a genuine Rule 30 simulation, and the hard branch may
occur at every sufficiently large scale.

## 7. The finite-exact hard-range scout

The companion performs one declared long computation.  It evolves

```text
x_0=1,
x_(t+1)=x_t xor ((2x_t) or (4x_t)) mod 2^68          (44)
```

for every

```text
0<=t<=2^27,                                          (45)
```

and records `nu_2(x_(2^m)-1)` for `0<=m<=27`.  The
`134217728` update universe is exact and uses constant memory.  Every
recorded valuation is below `68`; triangular compatibility of (1) therefore
proves that it is also the valuation of the untruncated packed row.

The exact list is

```text
(v_0,...,v_27)=
(1,3,4,6,7,9,15,16,24,25,27,29,34,36,
 37,39,41,43,48,49,51,54,55,58,60,63,64,66).         (46)
```

The first three blocks have wrap prefixes

```text
B_0: {1},
B_1: {2,3},
B_2: {4};                                            (47)
```

and `v_m<2^m` for every `3<=m<=27`.  Theorem 3.2 thus gives the
finite-exact conclusion

```text
boxed:
W intersect [1,2^28-1]={1,2,3,4},
H contains [5,2^28-1].                               (48)
```

Numerically, every depth

```text
5<=k<=268435455                                      (49)
```

is certified hard.  This is an exponentially amplified finite certificate,
not an asymptotic statement.

Smaller hostile controls in the companion check full boundary pairs through
time `256`, exact seed cycles through width `34`, the period floor,
no-endpoint and lift laws through width `33`, and both Mersenne functional
masks for `1<=m<=7` and `p/2^m in {2,4,8}`.  Every failure gate is an explicit
exception, never a Python assertion.

## 8. Hostiles and no-prize boundary

Two abstract scalar schedules show why the exact atlas alone cannot decide
wrap frequency.  They are not claimed to arise from the packed Rule 30 map.

1. If

   ```text
   kappa_r=2^r-1,                                    (50)
   ```

   then `v_m=2^(m+1)-1`, every dyadic block is fully wrapped, and its center
   word is `0^(2^m-1)1`.  The hard set is empty and the center one-density is
   zero.

2. If

   ```text
   kappa_r=2r-1,                                     (51)
   ```

   then `v_m=2m+1<2^m` eventually, so the hard set has density one.

Both scalar lift words avoid `111`, obey the strict valuation ceiling and
the dyadic period floor, and have unbounded power-of-two periods.  Additional
Rule 30 correlations across scales are therefore indispensable.

The staircase state count also cannot supply the missing fixed-seed lower
bound.  The same universal Rule 30 transducer contains the all-zero initial
row, whose staircase trajectory and center prediction are trivial.  Ambient
reachable-state richness does not imply richness along one distinguished
input.

Finally, (48) is finite.  THM-3458's periodic-ring construction gives a
genuine eventually periodic Rule 30 extension of every finite center prefix.
Neither a large hard interval nor a fitted density settles an infinite prize.

Accordingly this theorem proves no Rule 30 balance prize, nonperiodicity
prize, random-access prize, or general computational lower bound.

## 9. Verification and status

Reproduce the exact transcript with

```bash
python3 04-computation/rule30_dyadic_wrap_atlas_thm3493.py
python3 -O 04-computation/rule30_dyadic_wrap_atlas_thm3493.py
```

Both modes must agree byte-for-byte with

```text
05-knowledge/results/rule30_dyadic_wrap_atlas_thm3493.out.
```

The long gate makes exactly `2^27` packed updates modulo `2^68`; ordinary
and optimized runs execute the same explicit checks.  An independent
two-limb implementation reproduced the full valuation list, and both replay
modes match the stored transcript byte-for-byte.
