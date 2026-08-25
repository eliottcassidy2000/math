---
id: THM-4048
title: "Rule 30 seeded periodicity, dyadic balance, and query-model firewalls"
status: >
  PROVED REDUCTIONS + VERIFIED-EXACT + HOSTILE-AUDITED. For each fixed onset
  and period, center periodicity from that onset is exactly equivalent to
  vanishing of the corresponding seeded candidate's inverse left-boundary
  word; eventual periodicity is the existence of one such vanishing word.
  Its first nonzero bit otherwise occurs at the first failed period relation
  plus the period. THM-3493 wrap prefixes obey a
  new residue rigidity: under period T, late positive wrap lengths are bounded
  by T and constant on each eventual 2-power residue class modulo T. Prefix
  Cesaro balance is exactly equivalent to o(2^m) maximal discrepancy on
  dyadic shell prefixes and forces the dyadic wrap length to be sublinear.
  A target-blind low-owner residue observer gives a sufficient marked-current
  criterion after a wrap-controlled prefix. Exact center readout needs a
  base-addition carry, while fixed residue precision is not closed under Haar
  lifting because of a second quotient carry.
  For any observer, conditional variation equals block length minus twice
  optimal prediction error and bounds discrepancy; observer refinement
  improves prediction but weakens this balance certificate. An
  explicit balanced O(log n)-computable sequence has linear-scale
  right-special histories, refuting history depth as a query lower bound.
  Exact Rule 30 banks stop at 2^18 or 65535 as declared. The official Prize-3
  wording leaves non-equivalent model/lower-bound choices; no prize is solved.
source: codex/session-high-value-rule30-20260824, 2026-08-24
audit: >
  PASS. Two independent agents audited the first-defect index, wrap residue
  quantifiers, shell equivalence, observer mass identities, refinement
  direction, complexity-model fork, and easy right-special hostile. A final
  hostile audit found and retracted the false high-owner-only shell readout;
  it verified the repaired base carry, target-coordinate truncation,
  wrap-to-hard-density consequence, Haar quotient carry, and finite pairing
  tables. Packed prefixes are cross-checked against direct CA through time 512
  and THM-3468's frozen 2^18 prefix. Normal and optimized executions of all
  four companions match; every finite and all-scale scope is explicit.
depends_on:
  - THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary
  - THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
  - THM-3493-rule30-dyadic-wrap-atlas
  - THM-3507-rule30-normalized-dyadic-displacement-sibling-trace-and-assouad-spectrum
  - THM-3516-rule30-marked-van-der-put-carry-and-power-section-bridge
  - THM-4013-rule30-depth-fifteen-history-failure-and-adaptive-routed-repair
related:
  - THM-3468-rule30-radial-green-fold-innovation-discrepancy-and-fixed-seed-carrier-boundaries
  - THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum
  - THM-3500-rule30-dyadic-section-cut-defect-and-cross-depth-valuation-carrier
  - THM-4047-rule30-left-front-affine-monodromy-clock
  - THM-4050-rule30-half-arc-marked-cylinder-and-radius-nine-hostile
external:
  - Stephen Wolfram, "Announcing the Rule 30 Prizes", https://writings.stephenwolfram.com/2019/10/announcing-the-rule-30-prizes/ (2019; CITED for the problem statements and displayed Prize-3 formalization only)
  - Wolfram Rule 30 Prizes, https://rule30prize.org/ (CURRENT OFFICIAL LISTING checked 2026-08-24; problem statements and submission status only)
script: 04-computation/rule30_seeded_period_box_thm4048.py
output: 05-knowledge/results/rule30_seeded_period_box_thm4048.out
script_sha256: 32594d42fb0df7b0b8e3ceeb2efe4c7288a5b7161ca4ba011914a42a29b60d04
output_sha256: 74b344c9d7a0f3b52c13b8f1ccb76f1aa65500b1942d730a398ec0636c4bb432
secondary_script: 04-computation/rule30_shell_history_balance_thm4048.py
secondary_output: 05-knowledge/results/rule30_shell_history_balance_thm4048.out
secondary_script_sha256: a3a06c349b514614e861401586e76938d33b314ddd43a222148b7b1420c34980
secondary_output_sha256: 6536724391fa2b29c80f4a796b849e2660941b564d9b9289cb2f6b2a561d16a0
tertiary_script: 04-computation/rule30_observer_mass_reclassification_thm4048.py
tertiary_output: 05-knowledge/results/rule30_observer_mass_reclassification_thm4048.out
tertiary_script_sha256: fda5f43afb4fae516b0ccafe3ed8a1df26bf1dce9b9584d60c0a4d714b105c9f
tertiary_output_sha256: 25773ca2ea39fd855395021bf01edb50b620de06f1d11c6a3cf7d4728a205271
quaternary_script: 04-computation/rule30_low_owner_quotient_carry_thm4048.py
quaternary_output: 05-knowledge/results/rule30_low_owner_quotient_carry_thm4048.out
quaternary_script_sha256: 2f702ee9b122feb13828bf3d7395cd2fbdc066ec8de4d768205215e2c5b36032
quaternary_output_sha256: 8f93a0b1184541268c9a5f7a097f0b5f46566150d231a4dcf45bddcd8c63aefa
hash_basis: raw LF bytes
---

# THM-4048 -- exact interfaces for the three Rule 30 prizes

**PROVED REDUCTIONS + VERIFIED-EXACT + HOSTILE-AUDITED.** This theorem does
not promote finite irregularity to an asymptotic conclusion. It identifies
the proof obligations which survive the seed, balance, and machine-model
quotients, and gives cheap hostile sequences for the invalid implications.

Let `c_n=a_n(0)` be the center column of single-seed Rule 30.

## 1. Prize 1: a seeded inverse-boundary certificate

Fix `N>=0` and `T>=1`. Define the uniquely seeded `T`-periodic candidate

```text
y_n=c_n                                      (n<N),
y_n=c_(N+((n-N) mod T))                      (n>=N).  (1)
```

Thus `y` copies the whole physical block `c_N,...,c_(N+T-1)` before it
repeats. Apply THM-3456's inverse trace compiler with positive initial ray
`0^infinity`, and write its negative ray as

```text
beta^(N,T)=(beta_1,beta_2,...),       beta_d=x_(-d).  (2)
```

Then

```text
c is T-periodic from N    iff    beta^(N,T)=0^infinity.      (3)
```

If the relation fails, put

```text
e=min{n>=N:c_(n+T)!=c_n}.                              (4)
```

The first candidate/physical trace mismatch and the first nonzero inverse
bit occur exactly at

```text
d=e+T,
beta_1=...=beta_(d-1)=0,       beta_d=1.              (5)
```

### Proof

The copied seed block makes `y=c` through `N+T-1`. If equality holds for
every relation index `N<=k<e`, induction on `j=k+T` gives `y_j=c_j` through
`e+T-1`; at `e+T`, equation (1) gives `y_(e+T)=y_e=c_e!=c_(e+T)`.
THM-3456's triangular bijection therefore reconstructs the same zero left
bits as the physical seed before depth `d`. At depth `d`, left permutivity
forces the other binary input, namely one. Finally, `beta=0^infinity`
reconstructs the single-seed configuration itself, proving (3).

Consequently Prize 1 is exactly the obligation

```text
for every N>=0 and T>=1, some d>=N+T has beta^(N,T)_d=1.     (6)
```

Each witness is finite and triangular; no uniform method for producing it is
yet proved. The exact center prefix through `65535` proves (6) simultaneously
for every `0<=N<=32767`, `1<=T<=32768`. Across those shifts the longest run
of equal comparisons is `29`. This is **FINITE-EXACT**, not eventual
nonperiodicity.

## 2. Dyadic wrap-period rigidity

Use THM-3493's proved quantities

```text
q_m=2^m,
v_m=nu_2(R_(q_m)-1),
r_m=max(v_m-q_m+1,0).                                 (7)
```

Whenever `r_m>0`, its physical center block is exactly

```text
(c_(q_m),...,c_(q_m+r_m-1))=0^(r_m-1)1.              (8)
```

Assume hypothetically that `c_(n+T)=c_n` for every `n>=N`. Call a positive
wrap **late** when `q_m>=N`, and define its forced-one endpoint

```text
E_m=q_m+r_m-1.                                        (9)
```

Then all of the following are necessary.

1. `r_m<=T`.
2. No late endpoint `E_m` is congruent modulo `T` to a forced-zero position
   `z=q_j+k>=N`, where `0<=k<=r_j-2`.
3. If `q_m,q_j>=N`, `q_m=q_j mod T`, and both wrap lengths are positive,
   then `r_m=r_j`.

For (1), if `r_m>=T+1`, the forced one at `E_m` equals the forced zero at
`E_m-T`; both positions are at least `q_m>=N`. Statement (2) is the same
contradiction across two blocks. For (3),
if `0<r_m<r_j`, offset `r_m-1` is the shorter block's forced one and the
longer block's forced zero.

Write `T=2^a u` with `u` odd, and let `L=ord_u(2)`, taking `L=1` for `u=1`.
For `m>=a`,

```text
q_(m+L)=q_m mod T.                                    (10)
```

Hence on every sufficiently late class `m mod L`, all positive `r_m` have
one common value. Unbounded positive wraps, or recurrent unequal positive
lengths in one class, would prove nonperiodicity. What remains **OPEN** is an
autonomous cross-scale law for the physical `v_m` or `r_m`; THM-3500 retains
the full current sidecar rather than such a scalar recurrence.

## 3. Prize 2: exact dyadic shell criterion

Put

```text
s_n=2c_n-1,
D(N)=sum_(0<=n<N) s_n.                                (11)
```

The official prefix-Cesaro balance question is

```text
D(N)=o(N).                                            (12)
```

It is not uniform-in-start or Banach-density balance. For each dyadic shell,
define

```text
S_m(t)=sum_(0<=u<t) s_(2^m+u),       0<=t<=2^m,
A_m=max_(0<=t<=2^m)|S_m(t)|.                          (13)
```

Then the exact equivalence is

```text
D(N)=o(N)                 iff                 A_m=o(2^m).  (14)
```

If (12) holds, `S_m(t)=D(2^m+t)-D(2^m)` gives the forward implication
uniformly on the shell. Conversely,

```text
D(2^m)=s_0+sum_(j<m)S_j(2^j).                         (15)
```

A finite-head/geometric-tail split makes the right side `o(2^m)` when
`A_j=o(2^j)`. For `N=2^m+t`, use
`|D(N)|<=|D(2^m)|+A_m`.

The exact amplitudes through shell `m=17` are

```text
(1,1,2,2,3,3,8,16,12,52,36,25,65,183,249,180,367,889).  (16)
```

They are finite evidence only.

## 4. Observer error mass, not mismatch support

For a finite block `I`, an observer `pi:I->S`, and labels `c_i`, put

```text
n_(s,b)=#{i in I:pi(i)=s and c_i=b},
delta_s=n_(s,1)-n_(s,0),
E_pi=sum_s min(n_(s,0),n_(s,1)),
V_pi=sum_s |delta_s|.                                 (17)
```

Then exactly

```text
V_pi=|I|-2E_pi,
sum_(i in I)(2c_i-1)=sum_s delta_s,
|sum_(i in I)(2c_i-1)|<=V_pi.                         (18)
```

Here `E_pi` is the Hamming distance to the nearest `pi`-measurable labeling;
equivalently, it is the fewest records changed or deleted to make every fibre
label-compatible. Under the uniform law on `I`,

```text
V_pi/|I| = E | E[s_i | pi] |.                         (19)
```

If `pi'` refines `pi`, the triangle inequality gives

```text
E_(pi')<=E_pi,                    V_(pi')>=V_pi.       (20)
```

Thus added history improves deterministic prediction but weakens this direct
balance upper bound. The number of mixed fibres is only a support count: it
is positive iff `E_pi>0`, but is neither the error mass nor monotone under
refinement.

On shell prefixes `I_(m,t)=[2^m,2^m+t)`, let

```text
B_m(pi)=max_(0<=t<=2^m)(t-2E_(m,pi)(t)).              (21)
```

Then `A_m<=B_m(pi)`. Therefore `B_m(pi)=o(2^m)` proves Prize 2. For the
one-state observer, equality holds, so this is an exact reformulation rather
than a free simplification. The next proof-grade target is an explicitly
defined, all-scale, carry-blind coarse physical observer with a sign-reversing
pairing inside its fibres and only `o(2^m)` unpaired records, uniformly over
shell prefixes.

### 4.1 THM-4013 reclassification

THM-4013's projective histories are scale histories, not suffixes of center
bits. Its reported mismatch counts concern next-bank transition targets, not
center-only prediction, and its full carry state already determines center.
On the canonical records `1<=n<=65535`, deleting the whole carry state gives
optimal center errors at projective depths `0,...,6`

```text
(26149,20917,8661,1760,307,74,43).                    (22)
```

Keeping target-XOR but deleting the lower carry bit gives

```text
(23,18,5,2,0,0,0).                                   (23)
```

The first zero in (23) uses `65214` observer fibres for `65535` records and
is therefore nearly injective. Full carry has zero center error tautologically
at every tested history. None is an all-scale balance certificate.

### 4.2 A low-owner shell current and its two carry firewalls

Let `R_t` be the packed physical row from THM-3507/3516, and use their
ordinary-integer quantities

```text
J(t)=(R_t-1)/2,       q=q_m=2^m,       d_m=v_(m+1)-v_m,
U_m(r)=(R_(r+q)-R_r)/2^(v_m).                         (23a)
```

For every `0<=r<q`, the exact shell identity is

```text
J(q+r)=J(r)+2^(v_m-1)U_m(r).                          (23b)
```

THM-3516 gives `c_n=bit_(n-1)(J(n))`. Put `n=q+r`. If `n<v_m`, the shifted
term in (23b) lies above the target, while `J(r)` is supported through bit
`2r-1<n-1`; hence `c_n=0`. If `n>=v_m`, put

```text
L=n-v_m,              B_m(r)=floor(J(r)/2^(v_m-1)).
```

Splitting `J(r)` below bit `v_m-1` in (23b) gives the repaired exact readout

```text
c_n=bit_L(B_m(r)+U_m(r)).                              (23c)
```

The support bound gives `bit_L(B_m(r))=0`, but its lower bits can still carry
into `L`. Equivalently,

```text
c_n=bit_L(U_m(r)) xor gamma_m(r),
gamma_m(r)=1_[(B_m(r) mod 2^L)+(U_m(r) mod 2^L)>=2^L]. (23c')
```

For `L=0`, both residues modulo one vanish. The minimal physical carry is

```text
(m,q,v_m,r,n,L)=(2,4,4,2,6,2),
J(2)=12, B_2(2)=1, U_2(2)=399,
bit_2(U_2(2))=1, gamma_2(2)=1, c_6=bit_2(400)=0.       (23c'')
```

Thus the lower-shell term has no target bit of its own, but discarding its
base carry is unsound.

This exposes a direct bridge from Prize 1's wrap clock to Prize 2. By
THM-3493, `r_m<=q` and the shell begins with `0^(r_m-1)1` whenever `r_m>0`.
Taking the all-zero prefix gives

```text
A_m>=max(r_m-1,0).                                    (23d)
```

Consequently Prize-2 balance would force `r_m=o(q_m)`. If the wrap length is
not sublinear, the wrap prefix instead refutes balance. Moreover, THM-3493's
exact ledger

```text
|W intersect [1,2^M-1]|=sum_(m<M)r_m
```

and a finite-head/geometric-tail split then give
`|W intersect [1,N]|=o(N)` at arbitrary cutoffs. Thus balance would sharpen
THM-3493's necessary lower hard-depth density `1/2` to
`|H intersect [1,N]|/N->1`.

Fix `K>=1` and discard only the initial overflow length

```text
e_(m,K)=min(q,max(0,v_m+K-q))<=r_m+K-1.               (23e)
```

For `r>=e_(m,K)`, one has `L=q+r-v_m>=K`. The coarse physical observer

```text
pi_(m,K)(r)=U_m(r) mod 2^K
```

has at most `2^(K-1)` odd-residue fibres. It is determined by `R_(q+r)` and
`R_r` only modulo `2^(v_m+K)`; because `q+r>=v_m+K`, this truncation excludes
the center coordinate at bit `q+r` of `R_(q+r)`. It also forgets `B_m(r)`,
the moving high digit in (23c), and the full base carry in (23c'). Thus it is
genuinely target-coordinate-blind, not a target decoder. Put
`sigma_m(r)=2c_(q+r)-1` and, for odd `a mod 2^K`,

```text
C_(m,K,a)(t)=sum_(e_(m,K)<=r<t,
                   U_m(r)=a mod 2^K) sigma_m(r).       (23f)
```

The observer identity gives the exact conditional-variation bound

```text
A_m <= e_(m,K)
       +max_(e_(m,K)<=t<=q) sum_(a odd mod 2^K)|C_(m,K,a)(t)|.  (23g)
```

Thus it suffices, for one fixed `K`, to prove `r_m=o(q_m)` and an `o(q_m)`
bound for the marked current in (23g). For `K=1` there is only one fibre, so
this is the original signed discrepancy after deleting the controlled prefix,
not a simplification. There is a single causal partial involution realizing
the variation at every prefix: in each ordered fibre, greedily match an
arriving record to an earlier unmatched opposite sign. The remaining queue
has one sign and size `|C_(m,K,a)(t)|`. This involution reads the center labels;
an explicit label-blind operation pairing would be a stronger structural
proof mechanism, but no additional crossing hypothesis is needed for (23g).

The tempting fixed-precision scale induction is false. Exact Haar lifting is

```text
U_(m+1)(r)=(U_m(r)+U_m(r+q_m))/2^(d_m).               (23h)
```

The child modulo `2^K` needs the numerator modulo `2^(K+d_m)`, not merely the
two sibling residues modulo `2^K`. This loss occurs on the physical orbit in
the first wholly target-blind `K=2` shell: at `m=3`,

```text
q=8, v_3=6, d_3=1, q-v_3=K=2,

(U_3(0),U_3(8)) mod 8=(7,7),
(U_3(3),U_3(11)) mod 8=(7,3).                         (23i)
```

Both sibling pairs equal `(3,3)` modulo four, but after division in (23h)
their children are congruent respectively to `(7+7)/2=7=3 mod 4` and
`(7+3)/2=5=1 mod 4`. Thus (23h) is not closed on fixed low residues in
general; already the physical `K=2` quotient fails. An induction through
(23h) must retain the quotient carry (equivalently `K+d_m` parent precision)
unless it proves an additional physical relation. This is a stopping theorem
for that route, not an impossibility theorem for every coarse observer.

On the exact shells `m=12,...,17`, every declared
`K in {1,2,3,4,5,6,8}` is target-blind. The maximal prefix variations for
`K=(1,2,4,8)` are

```text
m=12: ( 65, 87, 135, 575)
m=13: (183,183, 257, 864)
m=14: (249,249, 279,1210)
m=15: (180,180, 545,1769)
m=16: (367,367, 500,2124)
m=17: (889,889,1104,3149).                            (23j)
```

Their square-root-like appearance is **FINITE-EXACT** only. At `m=17`, a
chronological label-blind within-fibre pairing has worst unmatched masses
`(65556,65488,65142,65075)`. Among the toggles `r->r xor 2^j` with
`j in {0,1,floor(m/2),m-2,m-1}`, the best whole-shell values are
`(65240,97812,122866,130566)`. The finite cancellation is therefore nonlocal
for these declared obvious operation pairings; it supplies no asymptotic
bound.

## 5. Prize 3: the model fork is load-bearing

Fix, for illustration, a deterministic finite multitape Turing machine,
binary input `n`, exact output for every `n`, a finite uniform program, no
advice/oracle/free preprocessing, and step runtime `T_M(n)`. Even then four
natural statements are not equivalent:

```text
(i)   every correct M has T_M not in o(n);
(ii)  every correct M has some c_M>0 with T_M(n)>=c_M n eventually;
(iii) every correct M has limsup_(n->infinity) T_M(n)/n=infinity;
(iv)  every correct M has max_(2^m<=n<2^(m+1))T_M(n)=Omega(2^m).  (24)
```

Statement (iii), the displayed 2019 formulation, says that no correct
`O(n)` algorithm exists; it is not the usual `Omega(n)` lower bound. At the
level of the inner runtime predicates, `T(n)=Theta(n)` meets the bound in
(ii) but refutes the predicate in (iii). Also, `T(n)=ceil(sqrt(n))` except
`T(2^m)=2^m` meets (i)'s non-`o(n)` predicate and (iv)'s shell predicate but
not (ii) or (iii); restricting those linear spikes to `n=2^(2^k)` destroys
(iv) while retaining (i); making the sparse spikes quadratic yields (iii)
without (ii) or (iv). Word-RAM, randomized, circuit, and nonuniform versions
require separate declarations. Thus a Prize-3 claim must name its model and
quantifiers before importing any observer statistic.

### 5.1 Linear history can be easy

Define the explicit infinite sequence

```text
x_n=(n mod 2) XOR 1_{n is a positive power of 2}.      (25)
```

It is balanced with prefix discrepancy `O(log N)`, because it changes the
alternating word only at powers of two. It is computable from the binary
digits of `n` in `O(log n)` bit time, and in constant standard word-RAM
operations.

Nevertheless, for every power `P>=8`, put `h=P/2-1`. The contexts of length
`h` immediately before positions `P` and `2P-2` are equal. Their intervals
contain no powers of two and begin with the same parity, so both are the same
alternating word. But

```text
x_P=1,                         x_(2P-2)=0.             (26)
```

Thus (25) has right-special histories whose length is linear in the follower
position, while remaining balanced and easy to evaluate. This universally
refutes

```text
large suffix-history depth  =>  linear query complexity.                (27)
```

For Rule 30 itself, all `2^h` contexts are observed and mixed through `h=12`
in the first `2^18` bits; a named length-`33` context has both followers, and
depth `34` gives exact lookup only because all observed contexts are unique.
These are finite table-memorization facts, not Prize-3 hardness.

## 6. Remaining proof obligations

The reductions generate three sharply typed next tasks.

1. **Prize 1:** produce nonzero inverse-boundary bits uniformly in `(N,T)`,
   or prove an all-scale wrap-length variation theorem retaining the current
   which scalar `r_m` forgets.
2. **Prize 2:** bound the target-blind current (23f) for one fixed `K>=2`, or
   construct a label-blind physical pairing that proves the same bound. Any
   scale recursion must retain the quotient carry in (23h); do not extrapolate
   (23j).
3. **Prize 3:** choose one model in (24), then seek a genuine locality,
   communication, perturbation, or crossing-sequence obstruction. State
   counts and history depth alone carry too little typed information.

The official prize page is an active listing accepting submissions. Per
MISTAKE-403, the repository's dated `OPEN` status is an inference from that
literal evidence, not a quotation that the page itself says "open."
