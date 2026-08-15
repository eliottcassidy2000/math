---
id: THM-3463
title: "Rule 30 Mealy sections, collision currents, period-lift restriction, and the universal-input complexity boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  A three-state
  invertible Mealy codec restores the moving center observer exactly and
  yields a suffix-parity plus intersection-current compiler.  A dual spacetime
  chart expresses every center departure from the Rule-150 baseline as an
  adjacent-11 collision current over a Green kernel.  The edge-period lift
  word has no 111 and obeys a
  sharper all-width period bound; innovation depths give an exact spatial
  Walsh cube.  Separately, a structured seed-neighborhood promise has exact
  decision-tree complexity n+1 and needs at least n binary AND gates.  These
  statements imply none of the three Rule 30 prizes.
source: root-rule30-diagonal-frontier-20260815
audit: >
  independent Mealy/section, Green-current, dyadic-baseline, period-cube,
  arrival-atlas, Boolean-derivative, decision-tree, AND-complexity,
  2-kernel, cost-model, scope, hash, and ordinary/optimized/stored replay
  audit; documentation clean
depends_on:
  - THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
related:
  - THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary
  - THM-3459-rule30-ternary-intersection-factorial-truth-lift-and-keller-boundaries
  - THM-2005-support-dirichlet-automatic-tournament-atlas
  - THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar
external:
  - Stephen Wolfram, "Announcing the Rule 30 Prizes", https://writings.stephenwolfram.com/2019/10/announcing-the-rule-30-prizes/ (2019; CITED for the problem statements and local rule only)
  - Wolfram Rule 30 Prizes, https://rule30prize.org/ (CURRENT OFFICIAL LISTING checked 2026-08-15; active listing and submission status only)
script: 04-computation/rule30_mealy_current_complexity_thm3463.py
output: 05-knowledge/results/rule30_mealy_current_complexity_thm3463.out
script_sha256: 73e5ab576dad91ae234337554dd5065cbbb922b0572ec31cc711318aa708da04
output_sha256: b69d2a3b3ef114bf1ab81cf574cfbac79309aaa656c2b26e13e651d50d8a64ae
hash_basis: raw bytes
---

# THM-3463 -- Rule 30 Mealy sections, collision currents, period-lift restriction, and the universal-input complexity boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

This theorem gives exact new structure on all three sides of the moving-
observer boundary in THM-3458.  It does **not** prove center nonperiodicity,
center density `1/2`, or a lower bound for random access to the distinguished
single-seed bit.

## 1. Conventions and the three-state edge automaton

Let `a_t(j)` be Rule 30 from the single seed,

```text
a_0(0)=1,       a_0(j)=0 for j!=0,
f(l,c,r)=l xor c xor r xor cr.                         (1)
```

As in THM-3458, read row `t` inward from its right edge:

```text
b_k(t)=a_t(t-k),
R_t=sum_(k>=0) b_k(t) 2^k,
R_(t+1)=Phi(R_t),
Phi(x)=x xor ((2x) or (4x)).                           (2)
```

Use the standing zero-extension convention `b_k(t)=0` for `k<0`.

Read one-sided binary words from low to high index.  For a Mealy state `g`,
write

```text
g=(g_0,g_1) sigma^e                                  (3)
```

when its root output is input xor `e` and its two sections are `g_0,g_1`.
Then the packed Rule 30 map is the state `A` in the three-state automaton

```text
A=(A,B),
B=(C,B) sigma,
C=(A,B) sigma.                                        (4)
```

Indeed, before input bit `x_k` is read, the transducer need only remember
`(x_(k-1),x_(k-2))`.  The states `10` and `11` have identical futures and
merge to `B`; `00` is `A` and `01` is `C`.  The output is

```text
x_k xor (x_(k-1) or x_(k-2)),                         (5)
```

which is exactly (2).  In particular,

```text
A(1w)=1 B(w),
R_n=A^n(10^infinity)=1 B^n(0^infinity).               (6)
```

The distinguished center therefore has the exact section formula

```text
c_n=a_n(0)
   = activity((B^n)|_(0^(n-1)))       (n>=1).          (7)
```

This faithfully restores THM-3458's moving observer inside
one fixed finite-state edge action.  The price is that the required section
depth grows with `n`.

## 2. The suffix-parity/intersection compiler

Fix `n>=1`.  Let `W_k` be the written product word representing
`(B^n)|_(0^k)`, and let

```text
p_k(i)=1  iff the i-th factor of W_k has root activity 1,
(Sp)(i)=xor_(j>i) p(j).                                (8)
```

Thus `S` is strict suffix parity.  In a product, the input reaching factor
`i` is `(Sp_k)(i)`.  From (4):

- the next factor is `B` iff that incoming bit is `1`;
- its next activity is the OR of that incoming bit and the indicator that
  the old factor was `B`.

Every section of a `B` factor has activity one, so `p_0=p_1=1^n`.  Eliminating
the `B` indicator gives the exact recurrence

```text
p_0=p_1=1^n,
p_(k+1)=S p_(k-1) or S p_k,
c_n=xor_i p_(n-1)(i).                                 (9)
```

Two suffix scans and one OR per level give an exact `O(n^2)` bit-operation,
`O(n)` memory compiler.  This is an upper-bound representation, not a
complexity lower bound.  Treating a whole suffix scan or a precomputed jump
as one primitive changes the cost model.

### 2.1 The missing scalar is an intersection correlation

Over `F_2`, OR has the typed decomposition

```text
u or v = u+v+uv.                                      (10)
```

For binary Hasse moments

```text
M_r(w)=xor_i binom(i,r) w_i,                           (11)
```

the hockey-stick identity gives

```text
M_r(Sw)=M_(r+1)(w).                                   (12)
```

But already at the scalar level,

```text
M_0(Sp or Sq)
 =M_1(p)+M_1(q)+M_0((Sp)(Sq)).                        (13)
```

The last term is not determined by the marginal moments.  An actual hostile
already occurs on the width-three section orbit.  In increasing-index order,

```text
S(111)=010,
S(010)=100.                                           (14)
```

The pairs `(111,111)` and `(111,010)` have the same binary `M_0,M_1`
marginals, and all four suffix images in question have parity one.  Their
suffix-image intersection parities are nevertheless `1` and `0`, so the
next OR parities differ.  These low Hasse marginals are therefore not a
closed center state; a labelled intersection mask is the missing sidecar.

### 2.2 Dyadic times are pure intersection current in this chart

Linearize (9) by deleting the product term in (10):

```text
u_0=u_1=1^n,
u_(k+1)=S(u_k+u_(k-1)).                               (15)
```

Write `u_k=H_k(S)1^n`, where

```text
H_0=H_1=1,
H_(k+1)(s)=s(H_k(s)+H_(k-1)(s)).                      (16)
```

For every vector `v` and `r>=0`,

```text
xor_i (S^r v)_i
 =xor_j binom(j,r) v_j.                               (17)
```

For `v=1^n`, this is `binom(n,r+1) mod 2`.  If `n=2^m>=2`, Lucas' theorem
and `deg H_(n-1)<=n-3` (with `n=2` direct) give

```text
xor_i u_(n-1)(i)=0.                                   (18)
```

Thus every dyadic center `c_(2^m)` with `m>=1` in the section chart is carried entirely
by propagated intersection masks from (10).  The linear Frobenius/Lucas
backbone contributes nothing.

## 3. The dual spacetime collision current

Put

```text
B_t(x)=sum_(k>=0) b_k(t)x^k,
L(x)=1+x+x^2,
d_k(t)=b_(k-1)(t)b_(k-2)(t),
D_t(x)=sum_k d_k(t)x^k
      =x^2 sum_i b_i(t)b_(i+1)(t)x^i.                 (19)
```

Equation (2) is equivalent over `F_2[x]` to

```text
B_(t+1)=L B_t+D_t,
B_t=L^t+sum_(s<t) L^(t-1-s)D_s.                       (20)
```

The central Rule-150 coefficient is always one:

```text
[x^t]L^t=1.                                           (21)
```

For even `t=2r`, Frobenius reduces (21) to the coefficient at `r`; for odd
`t=2r+1`, only the middle `x` term of `L L(x^2)^r` contributes.  Induction
starts at `t=0`.

Let

```text
K(n,q)=[x^q](1+x+x^2)^n.                              (22)
```

Taking the center coefficient in (20) gives the exact collision-current law

```text
c_t
 =1+sum_(s<t) sum_k d_k(s) K(t-1-s,t-k)   in F_2.     (23)
```

Every departure from Rule 150's constant-one center is therefore the parity
of actual adjacent-`11` events transported with their spacetime locations
through a ternary Green kernel.  Collision counts alone lose the kernel phase
and do not determine (23).

Frobenius gives the digit recursion

```text
K(2m,2r)     =K(m,r),
K(2m,2r+1)   =0,
K(2m+1,2r)   =K(m,r)+K(m,r-1),
K(2m+1,2r+1) =K(m,r).                                 (24)
```

Equivalently,

```text
L^n=product_(i:n_i=1)(1+x^(2^i)+x^(2^(i+1))).         (25)
```

The carry in choosing the three terms supplies a two-carry-state LSD-first
digit compiler for each Green weight.  A single `K(n,q)` is therefore
computable in `O(log(n+q))` Boolean steps.  The unresolved content is not the
linear transport; it is producing or compressing the nonlinear current
`(D_s)`.

For `h=2^m` and `0<=r<h`, (20) and
`L^h=1+x^h+x^(2h)` sharpen this to

```text
c_(h+r)+c_r
 =[x^(h+r)] sum_(j=0)^(h-1)L^(h-1-j)D_(r+j).          (26)
```

Thus every dyadic self-similarity defect is exactly the collision current
created during the intervening block.  Equations (18) and (26) are dual
charts of the same obstruction: sideways suffix intersections versus local
spacetime adjacent-`11` events.

## 4. The edge period tower cannot double three times in a row

Let `P_w` be the exact seed period of `Phi mod 2^w`, and define the lift bit

```text
epsilon_w=bit_w(R_(P_w)),
P_(w+1)=2^(epsilon_w)P_w.                              (27)
```

Directly, `epsilon_1=1`.  For `w>=2`, the bit recurrence in (2) gives

```text
epsilon_w
 =xor_(0<=t<P_w)(b_(w-1)(t) or b_(w-2)(t)).           (28)
```

For `w>=2`, `P_w` is even, so (28) is also the parity of the number of phases
where the displayed pair equals `00`.

If `epsilon_j=1`, lower bits repeat after `P_j` while

```text
b_j(t+P_j)=b_j(t)+1                                   (29)
```

for every `t`.  The difference of the two sides is constant in `t` because
both copies receive the same lower-bit cocycle, and it equals `epsilon_j` at
`t=0`.

Now suppose

```text
epsilon_(w-2)=epsilon_(w-1)=1.                        (30)
```

Put `p=P_(w-2)`.  Over the four consecutive `p`-blocks in the `P_w=4p`
cycle, (29) makes `(b_(w-1),b_(w-2))` assume each of `00,01,10,11` once for
every base phase.  Each pair therefore occurs exactly `p` times.  For
`w>=4`, `p` is even, so (28) yields

```text
epsilon_w=0.                                          (31)
```

The only earlier candidate is excluded by `epsilon_1 epsilon_2=10`.
Therefore the infinite lift word contains no `111`.

Since

```text
P_w=2^E_w,
E_w=sum_(j=1)^(w-1)epsilon_j,                         (32)
```

the initial `10` and the forbidden block give the all-width bound

```text
E_w<=ceil(2w/3)-1,
P_w<=2^(ceil(2w/3)-1).                                (33)
```

For every `w>=3`, this strictly improves THM-3458's generic `2^(w-1)`
ceiling; the two bounds agree at `w=1,2`.  It remains an upper bound on
edge-state periods, not a shortcut for the moving center.

A completed run ending at a zero parses uniquely into the codewords

```text
0, 10, 110,                                            (34)
```

with period multipliers `1,2,4`.  This is a faithful ternary run-length
address for the period tower.  It is not a Berggren arithmetic tree or a
three-way symmetry of Rule 30.

## 5. Exact spatial Walsh cube, marked-origin boundary

Let

```text
I_K={k<=K:epsilon_k=1},
m=|I_K|=log_2 P_(K+1),
T_h(k)=a_(h+k)(h)=b_k(h+k).                            (35)
```

Then the arrival-trace map

```text
Gamma_K: Z/P_(K+1)Z -> F_2^(I_K),
h |-> (T_h(k))_(k in I_K)                             (36)
```

is a bijection.

To prove this, order the innovation depths `k_1<...<k_m`.  The shift
`h |-> h+P_(k_r)` fixes every earlier coordinate because its bit period
divides `P_(k_r)`, while (29) complements coordinate `k_r`.  The phase group
has size `P_(K+1)=2^m`; successive innovation shifts therefore give a
triangular binary coordinate system, proving (36).

Consequently every innovation word occurs exactly once across one spatial
phase period, and every nonempty joint Walsh correlation of these coordinates
vanishes.  This is an exact balance/independence theorem, but it is spatial
and samples only innovation depths.  The prize center is the one marked
codeword `h=0`, and the omitted noninnovation coordinates are nontrivial
Boolean readouts of that phase.  THM-3458 proves that the depth-five bit word
has weight `3/8`; translating phase by five leaves its weight unchanged, so
the corresponding arrival readout is not balanced.  Uniformity of the
innovation address neither balances those readouts nor makes the marked
address temporally typical.

## 6. A linear lower bound on a seed-neighborhood promise

This section concerns arbitrary initial data, not just the single seed.  In
the Boolean ring

```text
F_2[x_i:i in Z]/(x_i^2-x_i),                           (37)
```

let `A_t(j)` be the time-`t` cell polynomial.  Define the finite Boolean
derivative `Delta_i F=F|_(x_i=0)+F|_(x_i=1)`.

Left permutivity gives the exact extreme-left form

```text
A_t(j)=x_(j-t)+H_t(j),                                (38)
```

where `H_t(j)` is independent of `x_(j-t)`.  Along the unique extreme-right
path, differentiating Rule 30 in its right input gives

```text
Delta_(x_(j+n)) A_n(j)
 =product_(s=0)^(n-1)(1+A_s(j+n-1-s)).                (39)
```

Set

```text
m_s=j+n-1-2s,
M_n(j)={m_0,...,m_(n-1)},
T_n(j)=M_n(j) union {j+n}.                             (40)
```

Factor `s` of (39) has extreme-left marker `m_s`, and by (38) its only
monomial containing that marker is the singleton `x_(m_s)`.  Peel markers
from smallest to largest.  To obtain all of `M_n(j)`, every factor is forced
to contribute its own singleton.  Hence the product in (39) has exactly one
monomial divisible by `M_n(j)`, namely `product_(i in M_n(j))x_i`.  Therefore

```text
Delta_(T_n(j)) A_n(j)=1,                              (41)
```

and `T_n(j)` is a maximal ANF monomial of degree `n+1`.

For the center,

```text
T_n={n,n-1,n-3,...,1-n}.                              (42)
```

Restrict the initial row to the affine seed-neighborhood cube

```text
delta_0 + sum_(i in T_n) u_i e_i.                     (43)
```

Equation (41) says that the XOR of the output over all `2^(n+1)` vertices is
one.  Thus the restricted function has degree exactly `n+1`.  Any exact
deterministic decision tree of depth below `n+1` would partition the cube
into constant leaves of even size, forcing total XOR zero.  Querying all
variables is an upper bound, so

```text
D_query((43) -> A_n(0))=n+1.                          (44)
```

Similarly, a Boolean circuit in which XOR and NOT are free and binary AND
gates are counted has degree at most one plus the number of AND gates in its
dependency subgraph.  To see this, replace a topologically first AND by a
fresh variable, apply induction to the remaining gates, and substitute back
the degree-two product; the degree increases by at most one.  Hence (43), and
also the universal-input cone function, need at least

```text
n binary AND gates.                                   (45)
```

This is a semantic linear lower bound robust on a structured neighborhood of
the seed.  It is **not** a lower bound for the single point `u=0`.

### 6.1 Why this does not settle Prize 3

- For fixed `n`, a nonuniform zero-input circuit can hardwire `c_n`.
- A meaningful random-access question needs one uniform machine taking `n`,
  its input representation, and explicit charges for advice, preprocessing,
  memory and bit operations.
- Unary input makes an `Omega(n)` reading bound vacuous; binary or random-
  access input is the relevant model.
- A radius-one causal circuit for **all** initial rows has depth `n`, but Rule
  60 has the same causal lower bound and a constant single-seed center.
- Literal iteration uses `n` calls to `Phi`; allowing jump maps `Phi^(2^k)`
  or hardwired advice destroys that syntactic count.

Thus (44)--(45) identify a robust-computation obstruction, while the reduction
from the binary time index to that promise problem remains missing.

There is also a source-level ambiguity in the 2019 formulation: its prose
gestures at excluding a sublinear shortcut, while its displayed `limsup`
predicate would exclude even `O(n)` algorithms.  Here “Prize 3” is only an
informal label.  No fixed-seed Prize-3 lower-bound statement is made until a
machine, index encoding, advice policy and bit-cost model are fixed.

## 7. A large finite-exact 2-kernel hostile

Let

```text
K_(e,r)=(c_(r+2^e n))_(n>=0).                         (46)
```

Exact enumeration gives three finite statements:

1. Among the `32,767` profiles with `0<=e<=14`, `0<=r<2^e`, the first
   sixteen bits take `25,830` distinct values.
2. After deleting the `n=0` coordinate, the fifteen positive-index bits
   `n=1,...,15` still take `20,794` distinct values in the same universe.
3. Among `0<=e<=15`, every one of the `4,096` binary words of length twelve
   occurs as the first twelve bits of some `K_(e,r)`.

Two independent row recurrences agree through time `131071`; the packed
recurrence supplies the longer hashed/enumerated prefix.  Therefore
any canonical LSD-first automaton whose states are the 2-kernel subsequences
would need at least `25,830` states.  For an `s`-state MSD-first binary DFAO on
standard positive binary representations, appending a fixed binary suffix
induces one of at most `2^s` Boolean output functions on its state set.  The
positive-index count `20,794>2^14` therefore forces `s>=15` without any
leading-zero convention.

The length-twelve statement also rules out any universal forbidden word or
nonzero Boolean relation on a twelve-bit initial window of every dyadic
kernel subsequence.  These are exact finite lower bounds and hostiles to a
small automatic explanation.  They do not prove that the 2-kernel is
infinite or that the center is nonautomatic.

## 8. Preservation, loss, and actionable frontier

| source -> target | preserves | destroys without sidecar | needed sidecar |
|---|---|---|---|
| packed edge -> three-state sections | exact seed center | fixed-depth state alone | section depth `n` and factor order |
| sections -> Hasse marginals | linear suffix transport | OR-intersection parity | labelled intersection masks |
| spacetime row -> collision count | number of adjacent `11` events | Green phase and cancellation | event locations and `K` weights |
| edge phase -> Walsh atlas | exact innovation-depth spatial joint law | marked origin and noninnovation readouts | address `h=0` plus all-depth Boolean readout functions |
| arbitrary cone -> seed cube | maximal ANF obstruction | the unperturbed single point | promise variables / reduction from `n` |

Three concrete next probes now replace vague randomness heuristics:

1. Find an infinite unique-channel family in (23), forcing uncancelled center
   changes and thereby approaching nonperiodicity.
2. Compile/control the noninnovation-coordinate functions on the innovation
   cube, then compare the marked address `h=0` with their spatial averages.
   Both steps are required before this can approach density `1/2`.
3. Either compress the intersection currents in (9)/(23), or prove a lower
   bound for them in a precise uniform binary-index machine model.

The empirically stronger formula `deg A_n(0)=2n-1` holds through the checked
range `2<=n<=8` (with degree `2` at `n=1`), but no all-`n` coefficient
proof is used here.  It remains finite-exact evidence only.

## 9. Exact verification and scope

Reproduce with

```bash
python3 04-computation/rule30_mealy_current_complexity_thm3463.py
python3 -O 04-computation/rule30_mealy_current_complexity_thm3463.py
```

The companion checks a direct local-rule/packed-row control through time
`128`, all input words through Mealy width `11`, and direct center agreement
through `n=512`.  It checks Duhamel reconstruction through `t=80`, the Green
digit identities through `m=32`, all `127` dyadic-current pairs with `h<=64`,
period lifts through width `36`, arrival atlases through depth `30`, exact
ANFs through `n=8`, and the stated 2-kernel universes.  The universal proofs are
the arguments above; the bounded checks audit their implementation and the
finite-exact claims.

The official prize page was checked on `2026-08-15`: it actively lists all
three questions and accepts submissions.  On that dated evidence the repo
treats them as open.  This theorem supplies no prize solution, no literature
novelty claim, and no transfer to LRC, the factorial conjecture, or the
Jacobian conjecture.
