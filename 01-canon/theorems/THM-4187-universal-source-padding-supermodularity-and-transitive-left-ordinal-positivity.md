---
id: THM-4187
title: "Universal-source padding supermodularity and transitive-left ordinal positivity"
status: >
  PROVED all-order universal-source padding supermodularity with its unique
  equality case, nonnegative transitive-left normalized local interaction,
  nonnegative transitive-left ordinal remainder with exact zero boundary,
  strict (OS+) for every transitive left factor against every no-sink right
  factor, the converse universal-sink/R_- dual, a sharp terminal-C3
  quantitative strengthening, and an exact two-parameter cycle-first
  negative family + VERIFIED-EXACT + INDEPENDENTLY AUDITED. General (OS+),
  arbitrary-left local positivity, the no-sink/no-source gate law, and the
  order-eleven asymmetric bank remain OPEN.
source: tournament-theta-sign-20260826
depends_on:
  - THM-4177-root-split-johnson-current-and-source-sink-boundary-census
  - THM-4181-ordinal-sum-capacity-transfer-and-parity-component-exchange
  - THM-4184-path-cover-parity-ordinal-cocycle-and-lollipop-positivity
related:
  - THM-4133-strong-cyclic-substitution-johnson-centrality-counterexample
  - THM-4169-prime-parent-one-vertex-augmentation-and-quartic-johnson-transfer
script: 04-computation/tournament_source_padding_supermodularity_thm4187.py
output: 05-knowledge/results/tournament_source_padding_supermodularity_thm4187.out
independent_audit_script: 04-computation/tournament_source_padding_supermodularity_thm4187_independent_audit.cpp
independent_audit_output: 05-knowledge/results/tournament_source_padding_supermodularity_thm4187_independent_audit.out
script_sha256: 541cb266af5f0cb36238f4b5267466563b492d4adaea15fc514424a70a9cfcc7
output_sha256: 2e1a67a0360d2415ffe80ba4048f711e715715307513608d99acf30f4813c01f
independent_audit_script_sha256: 0128797bed966f8f93b7163abe20302fddac156294c37a046fcf67fbba022dca
independent_audit_output_sha256: 1a812b853d7b6e8c134362c5b28e3ad2059d693c2dfef159cefab39b3f134947
hash_basis: raw LF bytes
primary_audit: >
  PASS. The inherited class-representative engine checks the literal padding
  identity on 76 factors, the exact nonnegative decomposition on all 5,776
  ordered factor-class pairs through order six, 46,208 transitive-left local
  interactions, both equality ledgers, the dual, the C3 sharpening, and the
  negative formulas through transitive-tail order 64.
independent_audit: >
  ACCEPT. A separate C++ literal-permutation referee constructs all 75
  labelled factors through order four without gentourng or subset DP, directly
  rebuilds 505 ordinal children, checks 5,625 labelled pair decompositions and
  45,000 transitive-left interactions, and verifies the dual, C3 formulas, and
  cycle-first hostiles. Clang O0/O3 and ASan/UBSan streams byte-match.
---

# THM-4187 -- universal-source padding supermodularity and transitive-left ordinal positivity

**PROVED all-order algebra and sign laws + VERIFIED-EXACT + INDEPENDENTLY
AUDITED.**

THM-4184 made the normalized ordinal remainder a cocycle, but its local
interaction still had no all-order sign.  The missing sign is now proved when
the first factor is transitive.  The mechanism is not an estimate on the
rank-two cross tensor.  It is the exact nonnegative padding defect of
THM-4177, applied before the ordinal cut is scalarized.

The central statement is that adjoining one universal source becomes no less
valuable after an arbitrary ordinal suffix is attached.  Its normalized
remainder is therefore monotone.  Iterating the source proves

```text
R_+(P_a,C)>=0
```

for every transitive tournament `P_a` and every nonempty tournament `C`, and
strictly proves `(OS+)` whenever `C` has no sink.  This is an all-order factor
family, not a finite census.

Factor order remains structural.  In the reverse direction
`R_+(C3,P_n)<0` for every `n>=1`, and the local interactions
`Theta(C3,P_b,P_c)` form an exact two-parameter negative family beginning at
THM-4184's hostile `Theta(C3,1,1)=-216`.

## 1. Conventions

All tournaments below are finite and nonempty.  Use THM-4177/4181/4184's
Hamilton count `H`, exposure capacity `c`, capacity mass `W`, gates `G_+` and
`G_-`, ordinal sum `triangleright`, and remainders

```text
R_eta(A,B)
 =G_eta(c(A triangleright B))
  -H(B)^2 G_eta(c(A))-H(A)^2 G_eta(c(B)),       eta in {+,-}.   (1)
```

For either sign, define the local interaction by

```text
Theta_eta(A,B,C)
 =R_eta(A triangleright B,C)-H(A)^2R_eta(B,C)
 =R_eta(A,B triangleright C)-H(C)^2R_eta(A,B).            (2)
```

The second equality is the normalized-coboundary identity of THM-4184,
which applies verbatim to either quadratic gate.  Write `Theta=Theta_+` and

```text
rho(A,B)=R_+(A,B)/(H(A)^2H(B)^2).                         (3)
```

Write `1=P_1` for the singleton tournament and `P_a` for the transitive
tournament of order `a>=1`.  The notation `P_0` is used only as an empty
prefix in telescoping sums; it is not inserted into a theorem quantifier.

For a tournament `Q`, put

```text
Delta(Q)=G_+(1 triangleright Q)-G_+(Q)=R_+(1,Q).          (4)
```

Let `x` be the new universal source.  For `v in Q`, set

```text
a_v(Q)=c_xv(1 triangleright Q),
F_v(Q)=sum_(u->v in Q) a_u(Q)-r_v(c(Q)).                  (5)
```

THM-4177's incoming-fan injection proves `F_v(Q)>=0`; its exact padding
identity is

```text
Delta(Q)
 =sum_v a_v(Q)[W(Q)-d_v(Q)]+4sum_v a_v(Q)F_v(Q).         (6)
```

Every `a_v` is positive, and every edge capacity is positive: the length-one
odd-path atom contributes `2H(Q-{u,v})>0`.

## 2. Universal-source padding is supermodular across an ordinal cut

> **Theorem 1 (exact source-padding supermodularity).** For all nonempty
> tournaments `B,C`,
>
> ```text
> Delta(B triangleright C)>=H(C)^2 Delta(B),              (7)
> ```
>
> with equality if and only if `B=C=1`.  Equivalently,
>
> ```text
> Theta(1,B,C)>=0,                                       (8)
> rho(1,B triangleright C)>=rho(1,B),                    (9)
> ```
>
> and equality in either statement occurs exactly at `(B,C)=(1,1)`.

### Proof and exact positive decomposition

Write

```text
h=H(B),          s=H(C),          T=B triangleright C.
```

For `i in B`, let `a_i=a_i(B)`.  In `T`, let `p_i` be the total cross
capacity incident with `i`, and put `Z=sum_i p_i`.  For `j in C`, let
`beta_j=a_j(T)` be its root capacity in `1 triangleright T`.

THM-4181's internal capacity transfer applied to

```text
1 triangleright T=(1 triangleright B) triangleright C
```

gives, on the `B` block,

```text
a_i(T)=s a_i,                                            (10)
c(T)|_B=s c(B).                                          (11)
```

No cross arc enters `B`.  Therefore

```text
r_i(T)=s r_i(B),
sum_(u->i in T)a_u(T)=s sum_(u->i in B)a_u(B),
F_i(T)=sF_i(B).                                          (12)
```

The total and incident masses are

```text
W(T)=sW(B)+hW(C)+Z,
d_i(T)=s d_i(B)+p_i.                                     (13)
```

Substitute `(10)--(13)` into the `B`-vertex terms of `(6)`.  They equal
`s^2Delta(B)` plus a residual.  Keeping the `C`-vertex terms literal gives
the exact identity

```text
Delta(T)-s^2Delta(B)
 =s sum_(i in B) a_i[hW(C)+Z-p_i]
  +sum_(j in C) beta_j[W(T)-d_j(T)]
  +4sum_(j in C) beta_j F_j(T).                          (14)
```

Every term in `(14)` is nonnegative.  The first bracket is the mass of all
internal-`C` and cross edges not incident with `i`; the second is the mass of
all edges not incident with `j`; the last uses the incoming-fan injection.

If `|C|>=2`, then `W(C)>0`, so the first sum is strictly positive.  If
`|C|=1` and `|B|>=2`, every `Z-p_i` contains the positive cross capacity of
another vertex of `B`, so it is again strictly positive.  The only remaining
case is `B=C=1`; directly, `Delta(P_2)=Delta(P_1)=0`.  This proves `(7)` and
its equality case.

Since `R_+(1,Q)=Delta(Q)`, the second form of `(2)` identifies the left side
of `(14)` with `Theta(1,B,C)`.  Division by `h^2s^2` proves `(9)`. QED.

The theorem is stronger than sign: `(14)` identifies the precise new atoms
created by the suffix.  No cancellation or asymptotic estimate is used.

## 3. Transitive-left interaction and `(OS+)`

> **Theorem 2 (transitive-left local interaction).** For every `a>=1` and
> all nonempty `B,C`,
>
> ```text
> Theta(P_a,B,C)>=0.                                    (15)
> ```
>
> Equality holds if and only if `a=1` and `B=C=1`.

### Proof

Because `H(P_k)=1`, repeated source insertion gives

```text
Theta(P_a,B,C)
 =sum_(k=0)^(a-1) Theta(1,P_k triangleright B,C).         (16)
```

Every summand is nonnegative by Theorem 1.  If `a=1`, its equality
classification is exactly that of Theorem 1.  If `a>=2`, the summand with
`k=1` has a middle factor of order at least two and is strict. QED.

> **Theorem 3 (transitive-left ordinal remainder).** For every `a>=1` and
> every nonempty tournament `C`,
>
> ```text
> R_+(P_a,C)>=R_+(1,C)=Delta(C)>=0.                      (17)
> ```
>
> The complete zero boundary is
>
> ```text
> (a,C)=(1,P_1), (1,P_2), (2,P_1).                      (18)
> ```
>
> In particular, if `C` has no sink, then
>
> ```text
> R_+(P_a,C)>0.                                         (19)
> ```

### Proof

Put `q_a(C)=R_+(P_a,C)`.  Equation `(2)` gives

```text
q_(a+1)(C)-q_a(C)=Theta(1,P_a,C)>=0.                    (20)
```

The initial value is `q_1(C)=Delta(C)`.  THM-4177 proves it strictly positive
for `|C|>=3`; direct evaluation gives `Delta(P_1)=Delta(P_2)=0`.  Theorem 1
makes the first increment zero only for `(P_a,C)=(P_1,P_1)` and every later
increment strict.  This proves `(17)--(18)`.

A nonempty no-sink tournament has order at least three, so its initial value
is already strict.  This proves `(19)`, hence `(OS+)` for every transitive
left factor at every order. QED.

### Transitive-prefix monotonicity and a source-free reduction

The same local theorem applies before an arbitrary left core.  For all
`a>=1` and nonempty `B,C`, equation `(2)` gives

```text
R_+(P_a triangleright B,C)
 =R_+(B,C)+Theta(P_a,B,C)>=R_+(B,C),                    (20a)
```

with equality exactly at `(a,B,C)=(1,1,1)`.  In particular, when `C` has no
sink the inequality is strict.  If general `(OS+)` fails, choose a
counterexample with the left factor of minimum order.  That factor cannot
have a source: deleting a universal source writes it as `1 triangleright B`,
and `(20a)` produces a smaller counterexample.  Equivalently, its first
strong component is not a singleton and therefore has order at least three.

This reduces the open problem to **source-free left factors**; it does not
reduce it to strong left factors.  The cycle-first hostiles in Section 5 show
that source-free left factors can already carry negative curvature when the
right factor has a sink.

### Converse dual

Converse preserves `H`, `D`, and the symmetric capacity coordinates, reverses
the sign of `C`, and satisfies

```text
(A triangleright B)^op=B^op triangleright A^op,
G_+(X^op)=G_-(X).                                       (21)
```

Thus, with

```text
Delta_-(X)=G_-(X triangleright 1)-G_-(X),                (22)
```

Theorem 1 dualizes to

```text
Delta_-(B triangleright C)>=H(B)^2Delta_-(C),            (23)
```

again with equality only for `B=C=1`.  Moreover

```text
R_-(C,P_a)=R_+(P_a,C^op)>=0,                             (24)
```

with the three zero cases dual to `(18)`, and it is strict when `C` has no
source.  The local dual is

```text
Theta_-(A,B,P_a)=Theta_+(P_a,B^op,A^op)>=0,              (25)
```

with equality only at `a=1`, `A=B=1`.

## 4. A sharp terminal-`C3` strengthening

Let `C3` be the directed triangle.  It has

```text
H(C3)=3,       G_+(C3)=0,
U_z(C3)=V_z(C3)=(1,2)                  for every z.       (26)
```

For a tournament `X`, write `h=H(X)`, `w=W(X)`, and let `d_i,o_i` be its
unsigned and outgoing capacity masses.  Define

```text
a_i=c_xi(1 triangleright X)=2V_i^1(X),
t_i=c_iz(X triangleright C3)=2U_i^0(X)+4U_i^1(X),        (27)
```

where `t_i` is independent of the cycle vertex `z`.  Put `A=sum_i a_i`.
THM-4184's parity balance gives

```text
A=w+2h,
sum_i t_i=3w+4h=3A-2h.                                  (28)
```

The `C3` specialization of THM-4181's block-gate formula is

```text
R_+(X,C3)
 =15sum_i t_i^2+9sum_i t_i(w-d_i+4o_i)
  -27w^2-108hw-120h^2.                                  (29)
```

Under universal-source padding `Y=1 triangleright X`, the old `t_i` stay
unchanged, the new source has `t_x=3A`, and

```text
w(Y)=w+A,
d_i(Y)=d_i+a_i,       o_i(Y)=o_i,
d_x(Y)=o_x(Y)=A.                                        (30)
```

Subtract `(29)` for `Y` and `X`.  Exact cancellation gives

```text
Theta(1,X,C3)
 =9[8A(3A-h)-sum_i a_i t_i].                            (31)
```

Since all coordinates are nonnegative,

```text
sum_i a_i t_i<=(sum_i a_i)(sum_i t_i)=A(3A-2h).
```

Therefore

```text
Theta(1,X,C3)>=27A(7A-2h)>=648h^2.                      (32)
```

The last inequality uses `A=w+2h>=2h`.  Equality forces `w=0`, hence
`X=1`; the singleton attains it.  Telescoping now proves the sharp all-order
statements

```text
Theta(P_a,X,C3)>=648a H(X)^2,                            (33)

R_+(P_a,X triangleright C3)
 =9R_+(P_a,X)+Theta(P_a,X,C3)
 >=648a H(X)^2.                                         (34)
```

Equality in either bound occurs exactly at `(a,X)=(1,1)`.  Equation `(34)`
extends THM-4184's transitive-prefix lollipops to an arbitrary middle
tournament while retaining a uniform quantitative margin.

## 5. The sharp cycle-first negative family

The transitive-left hypothesis cannot be replaced by a factor-order-blind
condition.  Put

```text
q_n=R_+(C3,P_n),          n>=1.                           (35)
```

In `C3 triangleright P_n`, the cycle capacities are `2`, the tail capacities
are `3c(P_n)`, and the capacity from any cycle vertex to tail vertex `j` is

```text
s_0=4,
s_j=3*2^j,                  1<=j<n.                      (36)
```

Indeed `V_0(P_n)=(0,1)`, while
`V_j(P_n)=(2^(j-1),2^(j-1))` for `j>=1`.  Consequently

```text
S=sum_j s_j=3*2^n-2,
Q=sum_j s_j^2=3*4^n+4.                                  (37)
```

For `n>=2`, using THM-4184's transitive degree formulas,

```text
sum_j s_j[W(P_n)-d_j(P_n)-4r_j(P_n)]
 =-2*4^n-(3n+14)2^(n-1)+24.                             (38)
```

Writing the left side of `(38)` as `L`, the exact block-gate formula reads

```text
q_n=18W(P_n)+30S+9L+9S^2-21Q.                           (38a)
```

Substitution of `(37)--(38)` into `(38a)` cancels the `4^n` terms and gives

```text
q_1=-72,
q_n=72-(27n+126)2^(n-1),              n>=2.              (39)
```

Thus `q_n<0` for every `n>=1`.  Since
`P_b triangleright P_c=P_(b+c)`, the second form of `(2)` gives

```text
Theta(C3,P_b,P_c)=q_(b+c)-q_b.                           (40)
```

Explicitly,

```text
Theta(C3,P_1,P_c)=144-(27c+153)2^c<0,                    (41)

Theta(C3,P_b,P_c)
 =-2^(b-1)[(27(b+c)+126)2^c-(27b+126)]<0,    b>=2.       (42)
```

The first value is

```text
Theta(C3,P_1,P_1)=-216,                                  (43)
```

exactly THM-4184's minimal hostile.  Hence that witness is the first member
of a two-parameter structural negative family, not an isolated accident.
By converse, `R_-(P_n,C3)=q_n<0` as well.

This does **not** refute `(OS+)`: every `P_n` has a sink.  It proves instead
that the no-sink condition and factor order are load-bearing, and that the
normalized cocycle has genuine signed curvature away from the transitive-left
cone.

## 6. Exact audits and type firewall

The primary certificate uses one `gentourng` representative of every
tournament class through order six:

```text
class counts 1,1,2,4,12,56,                 total 76.
```

It checks `(6)` on all `76` factors and `(14)` on all
`76^2=5,776` ordered factor-class presentations.  There is one zero row,
`(1,1)`, and the smallest positive row is `16`.  It also checks `46,208`
transitive-left local interactions, `608` remainders, the exact zero ledgers,
the converse dual, `(29)--(34)` on all `76` factors, `(39)` through `n=64`,
and `144` direct rows of `(41)--(42)`.  The ordered-pair stream digest is

```text
6a9d450d416f0b17606dbb703c0744fcb4882f03f30d5548fe901143ac691132.
```

The independent C++ referee constructs all `75` labelled tournaments of
orders one through four by literal permutations, not `gentourng` or subset
DP.  It directly rebuilds `505` ordinal children and checks `5,625` labelled
pair decompositions, `45,000` local interactions, the equality rows, dual,
terminal-`C3` identities, and negative family.  The two finite universes are
presentations: neither count is an unlabelled child-class census.  The
all-order theorems follow from `(14)`, telescoping, and the closed forms, not
from these finite checks.

The connection contract is

```text
source:       actual ordinal blocks and the universal-source capacity fan,
map:          attach the suffix before evaluating the padding defect,
preserved:    H, actual c, G_eta, R_eta, and the ordered factor cut,
new positive atoms:
              nonincident cross/internal-C mass and C-block fan defects,
lost by scalar rho alone:
              vertex ownership of those atoms and factor order,
needed sidecar:
              root capacities a_v, incoming mass r_v, and the cut rows p_i,
hostile:      C3 triangleright P_n and Theta(C3,P_b,P_c)<0.
```

Nothing here proves `(OS+)` for a nontransitive left factor, arbitrary-left
local positivity, THM-4177's all-order no-sink/no-source gate law, its strong
residual, the order-eleven asymmetric bank, exact Johnson cosets, or actual
response maximizers.  General `(OS+)` remains **OPEN**.

## 7. Replay

Primary exact audit:

```bash
python3 -B \
  04-computation/tournament_source_padding_supermodularity_thm4187.py
python3 -O -B \
  04-computation/tournament_source_padding_supermodularity_thm4187.py
PYTHONHASHSEED=4187 python3 -B \
  04-computation/tournament_source_padding_supermodularity_thm4187.py
```

Independent literal-permutation audit:

```bash
clang++ -std=c++20 -O0 -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_source_padding_supermodularity_thm4187_independent_audit.cpp \
  -o /tmp/thm4187-independent-O0
/tmp/thm4187-independent-O0

clang++ -std=c++20 -O3 -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_source_padding_supermodularity_thm4187_independent_audit.cpp \
  -o /tmp/thm4187-independent-O3
/tmp/thm4187-independent-O3

clang++ -std=c++20 -O1 -g -fsanitize=address,undefined \
  -fno-omit-frame-pointer -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_source_padding_supermodularity_thm4187_independent_audit.cpp \
  -o /tmp/thm4187-independent-san
/tmp/thm4187-independent-san
```

The three Python streams byte-match the frozen primary output.  Clang O0/O3
and ASan/UBSan byte-match the frozen independent output.  **QED for the
all-order source/sink supermodularity, transitive-factor signs, terminal-C3
sharpening, and cycle-first negative family only.**
