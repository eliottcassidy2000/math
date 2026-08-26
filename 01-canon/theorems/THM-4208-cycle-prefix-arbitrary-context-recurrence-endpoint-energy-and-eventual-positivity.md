---
id: THM-4208
title: "Cycle-prefix arbitrary-context recurrence, endpoint energy, and eventual positivity"
status: >
  PROVED exact fixed-context C-finite form, recurrence, and rational OGF for
  the C3-first transitive-tail contextual defect; exact endpoint-energy
  inequality and its converse dual, both with equality iff transitive, plus
  their exact ordinal-sum defect laws;
  exact coordinatewise rooted-chirality identities and a strictly positive
  unary endpoint correction with its universal-source equality boundary,
  including an exact collapse to oriented fans plus Hamilton endpoints;
  eventual strict positivity and increase for every fixed arbitrary nonempty
  context pair; exact 11-jet response factorization; and an exact C3/right-
  lollipop positive family + FINITE-EXACT A5 lower bound on all 54,937,744
  ordered factor-class pairs through order eight + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. The uniform arbitrary-context threshold and sharp
  n=5 lower bound are subsequently proved in THM-4212. General (OS+), the
  no-sink/no-source gate law, and the order-eleven asymmetric bank remain
  OPEN.
source: codex-tournament-normalized-holonomy-20260826
depends_on:
  - THM-4177-root-split-johnson-current-and-source-sink-boundary-census
  - THM-4181-ordinal-sum-capacity-transfer-and-parity-component-exchange
  - THM-4184-path-cover-parity-ordinal-cocycle-and-lollipop-positivity
  - THM-4187-universal-source-padding-supermodularity-and-transitive-left-ordinal-positivity
  - THM-4193-cycle-first-transitive-tail-crossing-and-transitive-context-positivity
related:
  - THM-4202-vertex-transitive-ordinal-remainder-positivity
  - THM-4212-cycle-prefix-uniform-arbitrary-context-threshold-and-tail-five-lower-bound
script: 04-computation/tournament_cycle_prefix_arbitrary_context_thm4208.py
output: 05-knowledge/results/tournament_cycle_prefix_arbitrary_context_thm4208.out
independent_audit_script: 04-computation/tournament_cycle_prefix_arbitrary_context_thm4208_independent_audit.cpp
independent_audit_output: 05-knowledge/results/tournament_cycle_prefix_arbitrary_context_thm4208_independent_audit.out
finite_census_emit_script: 04-computation/tournament_cycle_prefix_arbitrary_context_order8_emit_thm4208.py
finite_census_script: 04-computation/tournament_cycle_prefix_arbitrary_context_order8_census_thm4208.cpp
finite_census_output: 05-knowledge/results/tournament_cycle_prefix_arbitrary_context_order8_census_thm4208.out
script_sha256: 711d5e36637a45062025df44715c6fa3ce66379b5667f7d151ce0ce09197941b
output_sha256: daa81b52c0c5157d243e521419015b1da28dcead417f9e1629c2fa7aef39f1c7
independent_audit_script_sha256: 5edaa3b825291f11af05a9ac03b7a54a6e50e1c3c7d24cc6a259a5b00f884371
independent_audit_output_sha256: cc9a928e5e579d76d7b527694c7b22a42ef36c6e15d17d7abbb39dad0f6138f8
finite_census_emit_script_sha256: 59631d6dc58b294c2f2cd6bb2c0858db61bf3a4b2c7e534999e098ea30d7d9a4
finite_census_script_sha256: dd4f7ddf7023a2188de9abec7d1ed0d7e56946c19f1a9c7058ea9a3993345a50
finite_census_output_sha256: 87ad2ad762bc67698db556981896e6374d5dd78a9978fbe1a67cc9d8cf5b1453
gentourng_sha256: 89df605922cc574b28688248b7c256d24342cc615f887e89b2d096038970c110
hash_basis: raw LF bytes
primary_audit: >
  PASS. The inherited exact subset/exposed-word engine checks every 1,099
  labelled tournament through order five, 8,792 unary rows, 5,625 arbitrary
  context pairs through factor order four at n=0,...,8, 28,125 recurrence
  windows starting at n=0, 5,625 exact response-jet and endpoint-cocycle
  rows, both endpoint-energy equality boundaries and ordinal laws, 1,099
  literal coordinatewise rooted-chirality/fan-collapse rows and their
  universal-source boundary, 64 lollipop rows, and the rooted-chirality
  hostile. Normal, optimized, and fixed-hash streams byte-match.
independent_audit: >
  ACCEPT. A clean-room C++ referee imports no transfer code. It rebuilds
  capacities from literal odd directed paths, evaluates disjoint-edge gates
  directly, checks the endpoint injection object by object on all 1,099
  labelled tournaments through order five, directly reconstructs 139
  arbitrary-context recurrence windows including 121 beginning at n=0, and
  directly checks the lollipop through a=8. Clang O0/O3 and ASan/UBSan
  streams byte-match.
finite_census_audit: >
  PASS. One gentourng representative of each class through order eight gives
  7,412 factors and 54,937,744 ordered context pairs. An exact signed-128
  evaluator finds no negative response or normalized-bound failure and one
  equality, at the singleton pair. The order-eight sector has 47,334,400
  rows and no equality. This is FINITE-EXACT evidence, not an all-order proof.
---

# THM-4208 -- cycle-prefix arbitrary-context recurrence, endpoint energy, and eventual positivity

**PROVED all-order fixed-context recurrence and eventuality + exact endpoint
energy and response jets + exact lollipop family + FINITE-EXACT order-eight
context census + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-4193 proves a sharp, uniform statement in **transitive** contexts: the
cycle-first prefix `C3 triangleright P_n` becomes positive exactly at tail
length five. This theorem changes the quantifiers. The contexts are now
arbitrary, but fixed; the conclusion is eventual positivity with a threshold
that may depend on both contexts. The mechanism is an exact four-mode
recurrence whose leading coefficient is forced positive by a new endpoint-
energy injection.

This does not duplicate THM-4202. That theorem proves positivity for
vertex-transitive factor pairs and, when only the left factor is vertex-
transitive, gives an exact variance/covariance defect at every right order.
Here `A_n` is non-vertex-transitive for `n>=1`, both contexts are arbitrary,
and the retained rooted response jet records precisely the coordinates
killed by uniformity. Neither theorem proves the uniform arbitrary-context
`n=5` claim inside this file; THM-4212 subsequently proves that claim and the
sharp all-context threshold. General `(OS+)` remains open.

## 1. Conventions and exact scope

All quantified tournaments are finite and nonempty. Use THM-4181/4184's
Hamilton count `H`, exposure capacity `c`, capacity mass `W`, gate `G_+`,
ordinal sum `triangleright`, rooted parity vectors `U_i,V_i`, and remainder

```text
R_+(X,Y)
 =G_+(c(X triangleright Y))
  -H(Y)^2G_+(c(X))-H(X)^2G_+(c(Y)).                       (1)
```

Let `P_n` be the transitive tournament of order `n`; `P_0` is used only as
empty-tail notation. Put

```text
A_n=C3 triangleright P_n,              n>=0,              (2)

F_n(B,C)=R_+(A_n triangleright B,C)-R_+(B,C).             (3)
```

The proved statements are:

1. for every fixed `B,C`, `F_n` is an explicit exponential-polynomial and
   obeys one order-four C-finite recurrence from `n=0`;
2. endpoint defects `Delta_V(X)` and its converse-dual `Delta_U(X)` are
   nonnegative, with equality exactly for transitive tournaments, and obey
   dual weighted ordinal laws;
3. signed rooted coordinates count Hamilton starts/ends, making the unary
   endpoint correction strictly positive with an exact universal-source
   boundary, and collapse exactly to oriented fans plus those endpoint
   counts;
4. the leading two contextual coefficients have strict signs, forcing fixed-
   context eventual positivity and eventual strict increase;
5. `R_+(X,Y)` factors through explicit integer 11-jets;
6. `R_+(C3,P_a triangleright C3)` has one exact positive formula for every
   `a>=1`, with a separate `a=0` boundary;
7. a strong no-sink converse pair proves that unsigned coarse summaries do
   not determine arbitrary-context response.

## 2. Exact finite response jet

For a tournament `X`, abbreviate

```text
h_X=H(X),                 w_X=W(X),
s_X=(w_X/2,w_X/2+h_X).                                  (4)
```

THM-4184 proves that `s_X=sum_i U_i=sum_i V_i`. Let `d_i,o_i,r_i` be
capacity degree, outgoing capacity mass, and incoming capacity mass. Define

```text
L^+(X)=sum_i U_i(w_X-d_i+4o_i),
L^-(X)=sum_i V_i(w_X-d_i-4r_i),

Q_U(X)=sum_i U_i U_i^T,       Q_V(X)=sum_i V_i V_i^T.    (5)
```

> **Theorem 1 (exact response-jet identity).** For every nonempty `X,Y`,
>
> ```text
> R_+(X,Y)
>  =h_Xh_Yw_Xw_Y
>   +2h_Y <s_Y,L^+(X)>+2h_X <s_X,L^-(Y)>
>   +2<s_X,s_Y>^2+2<Q_U(X),Q_V(Y)>_F
>   +6s_Y^TQ_U(X)s_Y-10s_X^TQ_V(Y)s_X.                  (6)
> ```

### Proof

In THM-4181's ordinal block formula, the actual cross capacity is

```text
z_ij=2<U_i(X),V_j(Y)>.
```

Therefore its row sums, column sums, total, and square sum are

```text
p_i=2<U_i,s_Y>,                 q_j=2<s_X,V_j>,
Z=2<s_X,s_Y>,                   Z_2=4<Q_U,Q_V>_F.         (7)
```

Substitute `(7)` into THM-4181 equation (10). The two linear block terms give
the two `L` pairings. The cross-only quadratic term gives the last two lines
of `(6)`. Every internal gate cancels in the remainder `(1)`. QED.

Expanding the quadratic terms makes `(6)` an integer dot product

```text
R_+(X,Y)=lambda(X) dot mu(Y),                            (8)
```

where, writing `s=(s_0,s_1)` and symmetric Gram entries as `Q_00,Q_01,Q_11`,

```text
lambda(X)=(
 h w, 2L^+_0,2L^+_1,2hs_0,2hs_1,
 2s_0^2+6Q_U00,4s_0s_1+12Q_U01,2s_1^2+6Q_U11,
 2Q_U00-10s_0^2,4Q_U01-20s_0s_1,2Q_U11-10s_1^2),

mu(Y)=(
 h w,hs_0,hs_1,L^-_0,L^-_1,s_0^2,s_0s_1,s_1^2,
 Q_V00,Q_V01,Q_V11).                                    (9)
```

Consequently, for any fixed prefix `A`, define

```text
F_A(B,C):=R_+(A triangleright B,C)-R_+(B,C),

F_A(B,C)=[lambda(A triangleright B)-lambda(B)] dot mu(C). (10)
```

This is a sufficient finite response carrier. It is not a claim that any
listed coordinate can be discarded.

## 3. Exact arbitrary-context C-finite law

For a right tournament `X`, write `h=H(X),w=W(X)` and define

```text
ell_e=sum_i V_i^e(w-d_i-4r_i),
M_ef=sum_i V_i^eV_i^f.                                  (11)
```

Put `J_n(X)=R_+(A_n,X)`. For every `n>=1`, let `x=2^n`. Exact ordinal
transfer gives

```text
H(A_n)=3,                 W(A_n)=12x-6,
s_(A_n)=(6x-3,6x),

Q_U(A_n)=((15x^2/2-3,15x^2/2-3),
          (15x^2/2-3,15x^2/2+6)),

L^+(A_n)=(117x^2-(9n+72)x,
          117x^2-(9n+54)x-18).                          (12)
```

For example, each cycle vertex has rooted start state `(3x/2,3x/2)`. The
path vertices have states

```text
(3*2^(n-j-2),3*2^(n-j-2)), j<n-1;       (0,3), j=n-1.
```

The cycle-to-path capacities are `4` at the first path vertex and `3*2^j`
at vertex `j>=1`. Finite geometric summation gives `(12)`; no interpolation
is used.

Substitution of `(12)` into `(6)` gives

```text
J_n(X)=4^n alpha(X)+n2^n beta(X)+2^n gamma(X)+delta(X),  n>=1, (13)
```

where

```text
alpha=117(w+h)(w+3h)-345(M_00+2M_01+M_11),
beta=-18h(w+h),

gamma=36(ell_0+ell_1)+360(M_00+M_01)
      -36w^2-126hw-108h^2,

delta=12(M_11-M_01-8M_00)-18(hw+ell_0).                 (14)
```

### The exceptional unary endpoint and its cancellation

At `n=0`, the formal continuation of `(12)` has

```text
L^+=(45,45),              Q_U=((9/2,9/2),(9/2,27/2)),
```

whereas the actual cycle has

```text
L^+(C3)=(30,60),          Q_U(C3)=((3,6),(6,12)).        (15)
```

The resulting exact correction is

```text
J_0(X)=alpha(X)+gamma(X)+delta(X)+epsilon(X),

epsilon(X)=3[7H(X)^2-sum_i(V_i^0-V_i^1)^2].             (16)
```

### Coordinatewise rooted chirality and positivity of the correction

Let `Start_i(X)` and `End_i(X)` count Hamilton paths whose initial and
terminal vertices, respectively, are `i`. The total parity balance in
THM-4184 has the following coordinatewise strengthening:

```text
V_i^1-V_i^0=Start_i(X),          U_i^1-U_i^0=End_i(X).  (16a)
```

For the first identity, represent an object counted by `V_i^epsilon` as
`(P,Q)`, where the nonempty path `P` ends at `i`, `Q` is a Hamilton path of
the complement, and `epsilon=|P| mod 2`. When both paths are nonempty, write
`p,q` for their initial vertices. If `q->p`, move `q` from the beginning of
`Q` to the beginning of `P`. If `p->q` and `P` has at least two vertices,
move `p` from the beginning of `P` to the beginning of `Q`. If `Q` is empty
and `P` has at least two vertices, perform the latter move. Each move
preserves the endpoint `i`, flips the parity of `|P|`, and is reversed by
the next application of the same rule: after a `Q`-to-`P` move, the moved
vertex dominates the new first vertex of `Q`; after a `P`-to-`Q` move, the
moved vertex dominates the new first vertex of `P`.

The only unmatched objects have `P=(i)` and either `Q` empty or
`i->first(Q)`. Concatenation `iQ` bijects them with Hamilton paths starting
at `i`. They have odd `|P|`, proving the first identity in `(16a)`. Reversing
all paths identifies `V_i^epsilon` in the converse tournament with
`U_i^epsilon` in `X`, and Hamilton starts in the converse with Hamilton ends
in `X`; applying the proved identity there gives the dual `U`/terminal
identity.

The inherited source/sink padding formulas make `(16a)` still more concrete.
All `r_i,o_i,Start_i,End_i` in this paragraph are computed in the unpadded
parent `X`. For a universal-source extension `x triangleright X`, THM-4187
equation (27) and THM-4177 equation (12) give

```text
c_xi=2V_i^1=r_i+2Start_i.
```

Combining this with `Start_i=V_i^1-V_i^0`, and applying universal-sink
padding to `X triangleright z`, yields the exact coordinate collapse

```text
V_i=(r_i/2,r_i/2+Start_i),
U_i=(o_i/2,o_i/2+End_i).                               (16b)
```

Indeed, on the sink side the ordinal cross formula gives
`c_iz=2U_i^1`, while THM-4177 equation (12) gives
`c_iz=o_i+2End_i`; on the source side we use the capacity symmetry
`c_xi=c_ix`. In particular every oriented capacity fan is even, and

```text
v_i:=V_i^0+V_i^1=r_i+Start_i,
u_i:=U_i^0+U_i^1=o_i+End_i,

m_V=sum_i(r_i+Start_i)^2,
m_U=sum_i(o_i+End_i)^2.                                (16c)
```

Thus the rooted parity vectors contain exactly oriented fan data plus
Hamilton endpoint data; no further hidden rooted coordinate is needed.

Since the nonnegative `Start_i` sum to `H`,

```text
sum_i Start_i^2<=H^2.                                  (16d)
```

Equality in `(16d)` means that all Hamilton paths start at one vertex `u`.
This holds exactly when `u` is a universal source. The forward implication
needs one small argument: if `u` is not a universal source, take a Hamilton
path of `X-u` and insert `u` immediately after its last in-neighbor along
that path. There is such an in-neighbor, and every later vertex is dominated
by `u`, so the result is a Hamilton path not starting at `u`. Conversely, a
universal source has no possible predecessor and must start every Hamilton
path. Thus `(16)--(16d)` give the all-order bound

```text
epsilon(X)=3[7H(X)^2-sum_i Start_i(X)^2]>=18H(X)^2,

epsilon(X)=18H(X)^2  iff  X has a universal source.     (16e)
```

This positivity is stronger than what is required below: the recurrence
only needs the correction to cancel at the ordinal cut.

For `T=B triangleright C`, THM-4184's endpoint transfer scales every old
`B`-block difference `V_i^0-V_i^1` by `H(C)`. Every new `C`-block endpoint
state is balanced because it equals `(K(B),K(B))*V_j(C)`. Hence

```text
epsilon(B triangleright C)=H(C)^2epsilon(B).             (17)
```

THM-4184's weighted cocycle and `H(A_n)=3` give

```text
F_n(B,C)=8R_+(B,C)+J_n(B triangleright C)-H(C)^2J_n(B). (18)
```

Thus `(17)` cancels the exceptional unary endpoint. Define

```text
mathsf A=alpha(B triangleright C)-H(C)^2alpha(B),
mathsf B=beta(B triangleright C)-H(C)^2beta(B),
mathsf C=gamma(B triangleright C)-H(C)^2gamma(B),
mathsf D=delta(B triangleright C)-H(C)^2delta(B)+8R_+(B,C). (19)
```

> **Theorem 2 (arbitrary-context C-finite law).** For every fixed nonempty
> `B,C` and every `n>=0`,
>
> ```text
> F_n=mathsf A 4^n+mathsf B n2^n+mathsf C 2^n+mathsf D. (20)
> ```
>
> Consequently,
>
> ```text
> F_(n+4)=9F_(n+3)-28F_(n+2)+36F_(n+1)-16F_n,          n>=0. (21)
> ```

The characteristic polynomial is

```text
(t-4)(t-2)^2(t-1)=t^4-9t^3+28t^2-36t+16.
```

The ordinary generating function is rational:

```text
sum_(n>=0)F_n z^n
 =mathsf A/(1-4z)+2mathsf B z/(1-2z)^2
  +mathsf C/(1-2z)+mathsf D/(1-z).                      (22)
```

Thus `(1-z)(1-2z)^2(1-4z)` is a common denominator, not necessarily the
reduced denominator in every context. The sequence is C-finite, therefore
P-recursive and its OGF is D-finite. Once the finite context jets are known,
`F_n` is evaluable in `O(log n)` ring operations by binary powering. This
does not claim `O(log n)` bit complexity: the exact output has `Theta(n)`
bits. Nor is C-finiteness a finite-state assertion.

## 4. Endpoint energy, its dual, and the equality boundary

For a tournament `X`, put

```text
v_i=V_i^0+V_i^1,             a=sum_i v_i=W+H,
b=W+2H=a+H,                  m=sum_i v_i^2,

Delta_V(X)=a(a+2H)-3m=b^2-H^2-3m.                       (23)
```

The name `Delta_V` avoids collision with THM-4184's balanced parity vector
`E(X)`.

> **Theorem 3 (endpoint-energy inequality).** Every nonempty tournament
> satisfies
>
> ```text
> Delta_V(X)>=0,
> Delta_V(X)=0  iff  X is transitive.                   (24)
> ```

### Canonical insertion and the rooted injection

Fix a vertex `i`. Let `mathcal V_i` consist of pairs `(P,Q)` in which `P` is
a nonempty directed simple path ending at `i`, and `Q` is a Hamilton path of
the complementary subtournament. The empty `Q` is allowed with
`H(empty)=1`. By definition,

```text
|mathcal V_i|=v_i.                                      (25)
```

Use the following deterministic insertion. For a directed path
`Q=(q_1,...,q_k)` not containing `i`, scan from the left for the first `q_t`
such that `i->q_t`. Insert `i` immediately before that vertex; if there is no
such vertex, append `i`. Before the first such vertex every `q` dominates
`i`, so the result is directed. The empty path maps to `(i)`. Deleting `i`
recovers `Q`, making insertion injective. This is also the standard inductive
proof that every tournament has a Hamilton path.

Define the tagged map

```text
Phi_i: mathcal V_i -> mathcal H(X) disjoint_union_(j->i) mathcal V_j. (26)
```

If `P=(i)`, insert `i` into `Q` and use the Hamilton-path component. Otherwise
write `P=P'i`, let `j` be the terminal vertex of `P'`, and send `(P,Q)` to the
`j`-tagged object

```text
(P', insert_i(Q)) in mathcal V_j.                        (27)
```

Here `j->i`, and the new complementary path has gained precisely the removed
vertex `i`. The map is injective: in either component, deleting `i` from the
inserted path recovers `Q`; in the second component, appending `i` to `P'`
recovers `P`.

Consequently,

```text
v_i<=H+sum_(j->i)v_j.                                   (28)
```

Every `v_i` is strictly positive: pair the singleton path `(i)` with any
Hamilton path of `X-i` (or with the empty path when `X=(i)`).

Put

```text
sigma_i=H+sum_(j->i)v_j-v_i>=0.                          (29)
```

Every unordered vertex pair appears exactly once as an arc, so

```text
sum_i v_i sigma_i
 =Ha+sum_(j->i)v_iv_j-m
 =Ha+(a^2-m)/2-m
 =Delta_V(X)/2.                                         (30)
```

This proves nonnegativity. If equality holds, then every `v_i>0` and every
`sigma_i>=0`, so every `sigma_i=0`. For any arc `j->i`,

```text
v_i=H+sum_(k->i)v_k>v_j.                                (31)
```

A directed cycle would force a strict cyclic chain of scalar inequalities.
Hence the tournament is acyclic, and an acyclic orientation of a complete
graph is transitive.

Conversely, order the transitive tournament as
`0->1->...->r-1`. A rooted path ending at `i` is determined by an arbitrary
subset of its `i` predecessors, while the complement has its unique Hamilton
path. Thus

```text
v_i=2^i,          H=1,
a=2^r-1,          m=(4^r-1)/3,                           (32)
```

and `Delta_V=0`. This proves the equality classification without assuming
surjectivity of `(26)` or an unstated property of rooted objects. QED.

The energy is not uniquely `V`-sided. Put

```text
u_i=U_i^0+U_i^1,                 m_U=sum_i u_i^2,
Delta_U(X)=a(a+2H)-3m_U.                              (32a)
```

Apply Theorem 3 to the converse tournament and reverse every path. This
identifies its `V`-masses with the `U`-masses of `X`, while preserving
`H,W,a` and transitivity. Hence

```text
Delta_U(X)>=0,
Delta_U(X)=0  iff  X is transitive.                    (32b)
```

This dual controls the `U`-Gram appearing on the left-jet side of `(6)`.
Using the fan/endpoint collapse `(16c)`, the two energy inequalities become

```text
3sum_i(r_i+Start_i)^2<=(W+H)(W+3H),
3sum_i(o_i+End_i)^2  <=(W+H)(W+3H),                    (32c)
```

and equality in either line holds exactly for transitive tournaments. This
is the orientation-aware fan inequality available to the then-open uniform
`A_5` problem; THM-4212 subsequently uses it to close that problem.

## 5. Ordinal defect law and eventual positivity

THM-4184's endpoint transfer says that old `B`-block endpoint masses scale by
`H(C)`, while new `C`-block endpoint masses scale by `b_B`. Therefore

```text
a_(B triangleright C)=H(C)a_B+b_Ba_C,
b_(B triangleright C)=b_Bb_C,
m_(B triangleright C)=H(C)^2m_B+b_B^2m_C.               (33)
```

Substitution in `(23)` gives the exact weighted cocycle

```text
Delta_V(B triangleright C)
 =H(C)^2Delta_V(B)+b_B^2Delta_V(C).                      (34)
```

The start-state transfer has the reversed weights:

```text
m_U(B triangleright C)=b_C^2m_U(B)+H(B)^2m_U(C),       (34a)

Delta_U(B triangleright C)
 =b_C^2Delta_U(B)+H(B)^2Delta_U(C).                     (34b)
```

Thus the endpoint-energy mechanism is a converse-dual pair. The subsequent
leading-coefficient factorization uses `Delta_V` because `alpha` is a
right-factor endpoint statistic; `(32a)--(34b)` supplies the symmetric
left-jet control used by THM-4212 to close the uniform bound.

The first unary coefficient in `(14)` satisfies

```text
alpha(X)=117a(a+2H)-345m=6m+117Delta_V(X)>0.            (35)
```

Strictness follows from `m>0`. Equations `(33)--(35)` factor the leading two
contextual coefficients in `(20)`:

```text
mathsf A=b_B^2alpha(C)>0,
mathsf B=-18H(B)H(C)b_Ba_C<0.                            (36)
```

> **Corollary 4 (fixed-context eventual positivity).** For every fixed
> nonempty `B,C`, there exists `N(B,C)` such that `F_n(B,C)>0` and
> `F_(n+1)(B,C)>F_n(B,C)` for every `n>=N(B,C)`.

Indeed the positive `mathsf A4^n` term dominates the `O(n2^n)` remainder.
More explicitly,

```text
F_(n+1)-F_n
 =3mathsf A4^n+[mathsf B(n+2)+mathsf C]2^n,              (37)
```

which is eventually positive. This corollary alone gives a context-dependent
threshold and does not recover a uniform tail-five result outside transitive
contexts; THM-4212 subsequently supplies the missing lower-jet argument and
proves the sharp uniform threshold `n>=5`.

## 6. Exact C3/right-lollipop family

Let

```text
L_a=P_a triangleright C3,              a>=1.             (38)
```

Put `x=2^a` and `z_i=V_i^0+2V_i^1`. Endpoint transfer and THM-4184's exact
lollipop capacities give

```text
H(L_a)=3,                       W(L_a)=12x-6,

z_0=6,
z_i=9*2^(i-1),                  1<=i<a,
z_c=9x/2                        at each cycle vertex,

sum_i z_i^2=270*4^(a-1)+9.                              (39)
```

The capacity from spine vertex `i` to each cycle vertex is

```text
3*2^(a-i-1), i<a-1;             4, i=a-1.               (40)
```

Set

```text
Lambda_a=sum_i z_i(W-d_i-4r_i).                          (41)
```

For `a>=2`, split the directed edges into prefix, cross, and cycle arcs and
use

```text
Lambda_a=W sum_i z_i-sum_(u->v)c_uv(z_u+5z_v).
```

The three edge-category sums are

```text
prefix: (45/2)4^a+(27a-18)2^(a-2)-90,
cross:  (405/2)4^a+(81a-486)2^(a-2),
cycle:  162*2^a.
```

Finite geometric summation therefore gives

```text
Lambda_a=-9*4^a-(27a+180)2^a+108.                       (42)
```

Direct evaluation at `a=1` gives the same expression.

Specializing `(6)` to the left cycle, whose rooted vectors are all `(1,2)`,
gives for an arbitrary right factor `Y`

```text
R_+(C3,Y)
 =108H(Y)W(Y)+120H(Y)^2
  +18sum_i(V_i^0+2V_i^1)(W-d_i-4r_i)
  +9(3W(Y)+4H(Y))^2-84sum_i(V_i^0+2V_i^1)^2.            (43)
```

Substituting `(39)` and `(42)` proves

> **Theorem 5 (cycle-left lollipop boundary).** For every `a>=1`,
>
> ```text
> R_+(C3,P_a triangleright C3)
>  =162[36*4^a-(3a+20)2^a+4]>0.                         (44)
> ```

Positivity follows from `2^a>=a+1`, hence
`36*2^a>3a+20`. The degenerate endpoint must be evaluated separately:

```text
R_+(C3,C3)=3204,                                        (45)
```

whereas formal substitution `a=0` in `(44)` gives `3240`. The theorem proves
this family's positivity. It does not prove that these lollipops minimize
`R_+(C3,Y)` among all no-sink right factors.

## 7. Rooted chirality is a necessary sidecar

Use the lexicographic pair-bit convention. Let `X,Y` be the strong no-sink
order-six tournaments

```text
X=110011110111101,
Y=110001101110111.                                      (46)
```

Exact reconstruction gives

```text
(H,W,G_+)(X)=(H,W,G_+)(Y)=(43,338,22884).               (47)
```

The permutation `(2,0,5,4,3,1)` maps `X` to the converse of `Y` and maps
their symmetric capacity tensors exactly. Thus the weighted capacity graph,
capacity multiset, and every unsigned capacity summary coincide. Under
converse, however, `U` and `V` exchange, every arc reverses, and outgoing and
incoming capacity fans `o,r` exchange. Nevertheless,

```text
R_+(C3,X)=9,377,256,
R_+(C3,Y)=9,389,016.                                    (48)
```

Therefore order, `H`, `W`, `G_+`, the unsigned capacity tensor up to
isomorphism, and an unordered `{U,V}` quotient do not determine the one-sided
response. The destroyed information is ordered endpoint/fan chirality.
Formula `(6)` repairs the quotient by retaining the orientation-aware rooted
linear and Gram jets. By `(16b)`, these are equivalently the outgoing/incoming
capacity fans together with Hamilton end/start counts; there is no further
hidden rooted state. The lost sidecar is therefore not the `U/V` assignment
in isolation.

The connection contract is

```text
source:       actual ordered tournament factors with rooted U/V states,
map:          factor -> finite response jets lambda,mu,
preserved:    the exact ordinal remainder and contextual defect,
destroyed by unsigned quotient:
              which rooted state is a start/end and which capacity fan is
              outgoing/incoming,
repair:       orientation-aware rooted linear moments and 2x2 Gram matrices,
equivalently: outgoing/incoming fans plus Hamilton end/start counts (16b),
hostile:      the strong no-sink converse pair (46)--(48).                (49)
```

## 8. Complete finite factor-class census through order eight

**FINITE-EXACT, NON-LOAD-BEARING.** Take one `gentourng` representative of
every tournament class of orders one through eight:

```text
order:       1  2  3  4   5   6    7     8
classes:     1  1  2  4  12  56  456  6,880
total: 7,412.                                            (50)
```

For `A_5=C3 triangleright P_5`, the exact 11-jet evaluator checks every one
of the

```text
7,412^2=54,937,744                                      (51)
```

ordered context pairs. Every row satisfies

```text
F_5(B,C)>=10764H(B)^2H(C)^2,                            (52)
```

with unique equality at `B=C=P_1`. The next-smallest normalized value is
`68,364`, at `(P_2,P_1)`; equivalently, the minimum positive slack above
`(52)` is `57,600`.

The genuinely new exact-order-eight sector contains

```text
6,880^2=47,334,400                                      (53)
```

rows and no equality. Its minimum normalized value is

```text
158038868458804096/103829067,                            (54)
```

at

```text
B=1111110110011111011101111111, H(B)=477,
C=1111011110111111111111111111, H(C)=37.
```

The signed-128 evaluator certifies the exact bounds

```text
absolute dot bound =963776444032763562568,
ratio cross-product bound
 =183984884847066925795584715858888 <2^126.              (55)
```

Equations `(50)--(55)` extend THM-4193's finite census by one factor order.
They are evidence within this theorem, not its proof of the universal
arbitrary-context inequality `(52)`. THM-4212 subsequently proves `(52)`
analytically for all factors.

## 9. Scope firewall and honest frontier

The following remain **OPEN**:

```text
R_+(C3,Y)>0 for every no-sink Y,
general (OS+), the no-sink/no-source gate law,
the all-strong residual and order-eleven asymmetric bank.                (56)
```

THM-4212 subsequently proves

```text
F_n(B,C)>=10764H(B)^2H(C)^2 for all n>=5 and nonempty B,C,
F_n(B,C)>0 for every nonempty B,C iff n>=5,
```

with lower-bound equality exactly at `n=5,B=C=P_1`.

The endpoint-energy theorem controls the leading tail coefficient. It does
not sign the transient coefficients `gamma,delta` or the base remainder.
The rooted-chirality hostile shows why those lower jets cannot be inferred
from unsigned scalar data. Fixed-context eventuality by itself therefore
cannot be promoted to a uniform statement; THM-4212 closes the gap by using
the coordinatewise Start/End identities, dual endpoint energy, and mixed
same-vertex rooted moments.

## 10. Replay

Primary exact audit:

```bash
python3 -B \
  04-computation/tournament_cycle_prefix_arbitrary_context_thm4208.py
python3 -O -B \
  04-computation/tournament_cycle_prefix_arbitrary_context_thm4208.py
PYTHONHASHSEED=4208 python3 -B \
  04-computation/tournament_cycle_prefix_arbitrary_context_thm4208.py
```

Independent clean-room audit:

```bash
clang++ -std=c++17 -O0 -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_cycle_prefix_arbitrary_context_thm4208_independent_audit.cpp \
  -o /tmp/thm4208-independent-O0
/tmp/thm4208-independent-O0

clang++ -std=c++17 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_cycle_prefix_arbitrary_context_thm4208_independent_audit.cpp \
  -o /tmp/thm4208-independent-O3
/tmp/thm4208-independent-O3

clang++ -std=c++17 -O1 -g -fsanitize=address,undefined \
  -fno-omit-frame-pointer -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_cycle_prefix_arbitrary_context_thm4208_independent_audit.cpp \
  -o /tmp/thm4208-independent-san
ASAN_OPTIONS=detect_leaks=0 /tmp/thm4208-independent-san
```

Complete order-eight factor-class census:

```bash
clang++ -std=c++20 -O3 -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_cycle_prefix_arbitrary_context_order8_census_thm4208.cpp \
  -o /tmp/thm4208-order8-census
python3 -B \
  04-computation/tournament_cycle_prefix_arbitrary_context_order8_emit_thm4208.py \
  | /tmp/thm4208-order8-census
```

The three primary streams byte-match the frozen primary output. The C++ O0,
O3, and ASan/UBSan streams byte-match the frozen independent output. **QED for
the all-order symbolic statements `(6)`, `(16a)--(16e)`, `(20)--(24)`,
`(32a)--(32c)`, `(34)--(34b)`, `(36)--(45)`, and the hostile `(46)--(48)`;
FINITE-EXACT only for `(50)--(55)`.**
