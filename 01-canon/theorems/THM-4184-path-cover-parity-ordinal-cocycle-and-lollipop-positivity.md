---
id: THM-4184
title: "Path-cover parity, ordinal cocycle, and lollipop positivity"
status: >
  PROVED all-order two-path-cover parity balance, associative ordinal-chain
  capacity transfer, rank-one nonadjacent-block collapse, normalized R_+
  cocycle, and an infinite two-parameter transitive-spine/C3 ordinal-remainder
  positivity theorem + FINITE-EXACT strict positive local interaction on
  317,680 ordered factor-class presentations with the first two factors of
  orders one through six and the no-sink third factor of orders three through
  six + VERIFIED-EXACT + INDEPENDENTLY AUDITED. All-order no-sink local-
  interaction positivity, (OS+), the no-sink/no-source sign law, and the
  order-eleven asymmetric bank remain OPEN.
source: tournament-remainder-cocycle-20260826
depends_on:
  - THM-4177-root-split-johnson-current-and-source-sink-boundary-census
  - THM-4181-ordinal-sum-capacity-transfer-and-parity-component-exchange
related:
  - THM-4133-strong-cyclic-substitution-johnson-centrality-counterexample
  - THM-4169-prime-parent-one-vertex-augmentation-and-quartic-johnson-transfer
script: 04-computation/tournament_ordinal_cocycle_parity_thm4184.py
output: 05-knowledge/results/tournament_ordinal_cocycle_parity_thm4184.out
independent_audit_script: 04-computation/tournament_ordinal_cocycle_parity_thm4184_independent_audit.cpp
independent_audit_output: 05-knowledge/results/tournament_ordinal_cocycle_parity_thm4184_independent_audit.out
script_sha256: 9ab09bf8b70ee5dcf3d86698180a456d67f012655b49a16dfea9903caefbb39c
output_sha256: 785b67dc9e4dbb4075efff8152b2250df70277856e15337b32f09aefec97eb13
independent_audit_script_sha256: aff976d215ce1b70bfc23c6ace897cfd572ffd42191499204f771e9ab249cea9
independent_audit_output_sha256: 1c26d46a27188c47deed877e55a2c0f10d2a07d6c0d32ca3395666f5c807e8f2
gentourng_sha256: 89df605922cc574b28688248b7c256d24342cc615f887e89b2d096038970c110
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standalone subset/exposed-word engine checks the terminal-transfer
  involution on 15,810 two-path-cover objects in all 76 tournament-class
  representatives through order six, 5,776 pair sidecar products, 428 direct
  three-block transfers/cocycles, the complete 317,680-presentation declared
  local-interaction census, the minimal hostile, and the closed forms.
independent_audit: >
  ACCEPT. A separate C++ literal-permutation referee checks 1,099 labelled
  tournaments through order five and 78,418 two-path covers, 505 direct pairs,
  339 direct triples/cocycles, and 191,250 labelled no-sink-third-factor
  presentations. Clang O0/O3 and ASan/UBSan streams byte-match.
---

# THM-4184 -- path-cover parity, ordinal cocycle, and lollipop positivity

**PROVED algebra and infinite family + FINITE-EXACT local-interaction census
+ VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-4181 retained two rooted parity coordinates at an ordinal cut and proved
their total difference only at odd order. The missing even-order statement is
in fact true. A fixed-point-free terminal-transfer involution balances the
optional-path state in every nonempty tournament. Consequently an entire
intermediate strong component is a scalar parity filter: adjacent ordinal
blocks can still have genuine rank-two transfer, but every nonadjacent block
pair has rank-one cross capacity.

Associativity then gives an exact normalized remainder cocycle. The cocycle is
an identity, not a positivity proof. Its local interaction is positive in two
independent finite universes when the third factor has no sink, but unconditional
positivity fails minimally at `C3 triangleright 1 triangleright 1`. A separate
closed-form argument does prove `(OS+)` on the infinite two-parameter family

```text
R_+(P_a, P_b triangleright C3)>0,        a>=1, b>=0,
```

where `P_n` is the transitive tournament of order `n`.

## 1. Conventions and the parity algebra

All tournaments in the proved statements are finite and nonempty. Use the
Hamilton count `H`, capacity tensor `c`, mass

```text
W(X)=sum_(i<j)c_ij(X),
```

and gate `G_+=D+2C` of THM-4177/4181. Let `U_a^epsilon(X)` and
`V_a^epsilon(X)` be THM-4181's rooted start/end path-cover counts; singleton
zero-arc paths are included. Put

```text
P^epsilon(X)=sum_a U_a^epsilon(X)=sum_a V_a^epsilon(X),       (1)

E(X)=(E^0,E^1)=(H(X)+P^0(X), P^1(X)).                         (2)
```

The equality of the two sums in `(1)` is literal: every nonempty directed
simple path has one start and one end. The vector `E` counts ordered two-path
covers `(P,Q)` of `V(X)`, where `P` is allowed to be empty and the sign/parity
is that of `|P|`; `Q` is a chosen Hamilton path of the complement and may be
empty only when `P` spans.

For parity vectors use convolution in the group algebra of `C_2`:

```text
(x*y)^0=x^0y^0+x^1y^1,
(x*y)^1=x^0y^1+x^1y^0.                                    (3)
```

The empty product is `(1,0)`. Write `U_a^+=U_a^0+U_a^1` and
`V_a^+=V_a^0+V_a^1`.

## 2. The all-order terminal-transfer involution

> **Theorem 1 (all-order two-path-cover parity balance).** Every nonempty
> tournament `X` satisfies
>
> ```text
> P^1(X)-P^0(X)=H(X),                                      (4)
> E^0(X)=E^1(X)=K(X),
> K(X)=W(X)/2+H(X)=(W(X)+2H(X))/2.                          (5)
> ```

Thus THM-4181's odd-order identity extends to every order.

### Proof

Consider an ordered two-path cover `(P,Q)` counted by `E`. Since `X` is
nonempty, the two components are not both empty. Define a move as follows.

1. If `P` is empty, move the terminal vertex of `Q` across, making it the
   singleton path `P`.
2. If `Q` is empty, move the terminal vertex of `P` across, making it the
   singleton path `Q`.
3. If both are nonempty, let `p,q` be their terminal vertices. Completeness
   gives exactly one of `p->q` and `q->p`. In the first case delete `q` from
   the end of `Q` and append it to `P`; in the second delete `p` from the end
   of `P` and append it to `Q`.

The moved vertex is appended along a directed arc, so the result is again an
ordered two-path cover. It changes `|P|` by one. To see involutivity, suppose
for example that `q` moved from `Q` to `P`. If the shortened `Q` is nonempty,
its new terminal is the predecessor of `q`, hence dominates `q`; the comparison
rule moves `q` back. If it is empty, the empty-component rule moves `q` back.
The other direction and the singleton boundaries are identical. The move is
therefore fixed-point-free, sign reversing, and involutive.

Hence the signed cover sum is

```text
E^0-E^1=H+P^0-P^1=0,
```

which proves `(4)`. THM-4181 gives `P^0+P^1=W+H`; solving the two equations
gives `P^0=W/2`, `P^1=W/2+H`, and `(5)`. QED.

The nonempty hypothesis is sharp for this formulation: the empty tournament
has optional state `(1,0)`. Equation `(5)` balances only the **total** rooted
state. Individual vectors `U_a` and `V_a` need not have equal coordinates.

## 3. Associative sidecar transfer and rank-one separation

Let `T=A triangleright B`. Besides THM-4181's capacity formulas, the rooted
sidecars themselves obey the following exact transfer.

> **Theorem 2 (associative parity-state transfer).** With convolution `(3)`,
>
> ```text
> U_a(T)=U_a(A)*E(B),                   a in A,
> U_b(T)=H(A)U_b(B),                    b in B,
>
> V_a(T)=H(B)V_a(A),                    a in A,
> V_b(T)=E(A)*V_b(B),                   b in B,             (6)
>
> E(A triangleright B)=E(A)*E(B),
> K(A triangleright B)=2K(A)K(B).                            (7)
> ```

### Proof

A directed path starting in `A` either remains in `A` or crosses the ordinal
cut once and continues along a nonempty path of `B`. Treating the absent
`B` segment as an empty path with weight `H(B)` gives the first line of `(6)`.
A path starting in `B` cannot return to `A`, which gives the second formula.
The endpoint formulas are the exact reverse bookkeeping.

Summing the start formulas counts every optional path in `A triangleright B`.
Equivalently, an optional path decomposes uniquely into optional paths in the
two factors, not both absent unless the global path is absent. Complementary
Hamilton counts factor, proving the first equation in `(7)`. Substitution of
`E(A)=(K(A),K(A))` and `E(B)=(K(B),K(B))` proves the second. QED.

Now let

```text
T=X_1 triangleright X_2 triangleright ... triangleright X_s. (8)
```

Write `h_r=H(X_r)`, `K_r=K(X_r)`, and

```text
h_<i=product_(r<i)h_r,        h_>j=product_(r>j)h_r.      (9)
```

> **Theorem 3 (exact ordinal-chain capacity).** For vertices `u,v` in the
> indicated blocks,
>
> ```text
> c_uv(T)=(product_(r!=i)h_r)c_uv(X_i),       u,v in X_i, (10)
> ```
>
> and for `u in X_i`, `v in X_j`, `i<j`,
>
> ```text
> c_uv(T)=2h_<i h_>j
>   [U_u(X_i)*E(X_(i+1))*...*E(X_(j-1))*V_v(X_j)]^0.      (11)
> ```
>
> In particular, adjacent cross blocks have the rank-at-most-two matrix
>
> ```text
> 2h_<i h_>(i+1)[U_u^0V_v^0+U_u^1V_v^1],                 (12)
> ```
>
> while, if `m=j-i-1>=1`, the nonadjacent cross block is
>
> ```text
> 2^m h_<i h_>j (product_(i<r<j)K_r) U_u^+(X_i)V_v^+(X_j), (13)
> ```
>
> and therefore has ordinary rank at most one.

### Proof

An internal exposed path cannot leave and return to one block, giving `(10)`.
For `(11)`, a directed path from `u` to `v` has rooted endpoint segments in
`X_i,X_j`, an optional segment in each intermediate block, and no vertices in
the exterior blocks. This decomposition is unique; the exterior complement
contributes `h_<i h_>j`, and even total vertex parity is exactly the zero
coefficient of convolution `(3)`. THM-4177's odd-path formula supplies the
factor two. The empty intermediate product gives `(12)`. For `m>=1`,

```text
(K_(i+1),K_(i+1))*...*(K_(j-1),K_(j-1))
 =2^(m-1)(product K_r)(1,1),
```

which gives `(13)`. QED.

For three blocks, the new rank-one law is especially transparent:

```text
c_ac(A triangleright B triangleright C)
 =2K(B)U_a^+(A)V_c^+(C).                                (14)
```

This does not contradict THM-4181's genuine rank-two adjacent example: the
middle block in `(14)` is nonempty. Nor does it extend to arbitrary module
substitution, where a path can enter one module in several sections.

## 4. The normalized remainder cocycle

Recall

```text
R_+(A,B)=G_+(c(A triangleright B))
          -H(B)^2G_+(c(A))-H(A)^2G_+(c(B)).               (15)
```

Hamilton paths factor under ordinal sum. Since every tournament has positive
Hamilton count, define

```text
g(X)=G_+(c(X))/H(X)^2,
rho(A,B)=R_+(A,B)/(H(A)^2H(B)^2)
        =g(A triangleright B)-g(A)-g(B).                  (16)
```

> **Theorem 4 (remainder two-cocycle).** For all nonempty `A,B,C`,
>
> ```text
> rho(A,B)+rho(A triangleright B,C)
>  =rho(B,C)+rho(A,B triangleright C),                    (17)
> ```
>
> equivalently,
>
> ```text
> H(C)^2R_+(A,B)+R_+(A triangleright B,C)
>  =H(A)^2R_+(B,C)+R_+(A,B triangleright C).              (18)
> ```

### Proof

Expand either side of `(17)` using `(16)`. Both equal

```text
g(A triangleright B triangleright C)-g(A)-g(B)-g(C).
```

Multiplying by `H(A)^2H(B)^2H(C)^2` gives `(18)`. QED.

The mixed interval interaction is therefore well typed by

```text
kappa(A,B,C)
 =rho(A triangleright B,C)-rho(B,C)
 =rho(A,B triangleright C)-rho(A,B),                      (19)

Theta(A,B,C)
 =R_+(A triangleright B,C)-H(A)^2R_+(B,C)
 =R_+(A,B triangleright C)-H(C)^2R_+(A,B),                (20)

Theta=H(A)^2H(B)^2H(C)^2 kappa.                           (21)
```

Thus `kappa` is the second interval difference

```text
g(A triangleright B triangleright C)-g(A triangleright B)
 -g(B triangleright C)+g(B).                              (22)
```

The cocycle is formal once `(15)` is normalized. Its usefulness is that
Theorem 3 computes the actual three-block tensor without discarding the
middle parity state.

### Exact terminal-strong reformulation of `(OS+)`

Every no-sink tournament `B` is either strong, or has a unique decomposition

```text
B=P triangleright S,                                     (23)
```

where `S` is its terminal strong component, `|S|>=3`, and `P` is the nonempty
union of the preceding components. Hence `(OS+)` is exactly the conjunction
of

```text
rho(A,S)>0                                                (24)
```

for strong `S`, and the local three-block compensation inequality

```text
rho(A,P)+kappa(A,P,S)>0                                  (25)
```

for nonempty `A,P` and strong `S`. This is a typed reduction to a terminal
strong factor plus one local interaction, not a solution: `(25)` remains an
all-order inequality equivalent to the nonstrong part of `(OS+)`.

## 5. An infinite `(OS+)` family

Let `P_n` be the transitive tournament of order `n>=1`; use `P_0` only as the
empty prefix notation. Let `C3` be the directed triangle and

```text
L_n=P_n triangleright C3,        n>=0.                    (26)
```

Thus `L_0=C3`, and every `L_n` has no sink.

> **Theorem 5 (transitive-spine/C3 lollipop positivity).** For every
> `a>=1,b>=0`,
>
> ```text
> R_+(P_a,L_b)>0.                                         (27)
> ```
>
> More sharply,
>
> ```text
> R_+(P_a,L_b)>=R_+(P_a,C3)>=120,                         (28)
> ```
>
> with equality in the first inequality exactly when `b=0`, and global
> equality `120` exactly at `(a,b)=(1,0)`.

### Exact capacities and gates

Number the vertices of `P_n` increasingly. Its unique Hamilton path and the
odd-path formula give, for `i<j`,

```text
c_ij(P_n)=2                         if j=i+1,
c_ij(P_n)=2^(j-i-1)                 if j>=i+2.             (29)
```

The incoming and outgoing capacity masses at vertex `i` are

```text
r_i=0 (i=0), else 2^i;
o_i=2^(n-1-i) (i<n-1), else 0.                            (30)
```

Consequently `C(P_n)=sum_i(o_i+r_i)(o_i-r_i)=0` and, for `n>=2`,

```text
W(P_n)=2^n-2,
sum_e c_e^2=(4^n+24n-28)/9,
sum_i d_i^2=2(4^n-4)/3+(n-2)2^n.                         (31)
```

Using `D=(W^2+sum_e c_e^2-sum_i d_i^2)/2` gives

```text
G_+(P_1)=0,
G_+(P_n)=[4*4^n-9(n+2)2^n+24n+32]/18,       n>=2.       (32)
```

For `L_n`, the prefix capacities are `3c(P_n)`, the three cycle capacities
are `2`, and the capacity from prefix vertex `i` to each cycle vertex is

```text
t_i=3*2^(n-i-1)             for i<n-1,
t_(n-1)=4.                                             (33)
```

Indeed the rooted state at `i<n-1` in `P_n` is
`(2^(n-i-2),2^(n-i-2))`, while the last vertex has state `(0,1)`; pair it
with `(1,2)` at a cycle endpoint in THM-4181's cross formula.

Here

```text
sum_i t_i=3*2^n-2,             sum_i t_i^2=3*4^n+4.     (34)
```

Substitution of `(29)--(34)` into THM-4181's exact block-gate formula (with
the `n=1` endpoint evaluated directly) gives

```text
G_+(L_0)=0,
G_+(L_n)=2[37*4^n-9(n+4)2^n+6n-4],          n>=1.       (35)
```

For completeness, the only nontrivial left-linear sum in that substitution,
for `n>=2`, is

```text
sum_i t_i[W(P_n)-d_i(P_n)+4o_i(P_n)]
 =3W(P_n)^2+3(4^n-4)-3(n-2)2^(n-1)
  +4W(P_n)-2^(n+1).                                      (36)
```

Thus `(35)` is a finite geometric-sum identity, not an interpolated formula.

### Convexity proof of `(27)`

Put `F_n=G_+(L_n)` and `Delta_n=F_(n+1)-F_n`. Directly,

```text
Delta_0=120,
Delta_n=6[37*4^n-3(n+6)2^n+2],              n>=1,       (37)

Delta_1-Delta_0=528,
Delta_(n+1)-Delta_n=18*2^n[37*2^n-n-8]>0,   n>=1.       (38)
```

Hence the increments are positive and strictly increasing. Therefore

```text
F_(a+b)-F_b=sum_(j=b)^(a+b-1)Delta_j
            >=sum_(j=0)^(a-1)Delta_j=F_a.                (39)
```

Since `H(P_a)=1` and `H(L_b)=3`, equations `(15)` and `(39)` give

```text
R_+(P_a,L_b)=F_(a+b)-F_b-9G_+(P_a)
             >=F_a-9G_+(P_a)=R_+(P_a,C3).                (40)
```

For `a=1`, the last quantity is `120`. For `a>=2`, `(32)` and `(35)` give

```text
R_+(P_a,C3)
 =(3/2)[48*4^a-(9a+42)2^a-16]>0.                         (41)
```

After division by `2^a`, positivity follows at `a=2` from
`48*2^a>9a+42+16/2^a`; its left side doubles at the next step while the
nonconstant part of the right side grows by less than `9`. Strict convexity
in `(38)` gives the equality claims. QED.

This is an unbounded SCC-chain family: `L_b` has `b` singleton components
followed by one `C3` component. It does not classify arbitrary SCC chains.

## 6. Finite local-interaction census

**FINITE-EXACT.** The primary universe contains one `gentourng`
representative of every tournament class `A,B` of orders one through six and
every no-sink tournament class `C` of orders three through six:

```text
A and B class counts: 1,1,2,4,12,56,       total 76 each,
C no-sink counts:          1,2,8,44,        total 55.
```

Thus there are exactly

```text
76*76*55=317,680                                           (42)
```

ordered factor-class presentations. On every presentation,

```text
Theta(A,B,C)>0.                                           (43)
```

There are no zero rows. The exact normalized minimum is

```text
kappa(1,1,C3)=72,
Theta(1,1,C3)=648,       H(1)^2H(1)^2H(C3)^2=9.           (44)
```

The presentation-stream SHA-256 digest is

```text
e19aec23212629e20adda08c8ccf2a3386180b6259cb89defb55e7208491000e.
```

The primary engine also checks the explicit involution on `15,810` ordered
two-path covers in all `76` class representatives, `5,776` pair sidecar
products, `384` direct pair transfers, and `428` direct three-block transfers
and weighted cocycles.

The independent referee uses literal permutations rather than subset DP or
`gentourng`. It checks all `1,099` labelled tournaments through order five,
`78,418` explicit two-path-cover objects, `505` direct pair transfers, and
`339` direct three-block transfers/cocycles. Its separate labelled local
universe has all `75` labelled factors of orders one through four in the first
two positions and all `34` labelled no-sink tournaments of orders three and
four in the third, hence

```text
75*75*34=191,250                                          (45)
```

presentations. Again every `Theta` is positive and the normalized minimum is
`72` at `(1,1,C3)`.

These exact facts are evidence for, not a proof of, the all-order implication

```text
C has no sink => Theta(A,B,C)>0.                          (46)
```

Equation `(46)` remains **OPEN**.

## 7. Minimal hostile to unconditional local positivity

The no-sink condition in `(43)` cannot simply be deleted. In the complete
factor-class universe through order five, ordered by total factor order,
factor orders, and pair-bit labels, the first negative interaction is

```text
(A,B,C)=(C3,1,1),
Theta(C3,1,1)=-216.                                      (47)
```

It is already visible from

```text
R_+(C3,1)=-72,
R_+(C3 triangleright 1,1)=-216,
R_+(C3,1 triangleright 1)=-288.                          (48)
```

Both forms in `(20)` give `-216`. Total factor order three has the unique
all-singleton zero; there is no negative row through total factor order four.
The independent labelled referee finds the two labelled orientations of the
same `C3` witness at total order five.

Thus the overstrong conjecture `Theta>0` for all nonempty triples is
**REFUTED**. The first failed implication is treating a normalized coboundary
as an unconditional positive associator. The exact cocycle, parity balance,
rank-one separation, no-sink finite survivor, and lollipop theorem remain.

## 8. Type, loss, and consequence firewall

The connection contract is:

```text
source:       actual ordered tournament blocks with rooted U/V states,
map:          ordinal concatenation and C2 parity convolution,
target:       the actual capacity tensor and normalized gate g,
preserved:    H, every capacity coordinate, G_+, R_+, and associativity,
lost by bare SCC contraction:
              individual endpoint parity states and the chosen block cut,
needed sidecar:
              U/V at adjacent cuts; H and scalar K across intermediates,
cheapest test: direct exposed-word reconstruction and the C3,1,1 hostile.
```

Relabeling within each factor is harmless, but factor order is structural and
cannot be permuted. Formula `(13)` scalarizes a nonadjacent cross block; it does
not reconstruct adjacent rooted coordinates or arbitrary module sections.

Nothing here proves `(OS+)` beyond Theorem 5's family, the all-order local
inequality `(46)`, THM-4177's no-sink/no-source sign law, its strong residual,
the order-eleven asymmetric bank, exact Johnson cosets, or actual response
maximizers. The 317,680 and 191,250 rows are presentations, not child classes.

## 9. Replay

Primary exact audit:

```bash
python3 -B \
  04-computation/tournament_ordinal_cocycle_parity_thm4184.py
python3 -O -B \
  04-computation/tournament_ordinal_cocycle_parity_thm4184.py
PYTHONHASHSEED=4184 python3 -B \
  04-computation/tournament_ordinal_cocycle_parity_thm4184.py
```

Independent literal-permutation audit:

```bash
clang++ -std=c++20 -O0 -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_ordinal_cocycle_parity_thm4184_independent_audit.cpp \
  -o /tmp/thm4184-independent-O0
/tmp/thm4184-independent-O0

clang++ -std=c++20 -O3 -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_ordinal_cocycle_parity_thm4184_independent_audit.cpp \
  -o /tmp/thm4184-independent-O3
/tmp/thm4184-independent-O3

clang++ -std=c++20 -O1 -g -fsanitize=address,undefined \
  -fno-omit-frame-pointer -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_ordinal_cocycle_parity_thm4184_independent_audit.cpp \
  -o /tmp/thm4184-independent-san
/tmp/thm4184-independent-san
```

The three Python streams byte-match the frozen primary output. Clang O0/O3
and ASan/UBSan byte-match the frozen independent output. **QED for the
all-order parity identity, exact chain transfer, cocycle, lollipop family,
and stated finite universes only.**
