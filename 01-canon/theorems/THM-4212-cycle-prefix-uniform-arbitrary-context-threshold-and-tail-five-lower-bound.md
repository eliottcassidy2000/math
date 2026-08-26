---
id: THM-4212
title: "Cycle-prefix uniform arbitrary-context threshold and tail-five lower bound"
status: >
  PROVED exact eleven-coordinate A5 response difference; a stronger all-order
  right-context monotonicity with strict non-singleton boundary; the sharp
  uniform arbitrary-context inequality F_n(B,C)>=10764H(B)^2H(C)^2 for every
  n>=5; and the exact uniform positivity threshold n>=5, with lower-bound
  equality iff n=5 and B=C=P1 + VERIFIED-EXACT + INDEPENDENTLY AUDITED. This
  proves only the cycle-first transitive-tail prefix family; general (OS+),
  the no-sink/no-source gate law, and the order-eleven asymmetric bank remain
  OPEN.
source: codex-tournament-a5-analytic-20260826
depends_on:
  - THM-4187-universal-source-padding-supermodularity-and-transitive-left-ordinal-positivity
  - THM-4193-cycle-first-transitive-tail-crossing-and-transitive-context-positivity
  - THM-4208-cycle-prefix-arbitrary-context-recurrence-endpoint-energy-and-eventual-positivity
related:
  - THM-4177-root-split-johnson-current-and-source-sink-boundary-census
  - THM-4181-ordinal-sum-capacity-transfer-and-parity-component-exchange
  - THM-4184-path-cover-parity-ordinal-cocycle-and-lollipop-positivity
  - THM-4202-vertex-transitive-ordinal-remainder-positivity
script: 04-computation/tournament_a5_uniform_context_bound_thm4212.py
output: 05-knowledge/results/tournament_a5_uniform_context_bound_thm4212.out
independent_audit_script: 04-computation/tournament_a5_uniform_context_bound_independent_audit_thm4212.cpp
independent_audit_output: 05-knowledge/results/tournament_a5_uniform_context_bound_independent_audit_thm4212.out
script_sha256: 15aec170321b0038bc05fbc67ddc92fa37edcee4cf9b7702fd5405860535a3d3
output_sha256: 94f5631b67f82c083dbd2e37747c4f39206612b2a70e003b04a47dc14be4c7f6
independent_audit_script_sha256: 74c74a860988924b1f1052e55f4938dd7234678e0c626ad49e72af10a1d1820f
independent_audit_output_sha256: efab533a6676df847d928d51b48864db93afb8fe6f0dbeb8840a60aa81827984
hash_basis: raw LF bytes
primary_audit: >
  PASS. A clean exact SymPy derivation starts from THM-4208's A5 jet and
  ordinal-transfer identities, reconstructs all eleven coordinates of
  lambda(A5 triangleright B)-lambda(B), canonicalizes the right-context
  slack in endpoint variables, and checks every coefficient-debt and unary
  identity symbolically. Normal and python -O streams byte-match.
independent_audit: >
  ACCEPT. A standalone warning-clean C++17 referee imports no transfer or
  response-jet code. It rebuilds every composite tournament from labelled
  adjacency, Hamilton-subset DP, the literal odd-directed-path capacity
  formula, and the direct G_+ incidence kernel. It checks all 25 labelled
  pairs with factor orders at most three and total factor order at most four,
  including one equality and signed controls. Clang O0/O3 and ASan/UBSan
  streams byte-match.
---

# THM-4212 -- cycle-prefix uniform arbitrary-context threshold and tail-five lower bound

**PROVED sharp uniform arbitrary-context threshold + exact tail-five lower
bound + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-4208 proves an exact eleven-coordinate response jet for

```text
A_n=C3 triangleright P_n,
F_n(B,C)=R_+(A_n triangleright B,C)-R_+(B,C),
```

and supplies the endpoint-energy inequality and coordinate identities

```text
V_i^1-V_i^0=Start_i,             U_i^1-U_i^0=End_i.      (1)
```

It left the uniform threshold open. This theorem closes it. The mechanism
is not asymptotic dominance: at tail length five, the left response
difference lies in a cone dual to every non-singleton right context after
paying the exact endpoint-energy debt. THM-4187's transitive-left interaction
then propagates the bound to every later tail.

## 1. Statements

All tournaments below are finite and nonempty. Put

```text
A=A_5=C3 triangleright P_5.                              (2)
```

> **Theorem 1 (tail-five right-context monotonicity).** For all nonempty
> tournaments `B,C`,
>
> ```text
> F_5(B,C)>=H(C)^2F_5(B,P_1).                            (3)
> ```
>
> Equality holds in `(3)` if and only if `C=P_1`.

> **Theorem 2 (sharp uniform arbitrary-context bound and threshold).** For
> every `n>=5` and all nonempty tournaments `B,C`,
>
> ```text
> F_n(B,C)>=10764H(B)^2H(C)^2>0.                         (4)
> ```
>
> Equality in the first inequality holds if and only if
> `n=5` and `B=C=P_1`. Moreover,
>
> ```text
> F_n(B,C)>0 for every nonempty B,C    iff n>=5.         (5)
> ```

The exact equality value is THM-4193's crossing value
`F_5(P_1,P_1)=10764`.

## 2. Exact left response difference

Use all notation from THM-4208. For `B`, write

```text
h=H(B),                    w=W(B),
s=(s_0,s_1)=(w/2,w/2+h),
L_e=L_e^+(B),
q_ef=sum_i U_i^e(B)U_i^f(B),
n_ef=sum_i U_i^e(B)V_i^f(B).                            (6)
```

The mixed same-vertex moments `n_ef` are not an additional theorem
hypothesis. They arise when cross column sums in `A triangleright B` are
paired with the rooted `U` states of `B`.

Let `lambda` be THM-4208 equation (9), index its eleven coordinates from
zero, and define

```text
D(B)=lambda(A triangleright B)-lambda(B)=(D_0,...,D_10). (7)
```

> **Lemma 3 (exact A5 response difference).** The coordinates in `(7)` are
>
> ```text
> D_0=h(2286h+1151w),
>
> D_1=2[471168h^2+475182hw+8L_0
>       -1134n_00-1152n_01+119799w^2],
>
> D_2=2[473454h^2+476325hw+8L_1
>       -1134n_10-1152n_11+119799w^2],
>
> D_3=h(2286h+1151w),
> D_4=1151h(2h+w),
>
> D_5=474624h^2+476910hw+48q_00+119803w^2,
>
> D_6=2[476910h^2+478061hw+48q_01+119803w^2],
>
> D_7=479212h^2+479212hw+48q_11+119803w^2,
>
> D_8=-1390176h^2-1401606hw+16q_00-353279w^2,
>
> D_9=2[-1401606h^2-1407361hw+16q_01-353279w^2],
>
> D_10=-1413116h^2-1413116hw+16q_11-353279w^2.          (8)
> ```

### Derivation

THM-4208 gives the exact `A=A_5` data

```text
H(A)=3,             W(A)=378,             s_A=(189,192),

Q_U(A)=((7677,7677),(7677,7686)),
L^+(A)=(116064,116622).                                 (9)
```

Put `b=w+2h` and `T=A triangleright B`. Exact ordinal transfer gives

```text
H(T)=3h,                    W(T)=384w+762h,
s_T=(W(T)/2,W(T)/2+3h),

Q_U(T)=9Q_U(B)+(30717/4)b^2 ((1,1),(1,1)).             (10)
```

For completeness, keep the linear-moment transfer explicit. Set

```text
Xi=(b h/2)(116064+116622)
   +(b/2)384(w+h)(189+192)
   +3b(15354s_0+15363s_1).                              (11)
```

For an `A` vertex, its `U` state becomes balanced after convolution with
`E(B)`, while for a `B` vertex it scales by `H(A)=3`. The cross row and
column sums are

```text
p_i=2<U_i(A),s>,
q_j=378V_j^0(B)+384V_j^1(B).                            (12)
```

The two blockwise linear brackets therefore give

```text
L_0^+(T)=Xi+9L_0+3(381w+762h)s_0-3(378n_00+384n_01),
L_1^+(T)=Xi+9L_1+3(381w+762h)s_1-3(378n_10+384n_11).   (13)
```

Indeed, on the `A` block the bracket is the old bracket scaled by `h`, plus
`384(w+h)+3p_i`; on the `B` block it is three times the old bracket, plus
`381w+762h-q_j`. Equations `(9)--(13)` substituted into the eleven entries
of `lambda` give `(8)` by direct expansion. No interpolation or sign
estimate enters Lemma 3.

By THM-4208's exact response factorization,

```text
F_5(B,C)=D(B) dot mu(C).                                (14)
```

## 3. Exact right-context slack

For the right tournament `C`, write

```text
k=H(C),                 z=W(C),
v_i=V_i^0(C)+V_i^1(C),
t_i=V_i^1(C)-V_i^0(C)=Start_i(C),
g_i=z-d_i(C)-4r_i(C),

m=sum_i v_i^2,         p=sum_i v_i t_i,
tau=sum_i t_i^2.                                        (15)
```

The coordinatewise identity `(1)` is load-bearing: `v_i,t_i>=0`,
`sum_i v_i=z+k`, and `sum_i t_i=k`.

Subtracting the singleton right jet from `(14)` gives

```text
S(B,C):=F_5(B,C)-k^2F_5(B,P_1).                        (16)
```

Changing variables from `(V_i^0,V_i^1)` to `(v_i,t_i)` yields the exact
form

```text
S=K_0zk+E+K_2z^2+K_1zk+A_m m+A_p p+A_tau tau+A_k k^2, (17)
```

where

```text
K_0=2[473454h^2+476329hw+4L_0+4L_1
      -567(n_00+n_10)-576(n_01+n_11)+119799w^2],

E=sum_i [h(2294h+1151w)v_i+8h^2t_i]g_i,

K_2=476914h^2+478061hw+12q_00+24q_01+12q_11+119803w^2,

K_1=956122h^2+957273hw+48q_01+48q_11+239606w^2,

A_m=-1401626h^2-1407361hw+4q_00+8q_01+4q_11-353279w^2,

A_p=-11470h^2-5755hw-8q_00+8q_11,

A_tau=4[-5h^2+q_00-2q_01+q_11],

A_k=1413116h^2+1413116hw-16q_11+353279w^2.             (18)
```

Equation `(17)` is an identity, not a relaxation.

## 4. The left-factor coefficient debts

Put

```text
a=w+h,                b=w+2h,
N_0=n_00+n_10,        N_1=n_01+n_11.                   (19)
```

All rooted coordinates are nonnegative, and their totals are

```text
sum_i(U_i^0+U_i^1)=a,
sum_iV_i^0=w/2,       sum_iV_i^1=b/2.                  (20)
```

Therefore

```text
N_0<=aw/2,             N_1<=ab/2.                      (21)
```

Also `L_0,L_1>=0`, because `w-d_i+4o_i` is the nonnegative mass of
edges avoiding `i`, plus `4o_i`. The dual coordinate identity in `(1)`
gives two further exact positive quantities:

```text
q_11-q_00
 =sum_i End_i(B)[U_i^0(B)+U_i^1(B)]>=0,

q_00-2q_01+q_11=sum_i End_i(B)^2>=0.                   (22)
```

Finally,

```text
q_11<=[sum_i U_i^1]^2=b^2/4.                           (23)
```

The coefficient floors used below are

```text
K_0^-=945756h^2+950363hw+238455w^2,
K_2^-=476914h^2+478061hw+119803w^2,
K_1^-=956122h^2+957273hw+239606w^2,

A_m^-=-(1401626h^2+1407361hw+353279w^2),
A_p^-=-(11470h^2+5755hw),
A_tau^-=-20h^2,
A_k^-=1413100h^2+1413100hw+353275w^2.                 (24)
```

Every comparison has an explicit nonnegative debt:

```text
K_0-K_0^-
 =8(L_0+L_1)+1134(aw/2-N_0)+1152(ab/2-N_1),

K_2-K_2^-=12(q_00+2q_01+q_11),
K_1-K_1^-=48(q_01+q_11),

A_m-A_m^-=4(q_00+2q_01+q_11),
A_p-A_p^-=8(q_11-q_00),
A_tau-A_tau^-=4(q_00-2q_01+q_11),
A_k-A_k^-=4(b^2-4q_11).                               (25)
```

Thus no unsigned replacement of rooted chirality is being made: the
coordinatewise endpoint identity proves precisely the delicate differences
in `(22)` and `(25)`.

## 5. The right-factor endpoint bounds

THM-4208's endpoint energy applied to `C` gives

```text
m<=[(z+2k)^2-k^2]/3=(z^2+4zk+3k^2)/3.                 (26)
```

Nonnegativity and the totals in `(15)` give

```text
p<=(z+k)k,                 tau<=k^2.                   (27)
```

Moreover `d_i<=z` and `r_i<=z`, so

```text
g_i=z-d_i-4r_i>=-4z.                                   (28)
```

The weights in `E` are nonnegative. Hence

```text
E>=-4z[h(2294h+1151w)(z+k)+8h^2k].                    (29)
```

The signs in substituting `(26)--(27)` are explicit. For example,

```text
A_m m-A_m^-[(z^2+4zk+3k^2)/3]
 =4(q_00+2q_01+q_11)m
  +(-A_m^-)[(z^2+4zk+3k^2)/3-m]>=0,                   (30)
```

and the `p` and `tau` terms have the identical form with their nonnegative
debts in `(25)`. The `A_k` term uses its lower bound directly.

Substitution of `(24)--(29)` into `(17)` and exact collection gives

```text
S(B,C)>=L(h,w;k,z),                                    (31)

L= 2[794h^2+6505hw+3065w^2]z^2/3
  + [37096h^2+62387hw+21067w^2]zk/3
  - 4(2h+w)^2k^2.                                     (32)
```

If `C=P_1`, then `(16)` is exactly zero. If `C` is non-singleton, mark each
adjacency in every Hamilton path to obtain

```text
z=W(C)>=(|C|-1)H(C)>=k.                                (33)
```

The coefficients of `z^2` and `zk` in `(32)` are positive. Therefore
`z>=k` gives

```text
L>=k^2[38636h^2+75349hw+27185w^2]/3>0.                (34)
```

Equations `(16)` and `(31)--(34)` prove Theorem 1, including its equality
boundary.

## 6. The singleton-right estimate

The right response jet of `P_1` has nonzero coordinates `2,7,10`, each
equal to one. Thus `(8)` gives the exact identity

```text
F_5(B,P_1)-10764h^2
 =2240h^2+18746hw+6122w^2+16L_1+64q_11
  -2268n_10-2304n_11.                                 (35)
```

Because all coordinates are nonnegative,

```text
2268n_10+2304n_11
 <=(sum_i U_i^1)[2268sum_iV_i^0+2304sum_iV_i^1]
 =b(1143w+1152h).                                     (36)
```

Together with `L_1>=0`, equation `(35)` yields

```text
F_5(B,P_1)-10764h^2
 >=-64h^2+15308hw+4979w^2+64q_11.                     (37)
```

For `B=P_1`, one has

```text
(h,w,L_1,q_11,n_10,n_11)=(1,0,0,1,0,1),               (38)
```

and `(35)` is equality with value zero. If `B` is non-singleton, the same
marked-Hamilton-path argument gives `w>=h`; dropping `64q_11>=0` from
`(37)` then gives

```text
F_5(B,P_1)-10764h^2>=20223h^2>0.                       (39)
```

Combining `(3)` and `(35)--(39)` proves the `n=5` case of `(4)`, with
equality only at `B=C=P_1`.

## 7. Every later tail and sharpness of the threshold

Let `n=5+r` with `r>=1`. Associativity gives

```text
A_n=A_5 triangleright P_r.                              (40)
```

Set `B'=P_r triangleright B`. Since `H(P_r)=1`, Hamilton paths factor and
`H(B')=H(B)`. Add and subtract `R_+(B',C)` to obtain the exact decomposition

```text
F_n(B,C)
 =F_5(B',C)+[R_+(P_r triangleright B,C)-R_+(B,C)]
 =F_5(B',C)+Theta(P_r,B,C).                             (41)
```

The first term is at least `10764H(B)^2H(C)^2` by the proved tail-five
case. THM-4187 Theorem 2, equation (20a), gives

```text
Theta(P_r,B,C)>=0.                                      (42)
```

This proves `(4)` for every `n>=5`. For `r>=1`, the tournament
`B'=P_r triangleright B` is not a singleton, so the first term in `(41)` is
already strictly above the lower bound by the equality case at tail five.
Thus equality in `(4)` occurs exactly when `n=5` and `B=C=P_1`.

For the converse in `(5)`, THM-4193 gives the singleton-context rows

```text
F_n(P_1,P_1)=-216,-468,-900,-1332,-180
              for n=0,1,2,3,4, respectively.            (43)
```

Hence no `n<=4` is positive in every arbitrary context. Equations
`(4)` and `(43)` prove the sharp classification `(5)`. QED.

## 8. Connection contract and scope firewall

The proved connection is

```text
source:       actual ordered factors B,C with rooted U/V states,
map:          B -> D(B)=lambda(A5 triangleright B)-lambda(B),
              C -> mu(C)-H(C)^2mu(P1),
preserved:    the exact contextual response slack S(B,C),
controlled loss:
              endpoint distributions are replaced only after their exact
              Start/End differences and energy debts have been retained,
sidecars:      mixed same-vertex U/V moments, Start, End, oriented fans,
decisive test: symbolic coefficient debts (25) and the literal audit below.
```

This theorem proves the cycle-first transitive-tail family `A_n` exactly
from `n=5` onward. It does not prove general `(OS+)`, the no-sink or
no-source gate law, or the order-eleven asymmetric bank. The propagation in
`(41)` uses a transitive block `P_r`; it is not an arbitrary-prefix theorem.

## 9. Replay

Exact symbolic audit:

```bash
python3 -B 04-computation/tournament_a5_uniform_context_bound_thm4212.py
python3 -O -B 04-computation/tournament_a5_uniform_context_bound_thm4212.py
```

Both streams byte-match the frozen primary output.

Independent clean-room literal audit:

```bash
clang++ -std=c++17 -O0 -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_a5_uniform_context_bound_independent_audit_thm4212.cpp \
  -o /tmp/thm4212-independent-O0
/tmp/thm4212-independent-O0

clang++ -std=c++17 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_a5_uniform_context_bound_independent_audit_thm4212.cpp \
  -o /tmp/thm4212-independent-O3
/tmp/thm4212-independent-O3

clang++ -std=c++17 -O1 -g -fsanitize=address,undefined \
  -fno-omit-frame-pointer -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_a5_uniform_context_bound_independent_audit_thm4212.cpp \
  -o /tmp/thm4212-independent-san
ASAN_OPTIONS=detect_leaks=0 /tmp/thm4212-independent-san
```

The O0, O3, and ASan/UBSan streams byte-match the frozen independent output.
The literal universe is every labelled pair `B,C` with factor orders at most
three and `|B|+|C|<=4`, totaling 25 rows. It rebuilds composite tournaments
of order at most twelve without ordinal capacity transfer, finds the unique
equality `(P_1,P_1)`, and checks the signed controls

```text
F_5(P_1,P_1)=10764,
F_5(P_2,P_1)=68364,
F_5(P_1,P_2)=79128.
```

**QED for `(3)--(5)`, with equality exactly at the stated singleton
boundary.**
