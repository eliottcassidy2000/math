---
id: THM-4215
title: "Adjacent-cycle sharp uniform-context threshold seven"
status: >
  RESERVED / COMPLETE PROVISIONAL PROOF UNDER INDEPENDENT THEOREM AUDIT.
  No statement in this file is yet a proved dependency. The candidate proof
  gives the exact singleton formula for C3 triangleright C3 triangleright P_n;
  negative singleton witnesses for every n<=6; the sharp all-context bound
  F_Z7(B,C)>=967788H(B)^2H(C)^2 with equality only for singleton contexts;
  and, by the exact ordinal telescope and neutral transitive-tail propagation,
  uniform positivity exactly for n>=7. Primary symbolic and clean-room literal
  replay streams pass; status promotion awaits an external theorem referee.
source: adjacent-cycle-threshold-20260826
related:
  - THM-4187-universal-source-padding-supermodularity-and-transitive-left-ordinal-positivity
  - THM-4208-cycle-prefix-arbitrary-context-recurrence-endpoint-energy-and-eventual-positivity
  - THM-4212-cycle-prefix-uniform-arbitrary-context-threshold-and-tail-five-lower-bound
  - THM-4213-uniform-prefix-ordinal-semigroup-and-tail-five-cycle-language
script: 04-computation/tournament_adjacent_cycle_threshold_thm4215.py
output: 05-knowledge/results/tournament_adjacent_cycle_threshold_thm4215.out
independent_audit_script: 04-computation/tournament_adjacent_cycle_threshold_independent_audit_thm4215.cpp
independent_audit_output: 05-knowledge/results/tournament_adjacent_cycle_threshold_independent_audit_thm4215.out
script_sha256: 5b18e3cd48524df0342197dbfda43d1625a6255939a5f4cde7b0a316d19f7281
output_sha256: 3d0f912f778bfcf63acc3a3baa87adb9c9c16addd6cd9264250b886d5dc2d831
independent_audit_script_sha256: 8ffee22c9b1f395984ae75ebb85a23f6bb6a018862e60cd3d67c68160d1d4c90
independent_audit_output_sha256: d455e574d4cb2686baf690db189f974a2f83b4c7648069f44c23ae03db5de06b
hash_basis: raw LF bytes
primary_audit: >
  PASS. Exact ordinal reconstruction supplies the Z7 prefix coordinates and
  singleton boundary rows. A symbolic derivation then reconstructs all eleven
  coordinates of lambda(Z7 triangleright B)-lambda(B), verifies every
  coefficient debt, and proves the right-context and singleton-middle bounds.
  Normal and python -O streams byte-match.
independent_audit: >
  PASS COMPUTATIONAL REFEREE; EXTERNAL THEOREM AUDIT PENDING. A standalone
  warning-clean C++17 engine imports no response jet, endpoint inequality, or
  ordinal capacity transfer. From literal labelled adjacency and subset path
  DP it rebuilds Hamilton counts, odd-path capacities, G_+, all eight
  singleton boundary rows, the Z7 invariants, and five representative context
  rows including both one-sided C3 controls. Clang O0/O3 and ASan/UBSan streams
  byte-match.
---

# THM-4215 -- adjacent-cycle sharp uniform-context threshold seven

**RESERVED / COMPLETE PROVISIONAL PROOF UNDER INDEPENDENT THEOREM AUDIT.**

No statement below enters the proved dependency graph until the frontmatter
is promoted after external theorem audit.

THM-4213 found that a uniformly positive suffix does not by itself preserve
positivity:

```text
C3 triangleright C3 triangleright P_5
```

has singleton contextual defect `-338580`.  This proof candidate determines
the exact repair.  Two singleton tail vertices are both necessary and
sufficient: the sharp adjacent-cycle threshold is seven.

The mechanism is again finite but not finite-order.  At `Z_7`, the exact left
response difference lies in the dual cone cut out by the right endpoint-energy
inequality.  Its singleton-right face has one exact minimum.  The ordinal
telescope then propagates the certificate through every later transitive tail.

## 1. Statement

All quantified tournaments are finite and nonempty.  Use THM-4208/4213's
ordinal remainder `R_+` and contextual prefix defect

```text
F_A(B,C)=R_+(A triangleright B,C)-R_+(B,C).             (1)
```

For `n>=0`, put

```text
Z_n=C3 triangleright C3 triangleright P_n,              (2)
```

where `P_0` means that the final factor is omitted.

> **Candidate Theorem 1 (sharp adjacent-cycle threshold and floor).** For
> every `n>=7` and all nonempty tournaments `B,C`,
>
> ```text
> F_(Z_n)(B,C)>=967788H(B)^2H(C)^2>0.                   (3)
> ```
>
> Equality in the first inequality holds if and only if
>
> ```text
> n=7,                    B=C=P_1.                      (4)
> ```
>
> Moreover,
>
> ```text
> F_(Z_n)(B,C)>0 for every nonempty B,C    iff n>=7.    (5)
> ```

The coefficient `967788` is therefore the optimal Hamilton-normalized floor
at the crossing `n=7`.  For every `n>7`, inequality `(3)` is strict; no claim
is made here that `967788` is the optimal floor at those later tails.

As an immediate `(OS+)` consequence, `Z_n` is a left factor with positive
ordinal remainder against every no-sink tournament whenever `n>=8`.

## 2. Exact `Z_7` response difference

For the middle tournament `B`, use the notation

```text
h=H(B),                    w=W(B),
s=(s_0,s_1)=(w/2,w/2+h),
L_e=L_e^+(B),
q_ef=sum_i U_i^e(B)U_i^f(B),
n_ef=sum_i U_i^e(B)V_i^f(B).                           (6)
```

The mixed same-vertex moments `n_ef` are derived coordinates, not hypotheses.
Exact ordinal reconstruction of `A=Z_7` gives

```text
H(A)=9,                   W(A)=18414,
s_A=(9207,9216),

Q_U(A)=((17031141,17031141),
        (17031141,17031222)),

L^+(A)=(271372032,271454814).                          (7)
```

Let `lambda` be the eleven-coordinate left response jet in THM-4208 equation
(9), and put

```text
D(B)=lambda(A triangleright B)-lambda(B)=(D_0,...,D_10). (8)
```

> **Candidate Lemma 2 (exact adjacent-cycle response difference).** The
> coordinates in `(8)` are
>
> ```text
> D_0=h(331614h+165887w),
>
> D_1=2[1086773760h^2+1087499358hw+80L_0
>       -165726n_00-165888n_01+272056239w^2],
>
> D_2=2[1087105374h^2+1087665165hw+80L_1
>       -165726n_10-165888n_11+272056239w^2],
>
> D_3=h(331614h+165887w),
> D_4=165887h(2h+w),
>
> D_5=1087561728h^2+1087893342hw+480q_00+272056279w^2,
>
> D_6=2[1087893342h^2+1088059229hw+480q_01
>       +272056279w^2],
>
> D_7=1088225116h^2+1088225116hw+480q_11+272056279w^2,
>
> D_8=5[-651564000h^2-651895614hw+32q_00-163056847w^2],
>
> D_9=10[-651895614h^2-652061501hw+32q_01-163056847w^2],
>
> D_10=5[-652227388h^2-652227388hw+32q_11-163056847w^2]. (9)
> ```

### Derivation

Put `a=w+h` and `b=w+2h`, and let `T=A triangleright B`.  The exact ordinal
transfer identities give

```text
H(T)=9h,                  W(T)=18414b+18a,
s_T=(W(T)/2,W(T)/2+9h),

Q_U(T)=81Q_U(B)
       +(b^2/4)(q^A_00+2q^A_01+q^A_11)((1,1),(1,1)).   (10)
```

For an `A` vertex, its `U` state becomes balanced after convolution with the
optional state of `B`; a `B`-block `U` state scales by nine.  The corresponding
two blockwise linear moments are

```text
L_e^+(T)=C_A+81L_e+9[W(T)-9w]s_e
          -18(9207n_(e0)+9216n_(e1)),                 (11)

C_A=(b/2){h(L_0^+(A)+L_1^+(A))
     +(18414+18)a(9207+9216)
     +3X_A},

X_A=2[(q^A_00+q^A_01)s_0+(q^A_01+q^A_11)s_1].         (12)
```

Substitution of `(7)` and `(10)--(12)` into `lambda(T)-lambda(B)` expands to
`(9)`.  The primary audit checks all eleven polynomial identities exactly;
there is no interpolation.

THM-4208's response factorization now gives

```text
F_(Z_7)(B,C)=D(B) dot mu(C).                            (13)
```

## 3. Every right context is no worse than the singleton

For the right tournament `C`, write

```text
k=H(C),                    z=W(C),
v_i=V_i^0(C)+V_i^1(C),
t_i=V_i^1(C)-V_i^0(C)=Start_i(C),
g_i=z-d_i(C)-4r_i(C),

m=sum_i v_i^2,            p=sum_i v_it_i,
tau=sum_i t_i^2.                                          (14)
```

Thus `v_i,t_i>=0`, `sum v_i=z+k`, and `sum t_i=k`.  Define the
right-context slack

```text
S(B,C)=F_(Z_7)(B,C)-k^2F_(Z_7)(B,P_1).                 (15)
```

Changing from `(V_i^0,V_i^1)` to `(v_i,t_i)` in `(13)` gives the exact form

```text
S=K_0zk+E+K_2z^2+K_1zk+A_m m+A_p p+A_tau tau+A_k k^2, (16)

E=sum_i [h(331694h+165887w)v_i+80h^2t_i]g_i,           (17)
```

where

```text
K_0=2[40L_0+40L_1+1087105374h^2+1087665205hw
      -82863(n_00+n_10)-82944(n_01+n_11)+272056239w^2],

K_2=1087893382h^2+1088059229hw
    +120q_00+240q_01+120q_11+272056279w^2,

K_1=2176118458h^2+2176284345hw
    +480q_01+480q_11+544112558w^2,

A_m=5[-651895654h^2-652061501hw
      +8q_00+16q_01+8q_11-163056847w^2],

A_p=-5[331694h^2+165887hw+16q_00-16q_11],

A_tau=40[-5h^2+q_00-2q_01+q_11],

A_k=5[652227388h^2+652227388hw-32q_11+163056847w^2].  (18)
```

### Left-factor coefficient debts

Put

```text
N_0=n_00+n_10,            N_1=n_01+n_11.               (19)
```

Nonnegativity and the exact rooted totals give

```text
N_0<=a w/2,               N_1<=a b/2.                  (20)
```

THM-4208's coordinatewise endpoint identity gives

```text
q_11-q_00>=0,
q_00-2q_01+q_11>=0,
q_11<=b^2/4.                                            (21)
```

Use the coefficient floors

```text
K_0^-=2174044860h^2+2174998715hw+543946671w^2,
K_2^-=1087893382h^2+1088059229hw+272056279w^2,
K_1^-=(2h+w)(1088059229h+544112558w),

A_m^-=-5[651895654h^2+652061501hw+163056847w^2],
A_p^-=-5h(331694h+165887w),
A_tau^-=-200h^2,
A_k^-=815284195(2h+w)^2.                               (22)
```

Every comparison in `(22)` has an explicit nonnegative debt:

```text
K_0-K_0^-
 =80(L_0+L_1)+165726(aw/2-N_0)+165888(ab/2-N_1),

K_2-K_2^-=120(q_00+2q_01+q_11),
K_1-K_1^-=480(q_01+q_11),
A_m-A_m^-=40(q_00+2q_01+q_11),
A_p-A_p^-=80(q_11-q_00),
A_tau-A_tau^-=40(q_00-2q_01+q_11),
A_k-A_k^-=40(b^2-4q_11).                              (23)
```

Here `L_0,L_1>=0`: their weights are the nonnegative capacity mass avoiding
the marked vertex, plus four times its outgoing mass.  Thus `(20)--(23)`
justify every floor in the direction used below.

### Right-factor endpoint bounds

THM-4208's endpoint energy and elementary product bounds give

```text
m<=(z^2+4zk+3k^2)/3,
p<=(z+k)k,
tau<=k^2.                                               (24)
```

Also `g_i>=-4z`.  The weights in `(17)` are nonnegative, so

```text
E>=-4z{h(331694h+165887w)(z+k)+80h^2k}.                (25)
```

Substitute `(22)--(25)` into `(16)`.  The negative coefficients
`A_m^-,A_p^-,A_tau^-` multiply their upper bounds, while the other terms use
their lower bounds.  Exact collection gives

```text
S(B,C)>=L(h,w;k,z),                                    (26)

3L=[221548h^2+1879538hw+884602w^2]z^2
   +[3620176h^2+8140211hw+3040747w^2]zk
   -120(2h+w)^2k^2.                                   (27)
```

If `C=P_1`, equation `(15)` is exactly zero.  If `C` is non-singleton, marking
each adjacency in every Hamilton path gives

```text
z=W(C)>=(|C|-1)H(C)>=k.                                (28)
```

The `z^2` and `zk` coefficients in `(27)` are positive.  Setting `z=k` at
the lower boundary gives

```text
3L>=k^2[3841244h^2+10019269hw+3925229w^2]>0.           (29)
```

Equations `(15)` and `(26)--(29)` prove

```text
F_(Z_7)(B,C)>=H(C)^2F_(Z_7)(B,P_1),                    (30)
```

with equality if and only if `C=P_1`.

## 4. The sharp singleton-right face

The right response jet of `P_1` has nonzero coordinates `2,7,10`, each equal
to one.  Equation `(9)` therefore gives

```text
F_(Z_7)(B,P_1)
 =2[649462h^2+1209253hw+80L_1
    -165726n_10-165888n_11+320q_11+442261w^2].         (31)
```

The mixed same-vertex terms satisfy

```text
4*9*9207n_10+4*9*9216n_11
 <=9b(9207w+9216b).                                    (32)
```

Indeed the left side is a same-vertex product of nonnegative `U^1` and
weighted `V` coordinates; it is at most the product of their totals
`b/2`, `w/2`, and `b/2`.  Dropping `L_1>=0` from `(31)` and using `(32)` gives

```text
F_(Z_7)(B,P_1)
 >=967148h^2+1921004hw+718715w^2+640q_11.              (33)
```

For `B=P_1`, the exact coordinates are

```text
(h,w,L_1,q_11,n_10,n_11)=(1,0,0,1,0,1),               (34)
```

and `(31)` equals `967788`.  If `B` is non-singleton, the same marked-path
argument gives `w>=h`; hence `(33)` implies

```text
F_(Z_7)(B,P_1)-967788h^2
 >=-640h^2+1921004hw+718715w^2+640q_11
 >=2639079h^2>0.                                       (35)
```

Combining `(30)` and `(33)--(35)` proves the `n=7` case of `(3)`, including
the unique equality `B=C=P_1`.

## 5. Exact hostile rows and propagation to every later tail

Exact response-jet substitution gives the singleton row for every `n>=0`:

```text
F_(Z_n)(P_1,P_1)
 =108[2*4^n-(12n+102)2^n+1].                           (36)
```

At the eight boundary values this is

```text
n:       0       1       2        3        4        5        6       7
F:  -10692  -23652  -50868  -105300  -203796  -338580  -317844  967788.
                                                                    (37)
```

Thus every `n<=6` fails uniform positivity already at the singleton context
pair.  The nonmonotone hostile at `n=6` is load-bearing: extrapolation from
tail five would give the wrong crossing mechanism.

Now let `n=7+r` with `r>=1`.  Associativity gives

```text
Z_n=Z_7 triangleright P_r.                              (38)
```

THM-4213's exact prefix telescope, with `B'=P_r triangleright B`, gives

```text
F_(Z_n)(B,C)
 =F_(Z_7)(B',C)+F_(P_r)(B,C).                          (39)
```

Hamilton paths factor at the ordinal cut, so `H(B')=H(B)`.  The first term
in `(39)` obeys the proved-candidate `Z_7` floor.  THM-4187 gives

```text
F_(P_r)(B,C)>=0.                                       (40)
```

Therefore `(3)` holds for all `n>=7`.  If `r>=1`, the middle context `B'` is
not a singleton, so `(35)` and `(30)` make the first term strictly larger
than the floor.  Equality in `(3)` is consequently exactly `(4)`.  Together
with `(37)`, this proves the sharp equivalence `(5)`.

Finally, if `n>=8`, write `Z_n=Z_(n-1) triangleright P_1`.  For every no-sink
`C`, equations `(1)` and `(3)` give

```text
R_+(Z_n,C)
 =F_(Z_(n-1))(P_1,C)+R_+(P_1,C)>0,                    (41)
```

because THM-4187 makes the second term positive.  This proves the stated
adjacent-cycle `(OS+)` corollary.

## 6. Connection contract and scope firewall

The proved-candidate connection is

```text
source:       actual adjacent-cycle prefix and ordered context pair,
map:          B -> lambda(Z7 triangleright B)-lambda(B),
              C -> mu(C)-H(C)^2mu(P1),
preserved:    exact contextual response and Hamilton-normalized floor,
controlled loss:
              endpoint distributions are bounded only after retaining their
              exact Start/End differences and energy debts,
sidecars:      mixed same-vertex U/V moments, oriented capacity fans,
              Hamilton endpoints, and ordered strong-component positions,
decisive tests:
              the seven negative singleton rows, literal n=6/7 crossing,
              and one-sided nontransitive C3 context controls.
```

Once promoted, `(3)` says `Z_n` lies in THM-4213's uniformly positive ideal
for exactly the tails `n>=7`.  It can therefore serve as a new certified block
inside any ordinal word whose other factors lie in the nonnegative semigroup.

This candidate does not classify that semigroup or ideal; does not determine
the threshold for three or more consecutive cycle components; and does not
prove `(OS+)` for `Z_7`, arbitrary source-free tournaments, the all-strong
residual, the no-sink/no-source gate law, or the order-eleven asymmetric bank.
General `(OS+)` remains **OPEN**.

## 7. Replay

Primary exact symbolic audit:

```bash
python3 -B 04-computation/tournament_adjacent_cycle_threshold_thm4215.py
python3 -O -B 04-computation/tournament_adjacent_cycle_threshold_thm4215.py
```

Both streams byte-match the frozen primary output.

Independent clean-room literal audit:

```bash
clang++ -std=c++17 -O0 -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_adjacent_cycle_threshold_independent_audit_thm4215.cpp \
  -o /tmp/thm4215-independent-O0
/tmp/thm4215-independent-O0

clang++ -std=c++17 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_adjacent_cycle_threshold_independent_audit_thm4215.cpp \
  -o /tmp/thm4215-independent-O3
/tmp/thm4215-independent-O3

clang++ -std=c++17 -O1 -g -fsanitize=address,undefined \
  -fno-omit-frame-pointer -Wall -Wextra -Wpedantic -Werror \
  04-computation/tournament_adjacent_cycle_threshold_independent_audit_thm4215.cpp \
  -o /tmp/thm4215-independent-san
ASAN_OPTIONS=detect_leaks=0 /tmp/thm4215-independent-san
```

The O0, O3, and ASan/UBSan streams byte-match the frozen independent output.
The independent engine builds every composite from labelled adjacency and
recomputes Hamilton paths, odd-directed-path capacities, `G_+`, `R_+`, and
`F` without ordinal capacity transfer or response jets.

**Candidate QED for `(3)--(5)` and `(41)`; promotion awaits independent
theorem audit.**
