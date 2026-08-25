---
id: HYP-9080
title: "Tournament deletion-slack local unavoidability"
status: >
  OPEN CONJECTURE + PROVED REDUCTION IDENTITY + FINITE-EXACT EVIDENCE, NOT
  PROVED. For a tournament C and vertex x, put S=H-disc and
  chi_x(C)=2(S(C)-S(C-x)). THM-3729 and THM-4094 identify this exactly with
  a rooted boundary energy. OCF and principal-minor expansion give exact local
  and averaged support-moment formulas. Every nonconstant ear over every
  labelled order-6 base and the inherited 456-class order-7 bank has chi>=0,
  but a source over C3 has chi=-2.
  Existential/averaged unavoidability and H>=disc remain OPEN.
source: codex-snark-apex-260822870-20260825
depends_on:
  - THM-002-ocf
  - THM-3729-rooted-pfaffian-response-and-sign-root-deletion-average
  - THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness
related:
  - THM-1950-h-ge-disc-reduced-to-strongly-connected
  - THM-2290-context-selected-colored-pair-kernel-is-hafnian-complete
  - THM-4099-squarefree-gap-transfer-and-mixed-insertion-boundary
  - THM-4104-selected-order-eleven-strong-ear-solid-interval
  - THM-4111-uniform-ear-average-and-recursive-selected-bank-growth
  - THM-4113-maximal-noncrossing-half-kempe-atlas
  - THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature
  - THM-4115-uniform-ear-cut-walsh-variance-and-sharp-growth-refinement
  - THM-4118-ear-response-lattice-and-stateful-unit-component-intervals
script: 04-computation/tournament_deletion_slack_unavoidability_hyp9080.py
output: 05-knowledge/results/tournament_deletion_slack_unavoidability_hyp9080.out
script_sha256: 6ff91dacc06b63c7e7fcd34163ca647e35a814f683eb0aa9fdb5a8ebd54000e6
output_sha256: 0c41b86d3a445df8fcfcf2e2a49e7363f9410766de7c495a6ac02ac83610549b
auxiliary_script: 04-computation/tournament_alu_support_moment_ear_tariff_hyp9080.py
auxiliary_output: 05-knowledge/results/tournament_alu_support_moment_ear_tariff_hyp9080.out
auxiliary_script_sha256: 870403672f12654b6f1372af19978f21e95aa3a057f68ad48ca2716acdeff9f9
auxiliary_output_sha256: 41ffbeb17b83ff9ce087a7170d54c9ea44a292450e7b7c1c8bd458e51d7744a1
hash_basis: raw working-tree bytes; normalized-LF hash for the inherited n=7 representative bank
audit: >
  FINITE-EXACT PASS. The script exhausts all 33,866 labelled tournaments of
  orders 2 through 6 and the inherited 456-representative order-7 bank
  (353 strong classes). It computes H for every deletion by subset DP and
  separately evaluates disc and the rooted response with exact Bareiss
  determinants; these are formula cross-checks sharing one determinant
  kernel, not independent implementations. It also reconstructs full
  deletion-incidence fibers through order five. It checks
  410,408 rooted-response gates and 5,404 incidence gates. Normal and
  optimized runs reproduce the frozen transcript. Local unavoidability,
  its averaged strengthening, and H>=disc remain OPEN.
  The auxiliary companion reconstructs OCF support atoms by Boolean Mobius
  inversion and Pfaffian-square atoms by principal determinants on 34,322
  objects. It checks 273,848 moment gates, 1,096 direct ear formulas, and all
  2,031,616 nonconstant ears over labelled order-six bases plus 57,456 ears
  over the inherited order-seven representative bank. Its normal and
  optimized transcripts reproduce the frozen output. The combined ear tariff
  remains OPEN.
---

# HYP-9080 -- tournament deletion-slack local unavoidability

**OPEN + FINITE-EXACT.** The apex-cubic proof has two logically separate
halves: local reducibility and global unavoidability. The same shape gives a
precise proposal for the repo's open tournament inequality `H>=disc`, but no
planar theorem transfers.

## 1. Exact deletion charge

Let `C` be a finite tournament of order at least two, let `x` be a vertex,
and put

```text
T=C-x,                    S(W)=H(W)-disc(W).                 (1)
```

Let `u_x in {+1,-1}^V(T)` be the **actual incident sign root** of `x`; in
skew-matrix notation its coordinate at `v` is `K_C(v,x)`. Write

```text
Delta H_x=H(C)-H(T).                                        (2)
```

THM-4094 proves `Delta H_x>=0` and computes it from all insertion fibers and
orphan paths. THM-3729, equation (8), gives

```text
E_odd(T;u_x)=2 disc(C)-disc(T).                              (3)
```

Define the local deletion charge

```text
chi_x(C)=2 Delta H_x-E_odd(T;u_x)+disc(T).                  (4)
```

Substituting `(2)--(3)` yields the exact identity

```text
boxed:
chi_x(C)=2[(H-disc)(C)-(H-disc)(C-x)]
        =2[S(C)-S(C-x)].                                    (5)
```

Thus `chi_x` is an even integer. The sign condition

```text
chi_x(C)>=0  iff  S(C-x)<=S(C)                              (6)
```

is a genuine size-reducing certificate: deleting `x` does not increase the
slack which the open inequality seeks to keep nonnegative.

The rooted sign is load-bearing. Resetting it to the all-one vector changes
the energy in general by THM-3729. Likewise, a selected injection of old
Hamiltonian paths into new paths retains only `Delta H_x>=0`; THM-4094's
transitive-triple/`C_3` hostile proves that it loses the exact increment.

## 2. Local-unavoidability conjecture

> **Conjecture LU.** Every finite tournament `C` with `|C|>=2` has a vertex
> `x` such that
>
> ```text
> chi_x(C)>=0.                                               (7)
> ```

Conjecture LU implies the open inequality `H(C)>=disc(C)` for all finite
tournaments. The one-vertex base has `S=0`. Given `C`, choose `x` by `(7)`;
induction gives `S(C-x)>=0`, while `(6)` gives

```text
S(C)>=S(C-x)>=0.                                             (8)
```

Equivalently, a smallest counterexample to `H>=disc` would have

```text
chi_x(C)<0                 for every vertex x,               (9)
```

because all proper deletions would have nonnegative slack while `S(C)<0`.

LU is stronger than the desired inequality; no converse is known. A
nonnegative object may conceivably have every deletion with still larger
slack. THM-1950 reduces `H>=disc` to strong tournaments, but that theorem does
not prove `(7)` even on the strong stratum.

## 3. Averaged global-charge strengthening

Summing `(5)` over all vertices gives another exact identity:

```text
1/2 sum_x chi_x(C)
 = |C| S(C)-sum_x S(C-x).                                  (10)
```

This motivates the stronger conjecture

> **Conjecture ALU.** For every finite tournament `C`,
>
> ```text
> |C|S(C)-sum_xS(C-x)>=0.                                  (11)
> ```

ALU implies LU by averaging. It is the closest exact analogue to a
discharging total: the left side is a global charge, while one nonnegative
summand supplies a reducible deletion. It remains open. In particular, the
positive total charge of the apex-cubic paper supplies no inequality for
`(10)`.

### 3.1 Exact OCF/Pfaffian support-moment form

The averaged charge has a second exact expression which identifies where a
proof must act. Let `Gamma` range over collections of pairwise
vertex-disjoint directed odd cycles of `C`, including the empty collection,
and put

```text
U(Gamma)=the union of their vertex supports,
weight(Gamma)=2^|Gamma|.                                  (11a)
```

For every even vertex set `A`, put

```text
p_A=Pf(K_C[A])^2=det(K_C[A]),        p_empty=1.            (11b)
```

Then

```text
boxed:
|C|S(C)-sum_x S(C-x)
 =sum_Gamma |U(Gamma)| 2^|Gamma|
  -2^(1-|C|) sum_(A even)(2|A|-|C|)p_A.                   (11c)
```

Indeed, THM-002 gives `H(C)=sum_Gamma 2^|Gamma|`. A fixed collection
survives exactly `|C|-|U(Gamma)|` vertex deletions, so

```text
|C|H(C)-sum_xH(C-x)=sum_Gamma |U(Gamma)|2^|Gamma|.         (11d)
```

The principal-minor expansion gives

```text
disc(C)=2^(1-|C|)sum_(A even)p_A.                         (11e)
```

The same minor occurs in `disc(C-x)` exactly when `x notin A`, with twice
the normalization in `(11e)`. Hence

```text
|C|disc(C)-sum_xdisc(C-x)
 =2^(1-|C|)sum_(A even)(2|A|-|C|)p_A,                    (11f)
```

and subtracting `(11f)` from `(11d)` proves `(11c)`.

There is also an exact local form:

```text
chi_x(C)/2
 =sum_(Gamma:x in U(Gamma))2^|Gamma|
  -2^(1-|C|)(sum_(A even:x in A)p_A
             -sum_(A even:x notin A)p_A).                 (11g)
```

Thus ALU is a comparison between two **support moments**, not between the
unweighted numbers of cycles and minors. Only Pfaffian supports with
`|A|>|C|/2` have an adverse coefficient in `(11c)`; smaller supports help.
OCF/Mobius positivity proves that the first-family atoms are nonnegative, but
it gives no cross-family inequality against the large Pfaffian supports.
This is the first missing implication, so even the proved and independently
audited THM-4114 presence-cube positivity is not a proof of ALU: that theorem
controls the OCF atoms internally, not their comparison with the adverse
large-Pfaffian support moment.

### 3.2 Exact combined ear-cut tariff

The local charge also has a precise cut formulation. Let `T` be a tournament
with skew matrix `K`, put

```text
d=disc(T),                  R=(I-K^2)^(-1).                 (11h)
```

`I-K^2` is positive definite. Moreover `R` is the symmetric part of
`(I+K)^(-1)`. For `S subseteq V(T)`, adjoin an ear `x_S` with
`x_S->v` exactly when `v in S`; its incident root is

```text
u_S=1-2 1_S.                                               (11i)
```

THM-3729 gives the discriminant response

```text
g(S):=disc(T+x_S)
 =d/2(1+u_S^T R u_S)
 =g(empty)-2d cut_R(S).                                   (11j)
```

Use THM-4104's ear data

```text
w_ij=(Q(i,j)+Q(j,i))/2,
h_i=Start(i)-End(i)+(col_i(Q)-row_i(Q))/2,
sum_i h_i=0.                                               (11k)
```

Then

```text
H(T+x_S)-H(T)=cut_w(S)+h(S),       h(S)=sum_(i in S)h_i.   (11l)
```

Put

```text
a_0=g(empty)-d,                 wtilde=w+2dR.              (11m)
```

Combining `(11j)--(11m)` with `(5)` gives the exact paired formulas

```text
boxed:
chi_(x_S)(T+x_S)/2       =cut_wtilde(S)+h(S)-a_0,
chi_(x_(V-S))(T+x_(V-S))/2
                           =cut_wtilde(S)-h(S)-a_0.         (11n)
```

This isolates a stronger open tariff proposal:

> **Combined ear-tariff conjecture.** For every finite `T` and every
> nonempty proper `S subset V(T)`,
>
> ```text
> cut_wtilde(S)>=a_0+|h(S)|.                              (11o)
> ```

Equation `(11o)` would make the deletion charge nonnegative for both mixed
ear orientations. In particular, every vertex of a strong tournament has a
nonconstant incident cut, so `(11o)` would be much stronger than LU on the
strong stratum.

The signs cannot be read from THM-4104 alone. Although `w_ij>=0` and `R` is
positive definite, the off-diagonal entries of `R` have no general sign
control; positive definiteness does not imply `wtilde_ij>=0` or the tariff
`(11o)`. This is the exact presence-cube/ear-cube firewall left by THM-4114:
its proved directed-cut submodularity does not supply the missing
cross-family domination.

### 3.3 Exact Walsh lift of the combined tariff

THM-4115 restores the Walsh coordinates of the Hamiltonian ear cube. The
same transform applies exactly to the combined charge. Put

```text
Y_S=chi_(x_S)(T+x_S)/2,
epsilon_i(S)=+1 if i in S and -1 otherwise,
Wtilde=sum_(i<j)wtilde_ij.                                (11p)
```

Since `sum_i h_i=0`, equation `(11n)` becomes

```text
Y_S=Wtilde/2-a_0
    +(1/2)sum_i h_i epsilon_i
    -(1/2)sum_(i<j)wtilde_ij epsilon_i epsilon_j.          (11q)
```

Walsh orthogonality therefore gives the all-order identities

```text
E_S Y_S=Wtilde/2-a_0,
Var_S(Y_S)=1/4(sum_i h_i^2+sum_(i<j)wtilde_ij^2).          (11r)
```

This is a proved degree-two carrier, but it is not a tariff proof. Unlike
THM-4115's Hamiltonian weights, `wtilde=w+2dR` is not known coefficientwise
nonnegative at all orders, and mean plus variance does not control the
minimum of a labelled Boolean quadratic. Indeed

```text
min(Y_S,Y_(V-S))=cut_wtilde(S)-a_0-|h(S)|,                (11s)
```

so `(11o)` is exactly a nonconstant-cube support-floor problem. A cheap
necessary condition follows: if `(11o)` holds, then only the two constant
cuts can be negative, both equal to `-a_0`, and hence

```text
E_S Y_S >= -a_0/2^(|T|-1).                               (11t)
```

Neither `(11r)` nor `(11t)` is sufficient. The full labelled quadratic,
rather than its first two scalar moments, is the faithful next carrier.

There is nevertheless a sharp maximum consequence. Let `n=|T|`, put
`W=sum_(i<j)w_ij`, and use THM-4115's one-defect count `F_1(T)`. The diagonal
cofactor identity `disc(T-v)=2dR_vv` gives

```text
mu_Y:=E_S Y_S
 =W/2+d(1-tr R)/2
 =((n-1)H(T)+F_1(T))/4+d/2-(1/4)sum_v disc(T-v).          (11u)
```

Since `-K^2=K^T K` is positive semidefinite, `0<R<=I` in Loewner order.
Thus `u_S^T R u_S<=n`; Hamiltonian insertion never decreases `H`, and so

```text
Y_S>=L:=-d(n-1)/2.                                       (11v)
```

Writing `B=sum_i h_i^2+sum_(i<j)wtilde_ij^2` and
`M_Y=max_S Y_S`, average `(Y_S-L)(M_Y-Y_S)>=0`. Because
`mu_Y-L=(W+d(n-tr R))/2>0`, this proves

```text
M_Y>=mu_Y+B/(2(W+d(n-tr R))).                             (11w)
```

Equality holds for the `C3` base: `L=-1`, `mu_Y=5/4`, `Var(Y)=27/16`, and
`M_Y=2`. This can force a favorable **extension orientation**; it does not
control a fixed child's deletion charge, whose base changes with the deleted
vertex. Hence it proves neither LU nor ALU.

One scalar obstruction to the tariff also survives. Put

```text
Z_S=(Y_S+Y_(V-S))/2=cut_wtilde(S)-a_0,
Lh_S=(Y_S-Y_(V-S))/2=h(S).                                (11x)
```

The tariff implies `Z_S^2>=Lh_S^2` on every nonconstant cut, while the two
constant cuts contribute `a_0^2` each. Therefore it necessarily implies

```text
mu_Y^2+(sum_(i<j)wtilde_ij^2-sum_i h_i^2)/4
 >=2^(1-n)a_0^2.                                         (11y)
```

Violation of `(11y)` would refute the tariff. Satisfaction is not sufficient
because the scalar moment still forgets labelled incidence.

### 3.4 Stateful response lattice imported from THM-4118

THM-4118's proof is carrier-level and therefore applies to the combined
integer-valued quadratic `Y`, but its Hamiltonian parity conclusion does not.
For `n>=3`, define

```text
A_i^Y=Y_{ {i} }-Y_empty,
K_ij^Y=Y_{ {i} }+Y_{ {j} }-Y_{ {i,j} }-Y_empty,
delta_Y=gcd({A_i^Y} union {K_ij^Y}).                       (11z)
```

The degree-two identity `(11q)` gives

```text
Y_S=Y_empty+sum_(i in S)A_i^Y
              -sum_({i,j} subseteq S)K_ij^Y,
<Y_S-Y_empty: empty!=S!=V>_Z=delta_Y Z.                  (11aa)
```

If every coefficient in `(11z)` vanishes, `Y` is constant and the tariff is
decided directly. Otherwise make a graph on the nonconstant cuts, joining a
one-flip or one-in/one-out exchange exactly when its `Y` change is
`0` or `+/-delta_Y`. Every connected component maps onto a solid arithmetic
interval with step `delta_Y`, by the same discrete intermediate-value proof
as THM-4118.

This is the first stateful sidecar that genuinely survives the transfer. It
does not prove all component minima nonnegative, and THM-4118's Hamiltonian
fact `d_T=2` must not be copied to `delta_Y`. The decisive next experiment is
to classify the combined-charge unit components on the zero-gap order-seven
representatives and ask which labelled coefficient or Pfaffian sidecar
separates components that reach below zero from those that do not.

## 4. Pointwise monotonicity fails at order four

Use upper-pair bit order

```text
(01,02,03,12,13,23),       bit 1 meaning i->j.               (12)
```

At mask `2`, vertex `3` is a source over the directed triangle

```text
0->2->1->0.                                                   (13)
```

The full tournament and its source deletion have

```text
(H,disc,S)(C)   =(3,2,1),
(H,disc,S)(C-3) =(3,1,2),
chi(C)          =(2,2,2,-2).                                (14)
```

Hence the tempting statement `chi_x>=0` for every vertex is false. The
strongest survivor is existential unavoidability: the other three vertices
in `(14)` still have positive charge.

As an ear over the fixed `C3` base, the separation is even sharper. The two
constant source/sink cuts have `chi=-2`, while each of the six nonconstant
cuts has `chi=4`. Thus nonnegative OCF presence differences and the
submodular Hamiltonian ear response can coexist with a negative deletion
charge; the missing object is the combined Pfaffian tariff `(11o)`.

## 5. Finite-exact frontier

The exact scout finds:

```text
universe                         all-negative rows   minimum individual chi
labelled n=2                               0                    0
labelled n=3                               0                    0
labelled n=4                               0                   -2
labelled n=5                               0                   -2
labelled n=6                               0                   -6
all 456 n=7 isomorphism reps               0                   -8. (15)
```

The minimum averaged half-charge in every row of `(15)` is zero, attained on
transitive examples. Among the 353 strong order-seven classes, the smallest
best local certificate is already positive:

```text
mask 171: (H,disc,S)=(27,4,23),
chi=(44,34,20,44,20,20,34),       max chi=44.               (16)
```

These are `FINITE-EXACT` facts only. The order-seven bank is checked against
its normalized-LF SHA-256 before use. No extrapolation in `n` is licensed.

The auxiliary exact path adds:

```text
OCF/Pfaffian support-moment objects                       34,322
local and averaged support-moment gates                  273,848
direct ear-cut gates through base order four               1,096
nonconstant ears over all labelled order-six bases     2,031,616
nonconstant ears over 456 order-seven representatives     57,456
minimum combined off-diagonal weight                         1/2
bases with a negative combined off-diagonal weight              0
minimum nonconstant-ear chi                                     0
bases with a negative nonconstant-ear chi                        0
order-six bases attaining minimum chi zero                  10,368
order-seven representatives attaining zero                    101. (16a)
```

It also finds minimum combined tariff gap zero and no negative gap in either
universe. The response construction uses exact `Fraction` inversion of
`I-K^2` and a permutation reconstruction of `Start/End/Q`; direct
Hamiltonian-DP/determinant evaluation checks every ear through base order
four. Equation `(16a)` is evidence for `(11o)`, not a proof.

## 6. Typed relation to the snark proof

The source/target ledger is

```text
source:      a local island/configuration in an apex counterexample;
target:      the deletion packet (C-x,u_x,full insertion fibers);
source test: local configuration is reducible;
target test: chi_x(C)>=0;
source global step: discharging makes one configuration unavoidable;
target open step: prove LU or ALU;
destroyed by scalarization:
             the covariant root, insertion multiplicities, and orphans;
required sidecar:
             u_x plus THM-4094's full deletion-incidence response.          (17)
```

This is an architectural analogy with an exact target identity, not a graph
transformation between cubic graphs and tournaments. The apex theorem's
Euler charge, cyclic embeddings, Kempe moves, and computer-checked
configuration atlas do not exist on the tournament side.

THM-4104 and THM-4099 suggest the next native operation: expand `(4)` into
the full Boolean ear-cut response rather than its mean. THM-4111 is a
firewall: two parents can have the same aggregate ear mean and different ear
images, so average growth alone cannot prove `(7)`.

## 7. Reproduction and next decisive tests

Reproduce the current frontier with

```bash
python3 -B 04-computation/tournament_deletion_slack_unavoidability_hyp9080.py
python3 -B -O 04-computation/tournament_deletion_slack_unavoidability_hyp9080.py
python3 -B 04-computation/tournament_alu_support_moment_ear_tariff_hyp9080.py
python3 -B -O 04-computation/tournament_alu_support_moment_ear_tariff_hyp9080.py
```

High-value next tests are:

1. derive a local sign decomposition of `chi_x` from the complete
   `Start/End/Q` ear tensor of THM-4104;
2. seek a double-counting proof of ALU, retaining orphan paths and the actual
   rooted Pfaffian energy rather than separate upper bounds;
3. prove or refute the combined ear tariff `(11o)`; the order-seven
   representative bank still has minimum `wtilde_ij=1/2`, so replace the
   permutation reconstruction by subset DP and scan a complete order-eight
   representative bank;
4. restrict first to strong tournaments and classify equality or near-zero
   rows under deletion;
5. enumerate order eight by isomorphism classes before using random evidence;
6. construct a hostile with all `chi_x<0`, which would refute LU without
   refuting `H>=disc`, or prove a structural reason one cannot occur.

Neither LU nor ALU is a proved dependency of any current theorem. In
particular, **`H>=disc` remains OPEN.**
