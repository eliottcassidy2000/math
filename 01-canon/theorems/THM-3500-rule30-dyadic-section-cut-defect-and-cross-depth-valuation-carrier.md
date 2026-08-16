---
id: THM-3500
title: "Rule 30 dyadic section cut defect and cross-depth valuation carrier"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The Rule 30 Mealy-section
  masks admit an exact
  lossless dyadic split with a rank-one suffix cut term and a closed driven
  defect recurrence.  Reversal identifies those masks literally with finite
  windows of the physical edge current, so the defect is the time derivative
  of dyadic bit holonomy and its parity is the next-depth packed bit.  A
  first-carry Hasse criterion is exact, with a width-four cancellation as the
  cheapest hostile.  At every hard Mersenne endpoint the actual pointed
  terminal profile is a lossless two-diagonal Hasse transform of the actual
  phase current; its operator has a Frobenius cross-depth law, but the actual
  cross-depth current sidecar remains missing.  A bounded constant-memory
  replay gives v_28,...,v_32=(69,70,72,74,77).  No eventual no-wrap statement
  and no Rule 30 prize consequence are claimed.
source: root/rule30-sharp-unlocks/dyadic-section-defect/2026-08-16
audit: >
  PASS (2026-08-16), independent end-to-end proof and adversarial replay
  audit.  The auditor rederived the split/current/holonomy bridge, first-carry
  criterion, Mersenne Hasse carrier, operator-only boundary, and the finite
  wrap consequence.  The companion gates the section/current bridge for
  every 1<=n<=256, the dyadic split through m=12, packed/mask controls through
  m=18, every basis current for p<=256, and a constant-memory 2^32-update
  replay.  Ordinary, optimized, and stored transcripts are byte-identical.
depends_on:
  - THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
  - THM-3463-rule30-mealy-section-suffix-parity-current-and-complexity-boundary
  - THM-3481-rule30-cyclic-arc-norm-rank-and-marked-innovation-spectrum
  - THM-3489-rule30-packed-restart-and-pointed-pascal-face
  - THM-3493-rule30-dyadic-wrap-atlas
related:
  - THM-3442-hensel-lifting-finite-differences-and-valuation-tree
  - THM-3476-rule30-depth-four-transverse-jet-barrier-and-slack-pascal-atlas
  - THM-3488-rule30-inward-slack-monicity-and-parity-cartier-ramification
  - THM-3501-rule30-universal-cover-green-potential-and-slack-holonomy-seam
  - THM-3503-rule30-odometer-ultrametric-regrading-and-orbit-closure-dimensions
script: 04-computation/rule30_dyadic_section_cut_defect_thm3500.py
output: 05-knowledge/results/rule30_dyadic_section_cut_defect_thm3500.out
script_sha256: b65a043072aacece522f9c4a4ef0b539ef922b310a01712ac4dd97fa254b8694
output_sha256: 7b5a2128e1b8c28fb89ba1ef1b48146d38c31283acee9597c5fcb0f49a74bab2
hash_basis: raw bytes
---

# THM-3500 -- Rule 30 dyadic section cut defect and cross-depth valuation carrier

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

This theorem isolates two exact carriers which had previously appeared only
as marginal Pascal data.

1. Splitting the THM-3463 Mealy-section mask at a factor cut introduces one
   rank-one suffix-parity term.  Keeping that term gives a lossless nonlinear
   recurrence across dyadic widths.
2. At a hard Mersenne endpoint, recentering the actual terminal arc turns its
   operator into multiplication by `1+Y^n`.  This is a unit and therefore a
   lossless within-depth recoding of the actual current.

The first carrier is now also physical: after reversal, a section mask is
exactly a finite time window of the Rule 30 edge current.  Its dyadic defect
is the finite difference of a translated edge bit.  The second carrier has an
exact Frobenius law across Mersenne depths, but that law does not transport the
actual current from one depth to the next.  These two scopes must not be
blurred.

## 1. Inheritance, conventions, and live board

All vector algebra is over `F_2`; `or` is coordinatewise Boolean OR.  For
`x in F_2^n`, put

```text
e_n=(1,...,1),
(S_n x)_i=xor_(j>i)x_j,
M_r(x)=xor_(0<=i<n) binom(i,r)x_i.                    (1)
```

Thus `S_n` is strict suffix parity and

```text
M_r(S_nx)=M_(r+1)(x).                                (2)
```

The Hasse coordinate in (1) is indexed by the **Mealy factor position**.  It
is distinct from THM-3476's transport-slack Hasse jet and from the cyclic
phase Hasse coordinates used at the Mersenne endpoint below.  The three use
the same Pascal algebra, not the same coordinate.

The inheritance pass is:

1. closest proved mechanism: THM-3463's suffix-parity/intersection compiler;
2. canonical hostile: its width-three example showing that marginal Hasse
   moments do not close coordinatewise OR;
3. corrected near miss: retain the labelled cut defect, rather than only its
   parity; and
4. least-used sidecars: the factor-cut rank-one term, the time-translation
   holonomy, and the actual current under the Mersenne operator.

The live board is the section mask, factor cut, defect word, edge current,
time holonomy, valuation pivot, pointed Mersenne face, and missing cross-depth
current sidecar.

## 2. The universal section masks and their packed readout

For every width `n>=1`, let `p_k^(n) in F_2^n` be the activity mask of the
written product representing the section `(B^n)|_(0^k)` in THM-3463.  Its
proved recurrence is

```text
p_0^(n)=p_1^(n)=e_n,
p_(k+1)^(n)=S_np_(k-1)^(n) or S_np_k^(n),   k>=1.    (3)
```

Write the inward right-edge row as

```text
R_t=sum_(j>=0)b_j(t)2^j,
R_(t+1)=R_t xor ((2R_t) or (4R_t)).                  (4)
```

Since `R_n=1 B^n(0^infinity)`, the output activity of the section gives, for
every `n>=1` and `k>=0`,

```text
b_(k+1)(n)=M_0(p_k^(n)).                             (5)
```

Consequently, if `R_n!=1`,

```text
nu_2(R_n-1)=1+min{k>=0:M_0(p_k^(n))=1}.              (6)
```

For dyadic time `q=2^m`, retain the proved THM-3493 notation

```text
v_m=nu_2(R_q-1).                                    (7)
```

## 3. Exact factor-cut extension and the driven defect

Fix any `q>=1` and split a vector of length `2q` into increasing factor-index
blocks

```text
x=(L,U),
L=(x_0,...,x_(q-1)),
U=(x_q,...,x_(2q-1)).                               (8)
```

### Provisional Theorem 3.1 (universal split and lossless defect recurrence)

For every `(L,U) in F_2^q x F_2^q`,

```text
boxed:
S_(2q)(L,U)=(S_qL+M_0(U)e_q, S_qU).                 (A)
```

Split `p_k^(2q)=(L_k,U_k)`.  Then

```text
U_k=p_k^(q),
Delta_k=L_k+U_k,
z_k=M_0(p_k^(q)),
Delta_0=Delta_1=0.                                  (9)
```

Suppressing the superscript `(q)`, the exact OR form of the recurrence is

```text
boxed:
Delta_(k+1)
 =(S p_(k-1) or S p_k)
  +([S(p_(k-1)+Delta_(k-1))+z_(k-1)e]
     or [S(p_k+Delta_k)+z_ke]),             k>=1.   (B)
```

Equivalently, define

```text
a_j=Sp_j,
d_j=SDelta_j+z_je.                                  (10)
```

Using `u or v=u+v+uv`, where juxtaposition is coordinatewise intersection,
the same recurrence is

```text
boxed:
Delta_(k+1)
 =d_(k-1)+d_k
  +a_(k-1)d_k+a_kd_(k-1)+d_(k-1)d_k.                (C)
```

Thus the only nonlinear terms are three explicitly labelled intersections.
No marginal list of Hasse moments is silently treated as closed.

The map

```text
(L_k,U_k) -> (p_k^(q),Delta_k)=(U_k,L_k+U_k)         (11)
```

is invertible, with `L_k=p_k^(q)+Delta_k`.  Hence this is a lossless
within-depth carrier, not a quotient statistic.

### Proof

For a lower coordinate `i<q`, the strict suffix consists of the remaining
lower coordinates plus every upper coordinate.  Its parity is therefore

```text
(S_qL)_i+M_0(U).
```

For an upper coordinate only the strict upper suffix remains.  This proves
(A).  Applying (A) to the recurrence (3) gives

```text
U_(k+1)=S_qU_(k-1) or S_qU_k,
L_(k+1)
 =(S_qL_(k-1)+z_(k-1)e_q)
   or (S_qL_k+z_ke_q).                              (12)
```

The bases identify `U_k=p_k^(q)` by induction.  Substitute
`L_j=p_j^(q)+Delta_j` into (12), then add the upper recurrence.  This is (B).
Expanding both ORs as `u+v+uv` cancels their common `a_(k-1)a_k` term and
leaves exactly (C).  Invertibility of (11) is immediate.  `square`

In block-matrix notation, (A) is

```text
S_(2q)=[ S_q       e_q tensor M_0 ]
       [  0              S_q       ].               (13)
```

The off-diagonal map has rank one.  Calling it a cut extension is literal;
calling it a Bockstein would require an additional chain complex which is not
present here.

### Corollary 3.2 (cross-depth valuation carrier)

For every `q>=1`,

```text
boxed:
nu_2(R_(2q)-1)=1+min{k>=0:M_0(Delta_k)=1}.           (D)
```

In particular, for `q=2^m`, the left side is `v_(m+1)`.

Indeed,

```text
M_0(p_k^(2q))=M_0(L_k)+M_0(U_k)=M_0(Delta_k),        (14)
```

and (6) proves (D).  The full word `Delta_k`, rather than only its parity, is
needed to continue (B) or (C).

## 4. Literal edge-current and holonomy interpretation

The defect above is not merely analogous to a physical seam.  It is the
finite difference of an actual Rule 30 edge bit.

Put `b_j(t)=0` for `j<0`, and define the edge-time current and time difference

```text
Q_k^edge(t)=b_(k-1)(t) or b_(k-2)(t),
(Df)(t)=f(t+1)+f(t).                                 (15)
```

The packed update (4) says

```text
Q_k^edge=Db_k.                                      (16)
```

Let reversal be

```text
(rho_nx)(t)=x_(n-1-t),       0<=t<n.                 (17)
```

### Provisional Proposition 4.1 (Mealy masks are current windows)

For every `n,k>=1`,

```text
boxed:
(rho_np_(k-1)^(n))(t)=Q_k^edge(t),       0<=t<n.     (18)
```

For `n=2q`, define

```text
delta_k^(q)=rho_qDelta_(k-1),
H_k^(q)(t)=b_k(t+q)+b_k(t).                          (19)
```

Then, for `0<=t<q`,

```text
boxed:
delta_k^(q)(t)
 =Q_k^edge(t+q)+Q_k^edge(t)
 =H_k^(q)(t+1)+H_k^(q)(t),                          (E)

M_0(Delta_(k-1))=b_k(2q),       k>0.                (20)
```

Thus `Delta` is literally the time derivative of dyadic bit holonomy, and
its scalar parity is the endpoint of that holonomy.

### Proof

Reversal changes strict suffix parity into strict prefix parity:

```text
rho_nS_n=P_nrho_n,
(P_nf)(t)=xor_(0<=s<t)f(s).                          (21)
```

Both reversed mask rows at `k=1,2` are the all-one word.  Likewise
`Q_1^edge=Q_2^edge=1`.  For `j>=1`, telescoping (16) and `b_j(0)=0` give

```text
P Q_j^edge(t)=b_j(t).                                (22)
```

Hence

```text
Q_(k+1)^edge=P Q_k^edge or P Q_(k-1)^edge,           (23)
```

which is exactly the reversed recurrence (3) with the same two bases.  This
proves (18).

For the split at `q`, the two entries compared by
`rho_qDelta_(k-1)` are factor positions `q-1-t` and `2q-1-t` in the
length-`2q` mask.  Equation (18) identifies them with times `q+t` and `t`,
proving the first equality in (E).  The second follows because time
translation commutes with `D`.  Finally telescope over `0<=t<q`:

```text
xor_t delta_k^(q)(t)
 =H_k^(q)(q)+H_k^(q)(0)
 =b_k(2q)+b_k(q)+b_k(q)+b_k(0)
 =b_k(2q).                                          (24)
```

This is (20).  `square`

### 4.2 Exact overlap with the universal-cover slack holonomy

There is an exact scalar bridge to the THM-3501 object, but not
an identification of its full eventwise marker.  THM-3501 uses the phase
current

```text
Q_k^phase(h)=Q_k^edge(k+h)                           (25)
```

and a slack-polynomial potential `E_k(h;u)`.  THM-3501 needs only
`0<=a<k`; extend the same endpoint formula to every `a>=0` and denote it by
`widehat Omega_a(u)`.  At marker `u=1` and translation length `p=P_k`,

```text
widehat Omega_a(1)
 =E_k(-k+p+a;1)+E_k(-k+a;1)
 =b_k(p+a)+b_k(a)
 =H_k^(p)(a).                                       (26)
```

Therefore

```text
D_a(widehat Omega_a(1))=(rho_pDelta_(k-1))(a),
 0<=a<p,       with Delta formed at cut q=p.         (27)
```

This is a genuine commuting scalar specialization: both constructions share
the same translation finite-difference and cut-boundary object.  The packed
restart makes `widehat Omega_a(1)=epsilon_k` constant, so (27) vanishes at the seed
period.  However, the full polynomial `Omega_a(u)` retains Green-path slack
labels which the binary word `Delta` does not contain.  There is no proved
map recovering that marker polynomial, its event labels, or its low phase-
Hasse tomography from `Delta`.  The shared scalar cut-defect is therefore
more than analogy, but it is not a full marked chain equivalence and not a
literal Bockstein without extra chain data.

## 5. First carry, Hasse criterion, and the cheapest hostile

Let `q=2^m` with `m>=1`, and put

```text
t=v_m-1=min{k:M_0(p_k^(q))=1}.                       (28)
```

### Provisional Theorem 5.1 (first-carry criterion)

The defect has the exact initial shape

```text
Delta_0=...=Delta_t=0,
Delta_(t+1)=e_q+S_qp_(t-1)^(q).                     (29)
```

Consequently,

```text
boxed:
v_(m+1)=v_m+1
 iff M_1(p_(v_m-2)^(q))=1.                          (F)
```

If this Hasse bit is zero, then `v_(m+1)>v_m+1`; the recurrence (B) still
determines the later carry exactly.

### Proof

Before the pivot, every driver `z_j` is zero.  Starting from
`Delta_0=Delta_1=0`, either (B) or (C) therefore gives
`Delta_j=0` through `j=t`.  At the next step, `z_(t-1)=0` and `z_t=1`.
Coordinatewise, for bits `a,b`,

```text
(a or b)+(a or (b+1))=1+a.                           (30)
```

Apply this with `a=S_qp_(t-1)` to obtain (29).  Since `q` is even,
`M_0(e_q)=0`; (2) then gives

```text
M_0(Delta_(t+1))=M_1(p_(t-1)^(q)).                  (31)
```

All earlier defect parities vanish, so (D) proves (F).  `square`

The cheapest actual hostile is `q=4`.  In increasing factor-index order,

```text
v_2=4,
Delta_4=0011,       weight(Delta_4)=2,
Delta_5=0111,       weight(Delta_5)=3,
v_3=6.                                                   (32)
```

The first nonzero defect word can have even parity.  Therefore “the first
cut disturbance immediately carries” is false; (F), not nonzeroness of the
word, is the repaired statement.

## 6. The actual hard Mersenne endpoint operator

Now change from factor-index Hasse coordinates to cyclic **phase** Hasse
coordinates.  Let

```text
q=2^m,
n=2q-1,
p=P_n,
n<p.                                                 (33)
```

Thus `n` is a hard Mersenne endpoint.  In the phase module

```text
V_p=F_2[X]/(X^p-1)=F_2[Y]/(Y^p),       Y=X+1,        (34)
```

let `X` act as THM-3481's backward cyclic shift.  Let `Q_n` be the **actual**
Rule 30 terminal phase current and

```text
F_n=A_nQ_n,
N_n(X)=1+X+...+X^(n-1),
A_n=X^(-(n-1))N_n(X).                               (35)
```

Write `M_j` for phase Hasse coordinates at `X=1`, and recenter the physical
mark by

```text
widetilde F_n=X^nF_n.                               (36)
```

### Provisional Theorem 6.1 (lossless Mersenne carrier)

The actual recentered profile obeys

```text
boxed:
widetilde F_n=(1+Y^n)Q_n.                           (G)
```

Consequently its Hasse coordinates have the two-diagonal law

```text
boxed:
M_j(widetilde F_n)
 =M_j(Q_n)+1_[j>=n]M_(j-n)(Q_n),       0<=j<p,      (H)
```

and the triangular inverse

```text
boxed:
M_j(Q_n)
 =xor_(0<=r<=floor(j/n))M_(j-rn)(widetilde F_n).    (I)
```

The marked center has the equivalent current-tail and profile-face forms

```text
boxed:
c_n=widetilde F_n(0)
   =xor_(j=p-n)^(p-1)M_j(Q_n),                       (J)

boxed:
c_n=xor_(r=0)^(q-1)M_(p-n+2r)(F_n).                 (K)
```

Thus the current form uses all top `n` Hasse orders, while the pointed
profile form uses exactly the top `q` odd Hasse orders.

Finally, if

```text
C_m(Y)=1+Y^(2^(m+1)-1)=1+Y^n,                       (37)
```

then its exact Frobenius/Cartier section law is

```text
boxed:
C_(m+1)(Y)=1+Y+Y C_m(Y)^2,
C_m(Y)=1+Y(Y^2)^(q-1).                              (L)
```

### Proof

From (35),

```text
widetilde F_n=XN_n(X)Q_n=(N_(n+1)(X)+1)Q_n.         (38)
```

Here `n+1=2q` is a power of two.  Frobenius gives

```text
N_(2q)(X)=(X^(2q)+1)/(X+1)=Y^(2q-1)=Y^n,            (39)
```

which proves (G).  Multiplication by `1+Y^n` gives (H).  Its constant term is
one, so it is a unit in `F_2[Y]/(Y^p)` with inverse

```text
(1+Y^n)^(-1)=1+Y^n+Y^(2n)+...                       (40)
```

truncated below degree `p`; coefficient comparison proves (I).

The physical no-wrap mark is `F_n(p-n)`, which becomes
`widetilde F_n(0)`.  THM-3489's current-face formula simplifies because
`p-n-1=p-2q` has its low `m+1` bits zero: its Lucas coefficient is one
exactly on orders `p-n,...,p-1`.  This proves (J).  In the profile face,
`n-1=2q-2`; its submasks are exactly `2r`, `0<=r<q`.  THM-3489's Pascal
inversion gives (K).  Squaring `C_m` and using `2n+1=2^(m+2)-1` proves the
first identity in (L); the odd exponent `n=1+2(q-1)` proves the Cartier split.
`square`

The multiplier `C_m` is a lossless within-depth carrier by (I).  Equation
(L), however, relates only the operators.  It supplies no relation between
the actual currents `Q_n` and `Q_(2n+1)`.  The cheapest ambient hostile is
the pair `Q=0` and `Q=Y^(p-n)`: the same operator law applies to both, while
their marked functionals in (J) are opposite.  These are ambient currents,
not a claim about two realized Rule 30 depths.  An actual-current transport
sidecar is still required for a cross-depth center theorem.

## 7. Finite-exact valuation extension

The companion embeds a direct packed C gate with the frozen universe

```text
initial state:       R_0=1,
updates:             exactly 2^32,
retained bits:       0,...,119,
sample times:        2^m, 0<=m<=32,
reported statistic: nu_2(R_(2^m)-1).                (41)
```

Low packed bits evolve independently of discarded high bits.  Every reported
valuation is below `120`, so the truncation cannot change it.  The gate uses
constant memory; it is not the optional approximately `1.5 GiB` section-mask
route.  It reproduces the independently audited THM-3493 values through
`m=27` and gives the new finite-exact extension

```text
boxed:
(v_28,v_29,v_30,v_31,v_32)=(69,70,72,74,77).         (42)
```

By the proved THM-3493 dyadic wrap atlas, this finite transcript implies

```text
W intersect [1,2^33-1]={1,2,3,4},
[5,2^33-1] subset H.                                (43)
```

The exact upper endpoint is `2^33-1=8,589,934,591`.  This is
`FINITE-EXACT` pending independent replay.  It is not evidence for a stated
eventual bound `v_m<2^m`, and no asymptotic wrap-density conclusion is drawn.

The THM-3503 ultrametric regrading makes a future linear bound
`v_m=O(m)` the direct route to a positive lower Hausdorff-dimension statement
for the odometer orbit closure.  The literal section/current duality (18)--
(20) and the lossless defect recurrence (B)--(C) provide an exact carrier on
which such a bound could be attacked.  Nothing here proves that bound or any
dimension conclusion; THM-3503 is related context only, not a dependency.

## 8. Source/target/map ledger and failure boundaries

| Connection | Source -> target and map | Preserved | Lost | Required sidecar / cheapest hostile |
|---|---|---|---|---|
| factor cut | `p_k^(2q) -> (p_k^(q),Delta_k)` by `(L,U)->(U,L+U)` | every labelled mask bit; invertible | nothing within one row | recurrence needs both full words; scalar parity alone fails already at `q=4` |
| cut to valuation | `Delta_k -> M_0(Delta_k)` | the next-depth valuation pivot | all labelled intersections needed after the pivot | retain `Delta_k`; `0011` is nonzero but parity zero |
| Mealy to physical | `rho_np_(k-1)^(n) -> Q_k^edge|_[0,n)` | every bit on the finite window | times outside the window | extend `n`; checked for all `n<=256` in the companion |
| current to holonomy | `tau_qQ_k^edge=D(tau_qb_k)` | exact translated current defect and endpoint parity | the additive constant of `H` | one base value; at `q=P_k` restart supplies `epsilon_k` |
| eventwise holonomy | `Omega_a(u) -> Omega_a(1)=H_k^(P_k)(a)` | scalar physical bit holonomy | every slack/Green marker coefficient | full event-labelled potential; no recovery from `Delta` |
| Mersenne current | `Q_n -> widetilde F_n=C_mQ_n` | all phase-current information; explicit inverse | nothing within depth | marked phase calibration is retained in the recentering |
| Mersenne scale | `C_m -> C_(m+1)` by (L) | exact operator multiplier | which actual current occurs next | relation `Q_n -> Q_(2n+1)`; ambient `0` versus `Y^(p-n)` |

The equality/failure boundaries are therefore explicit:

1. (A)--(E) hold for every positive cut width; dyadicity enters the valuation
   notation and the first-carry simplification.
2. (F) requires even `q=2^m`, `m>=1`; at odd `q`, `M_0(e_q)=1` adds a term.
3. (G)--(L) require a power-of-two phase cycle and the hard Mersenne condition
   `n=2^(m+1)-1<p=P_n`.
4. (I) is lossless only within a fixed phase module.  It does not close the
   actual current across depths.
5. (42)--(43) stop at the displayed finite cutoff.

## 9. Reproduction and no-prize scope

Run

```bash
python3 04-computation/rule30_dyadic_section_cut_defect_thm3500.py
python3 -O 04-computation/rule30_dyadic_section_cut_defect_thm3500.py
```

and compare both byte-for-byte with

```text
05-knowledge/results/rule30_dyadic_section_cut_defect_thm3500.out
```

The default verifier freezes:

1. the Mealy/current bridge for every `1<=n<=256` and
   `1<=k<=min(n+4,96)`;
2. the split, defect, and holonomy identities for `q=2^m`, `0<=m<=12`, and
   `k<=min(2q+5,97)`;
3. independent packed/mask valuation controls through `m=18`;
4. the Mersenne carrier on every basis current for `p=2^d`, `2<=d<=8`, every
   Mersenne `n<p`, plus twelve operator scales; and
5. the constant-memory direct packed gate in (41).

Setting `RULE30_THM3500_LONG_SECTION=1` additionally runs a frozen independent
section-mask implementation through `m=32`; it is optional because its peak
mask memory is about `1.5 GiB` and it is not the sole validity gate.

This theorem proves neither eventual no-wrap, infinitely many hard
Mersenne endpoints, center nonperiodicity, center density `1/2`, nor a random-
access lower bound.
