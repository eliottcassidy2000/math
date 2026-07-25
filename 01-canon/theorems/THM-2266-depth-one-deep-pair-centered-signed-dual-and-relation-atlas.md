---
id: THM-2266
title: "Depth-one deep-pair centered signed dual and relation-atlas reduction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. In each strict
  first-depth-one profile (1,b,c), a centered product of the two deep
  negative-carrier literals can be subtracted from THM-2239's charge
  without changing its signed residual response. The resulting pointwise
  nonnegative dual has an exact adversarial three-core Bellman bound below
  961/6930 for exactly the 120 interior rows 3<=b<=c-2. It does not exclude
  those rows unconditionally: it forces one of six guard/shallow-core
  reduced products to be at most 757. The exact 150-profile table, 120/30
  split, 758 cutoff, and 77,646 rational primal-dual LP certificates are
  frozen and independently replayed.
source: codex-2026-07-25-depth-one-deep-pair
depends_on:
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-2080-unequal-comb-overlap-removes-depth-five
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
  - THM-2239-unrestricted-multicore-signed-dual-profile-exclusion
related:
  - THM-2239-unrestricted-multicore-signed-dual-profile-exclusion
  - THM-2246-depth-one-private-joint-two-step-fibre-cap
  - THM-2255-valuation-separated-pair-cap-and-exclusive-owner-mass
script: 04-computation/lrc14_depth_one_deep_pair_signed_dual_thm2266.py
output: 05-knowledge/results/lrc14_depth_one_deep_pair_signed_dual_thm2266.out
script_sha256: 96487ea44d94b6d395fbe47423bbcce98ee2fd182d9333e788de7dda8784d38e
output_sha256: 616686853c9ab477a0a24423a479906a93a59dd26feeeeb7a8a131945153cf4b
hash_basis: working-tree bytes (LF)
---

# THM-2266 -- depth-one deep-pair centered signed dual

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

This theorem repairs THM-2239's time-zero degeneration without pretending
that the remaining shallow response has disappeared. The repair retains the
first time-zero literal and centers the *intersection* of the two deep
literals. Its exact finite computation produces a sharp structural
alternative:

```text
120 interior first-depth-one rows
    -> one of six explicit reduced products is at most 757;

45 boundary rows
    -> not controlled by this packet.
```

The scalar ledger is not decremented. The result is a finite
relation-atlas reduction inside 120 rows, not a proof of LRC(14).

## 1. Scalar residual and the two checkpoint carriers

Retain the scalar five-unit/three-blocker branch of THM-2198 and THM-2239.
On `R/Z`, write

```text
D_a={x:||ax||<1/14},             C_H={x:||Hx||>1/7},

R=1_(C_H)-sum_(i=1)^5 1_(D_(q_i)),

R_+=max(R,0),                    R_-=max(-R,0).       (1)
```

The guard and all five `q_i` are thirteen-units. The three blocker speeds
have a strict first-depth-one profile

```text
c_1=13u_1,       c_2=13^b u_2,       c_3=13^c u_3,

2<=b<c,                         5<=c<=19,            (2)
```

where the `u_j` are thirteen-units. Put

```text
p=integral R_+=integral R_- >= delta_5:=961/6930.    (3)
```

Let `T(x)=13x mod 1` and

```text
X_(j,t)=1_(D_(u_j)) o T^t.                           (4)
```

The even checkpoint `0` and odd checkpoint `1` give

```text
supp(R_+) subset
 P_0={X_(1,1) OR X_(2,b) OR X_(3,c)},

supp(R_-) subset
 N_1={X_(1,0) OR X_(2,b-1) OR X_(3,c-1)}.           (5)
```

The problem at first depth one is exact: THM-2239's centered atom
`Z_(1,0)` is the constant one, so summing the three centered atoms forgets
whether the first literal in `N_1` actually fired.

## 2. Center the deep pair, not the shallow literal

Put `rho=-1/13` and retain THM-2239's nonnegative centered atoms

```text
Z_(j,t)
 =X_(j,t)-rho^t X_(j,0)+max(rho^t,0).                (6)
```

They satisfy

```text
integral R Z_(j,t)=0,

Z_(j,t)>=X_(j,t)>=0                                  (7)
```

pointwise.

The new object is the product of the two *deep* literals in `N_1`. Define

```text
sigma=(b-1) mod 2,                 sigma in {0,1},

d=b-1-sigma,                       d even,

F =X_(2,b-1) X_(3,c-1),

F_0=X_(2,sigma) X_(3,c-b+sigma),

epsilon=13^(-d).                                          (8)
```

The two products are exact translates:

```text
F=F_0 o T^d.                                           (9)
```

THM-2222 gives `L^dR=(-1)^dR` for the unnormalized transfer operator `L`
of `T`. Transfer duality therefore gives, for every bounded `G`,

```text
integral R(G o T^d)
 =13^(-d) integral (L^dR)G
 =13^(-d) integral RG                                (10)
```

because `d` is even. Consequently

```text
W:=F-epsilon F_0

satisfies       integral RW=0.                       (11)
```

Now set

```text
A
 =X_(1,0)+Z_(2,b-1)+Z_(3,c-1)-W,

q=1-A.                                               (12)
```

This is the key pointwise repair. If

```text
x_1=X_(1,0),       x_2=X_(2,b-1),       x_3=X_(3,c-1),
```

then (7), (8), and the nonnegativity of `epsilon F_0` give

```text
A
 >=x_1+x_2+x_3-x_2x_3
 >=0,                                                (13)
```

and the right side is at least one whenever
`x_1 OR x_2 OR x_3=1`. Thus

```text
A>=1 and q<=0 on N_1.                                (14)
```

Unlike `Z_(1,0)=1`, the literal `X_(1,0)` in (12) remains visible. Unlike a
raw subtraction of `F`, the correction `W` is exactly response-null.

## 3. The remaining response and the signed-dual inequality

Write

```text
s=integral R X_(1,0).                                (15)
```

Equations (7), (11), and (12) give

```text
integral RA=s,                 integral Rq=-s.        (16)
```

Using `R=R_+-R_-`, equation (16) rearranges exactly to

```text
p
 =integral R_+(1-q)+integral R_-q-s.                 (17)
```

The middle term is nonpositive by (5) and (14). Moreover `0<=R_+<=1`,
`A=1-q>=0`, and `R_+` is supported on `P_0`. Therefore

```text
p<=integral 1_(P_0) A-s.                             (18)
```

This is the honest stopping boundary. The new carrier improves the finite
upper bound for the first term in (18), but the shallow signed response `s`
must still be controlled.

## 4. Exact arbitrary-coupling Bellman majorant

For a future bit vector `xi in {0,1}^3`, exact root counting gives

```text
P(X_(j,t)=1 | X_(j,t+1)=xi_j)
 =(2-xi_j)/13.                                       (19)
```

As in THM-2239, no independence or joint Markov property is assumed. Over
each actual root fibre, enlarge the joint current law to the complete
coupling polytope

```text
K_xi={
 pi_eta>=0:
 sum_eta pi_eta=1,
 sum_(eta:eta_j=1) pi_eta=(2-xi_j)/13,  j=1,2,3
}.                                                   (20)
```

At the terminal end, enlarge the actual joint law to all laws on
`{0,1}^3` with each marginal equal to `1/7`.

The companion runs a backward Bellman recurrence whose state retains:

```text
the time;
whether the positive clause P_0 has fired;
the four labelled bits
  X_(2,b-1), X_(3,c-1),
  X_(2,sigma), X_(3,c-b+sigma);
the future three-bit vector.                         (21)
```

At the terminal state it returns zero if `P_0` did not fire, and otherwise
returns exactly

```text
X_(1,0)
+X_(2,b-1)+X_(3,c-1)
-rho^(b-1)X_(2,0)-rho^(c-1)X_(3,0)
+max(rho^(b-1),0)+max(rho^(c-1),0)
-X_(2,b-1)X_(3,c-1)
+epsilon X_(2,sigma)X_(3,c-b+sigma),                 (22)
```

which is `A`. Let the resulting rational maximum be `B(b,c)`. Fibrewise
root induction and the terminal Haar average prove

```text
integral 1_(P_0)A <= B(b,c).                         (23)
```

Combining (3), (18), and (23), every counterexample in profile `(1,b,c)`
must satisfy

```text
delta_5 <= B(b,c)-s.                                 (24)
```

### Exact LP audit

The eight columns of (20) are

```text
(1,eta_1,eta_2,eta_3)^T,       eta in {0,1}^3.       (25)
```

Exactly `58` four-column bases are invertible. The eight root polytopes
have

```text
(6,9,9,6,9,6,6,6)                                   (26)
```

vertices in lexicographic future-bit order; the terminal polytope has six.
For every Bellman objective the companion:

1. enumerates every exact rational primal vertex;
2. chooses an optimal basic distribution;
3. solves `A_I^T y=c_I` exactly;
4. checks all eight dual inequalities and exact primal-dual equality.

The full table makes `118,970` LP calls and contains `77,646` distinct
objective/right-side pairs. Every distinct pair receives a matching exact
dual certificate.

## 5. The complete 150-profile spectrum

The strict first-depth-one universe is

```text
S_strict
 ={(1,b,c): 5<=c<=19, 2<=b<c},

|S_strict|=150.                                      (27)
```

The exact computation gives

```text
B(b,c)<delta_5

iff             3<=b<=c-2.                          (28)
```

Thus the subtarget set has exactly `120` rows. The complementary `30` rows
are exactly

```text
{(1,2,c):5<=c<=19}

union

{(1,b,b+1):4<=b<=18}.                                (29)
```

The weakest strict inequality in (28) is the row `(1,3,5)`:

```text
B(3,5)=10146800138/74231495611,

delta_5-B(3,5)
 =145591760833/73489180654890
 =0.001981131909970645....                            (30)
```

Three frozen interior controls are

```text
B(4,7)
 =283293422086229/2120125746145771,

B(4,8)
 =47679834772180172/358301251098635299,

B(10,19)
 =152389911763114035508249271108093757321845
  /1150805888298461593768523506367492533368331.      (31)
```

The stored transcript prints the exact reduced fraction, target difference,
subtarget status, and profile-specific relation cutoff for all `150` rows.
Its complete-table digest is

```text
0f5e660eecf284406777d66781d04a2528fafac714b7ff99ca49b9462f890847.
                                                               (32)
```

The subtarget and boundary profile digests are respectively

```text
b5259a91c2a9fca337cce681e21455040a3c944fe8c63e86c38d2b803d899189,

c5e141aba821f6f285f1c9c8c44211bb2843e7d0c31687c2ea59ca07efdd077f.
                                                               (33)
```

These are Bellman inequalities, not profile exclusions: equation (24)
still contains `s`.

## 6. The shallow response is a six-pair relation spectrum

Put

```text
E_H={x:||Hx||<1/7},           C_H=(R/Z) minus E_H
```

up to null endpoints, and define

```text
epsilon_0
 =measure(E_H intersection D_(u_1))-2/49,

eta_i
 =measure(D_(u_1) intersection D_(q_i))-1/49.        (34)
```

Since `measure(D_(u_1))=1/7`, equation (15) becomes

```text
s
 =measure(C_H intersection D_(u_1))
  -sum_(i=1)^5 measure(D_(u_1) intersection D_(q_i))

 =-epsilon_0-sum_(i=1)^5 eta_i.                      (35)
```

Reduce the six ordered pairs by their gcds:

```text
(a_0,b_0)
 =(H/gcd(H,u_1), u_1/gcd(H,u_1)),

(a_i,b_i)
 =(u_1/gcd(u_1,q_i), q_i/gcd(u_1,q_i)),

K_i=a_i b_i,                         0<=i<=5.         (36)
```

THM-2080's exact fat-guard/danger formula writes

```text
epsilon_0=(2/K_0)F(x,y),

F(x,y)=min(x,y)+(x+y-1)_+-2xy.                       (37)
```

A four-region check on the two diagonals of the unit square gives
`|F(x,y)|<=1/8`, hence

```text
|epsilon_0|<=1/(4K_0).                               (38)
```

THM-1166's folded two-danger formula gives similarly

```text
|eta_i|<=1/(4K_i),                  1<=i<=5.          (39)
```

If all six reduced products satisfy `K_i>=758`, then (35), (38), and (39)
give

```text
s>=-6/(4*758)=-3/1516.                               (40)
```

The worst Bellman gap in (30) exceeds this possible debit by the exact
amount

```text
145591760833/73489180654890 - 3/1516

 =124783729079/55704798936406620
 >0.                                                 (41)
```

Therefore, for every interior profile in (28),

```text
B(b,c)-s
 <=B(3,5)+3/1516
 <delta_5,                                           (42)
```

contradicting (24). This proves the theorem's main conditional
conclusion:

```text
Every counterexample in one of the 120 interior profiles satisfies

min_(0<=i<=5) K_i <=757.                              (43)
```

The cutoff is sharp for this uniform six-defect estimate. Replacing `758`
by `757` changes the final difference to

```text
145591760833/73489180654890 - 3/1514

 =-10404015877/27815654877875865
 <0.                                                 (44)
```

For each profile the companion also freezes the least sufficient product
cutoff. Its histogram `(cutoff, number of rows)` is

```text
(240,28), (241,38), (264,11), (265,1), (268,12),
(297,1), (298,12), (299,1), (342,1), (467,11),
(468,2), (584,1), (758,1).                           (45)
```

The digest of all `120` profile-specific cutoffs is

```text
3874ce2ba7b2c26e873c4a028883f94673d58403f3f0bf401f5b161282472d6f.
                                                               (46)
```

Equation (43) is a finite primitive-ratio atlas. There are exactly `3,643`
ordered coprime pairs `(a,b)` with `ab<=757`, or `1,822` unordered shapes.
For whichever pair in (36) triggers (43), the corresponding two speeds
`v=ga`, `w=gb` satisfy the primitive bounded relation

```text
bv-aw=0,               max(a,b)<=757.                (47)
```

The content of (43) is the bounded height, not the tautological fact that
two integer speeds are rationally related.

## 7. Boundary anatomy and interaction with THM-2255

The full current first-depth-one ledger has `165` rows. This theorem
partitions it as follows:

```text
120 interior rows:
  3<=b<=c-2, reduced to the finite relation alternative (43);

15 repeated-first rows:
  (1,1,c), outside the strict packet;

15 next-first rows:
  (1,2,c), where d=0 and W=F-F=0;

15 adjacent-deep rows:
  (1,b,b+1), where the Bellman value is not subtarget. (48)
```

Thus exactly `45` rows remain outside the relation-atlas conclusion. The
last two classes in (48) are precisely the two geometric boundaries of the
strict interior: no room before the middle depth, or no room between the
two deep depths.

THM-2255 supplies a different, complementary observable. Its
valuation-separated pair cap forces quantitative exclusive-owner mass and
a post-expiration split on strict rows. The present packet instead detects
the signed covariance of the shallow core against the guard and five unit
combs. Neither theorem currently supplies the other's missing labelled
handoff. In particular:

```text
THM-2255 does not remove the six-pair response s;

THM-2266 does not turn an expiration spill into a successor collision.   (49)
```

No combination is claimed here to close any of the `165` profiles.

## 8. Connection and loss ledger

The proved connection is

```text
source:
  the positive checkpoint-0 carrier, negative checkpoint-1 carrier,
  transfer eigenfunction R, and the two deep negative literals;

target:
  an exact four-slot three-core Bellman capacity plus a six-pair
  reduced-product atlas;

map:
  translate the deep pair back by the largest even amount, subtract its
  response-null centered difference, and isolate the surviving shallow
  response as six explicit pair defects;

preserved:
  the time-zero first literal, the current deep pair, its parity-compatible
  reference pair, every exact profile time, single-core 1-versus-2 root
  counts, and all six reduced gcd products;

destroyed:
  the shared root digit, cross-core phases beyond the selected pair,
  guard/five-unit higher incidence, exact owner labels after expiration,
  and joint restrictions beyond one-coordinate root marginals;

cheapest hostile probes:
  b=2, c=b+1, the cutoff-757 sign in (44), and arbitrary couplings of all
  three root processes;

needed sidecar:
  enumerate or structurally exclude the 1,822 small primitive ratio shapes,
  or join the bounded shallow relation to THM-2255's labelled expiration
  handoff without forgetting owners.                               (50)
```

The useful conceptual move is portable: when a linear centered charge
degenerates at time zero, center a product of the remaining active literals
along an even transfer return. This preserves a nonlinear incidence
coordinate while keeping exact signed response zero.

## 9. Reproduction and independent verification

Run

```bash
python3 04-computation/lrc14_depth_one_deep_pair_signed_dual_thm2266.py
python3 -O 04-computation/lrc14_depth_one_deep_pair_signed_dual_thm2266.py
```

Both modes produce the stored transcript byte for byte. All load-bearing
checks use explicit exceptions, not Python assertions. The companion
independently freezes:

```text
the 150-row universe and 120/30 split;
the exact pointwise packet inequalities for every hostile bit assignment;
the even transfer shift and all four slot times;
the complete rational bound table and its digest;
the weakest row and exact margin;
the profile-specific cutoff spectrum and uniform cutoff 758;
the exact failure of cutoff 757;
the 3,643/1,822 finite-atlas counts;
all coupling-polytope vertices;
77,646 exact primal-dual LP certificates;
normal versus optimized transcript equality.                        (51)
```

Independent audit checked the unnormalized-transfer factor and all four
time shifts in (8)--(11), the pointwise implication (13)--(14), the signed
rearrangement (17)--(18), Bellman time orientation and terminal state,
every exact rational primal-dual certificate, the response sign in
(35), and both folded-defect bounds (38)--(39). Independent normal and
optimized replays match the stored LF-normalized transcript and the
frontmatter hashes. QED.
