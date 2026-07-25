# LRC14 depth one is an odd-transfer cone, not a missing coefficient

Source: codex-2026-07-25 depth-one Bellman audit, after THM-2234,
THM-2244, and incoming THM-2246.

Status split:

- **PROVED (analytic):** the nonnegative affine response-null cone is a
  bipartite transport cone between even and odd transfer times.  A shallow
  time-zero atom costs at least coefficient `13` on future atoms, with
  equality only for the bridge `X_0+13X_1`.
- **PROVED (exact hostile Bellman law):** the natural nonlinear owner gate
  fails uniformly on all `165` first-depth-one profiles.  There are two
  obstruction classes, according as one or two blockers have depth one.
  The obstruction can be supported entirely on the complement of the
  THM-2234 private stratum.
- **PROVED (exact local hostile family):** the single-owner sharp phase in
  THM-2246 extends to every finite fibre depth.  Its occupancy tends to
  `111/182`, and no depth choice can force more than the expansion ceiling
  `2197/1332` from the retained local data.
- **VERIFIED-EXACT (one best-row audit):** for `(1,10,19)`, the marginal
  private-mass Lagrange dual has its minimum at zero charge; exact primal
  active masses expose why.
- **FINITE-SCOUT:** explicitly retaining `P,T(P),T^2(P)` as conditional
  Bellman state still does not separate the active face in the relaxation
  using the separate root caps `10,12`.  THM-2246's sharper joint cap `112`
  changes the second-image floor but not the qualitative obstruction.

Nothing below excludes a valuation profile.  The current scalar ledger stays
`165` depth-one profiles plus the all-equal normalized-core branch of
`(3,4,5)` isolated by incoming THM-2250, and LRC(14) remains open.

## 1. The response quotient

Use the notation of THM-2244:

```text
rho=-1/13,
X_(j,t)=1_(D_(u_j)) o S^t,

integral R X_(j,t)=rho^t integral R X_(j,0).          (1)
```

For one labelled core, consider a finite nonnegative affine observable

```text
Y=sum_(t>=0) a_t X_t,                 a_t>=0.          (2)
```

Universal response-nullity is

```text
sum_t a_t rho^t=0.                                    (3)
```

Put `m_t=a_t 13^(-t)`.  Equation (3) says exactly

```text
sum_(t even)m_t=sum_(t odd)m_t.                        (4)
```

Choose any transport plan `gamma_(e,o)` between these two finite measures,
where `e` is even and `o` is odd.  Then

```text
sum_t a_t z^t
 =sum_(e even,o odd) gamma_(e,o)(13^e z^e+13^o z^o).  (5)
```

Thus the extreme rays are two-time parity bridges.  This is a literal
bipartite transport description of the positive response-null cone, not
just a coefficient estimate.

If `Y>=X_0` on the full independent Boolean box, then `a_0>=1` and every
other coefficient is nonnegative.  The odd bank in (4) has mass at least
one, so

```text
sum_(t odd)a_t>=13.                                   (6)
```

Equality in (6) forces all even mass to be the single unit at time zero and
all odd mass to be at time one:

```text
Y=X_0+13X_1.                                          (7)
```

For several labelled cores, an affine score has response

```text
sum_j (sum_t a_(j,t)rho^t) integral R X_(j,0).        (8)
```

The unrestricted transfer/Bellman model assumes no relation among the
labelled core responses.  Hence a universal affine certificate must cancel
core by core.  Cross-core cancellation is not a free coefficient trick; it
would be a new arithmetic sidecar.

This completely classifies the natural affine search that looked tempting
after THM-2244.  The coefficient `13` is optimal, not naive.

## 2. Nonlinear owner gating improves the score but preserves its kernel

Let `S` be the set of depth-one cores and `D` the deeper cores.  For a deep
core put

```text
s_j=lambda_j mod 2,
d_j=lambda_j-1-s_j,                  d_j positive odd,

Y_j=X_(j,lambda_j-1)+13^(-d_j)X_(j,s_j).              (9)
```

Define the shallow-only negative owner event and its shift

```text
H =
 (OR_(j in S) X_(j,0))
 AND (AND_(j in D) NOT X_(j,lambda_j-1)),

H^+=H o S
 =
 (OR_(j in S) X_(j,1))
 AND (AND_(j in D) NOT X_(j,lambda_j)).               (10)
```

For every measurable `H`,

```text
integral R(H+13 H o S)=0.                              (11)
```

Consequently

```text
Y=H+13H^+ + sum_(j in D)Y_j                            (12)
```

is response-null and dominates the selected negative carrier pointwise.
This is the clean nonlinear repair of the affine bridge: charge `13` only
when the positive clause is owned solely by shallow cores.

It still has an exact quotient kernel.  In THM-2233's four-bit arbitrary
coupling recursion, choose at every root step the product of the four exact
conditional Bernoulli marginals.  The resulting guard and danger chains are
independent and stationary, with marginals

```text
P(G_t=1)=5/7,                 P(X_(j,t)=1)=1/7.        (13)
```

This is a feasible Bellman law; it need not be an arithmetic realization.

If there is one shallow core, the event

```text
K={G_0=1, X_(shallow,1)=1,
              X_(deep 1,lambda_2)=X_(deep 2,lambda_3)=0}
```

has exact mass

```text
measure(K)=(5/7)(1/7)(6/7)^2=180/2401.                (14)
```

On `K`, the positive carrier is hit, `H^+=1`, and `Y>=13`.  Therefore the
THM-2233 four-bit arbitrary-coupling Bellman bound for this score is at least

```text
13 measure(K)=2340/2401
                    >961/6930                         (15)
```

by the exact margin `1986977/2376990`.  This covers all `151` profiles with
exactly one depth-one core, independently of the two deeper valuations.

For the `14` profiles `(1,1,c)`, require one fixed peeled shallow owner to
fire and the sole deep positive owner to be absent.  This event lies in the
peeled owner's nonprivate stratum and contributes

```text
13(5/7)(1/7)(6/7)=390/343
                    >961/6930                         (16)
```

by `339011/339570`.  If the full two-shallow union is used, the stronger
lower bound is `5070/2401`.

Equations (15)--(16) are the exact two-class obstruction census for all
`165` rows.

The optimized Bellman faces exhibit the same split more sharply.  The exact
best one-shallow row is

```text
(1,10,19):
125081347917950886017278696763997
/108520529618636773506886745738647
=1.152605395103693... .                              (17)
```

For exact representatives of the four `(lambda_2 mod 2,lambda_3 mod 2)`
classes, a terminal optimal law puts mass `5/7` on

```text
(guard,dangers)=(1,000)
```

and the remaining two atoms, each of mass `1/7`, on a complementary danger
pair with guard zero:

```text
(even,odd): {001,110},       (odd,odd):  {011,100},
(even,even):{000,111},       (odd,even): {010,101}.   (18)
```

The controls were `(1,2,5)/(1,10,19)`, `(1,3,5)/(1,11,19)`,
`(1,10,18)`, and `(1,9,18)`.  The two-shallow controls instead split by
the parity of the sole deep valuation.  The exact worst endpoint is

```text
(1,1,5): 6683153/2599051=2.571382015974292... .       (19)
```

Thus the active terminal coupling anti-aligns the guard with the danger
mass and uses parity only to decide which danger labels coalesce.  The
uniform rigorous classification remains the two feasible-law obstructions
(15)--(16); the parity table is an exact representative face audit, not a
claim that every optimizer is unique.

## 3. Why the private marginal is invisible

For one shallow core, THM-2234's private event is

```text
Q={X_(shallow,1)=0}.                                  (20)
```

The hostile event `K` in (14) is contained in `Q^c`.  Thus the bad
translated face and the known private mass can be disjoint.

This is visible exactly in the best numerical row `(1,10,19)`.  Let `E` be
the guard-positive carrier and set

```text
B(lambda,mu)
 =Bellman[1_E(Y-1+lambda+mu 1_Q)_+]
   -lambda(961/6930)-mu(2593/90090).                  (21)
```

The exact zero-charge value is

```text
B(0,0)
=88665671065548681775462541547897546079163
 /88523529869112430289886423566730194874487
=1.001605688302832... .                               (22)
```

There is a zero-charge optimal Bellman policy for which the mass of the
right-active face `E intersect {Y>=1}` and its private part are

```text
M=220749535108541786869/1461920290375446110677,

N=638256747361221833948/10233442032628122774739.       (23)
```

They satisfy, exactly,

```text
M-961/6930
 =124888879251390870641573
  /10131107612301841546991610 >0,

N-2593/90090
 =340277309661184095225023
  /10131107612301841546991610 >0.                     (24)
```

Since, pointwise for `lambda,mu>=0`,

```text
(f+lambda+mu 1_Q)_+
 >=f_+ +1_{f>=0}(lambda+mu 1_Q),                      (25)
```

the fixed policy in (23) proves

```text
B(lambda,mu)
 >=B(0,0)+(M-961/6930)lambda+(N-2593/90090)mu
 >=B(0,0)>0.                                         (26)
```

So the exact marginal Lagrange optimum is `(lambda,mu)=(0,0)`.  The private
invoice is not merely too small; it is attached to the wrong face.

A floating census of all `165` rows found no contrary slope.  Exact controls
at the two apparent slope minima were

```text
profile (1,14,15):
M=50450339183795443/358301251098635299
 >961/6930,

profile (1,2,5):
N=133307/2599051
 >2593/90090.                                        (27)
```

The uniform theorem is (15)--(16); the census and (27) are mechanism audits,
not an additional profile exclusion.

## 4. Conditioning on the private image bit

A second prototype retained three marked events in the Bellman state:

```text
J_0=P,                 J_1=T(P),                 J_2=T^2(P).
```

It enforced

```text
J_0 => guard AND Q AND a deep positive owner,
J_1 => Q AND a deep owner,
J_2 => a deep owner,                                  (28)
```

and conditioned the root coupling on the future marked bit.  A marked
`J_1` fibre had between `1` and `10` marked `J_0` roots; a marked `J_2`
fibre had between `1` and `12` marked `J_1` roots.  The mass floors used
were

```text
2593/90090,       2593/69300,       33709/776160,     (29)
```

where the last is THM-2246's incoming joint two-step improvement.

For `(1,10,19)`, zero charge again gave `1.0016056883028...`.  Every one of
the seven nonzero directions in the three nonnegative image charges had a
positive numerical right slope.  The coordinate slopes after subtracting
the three floors were approximately

```text
0.0335873741,          0.189988220,          0.221875643. (30)
```

The first agrees with the exact private slack in (24).  Charges
`0.1,0.25,0.5,1,2,4,8` along the common direction increased monotonically.

This is only a finite floating probe, and its transition model retains the
separate `10,12` caps rather than THM-2246's full co-adapted `112` table.
Its useful verdict is structural: a baseline-optimal policy can complete a
marked active private path through unmarked sibling roots.  Merely adding
the image-support bit does not force the private event onto the negative
defect face.

THM-2246 already proves the exact local stopping boundary: even the joint
floor `33709/776160` is below one danger-comb mass `1/7`.  The present audit
identifies the orthogonal stopping boundary: the active positive-part face
can live entirely on `Q^c`.

## 5. A persistent owner saturates arbitrarily deep local fibres

THM-2246 gives a locally sharp two-step phase with one labelled deep owner.
The same obstruction persists at every finite depth, so adding an owner
label without recording a switch or carry is still insufficient.

Fix `L>=2`, put

```text
N=13^L,                    M=13^(L-1),
H=1,                       a=2,
c_1=26,                    c_2=2N.                    (31)
```

Choose the terminal phase and five distinct unit coefficients by

```text
L even: z=1/2,   Q=2,
L odd:  z=55/56, Q=56,

q_i=1+QN i,                         1<=i<=5.           (32)
```

Both phases lie in `D_2`, `Qz` is integral, and every `q_i` is a
thirteen-unit.  For every root

```text
x=(z+k)/N,                          0<=k<N,
```

one has

```text
q_i x=x mod 1,                 c_2 x=2z mod 1.        (33)
```

Thus all five unit masks align with `D_1`, which is disjoint from the guard
`C_1`, while the single labelled owner `c_2` owns the entire `N`-root
fibre.

Count the locally private eligible roots

```text
E_L(z)={x:T^Lx=z, x in C_1, Tx notin D_2,
                      x notin union_i D_(q_i)}.       (34)
```

For `y=Tx`, exact guard-root counting gives

```text
#{x:Tx=y, x in C_1}=10-1_(C_1)(y).                   (35)
```

Among the `M` roots `y` of `z`, let

```text
A=#{y:y notin D_2},
I=#{y:y notin D_2 and y in C_1}.
```

Then `|E_L(z)|=10A-I`.  Since `13^2=1 mod 28`, direct interval counting
gives

```text
L even:
A=6(M+1)/7,               I=(9M-5)/14,
|E_L|=(111M+125)/14;

L odd:
A=6(M-1)/7,               I=9(M-1)/14,
|E_L|=111(M-1)/14.                                  (36)
```

No load-bearing endpoint occurs; the denominators `7,14,28` are coprime to
the `13`-power grids in (31)--(32).  At `L=2`, equation (36) recovers the
joint occupancy `112` from THM-2246.  The occupied proportion tends to

```text
111/182=0.609890109890... .                           (37)
```

Consequently a pure `L`-step fibre-cap argument using only the guard,
peeled complement, five unit labels, and one persistent deep owner cannot
force an expansion factor larger than

```text
L even: 182M/(111M+125),
L odd:  182M/(111(M-1)).                              (38)
```

Across `L>=2`, the largest of these ceilings is the `L=3` value

```text
2197/1332=1.649399399... .
```

Even multiplying THM-2234's private floor by that value gives only

```text
(2593/90090)(2197/1332)
 =438217/9230760
 <1/7.                                               (39)
```

This hostile family need not satisfy the global scalar cover.  Its precise
role is to delimit local information: arbitrary fibre depth plus a
persistent owner label still cannot pay one exhausted danger comb.  A
consumer must detect owner switching, shared ownership, or cross-fibre
winding.

## 6. The parity contrast with the successful depth-three score

Incoming THM-2243 first unions the three blocker atoms and then centers that
composite event across a two-step transfer interval.  The response factor is

```text
rho^2=1/169>0,
```

and two odd negative clauses make its hinge sign-safe.  At depth one there is
only one available odd clause and the bridge crosses an odd transfer
interval:

```text
rho=-1/13.
```

The factor `13` is therefore the parity boundary between the successful
depth-three mechanism and the failed depth-one mechanism.  Another pass over
affine coefficients, unions without owner labels, marginal private mass, or
unlabelled image support cannot cross it.

Incoming THM-2250 adds another useful hostile comparison.  Its same-time
distinct-core incidence cap

```text
integral X_(a,t)X_(b,t)<=1/14
```

closes every non-all-equal partition of `(3,4,5)`.  The same marginal tax is
invisible to the depth-one exclusive-owner kernel.  Indeed, under the
stationary product law used in (14),

```text
integral X_(a,t)X_(b,t)=1/49.
```

For any nonnegative packet of same-time pair penalties `theta_(t,ab)`, its
Lagrange bound is therefore at least

```text
hostile score
+sum_(t,ab) theta_(t,ab)(1/14-1/49)

=hostile score +(5/98)sum_(t,ab)theta_(t,ab).         (40)
```

Thus the THM-2250 tax strictly worsens the product lower witness in the
all-distinct Bellman partition.  This explains the mechanism difference:
the depth-three carrier forces repeated union incidence, whereas the
depth-one bad face is deliberately owned by exactly one shallow core.
Useful incidence must be conditional on an owner switch or cross-time
carry, not another unconditional same-time pair cap.

The next consumer must retain information that the hostile product law
forgets:

```text
which deep owner receives the private point;
whether that owner persists, switches, or ties on the next fibre;
the carry/root sheet on which the switch occurs; and
the co-adapted guard/peeled incidence from THM-2246.   (41)
```

In tournament language, the honest vertices are labelled owner-carry
obligations, not runners and not bare clauses.  Orient a pair by which owner
can absorb the other's marked subflow without paying a switch.  Ties must
retain their shared-root multiplicity.  A directed cycle would be a literal
owner-switch frustration cost; a transitive orientation would expose a
single persistent owner to peel.  That is the smallest new coordinate with
a plausible route around the odd-transfer cone.

The reusable meta-pattern is:

> When a transfer certificate fails at an odd boundary, classify the positive
> response-null cone before searching coefficients.  Then locate the exact
> carrier of the forced bridge and ask whether every proposed sidecar meets
> that carrier.  A large sidecar on a disjoint stratum is mathematically
> invisible.
