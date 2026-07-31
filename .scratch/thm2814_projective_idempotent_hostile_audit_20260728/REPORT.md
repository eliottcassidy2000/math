# Independent hostile audit of the proposed THM-2814

**Status:** `ACCEPT AFTER SCOPE REPAIR`.  The corrected two-branch
classification is mathematically sound and load-bearing.  The proposed
THM-2779 and THM-2807 applications also survive, with the qualifications
below.  This is a scratch audit, not a canon promotion or an LRC(14)
conclusion.

Audited sources:

- `.scratch/lrc_nonidempotent_twist_20260728/REPORT.md`;
- `.scratch/lrc_nonidempotent_twist_20260728/probe.py`;
- `01-canon/theorems/THM-2779-bockstein-symplectic-decoder-frame-torsor-and-heisenberg-root-degree-gate.md`;
- `01-canon/theorems/THM-2807-positive-graded-address-two-simplex-and-allocation-lift-boundary.md`.

The companion in this directory is independent: it imports none of the
proposed theorem code and exhausts fresh finite universes.

## 1. Exact row/column orbit theorem

Let `K` be a field and let

```text
v=(v00,v10,v01,v11) in (K^x)^4.
```

The row/column torus acts by

```text
vij |-> r_i c_j vij,       r_0,r_1,c_0,c_1 in K^x.       (1)
```

Its scalar kernel is

```text
(r0,r1,c0,c1)=(lambda,lambda,lambda^-1,lambda^-1).       (2)
```

The cross-ratio

```text
kappa(v)=v00 v11/(v10 v01)                               (3)
```

is invariant.  It is also complete.  One explicit normalization is

```text
(r0,r1,c0,c1)
 =(1,v00/v10,v00^-1,v01^-1),

v |-> (1,1,1,kappa(v)).                                  (4)
```

Consequently the effective action is free and its orbits are indexed by
`K^x`.

Over `F_q`, the unit locus has `(q-1)^4` points, split into `q-1` orbits of
size `(q-1)^3`; the stabilizer in the unreduced four-parameter torus has
size `q-1`.  The companion exhausts this classification over
`F_3,F_5,F_7,F_13`.

This proof works over a commutative ring if **all four entries are units**.
It does not work under the report's literal phrase “all four amplitudes are
nonzero” in the ambient commutative-ring scope:

- in `Z/6`, `(1,2,3,1)` has four nonzero entries but denominator
  `2*3=0`, so `(3)` is undefined;
- over `Z`, `(1,2,3,6)` and `(1,1,1,1)` both have formal cross-ratio one,
  but row/column unit gauges use only signs and cannot put them in one
  orbit.

**Required repair:** state the theorem over a field, or state that every
corner amplitude is a unit.

There is a second load-bearing hypothesis.  Four unrelated vertex-line
trivializations permit independent rescaling of the four coordinates and
destroy `(3)`: `(1,1,1,1)` can become `(1,1,1,2)`.  The line square must
come with the row/column factorization in `(1)` (equivalently, the tensor
identifications that make the cross-ratio line canonically trivial).

## 2. Why the two branches must not be merged

On a physically identified common coefficient line, independent
rank-one edge actions give the actual four elements

```text
w(1,alpha,beta,alpha beta),

mu=w(1-alpha)(1-beta),              kappa=1.              (5)
```

If the physical normalization permits only one common rescaling, the scalar
value changes covariantly and its nonvanishing is meaningful.  Over
`F_13`, the fixed normalization `w=1` has exactly

```text
(13-2)^2=121
```

ordered pairs with `alpha,beta` nonzero, nonidentity, and `mu!=0`.

Once independent row/column changes are gauge, additive `mu` is not
intrinsic.  The smallest hostile is over `F_5`:

```text
(1,2,2,4):       kappa=1, mu=1,
(1,1,1,1):       kappa=1, mu=0,                           (6)
```

and `(4)` carries the first vector to the second.

The warning is stronger than the `kappa=1` example.  Even the **vanishing**
of raw `mu` is not preserved when `kappa!=1`:

```text
(1,2,2,3):       kappa=2, mu=0,
(1,1,1,2):       kappa=2, mu=1.                           (7)
```

These are in the same row/column orbit.  Thus Branch B has the intrinsic
scalar `kappa`, while

```text
(kappa-1)w
```

is the additive defect in a chosen normal form (or a covariant section
after a base-corner trivialization).  It is not a gauge-invariant scalar.
The phrase “gauge-invariant same-atom face `(kappa-1)w`” in the scratch
report should be replaced accordingly.

There is a useful complete census inside the `kappa=1` orbit over `F_q`.
Write its elements uniquely as

```text
w(1,alpha,beta,alpha beta).
```

Then:

```text
number with mu=0              =(q-1)(2q-3);
number with each fixed mu!=0  =(q-2)^2;
total with mu!=0              =(q-1)(q-2)^2.             (8)
```

The exhaustions for `q=3,5,7,13` exactly reproduce `(8)`.  This is a
quantitative measure of how much additive information the projective
quotient discards.

## 3. Idempotent no-go and provenance

The field-idempotent statement passes.  Over a field, an idempotent scalar
is zero or one.  Fourfold co-support forces

```text
alpha=beta=1,
mu=w(1-alpha)(1-beta)=0.                                 (9)
```

The word “field” is essential.  In the product ring `F_2^4`, take

```text
alpha=(1,1,0,0),        beta=(1,0,1,0).
```

Then `1,alpha,beta,alpha beta` are all nonzero, while

```text
1-alpha-beta+alpha beta=(0,0,0,1) !=0.                  (10)
```

This does not contradict `(9)`; it pins its exact scope.

The linear provenance theorem also passes without positivity:

```text
Lambda(1)-Lambda(E)-Lambda(F)+Lambda(EF)
 =Lambda((1-E)(1-F)).                                   (11)
```

The companion checks `(11)` for all two Boolean masks on five atoms and all
`3^5` weights from `{-1,0,2}`, for `248,832` exact cases.  Hence a nonzero
coarse face can come from the joint-absent atom, but no linear
coarse-graining moves that contribution onto a raw fourfold-common atom.

## 4. THM-2779 application

THM-2779 fixes

```text
T e_r=e_(r+1),          M e_r=zeta^(-r)e_r.              (12)
```

With standard operator composition,

```text
(TM)e_r=zeta^(-r)e_(r+1),
(MT)e_r=zeta^(-(r+1))e_(r+1),

TM=zeta MT.                                               (13)
```

Thus the oriented path ratio `TM/MT` is `zeta`; reversing the square gives
`zeta^-1`.  Over `F_53`, the independent companion chooses the primitive
thirteenth root `zeta=10` and obtains:

```text
TM/MT=10,           reverse=16,
(1,1,1,zeta)=(1,1,1,10),
mu=zeta-1=9.                                             (14)
```

No nonzero one-dimensional scalar pair can satisfy `(13)`, since scalars
commute.  On projective basis lines, `M` is the identity and `T` is the
thirteen-cycle, so their two shadows commute.  More generally the scalar
commutator dies in `PGL`.

The application is therefore correct with this wording:

- `zeta` is the multiplier of two **lifted** paths in a projective action;
- `(1,1,1,zeta)` is its normalized scalar holonomy chart;
- it is not a literal four-corner square of one-dimensional commuting scalar
  toggles;
- the orientation must be stated;
- THM-2779 still supplies no same-ancestry physical coefficient fibre for
  the THM-2791/2806 atom.

Calling it a “projective commutator” is safe only if this means the central
multiplier of chosen linear lifts; the commutator of their images in `PGL`
is the identity.

## 5. THM-2807 application

The exact address data recheck as

```text
n0=3454614,       n+=3454627,       na=4143978,

n+-n0=13,
na-n+=169*4079=689351,
na-n0=689364,

(1,0)+(0,4079)=(1,4079).                                (15)
```

Hence every honest rank-one character has triangle boundary holonomy one.
For the explicit `F_53` character

```text
chi(Z1)=chi(Z2)=zeta=10,
```

the pure, vertical, and diagonal phases are

```text
(10,15,44),
10*15/44=1.                                              (16)
```

Here `15=zeta^4079=zeta^10`.  Completing the character algebraically to a
four-corner Segre square gives

```text
(1,10,15,44),          kappa=1,       mu=20!=0.           (17)
```

This is a useful hostile control: “no phase holonomy” does **not** mean
that every chart-dependent additive finite difference vanishes.  It means
that the honest translation boundary has `kappa=1`.  The fourth,
vertical-only vertex in `(17)` is an algebraic completion, not a claim that
THM-2807 contains a fourth selected physical cylinder.

Modulo `169`,

```text
(n0,n+,na)=(85,98,98).                                   (18)
```

The companion exhausts all `13` characters of `Z/169` valued in
`F_53^x`; every one assigns phase one to the vertical edge.  Therefore the
kernel residue

```text
4079=10 mod13
```

cannot be retained as a quotient character while simultaneously identifying
the two quotient vertices.  It requires the forgotten full-depth/kernel
coordinate.  The report's conclusion that THM-2807 supplies a flat honest
triangle but no physical allocation holonomy is correct.

## 6. Promotion verdict

The proposed result can be promoted after four wording repairs:

1. classify the field/nonzero or ring/unit locus, not arbitrary nonzero
   elements of a commutative ring;
2. require a row/column-factorized gauge, not four unrelated vertex gauges;
3. call only `kappa` an intrinsic scalar and call `(kappa-1)w` a normalized
   or covariant defect;
4. describe THM-2779's `zeta` as an oriented multiplier of linear lifts, and
   preserve the absence of a same-ancestry physical fibre.

With those repairs, the central conceptual distinction is exact:

```text
fixed common line:
    kappa=1 can coexist with physically meaningful mu!=0;

row/column projective gauge:
    raw mu is chart-dependent and kappa is the complete unit-locus orbit
    invariant.
```

Reproduction:

```bash
python3 .scratch/thm2814_projective_idempotent_hostile_audit_20260728/audit.py
python3 -O .scratch/thm2814_projective_idempotent_hostile_audit_20260728/audit.py
```

The two transcripts are byte-identical.
