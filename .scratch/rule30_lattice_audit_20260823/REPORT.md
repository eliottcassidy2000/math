# Rule 30 amplitude lattice: independent audit and all-period extension

**Status:** `PROVED (scratch, not canon) + VERIFIED-EXACT`.  The odd-period
sum-congruence claim passes an independent proof and exact replay.  It extends
to a complete Smith classification for every spatial period, and the even
projective boundary yields a separate nonclosing physical obstruction.

No canon file or theorem ID was edited or reserved.

## 1. Freshness and truth status

A fresh fetch during this audit put `origin/main` at `06734d90d`.  On that
live surface,
[THM-3778](../../01-canon/theorems/THM-3778-rule30-odd-period-finite-scale-cycle-projective-profile-no-go.md)
is still a **RESERVED / UNPROVED EMPTY STUB**.  The proof text visible in this
session worktree comes from the local ahead commit `1e140279a`; it is therefore
treated here as a candidate under audit, never as a proved dependency.

The elementary lattice theorem below is independent of THM-3778.  The
nonclosing physical corollary uses only current proved
[THM-3512](../../01-canon/theorems/THM-3512-rule30-van-der-put-haar-cocycle-and-profinite-automaton-boundary.md),
including its inherited no-consecutive-gap-one law.  The final join back to
the local THM-3778 proof is explicitly conditional on that proof receiving a
lawful status promotion.

## 2. Inheritance pass and concept board

- **Closest mechanism:** THM-3512's exact amplitude lift
  `U_m(s)+U_m(s+2^m)=2^(d_m)U_(m+1)(s)` and projective recurrence.
- **Canonical hostile:** at `n=3,r=2`, `A=(1,1,-2)` satisfies
  `T_3^2 A=A` and has zero total sum, but is not all odd.
- **Corrected near miss:** “even `n` is singular” is true but incomplete.
  The kernel has rank `n/2`, not merely the alternating line, and all powers
  have an exact free-plus-cyclic Smith classification.
- **Least-used sidecar:** the unprojectivized integral amplitude and its total
  sum.  Projectivization erases precisely the scalar on which this carry lives;
  THM-3511/3516's phase owner and ordinary carry are still needed by any
  operation-closed consumer.

The live concepts were:

| object | predicate | operation/invariant | lost coordinate | decisive test |
|---|---|---|---|---|
| `T_n^r Z^n` | exact integral image | Smith form and total sum | profile/owner | independent SNF |
| even `n=2^a m` | singular boundary | period embedding and pair-sum quotient | free differences | factor `T=EQ` |
| projective profile `g` | closing or nonclosing evolution | `R(g)_(j+n/2)=R(g)_j` | common amplitude | even-cycle hostile |
| physical translation ray | possible finite spatial period | no consecutive `d=1` | odd-core dynamics | constant-mode gap |
| prize interface | whether the result reaches a target | center-sequence semantics | center value/time cost | compare all three prizes |

The META-PATTERNS card used was “Audit saturation and basis covariance before
naming a lattice bridge”: the ambient module, image, kernel, cokernel,
saturation, and quotient loss are all explicit below.

## 3. Exact all-period theorem

For `n>=1`, label coordinates by `Z/nZ` and define

```text
(T_n A)_j=A_(2j)+A_(2j+1).                           (1)
```

Let `r>=0`, write

```text
n=2^a m,                   m odd,
b=min(a,r),                d=n/2^b,
s=r-b=max(0,r-a).                                    (2)
```

Let `Per_d(Z^n)` be the primitive sublattice of vectors with period dividing
`d`, so `B_(j+d)=B_j`.

### Theorem (image, cokernel, and Smith form)

One has exactly

```text
T_n^r Z^n
 = {B in Per_d(Z^n): sum_(j=0)^(d-1) B_j = 0 mod 2^s}.  (3)
```

Consequently

```text
rank(T_n^r)=d,
coker(T_n^r) = Z^(n-d) direct-sum Z/(2^s),             (4)
SNF(T_n^r)=diag(1^(d-1),2^s,0^(n-d)),                 (5)
```

where `Z/(1)=0`.  In particular, when `n` is odd, `a=b=0`, `d=n`,
and `s=r`, so

```text
T_n^r Z^n={B: sum B_j=0 mod 2^r},
coker(T_n^r)=Z/(2^r),
SNF(T_n^r)=diag(1,...,1,2^r).                         (6)
```

Thus the scout's odd-period statement and its claim that total coordinate sum
canonically generates the quotient are correct for every odd `n>=1` and
`r>=1`.

The kernel is equally explicit.  Put `q=2^b`.  Then

```text
ker(T_n^r)
 = {A: sum_(t=0)^(q-1) A_(qj+t)=0 for every 0<=j<d}, (7)
```

with basis

```text
e_(qj+t)-e_(qj),       0<=j<d, 1<=t<q.               (8)
```

It is a saturated lattice of rank `n-d`.  The Smith multiplicities can be
read without the unified notation as follows:

```text
0<=r<=a:
  rank=n/2^r,
  SNF=diag(1^(n/2^r),0^(n-n/2^r)),
  coker=Z^(n-n/2^r), with no torsion;

r>a:
  rank=m,
  SNF=diag(1^(m-1),2^(r-a),0^(n-m)),
  coker=Z^(n-m) direct-sum Z/(2^(r-a)).               (9)
```

### Proof: odd base

Let `S(A)_j=A_(j+1)` and `P_2(A)_j=A_(2j)`.  When `n` is odd, doubling is a
coordinate permutation.  Its permutation-cycle decomposition consists of
the fixed coordinate `0` and the cycles

```text
j -> 2j -> 2^2 j -> ... mod n,
```

whose length is the order of `2` modulo `n/gcd(j,n)`.  Thus `P_2` is
unimodular.  The shift `S` is one labeled `n`-cycle.  Consequently

```text
T_n=P_2(I+S),                    |det T_n|=2.          (10)
```

On the single shift cycle, the determinant is
`det(I+S)=1-(-1)^n=2`; equivalently, evaluate `x^n-1` at `x=-1`.
This permutation-cycle calculation is the exact odd/even boundary: for even
cycle length the same determinant is zero.  Also

```text
sum(T_n A)=2 sum(A).                                  (11)
```

Therefore `T_n^r Z^n` lies in the congruence lattice on the right of (6).
The two lattices have the same index `2^r`: the image index is
`|det(T_n^r)|`, while total sum modulo `2^r` is a surjection
`Z^n -> Z/(2^r)`.  Inclusion plus equal finite index proves equality.  The
quotient identification gives the cyclic cokernel and Smith form.

The constructive inverse in the scout is also correct.  For a target `B`, put

```text
D_k=B_(2^(-1)k),
A_k=(1/2) sum_(t=0)^(n-1) (-1)^t D_(k+t).             (12)
```

The finite alternating geometric series gives `(I+S)A=D`; every numerator
has parity `sum B`, so (12) is integral exactly when the one-step total sum is
even.  Repeating the unique rational inverse realizes (6) constructively.

The same index proof validates the scout's covariance extension: for odd `n`,
any product

```text
P_r(I+S) ... P_1(I+S)                                (13)
```

with arbitrary coordinate permutation matrices has the same image (6).
Each factor doubles total sum and has determinant of absolute value two.

### Proof: peel the two-part

For an even integer `k`, define the surjective pair-sum map and the period
embedding

```text
Q_k:Z^k -> Z^(k/2),      (Q_k A)_j=A_(2j)+A_(2j+1),
E_k:Z^(k/2) -> Z^k,      (E_k C)_j=C_(j mod k/2).     (14)
```

Then

```text
T_k=E_k Q_k,
T_k E_k=E_k T_(k/2),
T_k^r=E_k T_(k/2)^(r-1) Q_k.                         (15)
```

Every `Q_k` is onto: put the requested value in the even member of each pair
and zero in the odd member.  Iterating (15) peels one factor two from the
period at each scale until either the `r` operations are exhausted or the odd
core `m` is reached.  The remaining `s` operations act by `T_m^s`, whose image
is the odd-base congruence lattice already proved.  This is exactly (3).

More explicitly, let

```text
Q_(n,b):Z^n -> Z^d,
(Q_(n,b)A)_j=sum_(t=0)^(2^b-1) A_(2^b j+t),          (16)
```

and let `E_(n,d)` repeat a `d`-vector `n/d` times.  The iterated
intertwining gives the load-bearing exact factorization

```text
T_n^r=E_(n,d) T_d^s Q_(n,b).                         (17)
```

Here `Q_(n,b)` is onto and `E_(n,d)` is injective.  The middle operator is
the identity when `r<=a` and is nonsingular on the odd core when `r>a`.
Formula (17) proves the kernel formula (7)--(8), the rank, and the complete
free-rank behavior before any Smith computation is invoked.

The periodic sublattice `Per_d(Z^n)=E_(n,d)(Z^d)` is primitive: restriction
to the first `d` coordinates is a left inverse to repetition.  Hence
`Z^n/Per_d(Z^n)` is free of rank `n-d`, while the quotient inside the periodic
summand is `Z/(2^s)`.  This proves (4)--(5).

### Exact even boundary

At one even-period step,

```text
im(T_n)=Per_(n/2)(Z^n),
ker(T_n)=span_Z{e_(2j)-e_(2j+1):0<=j<n/2},           (18)
SNF(T_n)=diag(1^(n/2),0^(n/2)).                      (19)
```

Thus the earlier alternating kernel vector is valid but records only one of
`n/2` independent kernel directions.  There is no torsion at the first `a`
period-peeling steps.  Only after the odd core is reached does the single
cyclic `2`-carry appear, with exponent `2^(r-a)`.

## 4. Projective and nonclosing boundaries

For a labeled periodic profile on which the rational map is defined, write

```text
R(g)_j=-g_(2j)g_(2j+1)(1-g_(2j+2))/(1-g_(2j)).       (20)
```

If `n=2d`, direct substitution gives

```text
R(g)_(j+d)=R(g)_j.                                   (21)
```

This needs only the denominators in (20) to be nonzero.  It does not need the
stronger THM-3778 candidate saturation conditions `g_j!=0,1`.

Two consequences follow.

1. **Defined nonclosing evolution:** after each scale, one factor two is
   removed from the declared spatial period until the odd core is reached.
   Nothing forces subsequent odd-core evolution to close.
2. **Defined exact cycle:** if `R^q(g)=g` for some `q>=1`, then an even declared
   period was nonminimal.  Equation (21) makes `R^q(g)` period `n/2`, hence so
   is `g`; iterating reduces every finite cycle to its odd spatial core.

Therefore there are no genuinely minimal even-period projective cycles.  If
the local odd-period THM-3778 candidate is eventually proved and promoted,
its algebraic cycle classification automatically extends to arbitrary
declared `n` by this reduction.  Relative to current `origin/main`, that join
is **CONDITIONAL**, because THM-3778 is reserved.

### Physical nonclosing corollary: no dyadic spatial period

Retain THM-3512's proved physical profiles

```text
g^(k)_j=G_(M+k)(t+j2^(M+k)),
g^(k+1)=R(g^k),
d_l=nu_2(1-G_l(s)),                                  (22)
```

where every profile entry is an odd `2`-adic unit, each `d_l>=1` is independent
of phase, and consecutive gaps cannot both equal one.

Suppose `g^(0)` had period dividing `2^a`.  Repeated use of (21) would make
`g^a` constant, say `g^a_j=c`.  Then

```text
g^(a+1)_j=-c^2,
nu_2(1+c^2)=1                                        (23)
```

because every odd square is `1 mod 8`.  Applying the same calculation to the
next constant profile gives the next gap equal to one as well.  This
contradicts THM-3512's no-consecutive-gap-one law.

Hence:

```text
No physical Rule 30 translation-ray projective profile has spatial period
dividing a power of two.                              (24)
```

More generally, if a physical profile has any finite period `2^a m`, its
profile after `a` scales has a nonconstant odd period dividing `m`; in
particular `m>1`.  This is a nonclosing result: no scale return is assumed.
It does not exclude odd spatial periods `3,5,...` with nonclosing evolution.

## 5. Exact join to the local THM-3778 candidate

The lattice input used by the local candidate passes.  Moreover, its physical
invoice can be sharpened.  In its notation, projective return would give

```text
T_n^r A=2^D A',                 A'=uA,                (25)
```

with odd `n` and all coordinates of `A,A'` odd `2`-adic units.  Summing (25)
and using (11) gives

```text
2^r sum(A)=2^D u sum(A).                             (26)
```

Since an odd number of odd units has odd, hence nonzero, sum, (26) yields the
exact conclusions

```text
D=r,                         u=1.                    (27)
```

The candidate currently extracts only `nu_2(u)=0` and then `D=r`; the lattice
observer also fixes the projective return scalar.  It still does **not** force
`A` to be constant.  The candidate's separate spectral simplicity argument is
load-bearing for that step.

The stopping hostile is

```text
n=3, r=2, A=(1,1,-2),       T_3^2 A=A,  sum(A)=0.    (28)
```

It passes every total-sum carry test and lies on the algebraic fixed plane,
but its even coordinate violates the all-odd physical premise.  Conversely,
the ambient physical-typed one-step vector (not asserted to be realized by the
Rule 30 orbit)

```text
A=(1,5,9),        T_3 A=2(3,5,7),                    (29)
```

has a uniform positive normalizer, all-odd parent, and equal source/parent
sum, yet `(3,5,7)` is not proportional to `(1,5,9)`.  Thus the exact carry
does not manufacture projective closure.

The connection contract is:

```text
source:      labeled integral amplitude after r lifts
target:      free period-defect coordinates plus one cyclic 2-carry
map:         period projection; on the odd core, primitive-block sum mod 2^s
preserved:   integral lift existence, exactly
destroyed:   projective profile, scale, gap locations, phase owner, signed
             gauge, chronology, and center-cell values
sidecars:    lifted amplitude; THM-3511 owner; THM-3516 signed gauge/carry
test:        attach the quotient before projectivization and require closure
             under the next actual Rule 30 operation
```

## 6. Verification

The independent companion
[`rule30_lattice_all_period_audit.py`](rule30_lattice_all_period_audit.py)
imports none of the scout code and contains no Python `assert` statement.  It
checks:

- SymPy's independent exact integer Smith form for every `1<=n<=24` and
  `0<=r<=8` (SymPy `1.14.0`);
- the exact block factorization (17), rank, and all `812` displayed
  block-difference kernel basis vectors on that same universe;
- constructive preimages and free/torsion hostiles for every `1<=n<=32` and
  `0<=r<=10`;
- `T_(2m)=EQ`, the intertwining law, full pair-difference kernel, and exact
  even rank;
- arbitrary scale-by-scale coordinate permutations at odd periods through
  `n=17`, depth six;
- direct rational projective period-halving through `n=24`;
- exhaustive saturated finite-field transition hostiles for `p=5,7`,
  `n=2,4,6`, through depth three;
- the odd-square physical gap law at precisions `3..12`;
- the nonclosing control (29) and carry-only hostile (28).

Reproduce with:

```text
python3 -B .scratch/rule30_lattice_audit_20260823/rule30_lattice_all_period_audit.py
python3 -B -O .scratch/rule30_lattice_audit_20260823/rule30_lattice_all_period_audit.py
```

Both streams are identical to
[`rule30_lattice_all_period_audit.out`](rule30_lattice_all_period_audit.out).
An AST walk finds no Python `assert` node.  Raw-file SHA-256 values are:

```text
f8cada48d0b12b7d18a51b0621c23fd5af3508527dd89f938f18dc5bf298977d
  rule30_lattice_all_period_audit.py
b4e5d373680fc0cf93b9e430e7800f189428c86d7a62f0b2ddef4a45e59fc051
  rule30_lattice_all_period_audit.out
```

## 7. Prize scope and honest frontier

The official [Rule 30 prize page](https://rule30prize.org/) still lists the
three center-column problems: temporal nonperiodicity, asymptotic color
balance, and a linear computational lower bound.  This audit concerns a
labeled spatial profile of projective amplitude ratios along a translation
ray.  It proves none of those center-column predicates.

In particular, the results give no bounded innovation gaps, no finite Mealy or
signalizer closure, no center nonperiodicity, no density or balance, and no
query lower bound.  Every Rule 30 prize remains **OPEN** in the repository.
The honest next mathematical frontier is the odd-period **nonclosing** profile
with the unprojectivized amplitude, phase owner, signed gauge, and next-operation
closure retained.

## 8. Promotion recommendation

Two statements are promotion-ready after an independent prose audit:

1. the unconditional elementary all-period image/cokernel/Smith theorem
   (3)--(5), including the odd permutation-covariant corollary; and
2. the THM-3512-dependent physical exclusion of power-of-two spatial periods
   (24), which assumes no scale cycle.

The extension of the local THM-3778 cycle classification to all declared
periods should be promoted only together with, or after, lawful promotion of
THM-3778 itself.  It must remain conditional while the live theorem is a
reserved empty stub.
