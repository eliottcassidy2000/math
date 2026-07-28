# The double debt has a minimal quarter-turn local system

**Status: PROVED-ELEMENTARY ALGEBRA + VERIFIED-EXACT CONNECTION; PHYSICAL
APPLICATION OPEN.**  This note connects the coefficient-two debts in the
intrinsic `BABA` carrier and THM-2706's endpoint/ghost carrier with the
quarter torsion of the transverse `C_221` stalk.  It constructs a signed
local-coefficient cancellation, not a nonnegative semantic current, row
exclusion, or LRC(14) proof.

## 1. One incidence operator appears twice

Two independently derived recurrent four-cycles have the same integral
shadow:

```text
BABA guard debt:       (B_*+B_*^3)[A]=2[B];
THM-2706 ghost debt:   (B_*+B_*^3)[E]=2[G].              (1)
```

Before reflection quotient, both maps are the all-ones matrix

```text
[[1,1],
 [1,1]].                                                   (2)
```

Boolean support remembers only nonvanishing.  Reduction modulo two kills the
invariant quotient for the wrong reason.  Over the nonnegative integers the
two invoices remain load-bearing.

The useful reframe is to retain the two path orientations as a gain.  If `u`
is one step around the underlying four-cycle, then the symmetric debt
operator is

```text
u+u^(-1).                                                 (3)
```

## 2. Minimal integral cancellation lemma

Let `L` be a nonzero free abelian group and let `J` be an automorphism of
`L`.  Then

```text
J+J^(-1)=0       iff       J^2=-I.                        (4)
```

Indeed, multiply the first identity by `J`; the converse is immediate.
Thus a local system cancels `(3)` exactly when one step acts as an integral
complex structure.

This is impossible in rank one because

```text
GL_1(Z)={+1,-1},             J+J^(-1)=+2 or -2.           (5)
```

It is possible in rank two, minimally, with

```text
J=[[0,-1],
   [1, 0]],                 J^2=-I.                       (6)
```

More generally `(4)` makes `L tensor Q` a vector space over `Q(i)`, so the
integral rank is even.  The coefficient-two obstruction is therefore not
merely “parity”: its smallest faithful integral resolution is a rank-two
quarter-turn local system.

This is a signed-chain statement.  There is no `J`-equivariant augmentation
from `(6)` to the trivial nonnegative rank-one coefficient system, precisely
because `(5)` would restore the two units.  Hence `(4)` does not evade
THM-2644's nonnegativity and common-middle requirements.

There is a sharper pointed-cone obstruction.  Let `K` be a pointed cone in
`L tensor R`, so `K intersection (-K)={0}`, and let `J` be invertible.  There
are two useful forms:

```text
Jv,J^(-1)v in K and v!=0       => (J+J^(-1))v!=0;        (6a)

v in K, J(K) subset K, v!=0    => (J+J^(-1))v!=0.        (6b)
```

For `(6a)`, cancellation would put `Jv=-J^(-1)v` in both `K` and `-K`,
forcing `Jv=0`.  For `(6b)`, multiply a hypothetical cancellation by `J`:
`J^2v=-v`; forward invariance puts `J^2v` in `K`, so pointedness forces
`v=0`.  In particular `J^2=-I` cannot preserve a nonzero pointed cone, even
under the one-sided condition `J(K) subset K`.  If `J(K)=K`, both proofs
apply.  Thus quarter-turn cancellation is intrinsically signed; it cannot
itself be the one nonnegative transition required by THM-2644.

## 3. The delayed phase already carries the quarter turn

On the denominator-seventeen orbit, multiplication by `13` has exact order
four:

```text
4 -> 1 -> 13 -> 16 -> 4,
13^2=-1 mod17.                                            (7)
```

The two word endpoints are `{4,13}` and the two forbidden ghost midpoints are
`{1,16}`.  Applying `B` and `B^3` from either endpoint reaches both ghosts,
giving exactly `(2)`.  Moreover

```text
13+13^(-1)=13+4=0 mod17.                                 (8)
```

Thus the transverse phase is the mod-17 eigenline of the integral complex
structure `(6)`, with chosen eigenvalue `i_17=13`.

## 4. The conditional unit holonomy contains the dual quarter turn, but only after a lossy projection

The four normalized primitive-unit polynomials on the `BABA` carrier are

```text
11+6z^5,       5z^3,       11+z^4+2z^5,       8z^3.
```

Their conditional product in `F_13[z]/(Phi_7)` is

```text
U=9+2z+8z^2+7z^3+6z^4+9z^5,                              (9)
```

“Conditional” is load-bearing: the four displayed unit rows do not yet have
a proved common physical transport or gauge.  Equation `(9)` is their product
after identifying the coefficient algebras by hand; it is not yet the
holonomy of a constructed physical four-edge local system.

The product has exact order `168`.  The factorization

```text
Phi_7=(z^2+3z+1)(z^2+5z+1)(z^2+6z+1)       over F_13     (10)
```

is a product of three distinct irreducible quadratics.  Therefore the full
six-dimensional algebra has the CRT decomposition

```text
A=F_13[z]/(Phi_7)  ~=  K_3 x K_5 x K_6,
K_a=F_13[z]/(z^2+a z+1),             dim_F13(K_a)=2.    (10a)
```

The component orders of `U` are

```text
12, 168, 28.                                               (11)
```

On the middle component `z^2+5z+1`, one has

```text
U=(4+5z),                  ord(U)=168,
U^42=8,                    8^2=-1 mod13,
8+8^(-1)=8+5=0 mod13.                                  (12)
```

Hence the middle component contains an exact order-four quotient realizing
the same scalar cancellation `(4)`, with eigenvalue `i_13=8`.  The global
calculation is different.  In the polynomial basis `1,z,...,z^5`,

```text
U^42     =(8,0,10,8,8,10),
U^(-42) =(5,0,11,1,1,11),
S=U^42+U^(-42)=(0,0,8,9,9,8) != 0.                      (12a)
```

Under `(10a)`, this surviving symmetric debt is

```text
S=(-2,0,-2) in K_3 x K_5 x K_6.                         (12b)
```

The middle projection `pi_5:A->K_5` has rank two and the four-dimensional
kernel `K_3 x K_6`.  It sends the nonzero `S` in `(12a)` to zero, while both
discarded components remain `-2`.  Thus middle cancellation is a quotient
loss: the projection kills the whole displayed obstruction rather than
proving it vanishes.  No target-preserving, gauge-covariant physical
projector to `K_5` is known.

There is a second type loss.  `U^42` is the order-four holonomy extracted
from the conditional product of four rows; it is not an edgewise coefficient
automorphism `J`.  Identifying it with the one-step `J` of Section 2 would
require four actually transported edge gains, a common coefficient fibre,
and a proof that their cycle product is `(9)`.  None is constructed here.

## 5. Why the modulus is exactly `221`

The characteristic polynomial of `(6)` is `X^2+1`.  Its relevant roots are

```text
mod 13: {5,8},              choose 8 from (12);
mod 17: {4,13},             choose 13 from (7).           (13)
```

CRT combines the two chosen eigenvalues uniquely:

```text
i_221=47,
47=8 mod13,                 47=13 mod17,
47^2=-1 mod221.                                             (14)
```

All four square roots of `-1` modulo `221` are

```text
21,47,174,200.                                             (15)
```

Thus the already-discovered `C_221` stalk is not only the smallest state
remembering delayed mod-17 phase together with the old mod-13 carry/root
coordinate.  It is also the joint modulus on which the delayed and
coefficient quarter-turn eigenlines are reductions of one integral rank-two
local system.

This is an exact structural alignment, not yet a physical map.  The affine
`C_221` state update is not multiplication by the unit `47`; its multiplier
`13` is noninvertible modulo `221`, and its translation is state-dependent.
Equation `(14)` therefore supplies a candidate coefficient gain, not a new
chronology arrow.

## 6. The source-fibre wall is globally typed, not a failure of root search

The full source-fibre audit first establishes genuine strict scalar cospans
`E_3 -> D^6 Q_(3,{1,2})` on both the half and `C_221` fibres.  Their failure
to attach to a present packet has a uniform structural reason.  In the exact
present builder
`04-computation/lrc14_replica_dichotomy_typed_row_opus_20260727.py::build_F`,
every `F_(ell,s)` explicitly subtracts the unshifted `c_3`-danger comb:

```text
F_(ell,s) subset (D_(c_3))^c,
E_3 subset D_(c_3),
therefore E_3 intersection F_(ell,s)=empty.              (15a)
```

This holds for every `ell,s`, independently of phase, carry, private root,
or clock.  The finite census did not merely miss a favourable root; it
detected this global source/present type contradiction.

The exact finite ledgers show how the contradiction manifests:

```text
half carrier:
  each phase has 12,848 E_3 rows retaining shallow clock, rail, and private
  root when present is removed; the union of all 13 present labels is empty;

C_221 carrier:
  each full phase has 131,752 E_3 rows; the inherited root-6 carries 3/9
  contain none; the exhaustive all-carry re-root leaves only carries
  (0,6) over 4/17 and (6,12) over 13/17, with nonzero private roots 1/12;
  53,227 or 53,590 rows per surviving class pass clock/sector typing, but
  every one fails present and no reciprocal pair remains.                (15b)
```

## 7. The next decisive construction

The quarter-turn calculation says what a successful replacement should
retain:

1. an actual `E_3 -> Q_(3,{1,2})` source cospan;
2. the four-state endpoint/ghost or `BABA` incidence before reflection;
3. a lawful rank-two integral gain reducing to `i_17=13` on phase and
   `i_13=8` on the selected coefficient factor;
4. all three `Phi_7` components, or a proved target-preserving reason the
   two noncancelling components vanish;
5. a separate nonnegative/current argument if the goal is to invoke
   THM-2644 rather than merely build a signed mapping cone.

The all-carry, variable-root experiment is now saturated by `(15a)--(15b)`.
The cheapest positive experiment is a relative-present or present-free
mapping cone that attaches the scalar `E_3` current before rebuilding a
compatible unit sidecar.  It must alter the present object itself, not merely
its phase/carry/root labels.  In parallel, any use of the middle `Phi_7`
factor must construct a target-preserving projector and show why the
`K_3 x K_6` obstruction is physically null.

## 8. Exact reproduction

Run

```bash
python3 04-computation/lrc14_paired_quarter_turn_debt_local_system_20260728.py
python3 -O 04-computation/lrc14_paired_quarter_turn_debt_local_system_20260728.py
```

Both modes byte-match
`05-knowledge/results/lrc14_paired_quarter_turn_debt_local_system_20260728.out`.

SHA-256:

```text
script  abf625bd6ba7c874e4f8b0fcddc70c8ee269277116fcf5c1365aa0a4b3b21fb6
output  c8a8245272bc943dce12dcbb6483432d46d538e695f3f7e5e286a03171696bd1
```

No unit transport, middle-factor projector, semantic cospan attachment,
nonnegative transition, scalar-row exclusion, or LRC(14) conclusion is
proved.
