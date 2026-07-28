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

## 4. The conditional unit holonomy contains the dual quarter turn

The four normalized primitive-unit polynomials on the `BABA` carrier are

```text
11+6z^5,       5z^3,       11+z^4+2z^5,       8z^3.
```

Their conditional product in `F_13[z]/(Phi_7)` is

```text
U=9+2z+8z^2+7z^3+6z^4+9z^5,                              (9)
```

of exact order `168`.  The factorization

```text
Phi_7=(z^2+3z+1)(z^2+5z+1)(z^2+6z+1)       over F_13     (10)
```

gives component orders

```text
12, 168, 28.                                               (11)
```

On the middle component `z^2+5z+1`, one has

```text
U=(4+5z),                  ord(U)=168,
U^42=8,                    8^2=-1 mod13,
8+8^(-1)=8+5=0 mod13.                                  (12)
```

Hence the order-168 sidecar contains an exact order-four quotient realizing
the same cancellation `(4)`, now with coefficient eigenvalue `i_13=8`.

The loss ledger is essential.  On the other two factors, `U^42=-1`, so their
symmetric sums equal `-2`, not zero.  The whole `Phi_7` coefficient algebra
does **not** cancel.  One must justify projection to the middle factor, and no
transport identifying the four unit rows has yet been constructed.

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

## 6. The next decisive construction

The source-fibre audit gives the appropriate hostile constraints:

```text
half carrier:
  E_3 cospans exist, but dynamic present typing deletes them before the
  reciprocal cycle;

C_221 carrier:
  E_3 cospans exist, but the forced-carry/root-6 subgenerator deletes them
  before the reciprocal cycle.
```

The quarter-turn calculation says what a successful re-rooting should retain:

1. an actual `E_3 -> Q_(3,{1,2})` source cospan;
2. the four-state endpoint/ghost or `BABA` incidence before reflection;
3. a lawful rank-two integral gain reducing to `i_17=13` on phase and
   `i_13=8` on the selected coefficient factor;
4. all three `Phi_7` components, or a proved target-preserving reason the
   two noncancelling components vanish;
5. a separate nonnegative/current argument if the goal is to invoke
   THM-2644 rather than merely build a signed mapping cone.

The cheapest positive experiment is therefore an all-carry, variable-root
`C_221` re-root with the mod-221 gain `47` retained as a local coefficient.
The cheapest hostile is to prove that dynamic present typing or the two
discarded `Phi_7` factors make every such common-gauge square empty.

## 7. Exact reproduction

Run

```bash
python3 04-computation/lrc14_paired_quarter_turn_debt_local_system_20260728.py
python3 -O 04-computation/lrc14_paired_quarter_turn_debt_local_system_20260728.py
```

Both modes byte-match
`05-knowledge/results/lrc14_paired_quarter_turn_debt_local_system_20260728.out`.

SHA-256:

```text
script  f5669cdf93b18cfe36a1edf997fe2bd2577dd2b093df620bfa97126e16fe8d24
output  02cf323945444d7b49a42dc019e81a8f9d904ed36779400527be47ad5d338f64
```

No unit transport, middle-factor projector, semantic cospan attachment,
nonnegative transition, scalar-row exclusion, or LRC(14) conclusion is
proved.
