---
id: THM-2710
title: "Central-half literal-word nilpotence and prescribed-clock invisibility"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the
  canonical typed row, the changed delayed base B1(y)={13y+1/2} agrees
  identically with B0(y)={13y} at all three prescribed semantic clocks
  2,4,6.  Every six-state literal singleton-word cylinder in {Q_a,Q_b}
  under B1 is empty: an a state makes the unshifted source c1 dangerous
  two states later, while two initial b states make speed 14 dangerous at
  state five.  The exact 64-word partition is 32+16+8+4+4, and a strict
  bbba cylinder is nonempty.  The half and C_221 packet cycles both fail
  the literal endpoint at the unshifted-c1-safe bit and fail the required
  exclusive source leg.  The half phase does give the genuine predicate
  identities D_c1=B1^(-1)N_1 and D_c2=B1^(-3)N_1, moving both odd-delay
  debts to the radius-1/14 tooth about 1/2 inside guard safety.  It does
  not resolve the integral debt: on the doubled BABA reflection quotient
  the unsigned boundary is multiplication by 2 on each half-turn sheet.
  Thus unchanged literal endpoints are closed; an enriched semantic C2
  bibundle or direct macro cospan remains open.  No scalar row, endpoint
  current, row exclusion, or LRC(14) conclusion follows.
source: root/central-half-semantic-gate-audit-2026-07-28
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2461-temporal-blocker-word-cocycle-and-diagonal-polarized-repair-boundary
  - THM-2698-central-half-odometer-full-local-cycle-and-semantic-sidecar-boundary
  - THM-2701-literal-singleton-word-one-step-dilation-nilpotence
related:
  - THM-2623-guard-safe-danger-cospan-and-residual-unit-wall
  - THM-2644-odd-torsor-purity-return-gate-and-nonlinear-fixed-branch-decoder
  - THM-2693-odometer-skew-product-three-event-escape-and-uniform-delayed-depth-four-nilpotence
script: 04-computation/lrc14_phase_cycle_semantic_gate_probe_20260728.py
output: 05-knowledge/results/lrc14_phase_cycle_semantic_gate_probe_20260728.out
reflection: 07-reflections/lrc14-half-and-c221-cycles-miss-the-semantic-word-gate-20260728.md
script_sha256: a34618cf8d2d7266db44c750bb56ef4c40922888f92ec0f31a4619049b09872a
output_sha256: e6e5c6db41aad1e20c977eda3499875b34ede56c56a35b281eb79e2c068714d6
hash_basis: LF-normalized bytes
---

# THM-2710 -- central-half literal-word nilpotence and prescribed-clock invisibility

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The central half-turn in THM-2698 changes the one-step delayed base and
creates a strict full local packet cycle.  It does not change the literal
THM-2305 endpoint category.  This theorem proves both sides of that boundary:

```text
odd intermediate states:   the half phase moves narrow danger to a safe
                           half-tooth and can repair a local packet;

prescribed even endpoints: the half phase cancels identically, and every
                           six-state literal word still dies.             (1)
```

The first line is a real predicate-level bridge.  The second is the semantic
stopping theorem.  Neither line constructs the enriched `C_2` bibundle that
would be needed to turn THM-2698 into an endpoint-current transition.

## 1. Literal words and the two bases

Work on the canonical typed row

```text
(H,q1,...,q5,c1,c2,c3)
 =(1,14,27,40,53,66,13,13^3,2*13^5).                    (2)
```

For `x in T=R/Z`, use the representative `centered(x)` in
`[-1/2,1/2)`, and put

```text
D_v={x:centered(vx) in [-1/14,1/14)},
C_H={x:centered(x) notin [-1/7,1/7)},
A0=C_H intersection intersection_(i=1)^5 D_(q_i)^c.      (3)
```

With source owner `c1`, write `a=c2`, `b=c3`.  The literal singleton
terminal words of THM-2305 and THM-2701 are

```text
Q_a=A0 intersection D_(c2) intersection D_(c1)^c
       intersection D_(c3)^c,

Q_b=A0 intersection D_(c3) intersection D_(c1)^c
       intersection D_(c2)^c.                            (4)
```

In particular, the `c1` complement in (4) is unshifted and load-bearing.
Define

```text
B0(y)={13y},                 B1(y)={13y+1/2}.             (5)
```

For a word `w=w0...w_(d-1)` in `{a,b}^d`, its literal `B1` cylinder is

```text
Q_1[w]=intersection_(j=0)^(d-1) B1^(-j)Q_(w_j).          (6)
```

These terminal cylinders are already weaker than the full THM-2305 source
span.  Proving that they cannot recur is therefore a valid stopping result,
not a positive semantic construction.

## 2. The half phase is invisible at every prescribed semantic clock

Iteration of (5) gives

```text
B1^k(y)
 ={13^k y+(1/2)(1+13+...+13^(k-1))}
 ={13^k y+(13^k-1)/24}.                                 (7)
```

Since `13^2=1 mod24`, the affine constant in (7) is integral for every
even `k`.  The three blocker valuations are `1,3,5`, so their prescribed
THM-2305 clocks `k_j=lambda_j+1` are `2,4,6`, and

```text
(13^2-1)/24=7,
(13^4-1)/24=1190,
(13^6-1)/24=201117,

B1^k=B0^k                       for k in {2,4,6}.         (8)
```

Consequently, for unchanged source and terminal sets, every prescribed
atom-to-word entry

```text
mu(P_tau intersection E_j intersection B^(-k_j)Q_(j,sigma))             (9)
```

is identical for `B=B0` and `B=B1`.  This statement is deliberately about
the literal endpoint category.  An enlarged parity-labelled semantic fibre
could change the objects in (9), but no such fibre is supplied by (8).

## 3. An `a` state creates a two-step source debt

The even-clock identity already gives

```text
B1^2=B0^2,
c1 B1^2(y)=13*13^2 y=c2 y mod1.                          (10)
```

If `y in Q_a`, then `c2 y` is dangerous.  Equation (10) says that the
unshifted `c1` factor is dangerous at `B1^2y`, whereas every word in
`Q_a union Q_b` requires it safe.  Thus

```text
Q_a intersection B1^(-2)(Q_a union Q_b)=empty.           (11)
```

The conflict is exact at the half-open endpoints because the two phases in
(10) are equal on the circle; no interval relaxation is used.

## 4. Two initial `b` states create a state-five speed-14 failure

Suppose

```text
y in Q_b,                         B1y in Q_b,
z=B1^5y.                                                   (12)
```

For odd exponent five, (7) gives

```text
z={13^5y+1/2}.                                             (13)
```

Since `c3=2*13^5`, the two target-`b` conditions in (12) become

```text
centered(2z) in [-1/14,1/14),
centered(26z) in [-1/14,1/14).                            (14)
```

Indeed the half shifts disappear after multiplication by the even
coefficients `2` and `26`.  The first condition lets us write

```text
z=c+e mod1,
c in {0,1/2},                  -1/28<=e<1/28.             (15)
```

Now `26c` is integral and

```text
-13/14<=26e<13/14.                                      (16)
```

No nonzero wrap in (16) can enter the half-open interval in (14).  At the
negative extreme, adding one gives the excluded right endpoint `1/14`; the
positive extreme is itself excluded.  Hence

```text
-1/14<=26e<1/14,
-1/(28*13)<=e<1/(28*13).                                (17)
```

Since `14c` is integral,

```text
-1/26<=14e<1/26,                  1/26<1/14,             (18)
```

so speed `14` is dangerous at `z`.  Every literal word keeps speed `14`
safe, and therefore

```text
Q_b intersection B1^(-1)Q_b
    intersection B1^(-5)(Q_a union Q_b)=empty.            (19)
```

## 5. Every six-state literal word is empty

Let `w` have length six.  If the first `a` occurs at one of positions
`0,1,2,3`, equation (11), shifted to that position, kills the cylinder two
states later.  Otherwise its first four letters are `b`, and (19), using
only the first two, kills the sixth state.  This partitions all `64` words:

```text
first a at 0 -> c1 failure at 2:      32;
first a at 1 -> c1 failure at 3:      16;
first a at 2 -> c1 failure at 4:       8;
first a at 3 -> c1 failure at 5:       4;
first four letters b -> speed 14 at 5: 4.                (20)
```

Thus

```text
Q_1[w]=empty                    for every w in {a,b}^6,   (21)
```

and the literal `B1` language has no recurrent strongly connected
component.

The conclusion is not a one-state or two-state vacuity.  Exact direct
substitution shows that

```text
y0=12894291/80000000                                  (22)
```

lies in the strict cylinder `Q_1[bbba]`.  The minimum factor margins at its
four states are respectively

```text
10260037/560000000,
 1272987/560000000,
 6115129/280000000,
 3083243/560000000,                                    (23)
```

all positive.  No five-state witness is asserted, so (21) is an upper
nilpotence bound, not a claimed sharp index.

## 6. Exact typing of the half and `C_221` packet cycles

The delayed projections of the displayed packet cycles are

```text
THM-2698 half cycle:        11/24 -> 11/24,
transverse C_221 control:    4/17 <-> 13/17.              (24)
```

At all four listed occurrences the guard and all five ordinary speeds are
strictly safe, `c2` is dangerous, and `c3` is safe.  But `c1` is dangerous,
with distances

```text
1/24                         on the half cycle,
1/17                         on the C_221 control.        (25)
```

Both are strictly below `1/14`.  Thus these points fail only the unshifted
`c1`-safe bit of literal `Q_a`; they are the debt-enriched local packet, not
the endpoint in (4).

There is a second, independent source-side failure.  None of the four
physical event centres lies in any exclusive source `E_j`.  Ordinary
six-step dilation reaches `11/24`, `11/24`, `4/17`, `13/17`, respectively,
and each terminal point is in the deepest-owner fork `Q_(3,{1,2})`.  But its
starting centre is not in `E_3`.  Hence the terminal label does not supply
the source leg of the prescribed-clock span.

The `C_221` control is finite exact evidence, not a promoted packet theorem.
Its two phase labels `114,107` are inverse modulo `221`.  The projected
nonnegative weight has

```text
(M,E,delta,R,mu(0)^2)=(2,2,2,2,0),                       (26)
```

the equality/backtracking boundary of THM-2644, not its strict odd-torsor
gate.  The half cycle projects to the nontrivial element of `C_2`, where
reversal fixes that element.  This is THM-2644's sharp even-group hostile;
neither projection constructs a common physical semantic middle fibre.

## 7. The genuine `N_0 -> N_1` bridge

Let

```text
H(y)=y+1/2,
N_0={y:centered(y) in [-1/14,1/14)},
N_1=H(N_0)={y:centered(y-1/2) in [-1/14,1/14)}.           (27)
```

Thus `N_1=[3/7,4/7)` on the circle and lies strictly inside the broad
guard-safe set.  Since `13` is odd,

```text
B1^j=H^j B0^j.                                            (28)
```

For both odd delays relevant to source owner `c1`, (28) gives the global
half-open identities

```text
1_(D_c1)(y)=1_(N_0)(B0y)=1_(N_1)(B1y),
1_(D_c2)(y)=1_(N_0)(B0^3y)=1_(N_1)(B1^3y).               (29)
```

This is the exact predicate mechanism by which the half phase repairs the
raw guard conflict: it moves both narrow invoices from the tooth about zero
to the tooth about one half, which is guard-safe.  At `y=11/24`, `B1` fixes
the point, its distance from `1/2` is `1/24`, and (29) is realized strictly.
At the same point, however, the unshifted `c1` distance is still `1/24`, so
literal `Q_a` still fails exactly as in (25).  THM-2698 retains a nonzero
shifted source-clock complement instead.

Equation (29) preserves the two danger predicates.  It does not preserve the
future base point, the other guard/unit factors, a physical rail, a clock
label, or the terminal-word source leg.  It is therefore a predicate bridge,
not a semantic intertwiner.

## 8. The half sheet does not kill the integral debt coefficient

The distinction is already exact on the THM-2701 BABA debt orbit.  Write
`beta=B0` for the untwisted map, and label the four delayed points

```text
b0=7/170,       a1=91/170,       b2=163/170,       a3=79/170.             (30)
```

Reflection `J(y)=-y` swaps `b0<->b2` and `a1<->a3`.  Before quotienting,
the two odd-delay debt maps are

```text
a1 --beta--> b2,           a1 --beta^3--> b0,
a3 --beta--> b0,           a3 --beta^3--> b2.             (31)
```

Hence their unsigned matrix, with rows `(b0,b2)` and columns `(a1,a3)`, is

```text
beta_*+beta_*^3 = [[1,1],[1,1]].                          (32)
```

On the reflection coinvariant class, (32) is multiplication by `2`.

Now double the orbit by the independent half-turn sheet `epsilon in C_2`,
representing `(y,epsilon)` by `H^epsilon y`.  The two involutions `J` and
`H` commute but are distinct; generically they generate `V_4`.  Under the
changed base,

```text
B1(y,epsilon)=(beta y,epsilon+1).                         (33)
```

Both delay one and delay three are odd, so they land on the **same** toggled
half sheet.  Equations (31)--(33) give, after reflection quotient,

```text
partial(A,epsilon)=2(B,epsilon+1).                        (34)
```

Thus the central half sheet permutes the two `C_2` cokernel copies; it does
not cancel the coefficient `2`.  On the two characters of the half sheet,
the boundary is multiplication by `+2` or `-2`.  Reducing (34) modulo two
would erase precisely the nonnegative multiplicity that a lawful endpoint
current must retain.

This rules out the tempting identification of THM-2698's central parity with
an already constructed Bockstein resolving the BABA debt.  A genuine repair
still needs an `H`-stable free semantic fibre and physical intertwiners, or a
direct macro endpoint-current cospan which accounts for (34).

## 9. Consequence and scope

The proved boundary is

```text
unchanged literal THM-2305 endpoints under B1:
  prescribed clocks 2,4,6 unchanged;
  every six-state Qa/Qb cylinder empty;
  half/C221 centres fail source and terminal typing.       CLOSED

enriched half-phase category:
  N0 danger is retyped as guard-safe N1;
  local packet cycles exist;
  integral debt remains coefficient two per H sheet;
  no H-stable semantic bibundle or direct macro current.   OPEN          (35)
```

No scalar cover is asserted for the typed row.  No endpoint current, row
exclusion, or LRC(14) conclusion follows; the scalar ledger remains `165`.

## 10. Exact companion

Run

```bash
python3 04-computation/lrc14_phase_cycle_semantic_gate_probe_20260728.py
python3 -O 04-computation/lrc14_phase_cycle_semantic_gate_probe_20260728.py
```

Both executions byte-match

```text
05-knowledge/results/lrc14_phase_cycle_semantic_gate_probe_20260728.out
```

with LF-normalized hashes

```text
script: a34618cf8d2d7266db44c750bb56ef4c40922888f92ec0f31a4619049b09872a
output: e6e5c6db41aad1e20c977eda3499875b34ede56c56a35b281eb79e2c068714d6
```

The dependency-free companion uses only `Fraction` and finite exact group
arithmetic.  It checks the two affine cycles, every displayed factor bit,
all four physical source/terminal types, the even-clock constants, the exact
`64`-word certificate partition, the strict `bbba` witness and margins, and
the two THM-2644 hostile controls.  Independent normal and optimized runs
produce the stored output hash.

QED.
