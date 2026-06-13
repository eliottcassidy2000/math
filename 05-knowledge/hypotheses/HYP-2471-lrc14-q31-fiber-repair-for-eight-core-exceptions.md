# HYP-2471 - Q31 fiber repair for the HYP-2470 eight-core exceptions

**Status:** OPEN addendum with exact bounded certificates.

**Source:** codex-2026-06-13.  Addendum to HYP-2470; extends HYP-2465 and supports HYP-2469 / OPEN-Q-087.

**Computation:**

- `04-computation/lrc14_below_nine_core_q27_budget5_codex.py`
- `04-computation/lrc14_q31_exception_probe_codex.py`

Stored outputs:

- `05-knowledge/results/lrc14_below_nine_core_q27_budget5_codex.out`
- `05-knowledge/results/lrc14_q31_exception_probe_codex.out`

## Claim

HYP-2470 is the primary corrected eight-core theorem: in the carry window, a
primitive replacement row retaining at least eight speeds of

```text
CORE = 7*{1,...,12}
```

has either a Q27 witness or a plain-shell witness `q<=41`.

This addendum records an independent fibered repair for the same finite
exceptions.  Let

```text
Q31 = {d*m : d in {1,2,7,14}, 1 <= m <= 31} \ {1}.
```

Then every bounded replacement row retaining at least eight core speeds has a
Q31 witness.  Equivalently, the only four-deletion rows that evade Q27 in the
HYP-2470 census are not genuine no-Q31 rows.

The shell-41 repair is smaller and more transparent for the LRC predicate.  The
Q31 repair is useful because it preserves the divisor/fiber ladder grammar used
by HYP-2444, HYP-2465, and the Church-style retained-channel framing of
HYP-2469.

## Exact Evidence

HYP-2465 already proves the stronger Q27 statement for delete count `e<=3`:

```text
e=0:   1/1   infeasible
e=1:  12/12  infeasible
e=2:  66/66  infeasible
e=3: 220/220 infeasible
```

The new `e=4` budget-5 Q27 scan tests all `C(12,4)=495` deletion sets.  First
pass:

```text
Q27, e=4, budget=5:
  infeasible = 489
  feasible   = 2
  unknown    = 4  (20s limit)
  obligation range = 504..1274
```

The four time-limit cases resolve as Q27-infeasible with a longer limit:

```text
(14,21,28,70)  infeasible
(21,28,56,70)  infeasible
(21,28,70,84)  infeasible
(28,42,70,84)  infeasible
```

The only two Q27-feasible deletion shapes are exactly the HYP-2470 exceptions:

```text
D=(28,42,56,84)
D=(42,56,70,84)
```

Example Q27 blockers exist:

```text
D=(28,42,56,84), A=(91,322,350,504,832)
row=(7,14,21,35,49,63,70,77,91,322,350,504,832)
first_Q150=(31,10), min_q<=600=31, Bprime opens at 322.

D=(42,56,70,84), A=(91,119,156,700,1008)
row=(7,14,21,28,35,49,63,77,91,119,156,700,1008)
first_Q150=(31,5), min_q<=600=31, Bprime opens at 700.
```

But those two deletion shapes become infeasible when the obligation ladder is
widened from Q27 to Q31:

```text
D=(28,42,56,84): Q31 infeasible, 840 obligations.
D=(42,56,70,84): Q31 infeasible, 672 obligations.
```

Since Q27-infeasibility implies Q31-infeasibility, the Q31 ladder closes all
`e=4` deletion sets.

## Meaning

The naive hope "eight-core rows already force Q27" is false.  HYP-2470 found
the right corrected theorem by adding missing plain shells through `41`.

The Q31 computation says the same two failures are also explained by a retained
fiber channel: Q27 forgot the next few `m`-levels in the `{1,2,7,14}` ladder.
Once `m=28..31` is admitted, five added speeds cannot cover all obligations in
either exceptional deletion address.  The examples already expose plain
`q=31` witnesses and Bprime openings; Q31 infeasibility proves those openings
are structural rather than artifacts of a chosen MILP packet.

Thus the HYP-2469 proof grammar has two compatible repairs:

```text
scalar Q27 quotient passes in two shapes
+ plain shell<=41 side channel
=> HYP-2470 finite exception theorem

scalar Q27 quotient passes in two shapes
+ retained Q31 fiber ladder
=> HYP-2471 fibered repair
```

## Updated Descent Portal

Inside the carry window, a primitive row with no Q27 witness and no plain
`q<=41` witness must satisfy one of:

1. delete at least five speeds from `CORE`;
2. leave the replacement-row normal form;
3. trigger an owner/Bprime or low-clock exception;
4. descend to AP, Vstar, nonprimitive 2AP, or another named wall atom.

The next bounded computation is therefore delete count `e=5`, but it should not
be attempted as raw Q27/Q31 cardinality only.  Sample `e=5` Q27 MILPs already
show timeout pressure.  The right next carrier should keep typed budgets:

```text
13-clock spend,
divisor-fiber spend,
shell-27 class spend,
owner/Bprime opening,
support-load maximum,
low-clock escape.
```

## Tournament Analysis

The budget-5 script uses proof portals as vertices.  The observable is:

```text
(exactness, support retention, LRC leverage, side-channel, computability, descent)
```

The switch/gauge ranks portals by retained-channel strength, not scalar
denominator size.  The tie Hamiltonian path is declaration order.  Stored
fingerprint:

```text
score histogram = {0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed 3-cycles = 0
edge flips versus declaration = 1/21
leader = Q27_budget5_setcover
```

Assumption challenge: candidate vertices considered were runners, deleted core
speeds, added speeds, Q27 twists, Q31 twists, support-load obligations, Bprime
owner exits, outside-window carry fibers, and proof obligations.  The selected
quotient preserves the finite descent decision tree.  It destroys raw time
geometry and arbitrary speed magnitude data.
