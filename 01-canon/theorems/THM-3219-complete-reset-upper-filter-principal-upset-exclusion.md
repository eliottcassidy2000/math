---
id: THM-3219
title: "Complete reset upper-filter principal-upset exclusion"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  In the support-(1,3), bank-I2 selector model, every nonempty physical
  completion above THM-3209's complete quotient reset has strictly negative
  degree-five principal-upset response.  The exact value is
  -1440 times the fifth power sum of the added positive pole multiset.  Hence
  the whole 64-state principal filter, through depth sixteen, meets every
  degree-at-least-five selector cone only at the reset atom.
source: root/2026-08-02
audit: >
  The exact companion pins the transitive response helper and THM-3209,
  verifies the base complete/power functional rows, exhausts all 63 nonempty
  legal completions, and checks the formula, strict signs, depth census,
  sharp minimum, and the 25 multi-negative-alphabet cases where h_5[-tau]
  is nonzero.  An independent hostile audit rederived the common-prefix
  identities and every physical count.  Normal and `-O` replay byte-match the
  stored output.
depends_on:
  - THM-3209-depth-eight-complete-quotient-reset-and-negative-singleton-tangent
script: 04-computation/gmc_complete_reset_upper_filter_thm3219.py
output: 05-knowledge/results/gmc_complete_reset_upper_filter_thm3219.out
script_sha256: 57376e23c9774302c9b8bba0327b8edde2a00423b636c85c841b3d8a6d27a51d
output_sha256: efcf5353a06be238cabc1f788be36bb84a3e28e9e83f1758d92fad842be36c23
hash_basis: LF-normalized bytes
---

# THM-3219 -- complete reset upper-filter principal-upset exclusion

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3209 proves that the complete quotient reset is isolated along its four
one-extra-pole directions.  The isolation is much larger: one degree-five
facet excludes every nonempty completion of the reset inside the entire
physical prefix bank.  The mechanism is algebraic common-prefix invariance,
not a 63-case coincidence.

## 1. Bank and response

Retain THM-3209's reduced pole multiset and dominant quotient alphabet

```text
P=(8,7,6,5,5,4,4,3,3,2,2,2,1,1,1,1),
Q=(1,3,3,4,5,6,7,8).                                  (1)
```

Their multiset difference is

```text
R=P-Q=(1,1,1,2,2,2,4,5).                              (2)
```

For a physical prefix `sigma<=P`, write `Q^sigma=Q-sigma`.  The inherited
partition response at degree `N` is

```text
G_N^sigma(lambda)
 =Phi^sigma(h_N)m_lambda[Q^sigma]
  -Phi^sigma(m_lambda)h_N[Q^sigma].                    (3)
```

The principal coarsening upset at degree five is the singleton `{(5)}` and
`m_(5)=p_5`.

## 2. Common-prefix invariance at the first live degree

Before a common prefix is removed, the exact bank functional satisfies

```text
Phi(h_j)=0                 (0<=j<5),
Phi(h_5)=1440,
Phi(p_j)=0                 (1<=j<=5),
Phi(1)=0.                                                (4)
```

The first line is the order-five zero of the reduced numerator; the power-sum
line is the row-marginal cancellation already used by THM-3209.

For any finite virtual prefix `sigma`, the lambda-ring identity

```text
h_5[A-sigma]=sum_(j=0)^5 h_(5-j)[A] h_j[-sigma]         (5)
```

and `(4)` give

```text
Phi^sigma(h_5)=1440.                                    (6)
```

Likewise `p_5[A-sigma]=p_5[A]-p_5[sigma]`; because both
`Phi(p_5)` and `Phi(1)` vanish,

```text
Phi^sigma(p_5)=0.                                       (7)
```

Thus the numerator factor at the first live degree is invariant under every
common prefix, not only under prefixes in the physical bank.

## 3. Exact upper-filter formula

Let `tau` be any nonempty multiset of positive integers and put

```text
sigma=Q+tau.                                             (8)
```

Then `Q^sigma=-tau`.  Plethystic negation gives

```text
p_5[-tau]=-p_5[tau]=-sum_(r in tau) r^5.                (9)
```

Substituting `(6),(7),(9)` into `(3)` at `lambda=(5)` yields

```text
G_5^(Q+tau)((5))=-1440 sum_(r in tau) r^5<0.            (10)
```

Notice that no vanishing of `h_5[-tau]` is used.  Indeed
`h_5[-tau]=-e_5[tau]` can be nonzero when `tau` has at least five entries;
its coefficient in `(3)` vanishes because of `(7)`.  This is the precise
reason the singleton proof extends to the whole upper filter.

## 4. Physical filter and selector consequence

The legal completions are exactly the `64` submultisets `tau<=R`, including
the empty one.  Their nonempty depth census is

```text
|tau|:       1  2  3  4  5  6  7  8
#states:     4  8 12 14 12  8  4  1.                   (11)
```

At `tau=empty`, THM-3209 gives the all-degree zero response.  Every other
state is strictly negative on the same lawful degree-five upset by `(10)`.
The sharp smallest magnitude is `1440`, attained only by `tau={1}`; the full
completion has value

```text
-1440(3*1^5+3*2^5+4^5+5^5)=-6117120.                  (12)
```

Let `Delta_filter` be the probability simplex on these 64 states.  For every
selector cone imposing degrees through any `D>=5`, lawful upset inequalities
require the expectation of `(10)` to be nonnegative, whereas it is a convex
combination of zero and strictly negative numbers.  Therefore

```text
C_D^(<=16) intersect Delta_filter={delta_Q}       (D>=5). (13)
```

This closes every upper-filter direction from physical depths nine through
sixteen at once and supplies a one-sided local collar above the reset.  It
does not use the separately reserved global depth-nine certificate.

## 5. Scope and boundary

The source/target contract is

```text
source:      physical prefixes of the fixed support-(1,3), bank-I2 model;
map:         common-prefix lambda-ring subtraction;
preserved:   first live complete row and zero power-sum marginal;
destroyed:   information outside the principal filter above Q;
sidecar:     other upset rows for states not comparable with Q.             (14)
```

Equation `(10)` applies algebraically to any positive added multiset, but
only `tau<=R` is a legal state in this physical bank.  It does not classify
unrelated depth-ten-or-deeper states, prove strict selector feasibility,
construct a stopping rule, or imply `NC(2)`, `GMC(2)`, or `LRC(14)`.

## 6. Exact evidence

Run

```text
python 04-computation/gmc_complete_reset_upper_filter_thm3219.py
python -O 04-computation/gmc_complete_reset_upper_filter_thm3219.py
```

and compare LF-normalized bytes with the declared output.  The companion uses
only integer and `Fraction` arithmetic.  It reconstructs `(4)` directly,
enumerates all 63 nonempty physical completions, verifies `(6)--(12)`, and
explicitly includes all 25 cases with nonzero `h_5[-tau]`.  The all-prefix
quantifier in `(6),(7),(10)` comes from the displayed lambda-ring identities;
the finite run verifies the physical application.

An independent hostile audit checked the two distinct invariance mechanisms:
`Phi^sigma(h_5)=Phi(h_5)` uses the five vanished complete rows, whereas
`Phi^sigma(p_5)=Phi(p_5)` uses `Phi(1)=0`.  It then rederived all 64 physical
states, the depth census, the 25 nonzero-`h_5[-tau]` hostiles, both sharp
constants, and fresh normal/`-O` transcript equality.

QED.
