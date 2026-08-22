---
id: THM-3222
title: "Universal product-Gamma reset upper-filter collar"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On every
  one of the 115 maintained anchored supports and both THM-3110 product-Gamma
  banks, the distinguished Ferrers alphabet is a physical submultiset of the
  reduced pole alphabet.  Its first live complete-row coefficient is positive,
  and every nonempty physical completion has strictly negative degree-five
  principal-upset response.  Thus each of the 230 reset principal filters
  meets its selector cone only at the reset atom.  States outside those
  principal filters and arbitrary supports remain open.
source: root/2026-08-02
audit: >
  The assertion-independent exact companion pins THM-3110 and THM-3120 and
  their generating scripts, reconstructs all 230 reduced fractions, proves
  Q<=P in every case, checks the vanished complete and power rows, and runs
  5,442 direct virtual-alphabet response controls including every full collar.
  An independent hostile audit rederived the algebra and scope; fresh normal
  and optimized replays both byte-match the stored output and declared hashes.
depends_on:
  - THM-3110-arbitrary-anchored-product-gamma-dominant-tail-and-low-histogram-reduction
  - THM-3120-row-pole-prefix-newton-flag-positivity
related:
  - THM-3216-depth-nine-degree-fourteen-unique-reset-face-and-omega-cone-boundary
  - THM-3219-complete-reset-upper-filter-principal-upset-exclusion
script: 04-computation/gmc_universal_reset_upper_filter_thm3222.py
output: 05-knowledge/results/gmc_universal_reset_upper_filter_thm3222.out
script_sha256: 57265de7f1b163df60bb2c710d10d6ced1867de7542c300a1dd245ec32c4f4aa
output_sha256: 62ac94f3d495655188c0a5b4f3b619c5046a967b47f585f81d8f51e711c315b9
hash_basis: LF-normalized bytes
---

# THM-3222 -- universal product-Gamma reset upper-filter collar

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3219 isolates one complete reset in its whole physical upper filter.
The same mechanism is present on every maintained product-Gamma support and
in both signed response banks.  The nontrivial new fact is not merely that
the first live coefficient is positive--THM-3120 already proves that--but
that the distinguished quotient alphabet survives every denominator
cancellation and remains a physically available pole submultiset in all 230
cases.

## 1. The maintained finite universe

Let

```text
U={(a,b):1<=a<=10, a<b<=min(3a+4,21)}.                    (1)
```

Thus `|U|=115`.  Fix `(a,b) in U` and one THM-3110 bank `I_i`, `i=1,2`.
After the exact common-root deletion, write its signed residual-alphabet
functional as

```text
Phi_i(f)=sum_R epsilon_R f[S_R].                          (2)
```

THM-3120 reduces the complete-row generating function to

```text
sum_(n>=0) Phi_i(h_n)t^n=N_i(t)/D_i(t),
D_i(0)=1,
N_i(t)=t^5 P_i^num(t).                                   (3)
```

To avoid a clash between the numerator polynomial and the physical pole
alphabet, call the latter `P`.  It is the multiset of positive roots in the
reduced split denominator `D_i`.

THM-3110's unique positive Ferrers-dominant row is

```text
Q_1=(3b,2a,2a,a),
Q_2=(3b,a+b,2a,a,a).                                     (4)
```

Let `Q` be its residual root alphabet after the same common-root deletion.
Exact integer arithmetic in every case proves

```text
Q<=P.                                                     (5)
```

This containment is load bearing: it says the complete quotient reset is a
legal physical prefix, not merely a virtual alphabet.  Put

```text
R=P-Q.                                                    (6)
```

The complete census is

```text
cases                         230
|P|                           11..133
|Q|                            6..62
|R|                            5..71
number of submultisets of R   16..187406683791040512.     (7)
```

In particular every collar is nonempty.  The assertion is finite-exact on
`U`; no containment theorem for arbitrary supports is being inferred.

## 2. The invariant first live coefficient

Let

```text
c_(a,b,i)=Phi_i(h_5).                                     (8)
```

Equation `(3)` and THM-3120's positive pole flag give

```text
Phi_i(h_j)=0                    (0<=j<5),
c_(a,b,i)>0.                                             (9)
```

Across the 230 cases there are 225 distinct values, with sharp range

```text
32<=c_(a,b,i)<=410879700000.                             (10)
```

The lower endpoint is bank `I_2` at `(a,b)=(1,2)`; the upper endpoint is
bank `I_1` at `(10,21)`.

The signed row marginal of each THM-3110 bank vanishes.  Additivity of power
sums therefore gives

```text
Phi_i(1)=0,                    Phi_i(p_j)=0  (1<=j<=5).   (11)
```

For a finite common prefix `sigma`, lambda-ring subtraction yields

```text
h_5[A-sigma]=sum_(j=0)^5 h_(5-j)[A]h_j[-sigma],
p_5[A-sigma]=p_5[A]-p_5[sigma].                          (12)
```

Applying `(9)--(11)` to `(12)` proves the two distinct invariances

```text
Phi_i^sigma(h_5)=c_(a,b,i),
Phi_i^sigma(p_5)=0.                                      (13)
```

The first uses all five vanished complete rows; the second uses both the
power marginal and zero signed mass.  Neither is an asymptotic statement.

## 3. Universal upper-filter exclusion

For a physical prefix `sigma<=P`, write `Q^sigma=Q-sigma`.  The inherited
partition response at degree `N` is

```text
G_N^sigma(lambda)
 =Phi_i^sigma(h_N)m_lambda[Q^sigma]
  -Phi_i^sigma(m_lambda)h_N[Q^sigma].                    (14)
```

At degree five the principal coarsening upset is the singleton `{(5)}`, and
`m_(5)=p_5`.  Let `tau` be any nonempty submultiset of `R` and set

```text
sigma=Q+tau.                                              (15)
```

Then `sigma<=P`, while `Q^sigma=-tau`.  Hence

```text
p_5[-tau]=-sum_(r in tau)r^5.                            (16)
```

Substitution of `(13),(16)` into `(14)` gives the universal exact formula

```text
G_5^(Q+tau)((5))
 =-c_(a,b,i) sum_(r in tau) r^5 <0.                      (17)
```

No vanishing of `h_5[-tau]` is used.  In fact, because `|R|>=5`, every full
completion has

```text
h_5[-R]=-e_5[R]<0,                                       (18)
```

so all 230 maximal completions are explicit hostiles to that shortcut.
The globally sharp nonempty magnitude in `(17)` is `32`, attained by adding
one pole `1` to bank `I_2` at `(1,2)`.

At `tau=empty`, `Q^Q=0`, so every positive-degree partition response is zero.
Let `Delta_(a,b,i)` be the probability simplex on the submultiset states
`{Q+tau:tau<=R}`.  Every selector cone containing the lawful degree-five
upset inequality requires the expectation of `(17)` to be nonnegative.
It follows that, for every degree cutoff `D>=5`,

```text
C_D intersect Delta_(a,b,i)={delta_Q}.                   (19)
```

Thus the 230 reset atoms each have a strict one-sided collar throughout
their whole physical principal upper filter.

## 4. Relation to the global depth-nine dual

THM-3219 is the `(a,b,i)=(1,3,2)` specialization of `(17)`; its 64-state
census has `c=1440`.  THM-3216 is complementary rather than duplicated.
On that same face its first Farkas row is the degree-five local collar, but
the other seventeen rows stitch roughly 3,124 states lying outside the
principal filter.  The present theorem transports the local seed to 230
support-bank universes; it does **not** transport those other seventeen rows
or their positive dual weights.

The connection contract is therefore

```text
source:       230 reduced product-Gamma pole banks;
map:          distinguished reset plus common-prefix lambda subtraction;
preserved:    physical containment, first live coefficient, power marginal;
destroyed:    every state not comparable with the reset Q;
needed:       exterior upset rows and positive weights for global stitching. (20)
```

No new `NC(2)`, `GMC(2)`, arbitrary-support, stopping-rule, or `LRC(14)`
claim follows.  In particular, a union of local stars is not automatically
a global selector certificate.

## 5. Exact evidence

Run

```text
python 04-computation/gmc_universal_reset_upper_filter_thm3222.py
python -O 04-computation/gmc_universal_reset_upper_filter_thm3222.py
```

and compare LF-normalized bytes with the stored output.  The companion pins
both proved dependencies and both generating scripts.  It reconstructs all
230 reduced fractions and distinguished alphabets, checks `(5),(7),(9)--(11)`,
and performs 5,442 direct prefixed virtual-alphabet evaluations spanning every
case, every full collar, singleton extremes, and multi-letter controls.  It
also records a SHA-256 digest of every exact `P,Q,R,c` case record.  The
all-`tau` quantifier in `(17)` is the displayed algebraic identity, not an
attempt to enumerate up to `187406683791040512` states.

An independent hostile audit rechecked physical containment after the same
reduced-row cancellation, the two distinct common-prefix invariances, the
all-`tau` sign formula, the simplex implication, the full negative-alphabet
hostile, and the local-versus-global boundary.  Fresh normal and optimized
runs took independent paths through the assertion-free companion, matched
each other and the stored transcript byte-for-byte, and reproduced both
declared LF hashes.

QED.
