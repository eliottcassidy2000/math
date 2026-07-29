---
id: THM-2884
title: "Macro-semantic diagonal horn carrier and origin-even boundary"
status: >
  PROVED + VERIFIED-EXACT.  On every one of the 449 THM-2835 horn
  cylinders, adjoining the fully bare semantic word Q0 and imposing the
  honest Boolean relation E3_bit XOR outer_slot7 XOR outer_slot8 = 0
  produces a four-state V4-shaped carrier.  At the stepped endpoint
  origin the states are Q0(q3)=000, QA(q11)=110, QAB(q7)=011, and the
  fixed QB source=101.  Exactly 20 native cells retain
  I,J0,J3,J11,J7 and the other five factors, hence 20*449=8980 complete
  stepped packets; both q7 ancestry fillers survive.  The same
  origin-independent diagonal selector gives a nonzero signed q3
  coefficient in both exact fields at all 26 tested unit samples, but
  q11 and q7 are origin-even and cancel.  The V4 edge triangle is flat
  but does not extend to a nontrivial C13 action and is blind to the
  ancestry carry.  This is a physical carrier theorem, not a row
  exclusion or an LRC(14) proof.
source: root/lrc-macro-semantic-diagonal-2026-07-29
depends_on:
  - THM-2835-q11-semantic-word-horn-and-bockstein-blind-support-no-go
  - THM-2847-q3-q11-transverse-endpoint-edge-and-e3-realization-horn
  - THM-2868-two-origin-endpoint-projector-and-projective-kummer-carry-descent
  - THM-2882-event-twisted-all-q-coefficient-carry-lift
related:
  - THM-2788-physical-modular-odometer-versus-heisenberg-bockstein-extension
  - THM-2851-q3-q11-q7-bockstein-holonomy-and-realization-sidecar
  - THM-2878-endpoint-factor-exit-carry-transducer-and-flat-vertex-boundary
  - THM-2880-q0-q3-one-fibre-selector-provenance-obstruction
script: 04-computation/lrc14_macro_semantic_diagonal_horn_carrier_thm2884.py
output: 05-knowledge/results/lrc14_macro_semantic_diagonal_horn_carrier_thm2884.out
script_sha256: b739be20e741d5c061e0febcc8aef9b0f58f4ae8a648aa803610e0dad991929f
output_sha256: 8c3829b1052a641ca08a5e5bda86d9d5e8bd1584f5b2911c57c9fad9da41d4b6
hash_basis: LF-normalized bytes
---

# THM-2884 -- macro-semantic diagonal horn carrier and origin-even boundary

**PROVED + VERIFIED-EXACT.**

THM-2835 found a `449`-sheet semantic horn

```text
QB(source,a) -> QA(q11,a) -> QAB(q7,a+1),
```

but literal transport lost the macro `E3` factor at `q7`.  The missing
view was to stop asking `E3` and the outer semantic word to survive
separately.  They carry complementary parity defects.  Their diagonal
fibre product is an honest physical set.

The result is the first common physical carrier in this lane that contains
the fixed source, the `q3/q11/q7` horn, both `q7` carry fillers, and every
other native factor on a nonempty cell bank.  It also isolates what this
repair does **not** supply: the endpoint-origin character and the
Bockstein carry remain independent coordinates.

## 1. The fourth semantic word

Keep the THM-2835 scale

```text
p=13,                 D=13^5=371293,
T=297836897838480,    U=T/13,
I=[142004992589460,142005019034340),
J_q=I+431933040+qU.                                      (1)
```

Write `u,v in F2` for danger membership in outer slots `7,8`.  The three
previously named words are

```text
QA =(1,0),       QB =(0,1),       QAB=(1,1).              (2)
```

Adjoin the fully bare word

```text
Q0=(0,0),                                                (3)
```

with the same guard-good condition, all five units safe, owner safe, and
both outer slots safe.  This is not a formal complement: the companion
constructs it directly from the exact periodic combs.  It has

```text
132262 intervals,       measure numerator 45425355597700. (4)
```

Let `L` be the THM-2835 set of `449` ancestry labels.  The complete target
word on `L`, with no midpoint-only hits, is

```text
q       0    1 2 3 4 5    6    7     8 9 10 11 12
word    QB   Q0 Q0 Q0 Q0 Q0   QA   QAB   QA ... QA.      (5)
```

In particular,

```text
|Q0(q3)|=66104 globally,
L is contained in Q0(q3),                                (6)
```

and all `449` labels lie in `Q0` at each `q=1,...,5`.
At `q7`, both ancestry `a` and ancestry `a+1` lie in `QAB`, so the
zero-borrow and one-borrow fillers remain indistinguishable.

## 2. The macro-semantic diagonal

For endpoint origin `o in {0,1}`, let

```text
e_o(X)=1  iff the atom X lies in the o-chart E3 block.    (7)
```

The two exact endpoint truth banks are

```text
{q:e_0(J_q)=1}={0,3,11},
{q:e_1(J_q)=1}={0,11},                                   (8)
```

while `e_0(I)=e_1(I)=1`.

Define

```text
chi(e,u,v)=e+u+v in F2,
H=ker chi
 ={000,011,101,110}.                                     (9)
```

Equivalently, the physical diagonal set is

```text
D_o =
  (E3_o intersect {u XOR v=1})
  union
  (E3_o^c intersect {u XOR v=0}).                        (10)
```

Thus `(10)` is one origin-independent Boolean rule: select `E3` when the
outer parity is odd and its complement when the outer parity is even.
An atom lies in `D_o` exactly when its state `(e_o,u,v)` lies in `H`.

The fixed source and the three distinguished targets at the stepped
origin are

```text
source QB : 101,
q3 Q0     : 000,
q11 QA    : 110,
q7 QAB    : 011.                                         (11)
```

They are the four elements of `H`.  The target edge labels are

```text
q3 -> q11 : 110,
q11-> q7  : 101,
q3 -> q7  : 011,                                         (12)
```

and

```text
110 XOR 101 = 011.                                       (13)
```

So the local triangle is flat and the source state is exactly its middle
edge label.  This is the precise `V4`, not an analogy.

At the zero endpoint origin, `q3` instead has state `100`, in the opposite
coset of `H`.  Consequently the all-`q` diagonal supports are

```text
origin 0: {0,1,2,4,5,7,11},
origin 1: {0,1,2,3,4,5,7,11}.                            (14)
```

Their symmetric difference is exactly `{3}`.

## 3. Twenty complete physical cells

The six native factors are ordered

```text
(E3,clock,q1,q2,c2,c3).                                  (15)
```

Before replacing `E3` by the diagonal, the exact cell counts are

```text
full q3 cells                         56,
full q11 cells                        63,
full q3/q11 intersection              42,
q7 cells failing only E3              49,
q11 and q7-E3-only intersection       35.                (16)
```

Intersecting all three requirements leaves exactly

```text
clock=1,
s in {0,3,8,9,12},
t in {5,6,9,10},                                        (17)
```

hence `5*4=20` cells.  On each cell in `(17)` the five factors

```text
(clock,q1,q2,c2,c3)                                      (18)
```

contain all of

```text
I, J0, J3, J11, J7.                                      (19)
```

At the stepped endpoint origin, the diagonal `H` contains their
corresponding word/macro states:

```text
I and J0 : QB=101,
J3       : Q0=000,
J11      : QA=110,
J7       : QAB=011.                                      (20)
```

Therefore `(17)` and the `449` ancestry labels give

```text
20*449=8980                                               (21)
```

complete stepped-origin packets.  Since both `q7` ancestry fillers lie
in `QAB`, there are `17960` labelled filler occurrences if the two carry
choices are counted separately.

This repairs the `0/567` common-cell failure in THM-2835 only after
replacing the product demand “`E3` and word” by the correlated diagonal
`(10)`.  It does not assert that the original `E3` product cell was
nonempty.

## 4. Exact signed endpoint coefficient

For each target word, put

```text
p(q)=u(q) XOR v(q).                                      (22)
```

At both endpoint origins use the same selector:

```text
p(q)=1 -> choose E3,
p(q)=0 -> choose E3 complement.                          (23)
```

The selected endpoint atom is full exactly when `(e_o,u,v) in H`.
The companion evaluates `(23)` for the `26` inherited `91`-unit samples
in both exact cyclotomic finite-field embeddings

```text
352341050142921841,
956354278959359281.                                      (24)
```

For every sample and both fields,

```text
support(value at origin0 - value at origin1)={q3}.        (25)
```

At `q3`, the zero-origin complement is empty and the stepped-origin
complement is the whole atom, so

```text
signed q3 value = - base endpoint value !=0.              (26)
```

Multiplying by the inherited common source coefficient remains nonzero.
The `26` resulting `q3` currents are pairwise distinct in each field.

At `q11` and `q7`, both endpoint origins lie in `H`; their endpoint
values are equal and cancel.  Thus `(10)` supplies a canonical
origin-odd `q3` section but not a signed `q11` or `q7` coefficient.

## 5. Relation to the desired selector parity

On the four seam checkpoints

```text
(q0,q3,q11,q7),                                          (27)
```

the nonlinear indicator of a nonzero `H` state is

```text
1_{h !=0}=(1,0,1,1).                                     (28)
```

This equals the affine selector parity `1+delta_3` isolated in
THM-2882.  It gives a precise geometric meaning to that parity:
the three active words are the three nonzero points

```text
H\{0} = P^1(F2).                                         (29)
```

But `(29)` is a three-element set.  It is neither a subgroup, a coset,
nor a character fibre of `V4`; it is a nonlinear support condition.
Moreover it is only seam-local.  The exact all-`q` counterexamples are

```text
q in {1,2,4,5,6,8,9,10,12}.                              (30)
```

Most importantly, at `q11` and `q7` the two endpoint origins have the
same `H` state.  No origin-independent function of that state can output
different selectors at the two origins.  Equation `(28)` identifies the
needed selector **parity**, but realizing it still requires the separate
endpoint-origin character.

## 6. Sharp action and carry boundaries

The flat edge law `(13)` is local.  It is not a `C13` action on `H`.
Indeed a generator of such an action would be a permutation

```text
pi in S4,       pi^13=1.                                 (31)
```

Exact enumeration of all `24` permutations leaves only the identity.
Equivalently, neither `V4` nor `Aut(V4)=S3` contains an element of order
`13`.

Nor does `H` detect the Bockstein.  The two physical fillers

```text
(q7,a),       (q7,a+1)                                   (32)
```

both have state `011`.  The diagonal therefore preserves precisely the
support blindness already isolated by THM-2835.

The remaining invoice is now separated into two independent coordinates:

1. an endpoint-origin character or clutch that transports the unique
   stepped `q3` choice to `q11/q7` without reintroducing the cancelling
   origin;
2. an ancestry-carry sidecar distinguishing the two occurrences in
   `(32)`.

The event-twisted coefficient lift of THM-2882 addresses the second
coordinate algebraically, but does not manufacture the first physical
origin clutch.

## 7. Exact consequence and non-consequence

This theorem proves:

1. the fully bare `Q0` word and its exact all-`q` incidence on the
   `449` horn sheets;
2. an honest macro/word diagonal `H`, with a literal four-state
   `V4` geometry;
3. a `20`-cell, `8980`-packet physical repair of the macro/semantic
   mismatch;
4. a canonical stepped-only, nonzero signed `q3` coefficient;
5. the exact origin-even cancellation at `q11/q7`;
6. the seam-local nonlinear interpretation of the THM-2882 positive
   selector parity;
7. the global `C13`-action and ancestry-carry no-go boundaries.

It does **not** prove:

- a nonzero signed `q11` or `q7` current;
- an endpoint-origin clutch;
- a Bockstein-sensitive physical action;
- a row exclusion;
- the Lonely Runner Conjecture for `14` runners.

## 8. Reproduction

Run

```text
python3 04-computation/lrc14_macro_semantic_diagonal_horn_carrier_thm2884.py
python3 -O 04-computation/lrc14_macro_semantic_diagonal_horn_carrier_thm2884.py
```

Both executions are assertion-free and byte-match

```text
05-knowledge/results/lrc14_macro_semantic_diagonal_horn_carrier_thm2884.out.
```
