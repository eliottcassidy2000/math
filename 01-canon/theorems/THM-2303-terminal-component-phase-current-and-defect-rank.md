---
id: THM-2303
title: "Terminal-component phase current and its exact defect rank"
status: >
  PROVED + VERIFIED-EXACT. A rooted blocker handoff has an intrinsic directed
  component-current multigraph, and each component carries an exact U(1)
  Fourier current whose sum is the desired coefficient. Magnitudes and a
  spanning tree of relative phase transports determine vanishing; if the
  phase-transport graph has k nonzero blocks, exactly k-1 relative phases
  remain in the generic fiber. For two interval components, cancellation is
  equivalent to equal sinc amplitudes and antipodal relative base phase.
  On the multiplier-four THM-2299 carrier the rooted-energy/current-service
  quotient has exact linear continuation-defect rank one. An exact strict
  perturbation keeps the owner, target, clock, source and terminal root
  addresses, component masses, all rooted energies, pair relation, and
  anchored minor fixed while changing F-hat(4), E-hat(52), and W_4-hat(0)
  from zero to nonzero. The carrier is local, not a global scalar cover; no
  LRC(14) profile is excluded.
source: codex-2026-07-25-terminal-component-phase-current
depends_on:
  - THM-840-hamming-five-continuation-congruence-boundary
  - THM-2299-rooted-current-service-energy-and-base-phase-no-go
related:
  - THM-2075-safe-child-homeomorphism-and-wall-word-conjugacy
  - THM-2077-terminal-kakeya-needle-and-recursive-quarter-escape
  - THM-2276-shallow-owner-residue-aligned-crossing
  - THM-2296-prescribed-expiration-return-or-bounded-ancestry-resonance
script: 04-computation/lrc14_terminal_component_phase_current_thm2303.py
output: 05-knowledge/results/lrc14_terminal_component_phase_current_thm2303.out
script_sha256: 5d9671672d3af8dc8c269fe1b5df18e2276b7cda4f5d13fa6871f70269c52615
output_sha256: 3ceb60ce53bb53223787a4f5ff8579f71536e897cf23a3cb788a841a47de3d7b
hash_basis: working-tree bytes (LF)
---

# THM-2303 -- terminal-component phase current

**PROVED + VERIFIED-EXACT.**

THM-2299 shows that current service can fire every rooted character
pointwise while the prescribed pair coefficient cancels between two terminal
components. The missing object can be made exact. It is not a tournament of
components and it is not another nonnegative energy. It is a complex current
on a directed handoff multigraph, together with relative phase transport
between its parallel component edges.

## 1. The intrinsic component-current carrier

Let

```text
T(x)=13x mod 1.                                      (1)
```

Fix a source blocker label `j`, a terminal blocker label `t`, a time `k`,
and a measurable handoff set. Split it into finitely many pieces `E_e` so
that:

1. `T^k` is one-to-one on each `E_e`;
2. every point of `E_e` has source owner `j`;
3. every point of `J_e:=T^k(E_e)` has current target owner `t`; and
4. `E_e` and `J_e` lie in compatible proper circle intervals.

The labels `j,t` are vertices and each component `e` is a time-oriented edge

```text
j --(k, address a_e, component e)--> t.              (2)
```

Choose compatible real lifts. There is one address

```text
0<=a_e<13^k
```

such that the inverse branch on `J_e` is

```text
x=(y+a_e)/13^k.                                     (3)
```

For every integer frequency `n`, define the edge current

```text
I_e(n)=integral_(E_e) exp(-2*pi*i*n*x) dx.           (4)
```

Branchwise change of variables gives the exact phase-transport law

```text
I_e(n)
 =13^(-k) exp(-2*pi*i*n*a_e/13^k)
   integral_(J_e) exp(-2*pi*i*(n/13^k)*y) dy.        (5)
```

The separate factors in (5) depend on the compatible lift gauge, but their
product does not. At a divisible frequency `n=13^k h`, the address factor is
one and (5) is ordinary Fourier descent. At a nonzero rooted character, the
address factor is the discrete root character and the remaining fractional
exponential is the continuous base gauge. Thus (5) types the information
split in THM-2299:

```text
root address  x  continuous terminal base phase
                  = gauge-faithful edge phase.       (6)
```

If `E` is the disjoint union of the handoff components, then

```text
(1_E)_hat(n)=sum_e I_e(n).                           (7)
```

Rooted energy controls sizes of the currents, while (7) needs their complex
sum. Passing from the currents `I_e` to their magnitudes is therefore the
first information-losing map.

This is an intrinsic directed multigraph because time and the named
source/target owners orient every handoff. There is no intrinsic binary
comparison between two component edges. Exact target ties remain colored
ties or are made disjoint by the explicit priority partition of THM-2299.
A tournament would add a cosmetic orientation and forget the additive
current (7).

## 2. Phase transport is a spanning-tree problem

Fix one frequency and abbreviate the nonzero component currents by

```text
I_e=r_e xi_e,        r_e>0,   xi_e in U(1).          (8)
```

A **phase-transport edge** between components `e,f` is the exact ratio

```text
tau_(e,f)=xi_f/xi_e in U(1).                         (9)
```

It is not a service edge. It is a reversible edge in the component phase
groupoid, with

```text
tau_(f,e)=tau_(e,f)^(-1),   and
tau_(e,f)tau_(f,g)=tau_(e,g).                        (10)
```

> **Phase-tree lemma.** Component magnitudes and phase transports on a
> spanning tree determine `sum_e I_e` up to multiplication by one common
> element of `U(1)`. In particular, they determine whether the sum vanishes.
> If the available phase-transport graph has `q` connected components whose
> internal resultant currents are nonzero, then after quotienting the
> irrelevant common phase, `q-1` independent relative `U(1)` phases remain.

### Proof

Choose a root component. Along the unique path in a spanning tree, multiply
the ratios (9) to recover every `xi_e` relative to the root phase. Hence

```text
sum_e I_e=xi_root sum_e r_e (xi_e/xi_root),          (11)
```

and the factor `xi_root` cannot affect vanishing.

For a disconnected phase graph, perform this reconstruction inside each
connected component. Each nonzero block resultant is then known only up to
one block phase. One common simultaneous rotation is irrelevant, leaving
`q-1` independent relative phases. Conversely, specifying those `q-1`
ratios connects the blocks by a tree and reduces to (11). QED.

Thus the exact positive target after THM-2299 is concrete:

```text
construct enough lawful relative-phase identities to connect the
terminal component currents at the prescribed pair frequency.           (12)
```

Energy on every root character does not supply even one edge of this phase
tree.

## 3. Exact two-component cancellation

The smallest obstruction is already complete in two interval components.
Let

```text
I_s=(c_s-epsilon_s,c_s+epsilon_s),   s=1,2,

0<epsilon_s<1/(2m),                  m>=1,           (13)
```

using compatible proper real lifts. Direct integration gives

```text
(1_(I_1 union I_2))_hat(m)
 =A_1 exp(-2*pi*i*m*c_1)
  +A_2 exp(-2*pi*i*m*c_2),                           (14)

A_s=sin(2*pi*m*epsilon_s)/(pi*m)>0.                  (15)
```

Two nonzero complex vectors sum to zero exactly when they have equal
magnitudes and opposite directions. Therefore (14) vanishes if and only if

```text
A_1=A_2,
m(c_2-c_1)=1/2 mod 1.                               (16)
```

If both half-widths are below `1/(4m)`, sine is strictly increasing and the
first condition in (16) is simply

```text
epsilon_1=epsilon_2.                                (17)
```

Common translation rotates (14) and cannot affect its vanishing. The
gauge-invariant missing coordinate is exactly the relative base phase

```text
Phi_m=m(c_2-c_1) mod 1.                              (18)
```

Consequently, equal component widths together with any proof that
`Phi_m!=1/2` force the desired coefficient to be nonzero. This is the sharp
two-component positive lemma, not merely a sufficient semicircle estimate.

For more than two components, the same statement is the weighted phase
polygon

```text
sum_e I_e=0.                                         (19)
```

The phase-tree lemma says which coordinates reconstruct that polygon.

## 4. Exact defect rank of rooted energy

There are two distinct debts: component weights and component positions.
First freeze the antipodal phase directions in the two-component carrier
and allow its two component masses to vary as `(w_-,w_+)`.

On a one-sheet rooted handoff, every nonzero root character has the same
pointwise magnitude. After integration over the same named current-service
label, all twelve rooted-energy observations are scalar multiples of

```text
R(w_-,w_+)=w_-+w_+.                                  (20)
```

At the multiplier-four antipodal phases, after removing the common positive
sinc factor, the target current is

```text
N(w_-,w_+)=i(w_--w_+).                               (21)
```

Over the real field,

```text
ker R=span((1,-1)),
N(ker R)=i R,
dim_R N(ker R)=1.                                    (22)
```

THM-840's exact linear continuation-sidecar theorem therefore gives

```text
d_R(N)=1.                                            (23)
```

Equivalently, the single signed mass current

```text
S(w_-,w_+)=w_--w_+                                  (24)
```

is a minimum-rank repair when the phase directions are fixed. The companion
independently checks

```text
rank R=1,       rank(R,Re N,Im N)=2.                 (25)
```

This does not pay the positional debt. If the two component masses are
already retained, varying their relative placement leaves (20) unchanged
and changes (18). The next section proves that this second loss survives
every discrete label currently available.

## 5. Same labels and energies, opposite Fourier verdicts

Use the exact multiplier-four row from THM-2299:

```text
H=1,
(q_1,q_2,q_3,q_4,q_5)=(4,2,3,6,10),
(c_1,c_2,c_3)=(13,13^3,2*13^5).                    (26)
```

It has the primitive pair and independent relation

```text
13q_1-4c_1=0,
q_1-2q_2=0,                                         (27)
```

and their `(c_1,q_2)` anchored minor is

```text
8!=0 mod 13.                                        (28)
```

Put

```text
epsilon=10^(-12),       eta=10^(-12),

F_0=
 (-1/16-epsilon,-1/16+epsilon)
 union
 ( 1/16-epsilon, 1/16+epsilon),

F_1=
 (-1/16-epsilon,-1/16+epsilon)
 union
 ( 1/16+eta-epsilon,1/16+eta+epsilon),              (29)

E_s={(8+y)/13:y in F_s},        s=0,1.              (30)
```

The exact smallest source/terminal label margin in the `y` coordinate is

```text
3/540602608 > epsilon+eta.                           (31)
```

Thus both complete open carriers, not only their centers, have the same
strict labels:

```text
source:   c_1-private;
time 2:   c_2-only current service.                  (32)
```

Both have source root address `8`; their two terminal inverse-root addresses
are `(12,0)`. Their component half-widths and masses agree. On every terminal
root fiber the ancestry word has one sheet, so all twelve rooted characters
have the same pointwise magnitude and the same integrated energy in the two
carriers. The owner, target, clock, complete root-address multiset, component
masses, pair relation, second relation, and anchored minor are therefore
identical.

Their relative base phases are different:

```text
Phi_4(F_0)=1/2,
Phi_4(F_1)=1/2+4*10^(-12).                           (33)
```

The two-component criterion (16) gives the opposite exact verdicts

```text
(1_(F_0))_hat(4)=0,
(1_(F_1))_hat(4)!=0.                                (34)
```

The one-sheet lift identities give

```text
(1_(E_s))_hat(52)=(1/13)(1_(F_s))_hat(4).            (35)
```

Also `G_s=P1_(E_s)=(1/13)1_(F_s)`, and the character-four gauge of
THM-2299 gives

```text
(W_(s,4))_hat(0)=(1_(F_s))_hat(4).                  (36)
```

Hence every coefficient in (34)--(36) vanishes for `s=0` and is nonzero for
`s=1`.

This proves necessity, not only sufficiency, of an affine component-phase
sidecar. No combination of the fixed discrete data in (32), root energies,
component masses, or the anchored relation minor determines the prescribed
pair coefficient. A phase-tree edge between the two components would.

## 6. Connection and stopping boundary

The exact connection contract is

```text
source:
  THM-2299 rooted current-service word;

target:
  a complex current on the directed blocker-handoff component multigraph;

map:
  split by injective branches and named service, then integrate the
  gauge-faithful phase on every component;

preserved:
  owner, target, clock, branch address, component, and exact Fourier sum;

lost by energy:
  signed mass imbalance and relative affine base phase;

minimum fixed-phase linear sidecar:
  rank one on the canonical two-component carrier;

positional repair:
  one relative U(1) phase for two components, or a phase-transport spanning
  forest in general;

cheapest decisive next test:
  seek a global cover identity that supplies one lawful phase-tree edge
  between distinct prescribed-service components.                    (37)
```

The theorem does not claim that every genuine counterexample has only two
components or that the phase-transport graph is connected. It identifies
the exact additional coordinate and its reconstruction law. The hostile
carriers are local and do not satisfy a global scalar cover. No scalar
profile is excluded, and LRC(14) remains open.

## 7. Exact verification

The companion uses integer and `Fraction` arithmetic. It checks the scalar
row, both relations, anchored minor, exact margin (31), every source and
terminal label at every interval endpoint, source and terminal root
addresses, equality of the discrete/rooted-energy signatures, both rational
relative phases, and the exact ranks (25). The Fourier zero/nonzero verdicts
then follow from the proved formula (14)--(16), not from floating-point
sampling.

Reproduce with

```bash
python3 04-computation/lrc14_terminal_component_phase_current_thm2303.py
python3 -O 04-computation/lrc14_terminal_component_phase_current_thm2303.py
```

Normal and optimized transcripts are byte-identical to the stored output.
QED.
