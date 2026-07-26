---
id: THM-2396
title: "Common-core forty-nine-orbit word incompatibility"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. Assume
  THM-2393's sole no-clean residual M=1 and
  (C_1,C_2,c_1,c_2)=h(1,13,13,169). Bundle the seven
  THM-2394 root fibres into one generic Z/49Z orbit on which C_3 and
  c_3 are safe. THM-2391 controls the unique q_*-danger bin; THM-2394
  controls the other six. After an affine gauge fixes D_h, every
  physical orbit injects into a larger independent-translate word
  universe. Exact exhaustion finds 2,058 top-compatible four-q
  packings and rejects all of them: 14,386 B/packing cases have a
  forbidden nonzero hole, while the last 20 require D_(169h) to agree
  with D_h on four or five bins although any translate agrees on at
  most two. Thus the zero-clean common-core branch is empty and every
  remaining packet has positive clean-hole mass. No uniform mass floor,
  row exclusion, ledger decrement, target landing, or proof of LRC(14)
  is claimed.
source: codex-2026-07-26-common-core-forty-nine-orbit
depends_on:
  - THM-2391-blocker-caged-septimal-single-layer-address-reduction
  - THM-2393-c3-safe-double-fibre-capacity-and-common-core-residual
  - THM-2394-common-core-six-address-transversal-and-labelled-hole-trichotomy
related:
  - THM-2377-septimal-valuation-collision-and-bockstein-carry-gate
  - THM-2390-septimal-layer-kraft-peeling-and-heavy-word-reduction
  - THM-2395-common-core-successor-shell-and-forced-hole-escape-tax
script: 04-computation/lrc14_common_core_forty_nine_orbit_thm2396.py
output: 05-knowledge/results/lrc14_common_core_forty_nine_orbit_thm2396.out
script_sha256: 8a73dc2f5def873b8042149fb84d88dae8c6d0ee5e62e38be5c9f6b95f0d5112
output_sha256: 847d1f1c2132428bd24ce90efdc24f4ead5c55d9f2a86e94da3e4823ca68f08f
hash_basis: working-tree bytes (LF)
---

# THM-2396 -- common-core forty-nine-orbit word incompatibility

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2393 reduces every packet with no clean hole to one literal
common-core chain. THM-2394 then turns each high-safe seven-root fibre
into a six-address transversal. The missing move is not another
one-fibre estimate: the unique top `q_*` word couples seven such fibres
inside one orbit of length forty-nine.

Keeping that whole orbit gives a finite contradiction. The proof
deliberately enlarges the physical universe by letting all non-normalized
word translates vary independently. The enlarged universe is still empty.
Consequently the shared physical phase cannot be the reason for the
contradiction.

## 1. The generic forty-nine-root orbit

Retain THM-2393/2394's zero-clean boundary

```text
M=1,

(C_1,C_2,c_1,c_2)=h(1,13,13,169),

gcd(h,91)=1,                    delta=0.              (1)
```

Write

```text
q_*=7u,                 C_3=49V,                 c_3=13C_3.
                                                               (2)
```

The valuations in the last septimal lane justify all three integer
factorizations. The two high blockers have exact joint-safe base mass

```text
mu(D_V^c intersection D_(13V)^c)
 =1-2/7+mu(D_V intersection D_(13V))
 =1-2/7+1/91
 =66/91>0.                                             (3)
```

Choose a generic base `z` in this set, outside the finite endpoint sets
and the finite pullbacks of every almost-everywhere exceptional set.
On the orbit

```text
O_z={x_j=(z+j)/49:j in Z/49Z},                        (4)
```

both `C_3` and `c_3` are safe at every point. The top word `D_(q_*)`
depends only on `j mod 7` and, away from aligned endpoints, occupies
exactly one complete bin

```text
Q_(r_0)={j:j=r_0 mod 7}.                              (5)
```

For the other six bins, put

```text
y_r=(z+r)/7.
```

They lie in THM-2394's corrected divided-base set `Y_0`: `q_*/7=u` is
safe there, while `C_3/7` and `c_3/7` are safe because (3) is constant
on the seven base roots. Therefore:

- on the top bin (5), THM-2391 says that `D_h,D_(13h)` partition the
  two guard addresses and all four lower `q` addresses lie inside that
  guard pair;
- on every other bin, THM-2394 says that the guard and four lower `q`
  addresses are six distinct points, with unique hole

  ```text
  d=B,              or              d=A=C,            (6)
  ```

  where `A,B,C` are the addresses of
  `D_h,D_(13h),D_(169h)`.

The word constraints in (5)--(6) hold simultaneously on the same
forty-nine physical points. This simultaneity is the information lost by
one-fibre capacity and successor calculations.

## 2. A relaxed exact word universe

Write every orbit index uniquely as

```text
j=r+7s,                    r,s in F_7.                (7)
```

For a seven-unit `v mod 49` and a translate `t mod 49`, define the
ordinary and guard words

```text
O(v,t)
 ={j:vj-t in {0,...,6} mod 49},

G(v,t)
 ={j:vj-t in {0,...,13} mod 49}.                     (8)
```

Every ordinary word has one address `s` above each `r`; every guard
word has two. There are exactly

```text
1029 distinct ordinary words,
1029 distinct guard words.                           (9)
```

These banks contain every generic physical word in (4). Indeed, after
multiplying the strict danger inequality by forty-nine, a unit danger
arc selects seven consecutive integer residues and the guard selects
fourteen. The aligned six/thirteen-point cases were removed when `z`
was chosen generic.

Now let the physical `D_h` word be

```text
{j:hj-t_A in {0,...,6}}.
```

The affine permutation

```text
j -> hj-t_A                         mod 49            (10)
```

normalizes it to

```text
A=O(1,0),                 A(r)=0 for every r.         (11)
```

It sends the seven residue bins to residue bins, so the top bin remains
some `r_0`. The two correlated core words become

```text
B=O(13,t_B),

C=O(169,t_C)=O(22,t_C)              mod 49.          (12)
```

Every lower `q` word remains an arbitrary member of the ordinary bank
and the guard remains an arbitrary member of the guard bank.

The exact search makes one further relaxation: the starts `t_B,t_C`,
the four lower-`q` starts, and the guard start are allowed to vary
**independently**. In a physical orbit they arise from one base `z` and
are correlated. Thus every physical orbit maps into the searched
universe, but the converse is neither used nor asserted.

## 3. The finite incompatibility

The exact universe is filtered in the following order.

1. Choose the top residue `r_0` and a guard whose two top addresses
   include `A(r_0)=0`.
2. Choose a speed-thirteen translate `B` whose top address is the other
   guard address.
3. Choose four ordinary words. On the top row each must use one of the
   two guard addresses. On the other six rows their addresses must be
   pairwise distinct and avoid both guard addresses. Their unique
   missing address is the `K`-hole.
4. On each safe row require the labelled no-clean alternative (6):

   ```text
   hole=B(r),          or
   hole=0=C(r).                                      (13)
   ```

The complete exact ledger is

```text
relevant (top row,guard) pairs:              2058;

admissible four-q packings:                  2058;

packing/B tests:                            14406;

tests with a nonzero hole different from B: 14386;

remaining tests:                               20.   (14)
```

The last twenty tests are not numerical near misses. In eighteen of
them, (13) asks `C` to agree with `A=0` on four safe rows; in the other
two it asks for five. But the complete bank of forty-nine
speed-twenty-two translates has only thirteen distinct zero masks and

```text
#{r:C(r)=A(r)=0}<=2.                                (15)
```

Hence every last test fails. The enlarged universe has exactly zero
survivors.

For audit resolution, the intermediate distributions are

```text
number of lower-q candidates:
  22:588, 25:588, 26:294, 33:588;

number of four-q packings per relevant guard/top pair:
  0:1176, 1:294, 3:588;

required A=C rows in the final twenty tests:
  4:18, 5:2.                                         (16)
```

The `2058` intermediate packings and twenty final capacity boundaries
are positive controls: the contradiction does not come from separately
making the top-bin partition or the six safe-bin transversals
impossible. It appears only after the labelled common-core word `C` is
required to serve all six bins simultaneously.

An independent hostile audit reconstructed the physical embedding
before replaying the search. In particular, it checked directly that
`C_3=49V,c_3=13C_3` are safe on all forty-nine roots, that `q_*=7u`
selects exactly one residue row, and that the other six divided bases
lie in the corrected `Y_0`. It also verified analytically that

```text
22r mod 49 in {0,22,44,17,39,12,34},                (16a)
```

so no cyclic block of seven residues contains more than two of these
values. This independently proves the sharp cap (15), rather than
inferring it from the zero-survivor count.

## 4. Consequence and exact boundary

Sections 1--3 show that the assumptions in (1) cannot coexist. By
THM-2393, this was the only possible packet with

```text
delta=0.
```

Therefore every remaining scalar packet in the last septimal lane has

```text
delta>0.                                             (17)
```

THM-2392 then supplies a literal `q_*`-labelled charged word and a
same-parent `F_7 x F_13` tensor for each individual packet. The present
finite obstruction supplies no coefficient-independent lower bound for
`delta`: physical chamber widths can shrink with the coefficients, and
the relaxed word proof records only emptiness at equality.

This theorem does **not**:

- align the positive clean word with a terminal endpoint reference;
- force a target-thirteen phase after restoring the later factors;
- exclude a positive-clean scalar row;
- decrement the `165`-row ledger; or
- prove LRC(14).

The next sharp target is no longer a common-core zero-clean automaton. It
is a coefficient-dependent positive-clean landing theorem: preserve the
owner-resolved `F_7 x F_13` word from THM-2392 through the terminal
endpoint filtration without replacing its phase by unlabelled energy.

## 5. Exact companion

The dependency-free companion:

- constructs all `1029+1029` distinct word supports in (9);
- fixes `A`, exhausts all seven top rows, all top-compatible guards,
  all forty-nine `B` and `C` translates, and every four-q packing;
- verifies the full ledgers (14)--(16);
- checks all `980` final `C`-word tests directly;
- independently reconstructs the thirteen zero masks and the sharp
  two-row agreement cap (15); and
- retains the `2058` packing and twenty capacity-boundary positive
  controls.

Run

```bash
python3 04-computation/lrc14_common_core_forty_nine_orbit_thm2396.py
python3 -O 04-computation/lrc14_common_core_forty_nine_orbit_thm2396.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_common_core_forty_nine_orbit_thm2396.out
```

after LF normalization. Every executable assertion remains active under
optimized Python.
