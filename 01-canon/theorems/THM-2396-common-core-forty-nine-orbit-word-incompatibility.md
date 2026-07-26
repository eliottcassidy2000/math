---
id: THM-2396
title: "Common-core forty-nine-orbit word incompatibility"
status: >
  PROVED + VERIFIED-EXACT + TWICE INDEPENDENTLY HOSTILE-AUDITED. Assume
  THM-2393's residual speed/valuation hypotheses M=1 and
  (C_1,C_2,c_1,c_2)=h(1,13,13,169). Bundle the seven
  THM-2394 root fibres into one generic Z/49Z orbit on which C_3 and
  c_3 are safe. THM-2391 controls the unique q_*-danger bin; THM-2394
  controls the other six. After an affine gauge fixes D_h, every
  physical orbit injects into a larger independent-translate word
  universe. Exact exhaustion finds 2,058 top-compatible four-q
  packings and rejects all of them: 14,386 B/packing cases have a
  forbidden nonzero hole, while the last 20 require D_(169h) to agree
  with D_h on four or five bins although any translate agrees on at
  most two. Localizing the contrapositive shows that every high-safe
  orbit contains a clean root, so delta>=66/4459. Together with
  THM-2393 this gives the universal last-lane floor delta>=1/26754.
  The common-core charged-cell and transverse-tensor floors are
  33/115934 and 33/753571. No row exclusion, ledger decrement, target
  landing, or proof of LRC(14) is claimed.
source: codex-2026-07-26-common-core-forty-nine-orbit
depends_on:
  - THM-2391-blocker-caged-septimal-single-layer-address-reduction
  - THM-2392-clean-toothpick-or-bounded-cross-ancestor-cage
  - THM-2393-c3-safe-double-fibre-capacity-and-common-core-residual
  - THM-2394-common-core-six-address-transversal-and-labelled-hole-trichotomy
related:
  - THM-2377-septimal-valuation-collision-and-bockstein-carry-gate
  - THM-2390-septimal-layer-kraft-peeling-and-heavy-word-reduction
  - THM-2395-common-core-successor-shell-and-forced-hole-escape-tax
script: 04-computation/lrc14_common_core_forty_nine_orbit_thm2396.py
output: 05-knowledge/results/lrc14_common_core_forty_nine_orbit_thm2396.out
script_sha256: 927f6bb05dab4100415441392e3072b370c3a4837cfa35f5a76f576d2e0f8f85
output_sha256: c84eef5e9712edd12ec7e3fcbc4a01cefab28d2c82402a2079339dd922c6b7c0
hash_basis: working-tree bytes (LF)
---

# THM-2396 -- common-core forty-nine-orbit word incompatibility

**PROVED + VERIFIED-EXACT + TWICE INDEPENDENTLY HOSTILE-AUDITED.**

THM-2393 leaves one literal common-core chain as the only packet not
already assigned its universal clean-hole floor. THM-2394 turns every
clean-free high-safe seven-root fibre into a six-address transversal.
The missing move is not another one-fibre estimate: the unique top
`q_*` word couples seven such fibres inside one orbit of length
forty-nine.

Keeping that whole orbit gives a quantitative hitting theorem, not just
a contradiction at zero mass. The proof deliberately enlarges the
physical universe by letting all non-normalized word translates vary
independently. The enlarged universe is still empty. Consequently every
generic high-safe physical orbit contains a clean root.

## 1. The generic forty-nine-root orbit

Retain THM-2393/2394's residual speed and valuation data

```text
M=1,

(C_1,C_2,c_1,c_2)=h(1,13,13,169),

gcd(h,91)=1.                                           (1)
```

Global clean-set emptiness is **not** assumed. Use the actual clean set
and its arbitrary mass:

```text
K=1_(E_H)+sum_(q_i!=q_*) 1_(D_(q_i)),

S={K=0} intersection D_(q_*)^c intersection D_(c_3)^c
    intersection
    (D_(C_1) union D_(C_2) union D_(C_3))^c,

delta=mu(S)>=0.                                        (1a)
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
on the seven base roots. Fix such a generic `z` and suppose, temporarily,
that its forty-nine roots contain no point of `S`. Then:

- on the top bin (5), THM-2391 says that `D_h,D_(13h)` partition the
  two guard addresses and all four lower `q` addresses lie inside that
  guard pair;
- on every other bin, the proof of THM-2394 applies **locally**. Its
  scalar-cover inequality (9a) and collision-cage implication (9b) are
  pointwise. If (9c) failed at a `K`-hole, the safety of
  `q_*,c_3,C_3` on this bin and the absence of both `A,B` would make
  that root a point of `S`, contrary to the temporary assumption.
  Hence the guard and four lower `q` addresses are six distinct points,
  with unique hole

  ```text
  d=B,              or              d=A=C,            (6)
  ```

  where `A,B,C` are the addresses of
  `D_h,D_(13h),D_(169h)`.

Thus an `S`-free high-safe orbit satisfies the word constraints
(5)--(6) simultaneously on the same forty-nine physical points. This
simultaneity is the information lost by one-fibre capacity and successor
calculations.

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

## 4. Orbitwise localization and a uniform floor

Let

```text
H_0=D_V^c intersection D_(13V)^c,

N_S(z)=sum_(j in Z/49Z) 1_S((z+j)/49).              (17)
```

Sections 1--3 prove the orbitwise contrapositive

```text
N_S(z)>=1                         for almost every z in H_0.  (18)
```

Indeed, `N_S(z)=0` supplies precisely the temporary local hypothesis
used in Section 1 and would embed that physical orbit into the empty
relaxed universe. Endpoint sets and pullbacks of the pointwise null
exceptions remain null.

The forty-nine-root disintegration is exact:

```text
integral_T N_S(z) dz=49mu(S)=49delta.                (19)
```

Combining (3), (18), and (19) yields

```text
delta>=66/(91*49)
     =66/4459
     =396/26754.                                    (20)
```

This is coefficient-independent. THM-2392 therefore gives, within the
common-core packet, the explicit descendant floors

```text
fixed q_*-labelled charged cell:
  delta/52>=33/115934;

complete low-blocker charged cell:
  delta/78>=11/57967;

same-parent owner-resolved F_7 x F_13 tensor:
  delta/338>=33/753571.                              (21)
```

THM-2393 already proves `delta>=1/26754` outside the common-core chain.
Together with (20), every packet in the last septimal lane now obeys

```text
delta>=1/26754.                                     (22)
```

THM-2395's actual `delta=0` common-core branch is therefore superseded
as a live closure obligation. Its factor-free successor-role jet remains
a valid reusable counterfactual sidecar: the present theorem does not
identify any clean carrier with that role jet or with a terminal target.

This theorem does **not**:

- align the positive clean word with a terminal endpoint reference;
- force a target-thirteen phase after restoring the later factors;
- exclude a positive-clean scalar row;
- decrement the `165`-row ledger; or
- prove LRC(14).

The next sharp target is no longer a common-core zero-clean automaton or
a coefficient-dependent positivity argument. It is a uniform
positive-clean landing theorem: preserve the owner-resolved
`F_7 x F_13` word from THM-2392 through the terminal endpoint filtration
without replacing its phase by unlabelled energy.

## 5. Exact companion

The dependency-free companion:

- constructs all `1029+1029` distinct word supports in (9);
- fixes `A`, exhausts all seven top rows, all top-compatible guards,
  all forty-nine `B` and `C` translates, and every four-q packing;
- verifies the full ledgers (14)--(16);
- checks all `980` final `C`-word tests directly;
- independently reconstructs the thirteen zero masks and the sharp
  two-row agreement cap (15);
- retains the `2058` packing and twenty capacity-boundary positive
  controls;
- checks the exact disintegration floor (20), all three descendant
  cell floors (21), and the universal minimum (22).

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
