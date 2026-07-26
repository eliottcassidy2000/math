---
id: THM-2381
title: "One-top-one-blocker septimal root-stalk closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. In
  THM-2367's only-c_3-dominant regime M>0, nu_7(H)<M, the k=2,
  (t,b)=(1,1) alternative is empty. The guard, four lower unit masks,
  and the unique strictly lower blocker have total seven-bin weight
  exactly seven, so they form a one-fold partition off the unique top
  unit, the equal-depth blocker, and c_3. On a thirteen-root stalk, a
  low-danger/top-blocker-safe/c_3-safe quotient phase would force the
  guard word of size three or four inside the top-unit word of size one
  or two. Such phases have exact mass 36/343: on a 7^(M+1)-root fibre,
  the low and equal-depth blocker words have sizes N/7 and intersection
  N/49, while the higher blocker is constant and safe on base mass 6/7.
  This empties one of the four k=2 septimal alternatives uniformly, but
  removes no thirteen-adic row, leaves the ledger at 165, and does not
  prove LRC(14).
source: codex-2026-07-25-one-top-one-blocker-root-stalk
depends_on:
  - THM-2367-septimal-root-averaging-graft-and-cover-alignment
related:
  - THM-2378-hard-septimal-root-stalk-closure
  - THM-2377-septimal-valuation-collision-and-bockstein-carry-gate
  - THM-2382-saturated-septimal-seven-bin-root-fibre-closure
script: 04-computation/lrc14_one_top_one_blocker_root_stalk_thm2381.py
output: 05-knowledge/results/lrc14_one_top_one_blocker_root_stalk_thm2381.out
script_sha256: 36632876efc5d32ae376edd731df66a801ba906a6cf1636237da3bf7fa32ca50
output_sha256: 4c1b23909d5cf162c5194fbaee47a569a76c5f15e87ae097346c376cfd35afa2
hash_basis: working-tree bytes (LF)
---

# THM-2381 -- the one-top-one-blocker septimal branch is empty

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2378 empties the hard `W=k=1` lane by retaining the pointwise
thirteen-root word which signed currents and chamber integrals lose. The same
coordinate also reaches one of the three nonsaturated `k=2` alternatives.
The mechanism is

```text
(t,b)=(1,1)
  -> seven units of lower capacity in every unoccupied top bin
  -> exact lower partition off three absorbers
  -> guard root word 3/4 inside top-unit root word 1/2
  -> contradiction on a mixed-depth anti-shield set of mass 36/343.  (1)
```

No drift promotion, target landing, or tournament orientation is needed.

## 1. The live `k=2`, `(t,b)=(1,1)` notation

Retain THM-2367 Section 6:

```text
M=max(nu_7(H),nu_7(q_1),...,nu_7(q_5))>0,

nu_7(H)<M,

nu_7(c_3)>M,

k=#{j in {1,2}:nu_7(c_j)<=M}=2.                    (2)
```

Suppose the surviving profile is

```text
(t,b)=(1,1),                                       (3)
```

where `t` counts ordinary unit labels at depth `M` and `b` counts
`c_1,c_2` at depth `M`. Name the roles

```text
nu_7(q_*)=M,              nu_7(q_i)<M for q_i!=q_*,

nu_7(c_low)<M,            nu_7(c_top)=M,

{c_low,c_top}={c_1,c_2}.                            (4)
```

As in the canonical first-depth-one scalar model, all three blockers are
divisible by thirteen and the six labels `H,q_1,...,q_5` are units modulo
thirteen. Put

```text
mathcal A=D_(q_*) union D_(c_top) union D_(c_3),

L
 =1_(E_H)
  +sum_(q_i!=q_*)1_(D_(q_i))
  +1_(D_(c_low)).                                   (5)
```

The widths in `L`, measured in `1/7`-arc units, sum to

```text
2+4+1=7.                                            (6)
```

## 2. Seven-bin exactness off the three absorbers

Let

```text
N=7^(M+1).
```

Consider a generic orbit

```text
O_x={x+j/N:j in Z/NZ}.                              (7)
```

Because `N|c_3`, the `c_3` mask is constant on this orbit. A factor
`v` of septimal depth exactly `M` occupies one complete residue bin
`j mod 7`. A factor of depth `r<M` and width `a/7`, with
`a in {1,2}`, occupies exactly

```text
aN/49
```

points in each of the seven bins. Indeed, after fixing `j mod 7`, its
phase runs over a `7^(M-r)`-grid, each point with multiplicity `7^r`.
Genericity excludes the finitely many aligned endpoints.

It follows from (6) that in every bin the total incidence of the lower
load `L` is exactly

```text
N/7.                                                (8)
```

On a `c_3`-safe orbit, take a bin occupied by neither `q_*` nor
`c_top`. The scalar cover forces its `N/7` points to be covered by the
lower masks. Since their total incidence is also `N/7`, they form an
exact one-fold partition there. Therefore

```text
L=1                         almost everywhere on mathcal A^c. (9)
```

This is the precise analogue of THM-2378's hard-lane identity. The
identity is not asserted in the other `k=2` profiles: their lower
weights or absorber words differ.

## 3. The thirteen-root word squeeze

Write

```text
c_low=13C,             c_top=13A,             c_3=13B.       (10)
```

Choose a quotient phase `y` away from the finite set of blocker
endpoints and all endpoint pullbacks under

```text
x_h=(y+h)/13,                         h in F_13.      (11)
```

Every blocker mask is constant on this stalk:

```text
1_(D_(c_low))(x_h)=1_(D_C)(y),

1_(D_(c_top))(x_h)=1_(D_A)(y),

1_(D_(c_3))(x_h)=1_(D_B)(y).                         (12)
```

Suppose

```text
y in D_C minus (closure(D_A) union closure(D_B)).    (13)
```

Then `c_low` contributes one to `L(x_h)` on all thirteen roots, while
both blocker absorbers are absent. Whenever `x_h` is also `q_*`-safe,
it lies outside `mathcal A`; (9) says `L(x_h)=1`. Hence every other
lower mask, including the guard, is absent at that root. Thus

```text
{h:x_h in E_H} subset {h:x_h in D_(q_*)}.            (14)
```

Multiplication by a thirteen-unit permutes the root grid. A centered
arc of length `2/7` meets it in three or four points, while a centered
arc of length `1/7` meets it in one or two. Therefore (14) would give

```text
3<=#{h:x_h in E_H}
  <=#{h:x_h in D_(q_*)}<=2,                          (15)
```

which is impossible. Consequently a scalar cover in this branch would
force the containment

```text
D_C subset closure(D_A) union closure(D_B).          (16)
```

## 4. The mixed-depth anti-shield has exact mass `36/343`

The required anti-shield is slightly different from THM-2378's
`5/49` lemma because one target now lies at depth exactly `M`.

> **Mixed-depth anti-shield lemma.** If positive integers `C,A,B`
> satisfy
>
> ```text
> nu_7(C)<M=nu_7(A)<nu_7(B),                         (17)
> ```
>
> then
>
> ```text
> measure(
>   D_C minus (closure(D_A) union closure(D_B))
> )=36/343.                                          (18)
> ```

**Proof.** Again put `N=7^(M+1)` and use the inverse fibre

```text
y_j=(z+j)/N,                         j in Z/NZ.       (19)
```

Since `N|B`, the `B`-mask is constant on the fibre. Its safe base set

```text
G=T minus closure(D_(B/N))
```

has measure `6/7`. For every generic `z`, the low and equal-depth
masks obey

```text
#{j:y_j in D_C}=N/7,

#{j:y_j in D_A}=N/7,

#{j:y_j in D_C intersection D_A}=N/49.              (20)
```

For the last count, `D_A` first selects one residue class modulo seven.
Inside that class the `C`-phase traverses a
`7^(M-nu_7(C))`-grid, so exactly one seventh of its `N/7` points are
`C`-dangerous.

Thus every `z in G` has exactly

```text
N/7-N/49=6N/49                                      (21)
```

roots in `D_C` but outside both target masks. Disintegration under
`y -> Ny` gives

```text
(6/7)(6/49)=36/343,                                  (22)
```

proving (18). Boundary fibres form a null set. QED.

Equations (16) and (18) contradict one another. Hence:

> **Branch closure.** No scalar cover satisfies the only-`c_3`-dominant
> hypotheses (2) together with `k=2`, `(t,b)=(1,1)`.

## 5. Exact scope and the new residual

This theorem and THM-2378 refine THM-2367's list as follows:

```text
k=1, (t,b)=(1,0)                         empty [THM-2378],

k=2, (t,b)=(1,1)                         empty [this theorem],

k=2, (t,b) in {(1,0),(2,0),(5,2)}        still open. (23)
```

The proof does not silently cover the remaining cases:

- for `(1,0)`, the lower load has weight eight, so the exact
  one-fold identity (9) is unavailable;
- for `(2,0)`, exactness survives, but two top-unit root words can have
  total size four, so the guard support squeeze alone does not
  contradict them; simultaneous low-blocker danger instead forces a
  new two-source/one-target containment;
- for saturated `(5,2)`, the seven top masks themselves partition the
  septimal bins, and the faithful object is their moving permutation
  word.

The statement is uniform in the thirteen-adic roles of `c_1,c_2`; it
does not remove a thirteen-adic valuation row. The ledger remains
`165`, and LRC(14) remains open.

The formerly intended THM-2381 two-high drift calculation concerned
THM-2378's now-empty hard lane. Its exact reduced carrier and residue
sweep remain valid as a reusable abstract sidecar, but are not needed
for this live branch closure and are not asserted here.

## 6. Exact companion

The dependency-free companion:

- checks the generic thirteen-root word sizes `1/2` and `3/4`;
- checks the per-bin `aN/49` count for lower factors and the one-bin
  support of depth-`M` factors through three septimal layers;
- checks the mixed-depth counts `(N/7,N/7,N/49)` on exact fibres;
- independently measures
  `D_1 minus (D_7 union D_49)` by its rational wall cells and obtains
  `36/343`; and
- records the precise branch residual (23).

Run

```bash
python3 04-computation/lrc14_one_top_one_blocker_root_stalk_thm2381.py
python3 -O 04-computation/lrc14_one_top_one_blocker_root_stalk_thm2381.py
```

Both transcripts must byte-match

```text
05-knowledge/results/lrc14_one_top_one_blocker_root_stalk_thm2381.out
```

after LF normalization. Every executable check raises explicitly under
optimized Python.

## 7. Independent hostile audit

An independent read-only audit reconstructed the proof from the scalar-cover
orbit counts rather than trusting the companion. It checked:

- every subtop factor contributes exactly its width times `N/49` in a
  fixed top residue bin;
- lower weight seven gives (9) pointwise, not merely after integration;
- `c_top` is not treated as fibre-constant: its word has size `N/7`,
  meeting the low word in exactly `N/49`;
- the higher `c_3` word is constant and its quotient-safe base has mass
  exactly `6/7`;
- the thirteen-root squeeze, endpoint qualifications, branch typing,
  and zero-row consequence.

The audit supplied sharp failure controls. If the low depth is raised to
`M`, take `C=A` and the anti-shield can vanish. If the alleged high depth
is lowered, take `B=C`. If the lower total width differs from seven, the
one-fold inference fails. It also caught that Python `assert` would vanish
under optimization; all checks now raise explicit `RuntimeError`s. Normal and
optimized transcripts byte-match the stored output. QED.
