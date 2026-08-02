# LRC14 reflected tails: orientation correction and repaired `3m >= 2D` cone

**Correction + frozen exact referee, 2026-08-01.**  The conclusion of
[`lrc14-reflected-two-thirds-cone-proof-codex-20260801.md`](lrc14-reflected-two-thirds-cone-proof-codex-20260801.md)
was not justified by its published tail policies: two split policies confused
an oriented ratio chart with a physical ordering of packet levels.  This note
retracts that proof route, preserves the exact finite data that survives, and
gives a replacement proof of the same scoped conclusion

```text
D >= 6,                       3m >= 2D,
m=min q_e,                    D=max q_e-min q_e.             (1)
```

The repair applies only inside the sufficient reflected `k=1` family of
[THM-2941](../01-canon/theorems/THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary.md).
It is not a physical-survivor census and not a proof of LRC(14).  In
particular, it does **not** repair or promote the cap-`3` / half-cone artifact.

## 1. The precise quantifier failure

Write an ordered physical pair as `(i,j)` and its ratio as

```text
r=q_j/q_i.
```

A `low` lane means `r<1`, hence `q_j<q_i`; a `high` lane means `r>1`,
hence `q_i<q_j`.  Reversing both the pair and the side does not reverse the
physical event.  On

```text
H=(1,2,3,4,6,12),
```

the old lanes were `low(1,0)` and `high(0,1)`.  Both say `q_0<q_1`.
The admissible assignment `(q_0,q_1)=(2,1)` satisfies neither.  Thus the two
channel tables were individually correct but did not cover all assignments.

On

```text
H2=(1,3,4,6,8,12),
```

the old lanes were `low(2,1)` and `high(0,1)`.  Their simultaneous failure is
the consistent order `q_2<q_1<q_0`; the cap-admissible assignment
`(q_0,q_1,q_2)=(5,4,2)` is an exact witness.  In tournament language, negate
each available weak-order predicate.  The failure edges here are acyclic, so
a hostile level order exists.  This is a useful diagnostic, not an equivalence
between the reflected certificate and a tournament.

The first failed implication was therefore

```text
separate oriented lane checks  =>  assignment cover.        (2)  FALSE
```

The CSP census, located margins, primitive overlaps, and each lane conditional
on its declared orientation remain exact data.  The predecessor's overall
conclusion is nevertheless retracted and is not used as a theorem below.

## 2. Independent reconstruction of the surviving proof graph

The repair referee imports pinned predecessor code and data, but reruns every
surviving downstream obligation:

- it regenerates all `561` residual bodies;
- it rechecks the independent universal repeated-level theorem on each body;
  its two exception bodies have empty intersection with this universe, so
  `P=Q` may be omitted from every later channel bank;
- it regenerates `6,358` constrained edges and the `15,995`-element primitive
  bank at cap `5/2`;
- it recomputes the CSP census `389` closed / `172` traps in both search
  orders, with trap digest
  `b8e1f964d610af641a229d06951c4805167af4310429740501809034e9b2a716`;
- after removing all four split rows on `H,H2`, it recomputes all `160`
  constrained located policies and all `1,566` physical orientation controls,
  with located digest
  `e0a6c267f4f77017ea56229f3f2c80ba1b1a9c00b065ad9ab0540397ce238a75`;
- it identifies exactly ten remaining unconstrained bodies with genuine
  `all`-interval policies and rechecks their analytic tails and all `300`
  finite heads.

Thus no statement below inherits the predecessor's retracted conclusion.
Only independently regenerated components enter the repaired assembly.

## 3. A same-pair closed Farey cover

For `H`, freeze the physical ordered pair `(0,1)` once and for all.  The cone
implies `r=q_1/q_0 in [2/5,5/2]`.  Five closed intervals cover it without a
change of owner:

| interval for `q_1/q_0` | tail start | endpoint numerator |
|---|---:|---:|
| `[2/5,3/7]` | 48 | `16/5` |
| `[3/7,1/2]` | 38 | `22/7` |
| `[1/2,2/3]` | 26 | `3` |
| `[2/3,11/9]` | 16 | `8/3` |
| `[11/9,5/2]` | 7 | `14/9` |

Every right endpoint is exactly the next left endpoint.  The intervals are
closed, so the overlaps are harmless; the exact assignment controls
`(2,1)` and `(1,2)` land in lanes `(2,3)` and `5`, respectively.  For `H2`,
one whole interval `[2/5,5/2]` on the same fixed pair `(0,1)` works from
scale `14`, with endpoint numerator `26/5`.

This fixed-owner cover is analogous only at the level of information type to
the joint-label lift in [THM-3072](../01-canon/theorems/THM-3072-a4-flag-three-c2-tomography-and-edge-cycle-cospan.md):
the old construction kept a ratio chart and a renamed pair as separate
marginals, while the repair carries their coupling.  There is no claimed
`A4`/LRC map, action, or intertwiner.

## 4. Coupling packet scale to singleton debt

From `(1)`,

```text
q_max=m+D <= 5m/2.                                        (3)
```

For a reduced channel `(P,Q)` at common scale `g`, the selected pair has
maximum level `g max(P,Q)`.  Hence

```text
m >= ceil(2g max(P,Q)/5).                                 (4)
```

The six-singleton debt at a common lower level `ell` is

```text
Debt_E(ell)=sum_(e in E) e/[7(ell L-e)],
L=14 lcm(E).                                               (5)
```

It is strictly decreasing in `ell`.  Every replacement head therefore uses

```text
ell(g)=max(4,ceil(2g max(P,Q)/5)).                         (6)
```

The `1,566` regenerated constrained controls and the `300` unaffected
whole-tail heads use the independently forced integer-corner bound `m>=4`.
The five analytic `H` tails use the weaker channel-free consequence
`ell(s)=max(4,ceil(2s/5))`, where every tail channel has
`min(P,Q)>=s`.  The `H2` whole tail retains the independently safe base level
`4`; its finite heads still use `(6)`.

The six replacement lanes contain exactly `605` reduced finite-head channels.
Every one is checked directly until its first uniform transported scale.  The
longest direct prefix has length four, for channel `(7,3)` on `H`, from
`g_0=2` through the tail start `g=6`.  The weakest first-tail margin is

```text
4917100104747739 / 477807176182256408075                 (7)
```

on channel `(42,17)` of `H`.

## 5. Why one first-tail check proves every later scale

Let `a,b` be the two labels.  At a body-safe cell the primitive skeleton is
independent of `g`, while the homotopy displacement is

```text
eta_g=g(Qa-Pb)/(PgL-a).
```

For every later integer `g`, exact subtraction gives

```text
|eta_g|-|eta_(g+1)|
 = |Qa-Pb|a/[(PgL-a)(P(g+1)L-a)] >= 0,                   (8)
```

strictly unless `Qa=Pb`.  The transport displacement also falls strictly:

```text
a/(PgL-a)-a/(P(g+1)L-a)
 = aPL/[(PgL-a)(P(g+1)L-a)] > 0.                         (9)
```

The favourable factor `c_g^{-1}=1+a/(PgL-a)` is discarded only after the
primitive bracket is positive.  Meanwhile `(6)` is nondecreasing and `(5)`
is nonincreasing.  Therefore, once

```text
skeleton-2|eta_g| > Debt_E(ell(g))                        (10)
```

holds at the first tail scale, it holds at every later scale.  The referee
checks the direct reflected overlap against the transport lemma at that first
scale, verifies its slope gate `|eta|<1`, and uses `(8)`--`(10)` for the
infinite suffix.  Earlier scales are direct exact computations.  No sampled
finite suffix is promoted to a tail.

For the analytic lane scale `s`, the exact one-step gain is the sum of

```text
(12/49)(1/s^2-1/(s+1)^2) > 0,
Na/[(Ls-a)(L(s+1)-a)] > 0,
Debt_E(ell(s))-Debt_E(ell(s+1)) >= 0.                     (11)
```

Thus each positive lane threshold is uniform for all larger `s` as well.

## 6. Repaired conclusion and boundary

The regenerated `389` coarse closures, `160` two-orientation located
policies, ten unaffected whole tails, and two repaired bodies exhaust the
`561` residual bodies.  Consequently the reflected THM-2941 sufficient
family closes throughout `(1)`.  Its possible certificate-failure region is
again confined to

```text
D >= 6,                       1 <= m < 2D/3.              (12)
```

This restores the cap-`5/2` result only.  Any stronger cone must repeat the
physical assignment-cover audit; the former cap-`3` conclusion remains
retracted pending its own repair.

## 7. Reproduction and seals

```bash
python3 04-computation/lrc14_j7_reflected_two_thirds_tail_orientation_repair_thm2941.py
python3 -O 04-computation/lrc14_j7_reflected_two_thirds_tail_orientation_repair_thm2941.py --output /tmp/lrc14_two_thirds_repair_O.out
cmp 05-knowledge/results/lrc14_j7_reflected_two_thirds_tail_orientation_repair_thm2941.out /tmp/lrc14_two_thirds_repair_O.out
```

Ordinary Python, optimized Python, and a fresh stored-output replay are
byte-identical.  Final seals:

```text
source_sha256=420d21d0c65cfbd57d8cc926c04993b14e61d2554fb291a836af949d1fd664ed
output_sha256=87db762ba3ff6794748a2a19a4bf192a8cef8d6cebe7c619c871579da669864f
semantic_sha256=525c8d2d89cf11c9fcc6260dbf81e57d0e47c999049c89b9148ae731dbac4730
lane_digest=1d90551e665eb8f4868affb7650fdd4e338095812390857525af7a282e87acd3
replacement_head_digest=0253483de002db1e2fa32a3feca1811f28b355b43d6f2426d19ba6fdfe58ffdc
unaffected_whole_head_digest=4787fd19180878e0b5d502609121875eb8a4fc1f45e64a15b8be533fa5547dd7
unaffected_monotonicity_digest=e23f2a4c5e572391d03ec3019239f1818ff0cb92740ec549b8a107a6b96622eb
```

The result file records the independent CSP, located-policy, whole-tail,
endpoint, assignment, singleton-debt, and infinite-monotonicity controls.
