---
id: THM-1153
title: THE MULTI-DELETION HARMONIC CROWN — every actual LRC counterexample obeys an all-N reciprocal deletion law; at N=14 it compresses the top seven speeds by an exact triangular recurrence and exposes r=7 as the zero-coefficient wall
status: PROVED analytic implication, conditional only where the displayed lower-case LRC bound is invoked. For N=14, r=1,...,6 use the repository's named LRC(13),...,LRC(8) citation inputs. The all-N coefficient, scalar corollaries, ordered recurrence, and tournament audit are finite-exact replayed. This does NOT cross r=7 or prove LRC(14)
source: codex-2026-07-18-S75
depends_on:
  - THM-883
  - THM-1015
  - THM-1008
  - lower-case LRC citation
related: [THM-1149, THM-1152, THM-1025, HYP-7678]
script: 04-computation/lrc14_multideletion_harmonic_crown_codex_20260718.py
output: 05-knowledge/results/lrc14_multideletion_harmonic_crown_codex_20260718.out
lean:
  - 04-computation/lean/TournamentH7/TournamentH7/FragmentationLemma.lean
  - 04-computation/lean/TournamentH7/TournamentH7/LRCMultiDeletionCrown.lean
---

# THM-1153 -- the multi-deletion harmonic crown

Use the standard relative-speed convention: an `N`-runner instance is a set
`V` of `N-1` distinct positive integer speeds and

```text
M(V)=max_t min_(v in V) ||vt||.
```

The new object is not a single deletion's essential region. It is the
simultaneous reciprocal obligation carried by every deletion cut.

## 1. All-N harmonic-crown theorem

> **Theorem A (multi-deletion harmonic crown).** Let `M(V)<1/N`. Choose
> `S subset V` with `1<=|S|=r<N/2`, put `W=V\S`, and let `m=max W`. Then
>
> ```text
> sum_(s in S) 1/s >= (N-2r)(M(W)-1/N)/m.            (1)
> ```
>
> If the lower instance supplies `M(W)>=1/(N-r)`, then
>
> ```text
> sum_(s in S) 1/s >= r(N-2r)/(N(N-r)m).             (2)
> ```

**Proof.** Write `mu=M(W)`. If `mu<=1/N`, the right side of (1) is
nonpositive and (1) is immediate. Assume henceforth that `mu>1/N`, and choose
a maximizing real lift `t0`. The function

```text
f_W(t)=min_(w in W) ||wt||
```

is `m`-Lipschitz. Therefore the closed interval

```text
I=[t0-(mu-1/N)/m, t0+(mu-1/N)/m]
```

has length `L=2(mu-1/N)/m`, and every core speed is safe at radius `1/N`
throughout `I`. There is no circle-seam loss: work on the real lift, where
all danger combs are periodic.

Since `M(V)<1/N`, every point of `I`, including its endpoints, lies in the
strict-open danger comb of at least one deleted speed. For one deleted speed
`s`, fragmentation gives

```text
measure(I intersect D_s)<=2L/N+2/(Ns).
```

Covering `I` by all `r` deleted combs and summing yields

```text
L(1-2r/N)<= (2/N) sum_(s in S) 1/s.                  (3)
```

Substitution proves (1). Finally

```text
M(W)-1/N >= 1/(N-r)-1/N = r/(N(N-r)),
```

which proves (2). This is exactly the measure-theoretic content of the
kernel-checked `FragmentationLemma.killer_budget`; the new step is feeding it
the cardinality-sensitive lower-LRC fattening. ∎

The same non-strict conclusion holds under `M(V)<=1/N`. If `M(W)<=1/N`, it
is trivial. Otherwise choose
`alpha=1/N+epsilon<min(M(W),1/(2r))`, apply fragmentation at radius `alpha`,
and obtain

```text
sum_(s in S)1/s >= (M(W)-alpha)(1-2r alpha)/(alpha m).
```

Letting `epsilon` decrease to zero gives (1).

If `s0=min S`, then `sum 1/s<=r/s0`. Consequently (2) gives

```text
s0 <= N(N-r)m/(N-2r).                                (4)
```

For `r>=2`, distinctness makes the reciprocal comparison strict, hence (4)
is strict. At `r=1` it is only non-strict.

## 2. LRC(14): the exact six-rung crown

For `N=14`, define

```text
c_r=r(7-r)/(7(14-r)).
```

Every actual LRC(14) counterexample satisfies

```text
m sum_(s in S) 1/s >= c_r.                            (5)
```

| `r` | `c_r` | `min(S)/m` consequence |
|---:|---:|---:|
| 1 | `6/91` | `<=91/6` |
| 2 | `5/42` | `<84/5` |
| 3 | `12/77` | `<77/4` |
| 4 | `6/35` | `<70/3` |
| 5 | `10/63` | `<63/2` |
| 6 | `3/28` | `<56` |
| 7 | `0` | no information |

The theorem's informative range is `r<7`. At `r=7`, the same budget has the
zero boundary continuation displayed in the table; for `r=8` its algebraic
coefficient is negative and the continuation is trivial. The vanishing at
seven is exact, not numerical looseness. Seven danger combs
each have bulk duty `1/7`, so first-order density can tile the protected
needle with no deficit. This derives the apex-7 wall for the
bulk-plus-endpoint fragmentation/union-bound scheme directly from lower-case
fattening.

## 3. Retaining harmonic mass sharply compresses the top seven

Order the speeds `v1<...<v13` and normalize `x_i=v13/v_i`. Apply (5) to the
top `r` speeds. The retained core has maximum `v_(13-r)`, so

```text
c_r x_(13-r) <= sum_(i=14-r)^13 x_i.                 (6)
```

The compact/dominance reduction gives `x12<13`, while `x13=1`. Using the
full sums in (6), instead of replacing every reciprocal by the smallest one,
gives

| `r` | conclusion |
|---:|---:|
| 2 | `x11<588/5` |
| 3 | `x10<25333/30` |
| 4 | `x9<204967/36` |
| 5 | `x8<8403647/200` |
| 6 | `x7<613466231/1350` |

Thus every compact actual counterexample obeys

```text
v13/v7 < 613466231/1350 = 454419.430370... .          (7)
```

Multiplying only the five scalar bounds after the compact seed gives
`173044872`; retaining harmonic mass improves this by a factor
`380.804297...`. The cruder auxiliary-`1/13` product calculation gives
`20388441216/7`; explicitly it is

```text
13 product_(r=2)^6 [13r(14-r)/((r-1)(13-2r))].
```

Thus (7) is better by a factor `6409.572885...`.

Without THM-1008's compact seed, the `r=1` crown gives `x12<=91/6`, and the
same recursion still gives

```text
x7<=8500889201/16200 = 524746.246975... .              (8)
```

## 4. The ladder's toothpick functional form

Let `T_r=sum_(i=14-r)^13 x_i`. At the extremal triangular relaxation, each
rung appends the new tooth `T_r/c_r`, hence

```text
T_(r+1)=T_r(1+1/c_r),
1+1/c_r=(98-r^2)/(r(7-r)).                            (9)
```

This is the clean H-drift/toothpick law that scalar ratios hide. The state is
cumulative reciprocal mass; one new deletion copies and magnifies the entire
preceding stalk. The pole at `r=7` is the same obstruction as the zero in
(5), seen in the dual recurrence coordinate.

## 5. Relation to the previous frontier

- THM-883 already proves fragmentation and an auxiliary `1/13` killer box.
- THM-1015 says that, on a protected interval at radius `1/14`,
  `sum 1/s < L(7-r)` forces survival. Theorem A is its contrapositive after
  inserting the exact lower-case length `L>=r/(7(14-r)m)`.
- THM-1152 uses the same density engine for two replacements around the AP
  family, with exact essential-region component lengths. Theorem A is the
  family-independent recursive version.
- THM-1149's essential crown lives at the stronger auxiliary radius `1/13`.
  The present theorem instead works at the actual LRC radius `1/14`; this is
  why its seven-wall coefficient is exactly zero.

The lower bounds used in (2) remain explicit inputs. In the LRC(14)
specialization they are the repository's named LRC(13) through LRC(8)
citation hypotheses; this theorem does not independently prove those cases.

## 6. Kakeya needle, Fano/chi7, and tournament audit

The protected interval `I` is the correct Kakeya needle: a one-dimensional
carrier that the **union** of the deleted combs must cover, with endpoint
leakage charged by `1/s`. An individual deleted speed may be redundant on
this chosen needle. For `r<=6`, bulk duty leaves positive length and harmonic endpoint
debt controls scale. At `r=7`, only overlap, ownership, and boundary incidence
can distinguish a cover from a partition.

The runner-order tournament is transitive (score histogram `0,...,12`, no
cycles, singleton SCCs, one Hamiltonian path). It remembers order and loses
the theorem's reciprocal cut mass. Taking the seven deletion obligations as
vertices is also transitive (scores `0,...,6`): it records dependency rank and
the wall, but still loses tooth phase, pair overlap, triple incidence, and
boundary excess. The faithful current carrier is

```text
weighted nested-deletion hypergraph
  + protected-needle endpoint/owner sidecar.
```

This challenges the runner-vertex assumption explicitly. At the wall, a
Fano/`chi_7` probe must be genuinely hypergraphic: seven comb families have
total bulk duty exactly one, while THM-1149's private/triple balance warns
that pairwise orientations cannot see the compensating triple-overlap debt.
That is a precise next target, not a claimed closure.

## 7. Verification and honest frontier

The new Lean module formalizes the actual-radius needle-budget consumer and
the exact triangular ceiling; production of the lower-core interval remains
an explicit analytic input. The companion exact script checks the all-`N` algebra for every `3<=N<=100`,
all legal `r`, every displayed fraction, both recurrence seeds, the scalar
comparison factors, and the deletion-obligation tournament fingerprint.
Normal and optimized runs are byte-identical. Frozen hashes are

```text
source  a902bbcca9bd82b049c59359e185ba2206e7adf5e18984ec700225bc6d2de252
output  957c9b47cf23b0749b010d9385d41cf5a83b0e9ba540a391810dae03211b0de5
```

What is proved is a universal recursive compression through the sixth
deletion and an exact explanation of where first-order density stops. The
remaining frontier includes the seventh-comb overlap/Fano debt, crown
collapse, the uniform `r=5`/`r=6` tails, the twelve-speed equality
classification, and LRC(14).
