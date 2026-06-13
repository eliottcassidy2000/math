---
id: HYP-2239
status: OPEN method hypothesis with S663 LRC perturbation evidence and S664 cycle atlas
source: codex-2026-06-05-S663/S664
related:
  - HYP-2238
  - HYP-2237
  - HYP-2236
  - HYP-2235
  - HYP-2231
  - HYP-2230
  - HYP-2171
  - HYP-2167
  - THM-401
---

# HYP-2239: Triune Carrier Applications to LRC and Sibling Problems

## Claim

HYP-2238's sum/product/fraction trinity transfers to LRC as:

```text
sum face      -> additive wall/pair packets
product face  -> local product/gcd/obstruction shell
fraction face -> recursive carry/owner/continuant boundary state
```

For LRC `n=14`, the fraction face is not optional.  The carry vector

```text
v = r + 27k
```

is the continued-fraction-like boundary state: it records how a `Res_27`
shadow is lifted into the `n=14` clock, parity, apex, and owner/certificate
world.  If `k` is forgotten, floor and strict rows collide.

## S663 LRC Evidence

`04-computation/triune_carrier_applications_s663.py` builds triune records for
the three known least-positive floor shadows:

```text
AP, Vstar, 2AP.
```

All three have exact floor value:

```text
M = 1/14.
```

They also share the same product face:

```text
C=27 gcd shells = ((1,9), (3,3), (9,1))
gcd-shell mass = 27.
```

The base carry support is empty for all three.  This is the zero-carry floor
section.

The script then adds local `+27` carry perturbations of Hamming weight one and
two over each floor atom.  Every checked perturbation is strict:

```text
AP:    weight 1: 13 strict, min M=1/13
AP:    weight 2: 78 strict, min M=1/12
Vstar: weight 1: 13 strict, min M=2/25
Vstar: weight 2: 78 strict, min M=1/12
2AP:   weight 1: 13 strict, min M=1/13
2AP:   weight 2: 78 strict, min M=1/12
```

No local perturbation is below floor, and none remains exactly floor.

## Projection-Collision Test

The decisive test groups all base and perturbation rows in two ways.

First, keep only the sum/product shadow:

```text
Res_27 shadow additive packets + C=27 gcd/product shell.
```

This leaves three mixed groups:

```text
AP group:    1 floor + 91 strict
Vstar group: 1 floor + 91 strict
2AP group:   1 floor + 91 strict
```

Then add the fraction face:

```text
carry word k, carry support, mod-14 word, parity word,
apex/pair-apex data, and continuant.
```

The mixed groups disappear:

```text
full triune key: 0 mixed floor/strict groups.
```

So, in this bounded atlas, the carry/fraction face exactly repairs the leak
created by the scalar sum/product quotient.

## Cross-Problem Transfer

The same pattern suggests the following table.

| Problem | Sum Face | Product Face | Fraction / Boundary Face |
|---|---|---|---|
| LRC `n=14` | odd-wall and pair-sum packets | `C=27` gcd shells, Pillai mass, local obstructions | carry word `k`, owner route, continuant |
| OCF / `H(T)` | odd-cycle packet coefficients | strong-component product and forbidden `H` packets | deletion/substitution boundary continuant for Hamiltonian paths |
| Tournament decks | deleted-card loss sums | deck product/scissors component packets | paired card-owner derivative |
| Unit distance | unit-edge spine plus tile/bulk packets | direction support, norm/unit-shell factors | point-deletion frontier or ear owner |
| Finite-field Kakeya/Falconer | distance and pinned-distance packets | line-direction products, concurrency factors | owner of pins and line-choice recursion |
| Goldbach/Lemoine | prime-pair sums `E=p+q`, `O=p+2q` | singular-series local obstruction products | ordered-pair reconstruction `q=O-E`, `p=2E-O` |
| CH / forcing | cardinal/equinumerosity sums | model-local consistency factors | generic extension boundary state |
| pi/e trace-norm | trace sums and power sums | norm/discriminant products | branch sheet of `T^2-S*T+P` |

S663 ranks these routes by triune fit, finite testability, and proof leverage.
The application-route tournament is transitive:

```text
LRC14 carry-owner theorem
> OCF / H(T)
> Tournament decks
> Unit distance
> Finite-field Kakeya/Falconer
> pi/e trace-norm
> Goldbach/Lemoine
> CH / forcing.
```

The ranking is not a claim that the later routes are unimportant.  It says the
first routes have the best combination of exact finite tests and immediate
repo payoff.

## S664 Everywhere-Cycle Upgrade

`04-computation/triune_cycle_everywhere_s664.py` upgrades the trinity from a
static record to a projection-repair loop:

```text
sum -> product -> fraction -> sum.
```

Here:

```text
sum      = additive aggregation: packets, traces, walls, counts, moments
product  = local factorization: gcd/norm shells, components, sieves
fraction = boundary memory: carry, owner, branch, lift, deletion, recursion
```

The three arrows mean:

```text
sum -> product
  additive packets expose or localize factor data;

product -> fraction
  local factors are still a quotient and need a lift/owner/branch;

fraction -> sum
  boundary state re-expands into the next additive packet ledger.
```

The S664 face-dependency tournament uses vertices

```text
sum, product, fraction, raw_scalar
```

and domain-weighted dependency support as the observable.  It recovers the
same nontransitive core as S662:

```text
score_hist={0: 1, 2: 3}
directed_3cycles=1 [('sum', 'product', 'fraction')]
sccs=[['fraction', 'product', 'sum'], ['raw_scalar']]
hamiltonian_paths=3
```

with weights:

```text
sum -> product       support 54
product -> fraction  support 52
fraction -> sum      support 51
```

All three faces beat `raw_scalar`.

S664's strongest exact projection-repair instances are:

| Domain | Projection Leak | Fraction Repair |
|---|---|---|
| LRC14 | `Res_27` sum/product shadow mixes floor and strict rows | carry-continuant / owner state |
| Tournament decks | full deck plus global scalar repairs collide through `n=6` | paired `(card, deleted score)` |
| Pi/e trace-norm | trace/norm leaves the two-root branch unresolved | discriminant or branch sheet |

It also adds weaker but actionable labs:

| Domain | Next finite probe |
|---|---|
| OCF / `H(T)` | encode deletion/substitution macro-words as continuants |
| Unit distance | compare edge sums + direction/norm products + point-deletion owners on `n=21/22` cores |
| Goldbach/Lemoine | track shared prime-pair columns under `E=p+q`, `O=p+2q`, `q=O-E`, `p=2E-O` |
| Perfect/aliquot | attach aliquot orbit/preperiod owner to divisor-sum/product shadows |
| A000568/LRC quotient | treat every quotient count as needing its marked observer fiber |

## Assumption Challenge

Candidate vertices considered for Tournament Analysis:

```text
runners, residues, wall pairs, gcd shells, carry words, proof atoms,
owner routes, problem domains, application routes, representation faces,
files, deleted vertices, unit-distance points, prime pairs, models, generics,
and proof obligations.
```

S663 uses application routes as vertices for the cross-problem tournament, but
uses LRC rows as objects in the projection-collision audit.  Raw runners are
not the right vertex set for the proof theorem; the load-bearing object is the
triune record.

S664 uses faces as vertices for the cycle tournament and domains as vertices
for the route tournament.  The preserved predicate is not raw truth of each
problem.  It is:

```text
which carrier repairs a projection collision?
```

The bounded test destroys full owner/certificate detail: it keeps only carry,
apex, pair-apex, parity, and continuant data.  Future work should attach the
HYP-2165 owner route explicitly.

## Next Experiments

- Build `lrc14_triune_owner_derivative_s665.py`: add HYP-2165 owner routes to
  the carry-continuant state and test whether any mixed floor/strict groups
  survive deeper perturbations.
- Build an OCF continuant toy: encode deletion or substitution macro-words as
  continuants and compare to direct Hamiltonian-path DP.
- Apply the triune ledger to unit-distance `n=21/22` cores: edge sums,
  direction/norm product support, and point-deletion frontier owners.
- Revisit pi/e trace-norm shadows with a branch-sheet/fraction face rather
  than only the Vieta trace/norm pair.
- Build a unit-distance `n=21/22` triune impairment lab: hold edge counts and
  direction/norm support fixed while changing point-deletion owner state.

**See:** `04-computation/triune_carrier_applications_s663.py`;
`05-knowledge/results/triune_carrier_applications_s663.out`;
`07-reflections/triune-carrier-applications-lrc-s663.md`;
`04-computation/triune_cycle_everywhere_s664.py`;
`05-knowledge/results/triune_cycle_everywhere_s664.out`;
`07-reflections/triune-cycle-everywhere-s664.md`; HYP-2238, HYP-2237,
HYP-2236, HYP-2235, HYP-2231, HYP-2230, HYP-2171, HYP-2167, HYP-2165,
HYP-2164, THM-401.
