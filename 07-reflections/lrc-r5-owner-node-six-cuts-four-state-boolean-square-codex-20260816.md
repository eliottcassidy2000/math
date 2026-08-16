# A five-state source Gray path restricts to a four-state Boolean square

**Status: FINITE-EXACT CANDIDATE, INDEPENDENT AUDIT PENDING.**  Replacing the
collapsed seven-cell coordinate of MISTAKE-417 by the actual source-support
state gives a nonseparable `4 x 13` table of full rank four.  It is one exact
owner-node observable on the canonical `r=5` host, not an exact relation
address, a `U_clock` chronology, a uniform-row theorem, or LRC(14).

## Inheritance pass

The closest proved mechanism is THM-2471's root-service disintegration.  Its
exact `r=5` source profiles satisfy `U_uV_u=0` pointwise, and the corrected
common-owner computation at commit `5a987bef1` multiplies those profiles by
the actual `U_full` endpoint indicators before integration.  Its nonzero
thirteen-class bridge survives, but MISTAKE-417 shows why its seven-cell
support cannot be used: the endpoint owner factor forces `ell=0`.

The canonical hostile is therefore any refiner implied by `OWNER`; Fourier
transforming a boundary delta populates all formal cell frequencies without
creating tensor rank.  The least-used surviving sidecar is the support state
of the source service itself.  It varies inside owner cell zero and is not a
function of the old seven-cell label.

## The six-cut completion and the five-state source path

On every exact source segment, `supp(V)` is the complement of `supp(U)`.  The
only roots ever entering `supp(U)` are

```text
R = {0,6,12} subset F_13.
```

For any nonempty proper subset `S` of `R`, orient every cross-cut pair from
`S` to `R\S`.  Among the three unordered root pairs, each such state has

```text
one missing edge, two one-way edges, zero two-way edges.
```

Thus the natural *combinatorial completion* has size six.  It is not a
six-vertex tournament, but the family of all six oriented cuts of a
three-vertex graph, equivalently six three-vertex tournaments with one edge
deleted.  Complement reverses every arc and pairs these abstract states into
three antipodal pairs:

```text
{0}       <-> {6,12},
{0,6}     <-> {12},
{6}       <-> {0,12}.
```

This gives a precise interpretation of the earlier “tournament with missing
or both-way edges” prompt.  There are no both-way edges here.  The binary
observable is membership across a source-service cut; missing edges are
within one side of that cut.

The source does **not** realize all six cuts.  Its compressed support sequence
over the full common base is the exact Gray path

```text
{0} -> {0,6} -> {6} -> {6,12} -> {12},
```

whose single-root toggle word is

```text
(6,0,12,6).
```

The sixth cut `{0,12}` is globally absent.  In units of the whole circle, the
five path-state measures are

```text
(1,12,2,12,1)/28.
```

This distinguishes three levels that must not be conflated: six abstract
nontrivial cuts, five physically realized source states, and the four states
visible after the owner restriction.

## Why the owner sees exactly four states

Inside `U_full` owner cell zero, the physically realized centre `{6}` is
absent.  The four visible states are

```text
state  bits (component,multiplicity)  source subset
  0              (0,0)                  {0}
  1              (0,1)                  {0,6}
  2              (1,0)                  {12}
  3              (1,1)                  {6,12}.
```

The first bit records the left/right connected component of the owner set;
the second records singleton/doubleton source support.  Complement acts as

```text
state -> state XOR 3,
```

so the two surviving antipodal pairs form a genuine Boolean square.  This is
not a cosmetic `V_4` label: both bits are defined pointwise by intrinsic
Boolean observables before integration.

The owner-excluded realized state is exactly

```text
{6},
```

which lives in the owner-inactive middle region.  Its complement `{0,12}` was
already absent on the full source circle; the owner does not remove it.
Hence the passage is `six-cut completion -> five-state physical path ->
four-state owner square`, not deletion of a complement pair chosen to make a
square.

The exact owner-cell measures are strikingly balanced:

```text
mu(state_i)=1/28,  i=0,1,2,3,
sum_i mu(state_i)=1/7.
```

On the common integer coordinate, each numerator is
`345867129050087785140`.

## The common-base table

For every `U_full` character and endpoint guard shift, the parent one-base
integrand is split by the four source states before integration:

```text
Gamma_(alpha,beta,tau)(s)
 = integral 1_(state=s)(y) Q(13y) exp(2*pi*i*57122y)
     [sum_u U_u(y)E_u(y)] [sum_q V_q(y)E_q(y)] dy.
```

Summing `s` recovers every row of the corrected parent delta-cell current
exactly.  Five independently built literal guarded endpoint sets reproduce
the delete-and-restore guard rows.  The same-root sector remains pointwise
zero before integration.

After inversion to the thirteen `(1,0,t)` residue classes, the coordinate
ranks are

```text
coupled table       4,
source-weight erasure 4,
doubly centred      3.
```

Rank three is maximal after removing the constant state direction.  In
particular this cannot be another rank-one boundary-delta artifact.

## Full Walsh x F_13 closure

Use the intrinsic Walsh characters of the two state bits and the inherited
`F_13` relation character.  The exact support census is

```text
(total,DC,V4-axis,F13-axis,mixed)
  = (52,1,3,12,36).
```

After output ANOVA, exactly all `36=3*12` genuinely mixed modes survive:

```text
(36,0,0,0,36).
```

Both one-bit marginals remain nontrivial.  The chamber bit and multiplicity
bit each give a rank-two `2 x 13` table with full census

```text
(26,1,1,12,12).
```

Thus neither bit was merely carried by the other.

At the fixed relation class `(1,0,6)`, all four state values are nonzero:

```text
(79866518267205440406168,
 310652419985144092895775,
 367085997033220164935953,
 315468006625786970755333).
```

All four Walsh channels are also nonzero.  Their trivial Walsh sum is
`317699132065964946247468`, exactly the corrected parent fixed-class value.
Nonzero reduction in the Lucas-certified split field proves
characteristic-zero nonvanishing of each corresponding cyclotomic quantity.

## What the source-erasure says

Replacing `U_u,V_q` by unit weights while retaining the source-defined state
partition still gives rank four and full `52/52` support.  This is an honest
hostile to a causal overread: the nonseparability is not proved to be caused
by the numerical source weights.  What is proved is that the *typed support
partition* supplied by the source service is load-bearing and compatible with
the actual endpoint integrand.  Erasing the partition itself, or either
coordinate after centering, removes the corresponding Walsh axis.

## Connection contract

| field | exact answer |
|---|---|
| source | THM-2471 unsplit root-service support plus actual `U_full` owner-node integrand |
| target | four-state Boolean square times the thirteen refined relation residues |
| map | `(y,u-support) -> (left/right bit, singleton/doubleton bit)` before integration |
| preserved | common base, source roots and weights, endpoint sheets/guard, complement involution, residue phase |
| destroyed | finite ancestry sheets `a,e',b`, exact address orbit, source/arrival time, U_clock |
| hostile | MISTAKE-417 delta-cell rank one; source-weight erasure; each one-bit marginal |
| decisive positive | rank `4`, centred rank `3`, all `36` mixed Walsh-residue modes |

The object should not be called a tournament of size four.  Its four vertices
are *states of partial tournaments* on three roots.  The correct algebra on
the four states is the Boolean square; the six-object extension is a
combinatorial cut completion, while the actual unrestricted source object is
the five-vertex Gray path above.

## Remaining frontier

This closes the immediate replacement demanded by MISTAKE-417: there is a
genuine rank-greater-than-one Boolean ancestry coordinate on the same
owner-node integrand.  It does not yet provide the semantic map needed for a
row exclusion.  The next minimal tests are:

1. retain one inverse-ancestry sheet or exact-address sidecar together with
   the four states and test whether rank three and the fixed class survive;
2. determine whether the realized centre `{6}` can be transported by a later
   `U_clock` owner without changing temporal copy, and separately whether the
   missing cut `{0,12}` can ever be generated; and
3. only then ask whether the toggle palindrome `(6,0,12,6)` has a lawful
   recurrence or tree action.  No Fibonacci or ternary-tree identification is
   currently proved.

## Reproduction

```text
python -B 04-computation/lrc_r5_ufull_owner_node_boolean_square_refiner_probe_20260816.py
python -B -O 04-computation/lrc_r5_ufull_owner_node_boolean_square_refiner_probe_20260816.py
```

Normal and optimized transcripts are byte-identical.  The semantic digest is
`bae28345b0b1aea35b244bfbf04123414f0c8fbf9eeca98e39d2b94dd6d107ec`.
