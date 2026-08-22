# Six pointed cuts form a full-rank ancestry path over the owner square

**Status: FINITE-EXACT CANDIDATE, INDEPENDENT AUDIT PENDING.  LRC(14) remains
OPEN.**  Marking the active source tail root refines the four owner-visible
cut states to exactly six realized pointed states.  Their tensor with source
root difference and endpoint relation has maximal ranks `(6,12,13)` and
maximal three-way ANOVA ranks `(5,12,12)`.  An equal-tail lift with the same
four-state marginal has ranks only `(4,12,13)` and `(3,12,12)`, proving that
the two added tail directions are genuine.

## Inheritance and the six-object correction

The independently audited Boolean square has states

```text
0:{0},   1:{0,6},   2:{12},   3:{6,12}.
```

The root-difference refinement then keeps the colour `s=u-q` of every directed
cut arc.  It still sums the absolute tail root `u`.  The independent
connection obstruction identified this forgotten absolute root as part of
the smallest surviving source record.

Marking `u` produces exactly six owner-visible pairs, in geometric order:

```text
(0,0) -> (1,0) -> (1,6) -> (3,6) -> (3,12) -> (2,12).
```

This is the natural size-six carrier.  It is not the abstract family of six
nontrivial cuts of `{0,6,12}`: one of those cuts, `{0,12}`, never occurs.
Nor is it a tournament on six vertices.  Each vertex is a **pointed directed
cut state** on the thirteen-root fibre.

The five path edges have exact meanings:

```text
enter the left doubleton while retaining tail 0,
switch marked tail 0 -> 6 inside that doubleton,
cross the owner-inactive centre while retaining tail 6,
switch marked tail 6 -> 12 inside the right doubleton,
leave the right doubleton while retaining tail 12.
```

The middle edge is an owner-gap bridge, not a direct Gray edge: the omitted
physical support state `{6}` lies between its endpoints.

## Construction and exact marginal

For each actual endpoint-selected source pair `(u,q)`, the integrand retains

```text
(pointed_state=(state,u), difference=u-q)
```

before source weighting, endpoint multiplication, word jump, harmonic phase,
and integration.  Inverting the full endpoint character bank gives

```text
P(pointed_state,s,t),

pointed_state in the six-state path,
s in F_13,
t in F_13 for the endpoint relation (1,0,t).
```

Summing pointed states over each four-state fibre recovers the previously
pinned `V_4 x F_13 x F_13` tensor exactly, both for the weighted and
support-only banks.  The same-root slice `s=0` remains pointwise zero.  Three
literal guarded controls agree with independently rebuilt rows.

The segment-count vector is palindromic:

```text
(41935,51187,51187,51187,51187,41935).
```

These are integrator segment counts, not probability masses; their role is an
exact census and reflection control.

## The marked tail adds exactly two directions

The raw flattening ranks `(pointed,difference,relation)` are

```text
weighted       (6,12,13),
support-only   (6,12,13),
flat-tail      (4,12,13).
```

Here `flat-tail` divides each four-state parent entry equally among its one or
two allowed marked roots.  It has exactly the same parent marginal.  After
centering all three axes, the ranks are

```text
weighted       (5,12,12),
support-only   (5,12,12),
flat-tail      (3,12,12).
```

Thus the actual and support-only pointed carriers are maximal in every axis,
while forgetting which tail in a doubleton carries the arc loses exactly two
state directions.  Numerical source weights again are not needed for the rank
claim: the typed support incidence is load-bearing.

At fixed relation `(1,0,6)`, the weighted and support-only `6 x 13` matrices
have rank six and centred rank five.  Each of the two within-doubleton tail
contrasts is nonzero at all twelve nonzero differences.  Both contrasts vanish
identically in the flat-tail hostile.

## The five path-edge spectra

There is no canonical `C_6` group law on the six states.  A six-point DFT would
therefore be cosmetic.  Instead, use the geometrically ordered path incidence
map, whose five edge contrasts span the entire nonconstant state space, and
Fourier-transform only the two genuine `F_13` coordinates after centering
them.

For both weighted and support-only banks, the mixed support counts are

```text
(120,144,120,144,120).
```

The missing modes are complete relation-frequency lines:

```text
edge 0, enter doubleton:        t_hat = +/-1,
edge 1, switch tail 0 -> 6:     none,
edge 2, owner-gap bridge:       t_hat = +/-6,
edge 3, switch tail 6 -> 12:    none,
edge 4, leave doubleton:        t_hat = +/-1.
```

Every marked-tail switch mode survives (`144/144`).  The flat-tail hostile
kills those two edge spectra completely, while its three state-transition
edges have `144/144`.  This separates two kinds of transport that the
four-state marginal conflates:

```text
state transport across support chambers,
tail transport inside a fixed doubleton.
```

The palindromes in path order, segment counts, edge support, and missing lines
are exact reflection structure.  They are not yet a recurrence theorem.

## Tournament, tree, and recurrence boundary

The pointed path gives a lawful interpretation of the user's sizes four and
six:

```text
four = unpointed owner-visible cut states,
six  = their realized marked-tail fibres.
```

The source arcs are directed from `S` to `F_13 minus S`; within-side pairs are
missing, and there are no both-way arcs.  Root difference is the translation-
invariant arc colour.  The six path vertices are states of these partial
tournaments, not contestants in another tournament.

The edge-type word

```text
state, tail, gap, tail, state
```

and spectral word

```text
120,144,120,144,120
```

are finite palindromes.  Calling either Fibonacci, ternary, or substitutive
would require a second scale and a proved transition operator carrying the
same typed records.  The next recurrence test should transport this pointed
carrier to `U_clock` or another lawful owner depth and compare edge operators,
not merely compare counts.

Encoding these six states as a subset of the natural numbers, or as terms in
the harmonic series, is always possible but loses the cyclic difference and
owner geometry unless a preservation theorem is supplied.  The useful
harmonic object here is the finite-character response of arc-coloured subset
indicators.

## Connection contract

| field | exact content |
|---|---|
| source | actual THM-2471 source roots and actual THM-3514 endpoint indicators on one owner base |
| target | six pointed cuts times source difference times endpoint relation residue |
| map | `(state,u,q) -> ((state,u),u-q)` before multiplication/integration |
| preserved | absolute source tail, cut state, root difference, source weights/support, endpoint sheets/guards, harmonic phase |
| destroyed | deep label `theta`, inverse ancestry sheets, horizons, exact THM-2334 address, chronology |
| hostile | exact four-state marginal, support-only weights, equal-tail lift, same-root zero, literal guards |
| decisive positive | rank gain `4->6`, centred gain `3->5`, two full `144/144` tail-switch edge spectra |

This closes the cheap absolute-source-root refinement on one host.  It does
not yet retain the full `(u,s,ell,theta;a,d,C,D)` relation isolated by the
connection obstruction.  The next serious sidecar is an independently
defined deep/horizon or inverse-ancestry label on these same pointed records.
Only after that should one attempt an `ell`-independent THM-2334 address map.

No grouped `C(a;X,m)`, `U_clock` chronology, physical current, uniform-row
statement, row exclusion, or LRC(14) conclusion follows.

## Reproduction

```text
python -B 04-computation/lrc_r5_ufull_owner_node_pointed_six_state_root_difference_probe_20260816.py
python -B -O 04-computation/lrc_r5_ufull_owner_node_pointed_six_state_root_difference_probe_20260816.py
```

Normal, optimized, and stored transcripts are byte-identical.  Script,
output, and semantic LF SHA-256 are

```text
52b5af635b394e8f6dda59746d369b4b62a73da5ee89c6ca8e758426a7e81b76
40ca4ce5e6e4428c833e41e412e8acc58a1a1bb25bbfb542bc43a7319fe68e7e
5001bed534a9b8f953529101f0b7a51cf6994c7dcb29b7aed576aec239078384
```
