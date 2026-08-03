# One normal jet is first-order; the endpoint `K4` needs mixed curvature

**Status: FINITE-EXACT PARTIAL SCOUT / NON-CANONICAL SYNTHESIS.**  The exact
companion is
[`lrc14_horn_first_normal_lift_partial_scout_20260803.py`](../04-computation/lrc14_horn_first_normal_lift_partial_scout_20260803.py),
with its frozen
[`output`](../05-knowledge/results/lrc14_horn_first_normal_lift_partial_scout_20260803.out).
It exhausts all thirteen inherited first-normal labels on the canonical
THM-3285 horn unit subatoms.  No label passes the prescribed pre-Fourier
gates, so no new endpoint-origin covariance scan is entered.  This is not a
theorem against semantic-cell reselection, a new two-object coupling, or an
arbitrary endpoint correspondence, and it has no theorem ID.

The LF-normalized hashes are

```text
script  02f6421702293cb71adf6deef3a2381da1db6c3afeb110f3ff56a05010cccf20
output  872b3ee2e6f40b469aa6660f5c56c4c3a0de88b5a86994d6ef1fe9de99889d10
```

## Inheritance pass

The closest proved positive mechanism is
[THM-3285](../01-canon/theorems/THM-3285-same-ancestry-allocation-switching-horn.md):
on the fixed rail-eight label `(clock,sigma,tau)=(1,0,3)`, the three address
cylinders are literal same-ancestry whole cylinders with allocation word

```text
R -> M -> R.
```

The corrected near miss is the frozen
[endpoint-origin scout](lrc14-endpoint-origin-k4-horn-partial-scout-codex-20260803.md).
It proves that the middle has all `169` origins in two exact fields, but the
outer typed source companions and the full-depth endpoint covariance fail.
Most importantly, `R` means loss of a fully typed source section; it is not a
source-allocation-absence bit.

The least-used positive sidecar is
[THM-2813](../01-canon/theorems/THM-2813-affine-lift-transvection-and-projective-horn-decoder.md).
Its relative lifts are

```text
A_t(y)=y+t*13^5*(y-7 mod 13),             t in F_13.
```

They fix the residue-seven sheet and recover `t` sharply on the adjacent
residue-eight sheet.  The canonical hostile is
[THM-2820](../01-canon/theorems/THM-2820-boolean-idempotent-rigidity-and-norm-top-cotangent-jet-no-go.md):
the target-only successor commutator has a nonzero normal signal, but its
common rechart is pure gauge and its nonzero raw mixed aggregate lives on
joint-absent cells.  Finally,
[THM-2772](../01-canon/theorems/THM-2772-carrier-allocation-pullback-k4-segre-and-mixed-face-obstruction.md)
says that the desired physical endpoint certificate is a **mixed** face on
one common endpoint atom, produced conditionally by two independent endpoint
translations.

The live concept board was therefore:

1. the `R-M-R` same-ancestry horn;
2. the successor/transvection commutator;
3. all six semantic gates `(E3,clock,q1,q2,c2,c3)`;
4. the raw pre-DFT `(bare,source,target,both)` square;
5. the `169`-origin endpoint sidecar; and
6. first-normal displacement versus mixed endpoint curvature.

## The lawful first-normal carrier

Write the three source addresses as

```text
n0=3454614,   n+=3454627,   na=4143978.
```

All satisfy `n=7 mod 13`.  Their successor target addresses `d(n)=n+1`
are residue eight.  Consequently the actual square is

```text
A_t(n)=n,
A_t(d(n))=d(n)+t*13^5                         mod 13^6.
```

The physical address scale turns the target displacement into

```text
q_alloc=7t mod 13,
target unit atom |-> target unit atom + q_alloc*T/13.
```

The scout constructs every shifted target atom by precisely this inherited
map and verifies that inverse translation returns the frozen atom.  Thus the
raw target-carrier copy is transported lawfully.  This does **not** reprove
literal equality of the full THM-2791 contributor set after the fixed
semantic section; that stronger ancestry statement is deliberately not
claimed off sheet.

The normal label is fully visible:

```text
t=2 q_alloc mod 13.
```

So the experiment is not failing because the first jet collapsed.  It fails
after the jet has done exactly what THM-2813 promises.

## Gate 1: the outer typed source never moves

The source atoms remain on the fixed residue-seven sheet for every `t`.
Their six-factor signatures are therefore invariant:

```text
n0: (0,1,1,1,0,1),
n+: (1,1,1,1,1,1),
na: (0,1,1,1,1,1).
```

At `n0`, `E3` and `c2` are absent.  At `na`, `E3` is absent.  Hence

```text
outer typed co-support survivors = 0/13.
```

This is the first horn-global failure.  A target-relative normal motion
cannot repair a source semantic gate on a point that the action fixes.  The
failure is stronger and more precisely typed than “the covariance scan found
nothing”: no normal label is eligible to enter such a scan.

The target section itself is complete at all three vertices only for

```text
(t,q_alloc)=(0,0),(6,3).
```

Thus the middle has exactly two local six-gate candidates.  Label zero is the
frozen control.  The only nonzero candidate is `t=6`, not `t=1`.

## Gate 2: the raw common face remains flat

For every vertex and normal label, the exact raw carrier masks are

```text
source mask = delta_0,
target mask = delta_q.
```

Before endpoint or carrier Fourier transform, the inherited Boolean square at
harmonics `(k,l)` is

```text
(B,S,T,H)
 =(1, delta_0(k), delta_q(l), delta_0(k)delta_q(l)).
```

Its unique co-support point is `(k,l)=(0,q)`.  There

```text
(B,S,T,H)=(1,1,1,1),
B-S-T+H=0.
```

The full raw mixed aggregate is nonzero on exactly the other
`12^2=144` joint-absent cells:

```text
Omega=(1-delta_0(k))(1-delta_q(l)).
```

It is disjoint from the common support.  The target-axis normal action merely
reindexes `delta_0` to `delta_q`; it cannot change this factorization.  The
scout checks all

```text
13 labels * 3 vertices * 169 harmonic pairs = 6,591
```

raw rows, and retains their cross-product with all `169` endpoint-origin
labels.  Flatness at co-support is uniform before DFT, so no origin or field
choice can repair it at this stage.

This is the first strong **local** failure.  Even the middle-complete labels
`t=0,6` have zero raw mixed face.

## Strongest survivor: `t=6`, `q_alloc=3`

The nonzero label `t=6` is informative:

```text
normal decoder:              t=2*3=6,
middle source gates:         111111,
middle target gates:         111111,
raw carrier co-support:      (0,3),
raw K4 there:                (1,1,1,1),
raw mixed face there:        0.
```

This recovers the factor-side positive signal isolated by THM-2820: `q=3`
is the unique nonzero complete-section successor displacement in its fixed
cell.  It does **not** recover THM-2820's missing common raw face.  It also
does not inherit that theorem's endpoint statement, because THM-2820 used a
different `(sigma,tau)` cell.

The distinction is sharp.  In THM-2820's external control cell,

```text
q=3  keeps the complete fixed section but loses the translated endpoint companion;
q=11 keeps the two-point endpoint edge but loses q2.
```

In the present canonical `tau=3` horn, `q=11` is weaker still: at the middle
it loses both `E3` and `q2`.  Replacing the semantic cell to improve that
signature would be a new typed branch.  It cannot be silently folded into
the THM-3285 horn because carrier weight, literal ancestry, address incidence,
and both endpoint origins would all need fresh verification.

## Named hostile controls

### Normal label zero

At `t=0`, the atoms, masks, and six-factor signatures exactly reproduce the
pre-DFT frozen scout.  Its pinned downstream computation remains the positive
control:

```text
middle endpoint support: 169/169 in both fields,
raw common twist:         flat,
vertical affine covariance: absent.
```

The new gate order does not rerun those expensive rows because the raw face
has already failed.

### `tau=12`

The independently audited target-only simplex is reconstructed before any
normal label:

```text
(A,B,M,L,R) piece counts=(240,241,0,240,241).
```

Thus `M` is empty before address restriction.  A normal target shift cannot
turn this hostile into a common atom.

### THM-2820 exclusive cells

The pinned `q=3` versus `q=11` either-or is used only as a comparison.  It
confirms that a complete factor section and a two-origin endpoint edge are
distinct resources.  No endpoint value from that different cell is imported
into the canonical horn.

## Why one normal jet cannot be the missing map

There are two coordinates both informally called `q`, and the computation
shows why they must not be identified:

```text
q_alloc = target-axis carrier/Hasse displacement (0,7t),
q_E     = endpoint difference L-R in THM-2772's endpoint plane.
```

THM-2813 and THM-2820 construct `q_alloc`.  The desired pullback needs

```text
(L,R,q_E,Delta),             q_E=L-R,   Delta=det(L,R).
```

No inherited map sends the transported raw carrier copy to such an endpoint
pair.  More conceptually, the normal jet is a one-axis **first difference**.
THM-2772's physical escape is a two-axis **mixed second difference**:

```text
F00-F10-F01+F11=det(s,t).
```

Translating only the target axis can move the location of a delta mask, but
cannot manufacture the independent source increment needed for nonzero
`det(s,t)`.  The exact raw factorization above is the finite manifestation of
that derivative-order mismatch.

## Connection contract and next test

```text
source:
  the canonical THM-3285 integral unit subatoms and all thirteen THM-2813
  relative lifts;

target tested:
  the raw Boolean K4 crossed with the 169 endpoint-origin labels, before DFT;

map:
  source address n fixed, successor target d(n) moved by A_t, then exact
  physical target shift q_alloc=7t;

preserved:
  normal label, source address, successor relation, transported raw carrier
  copy, clock, sigma, tau, all six separately tested semantic gates, open
  endpoints, and the full endpoint-origin label universe;

destroyed or absent:
  outer source E3 (and n0 c2), a nonflat raw face on common support, a literal
  off-sheet full-section ancestry audit, endpoint difference q_E, endpoint
  pair (L,R), determinant Delta, and vertical covariance;

verdict:
  one normal jet faithfully decodes t but does not supply the missing typed
  ancestry-to-(L,R,q_E,Delta) map.
```

The cheapest genuinely new positive test is no longer another origin scan or
another single normal label.  It is one of two explicitly separate branches:

1. construct a lawful **two-axis** physical square on one common raw atom,
   with independent source and target endpoint increments and all six gates;
   test the mixed face before DFT; or
2. reselect a semantic cell that repairs the factor/endpoint `q=3` versus
   `q=11` tradeoff, then independently rebuild its carrier weight, ancestry,
   address horn, and endpoint pair before comparing it with THM-3285.

The first branch is mathematically cleaner: it targets the missing derivative
order directly.  But canon currently supplies no lawful source-normal
connection paired with the target successor commutator.  Until such a
two-cell exists, assigning `q_alloc` to `q_E` would be exactly the unproved
syntax bridge this experiment was designed to avoid.

No root/Čech map, current cycle, relation-row exclusion, or `LRC(14)`
conclusion follows.
