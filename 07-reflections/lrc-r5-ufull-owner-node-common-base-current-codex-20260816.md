# The actual source service and U_full endpoint factors now coexist on one owner-node base

**Status: FINITE-EXACT CANDIDATE, NOT A THEOREM.**  This experiment moves the
actual `U_full` endpoint Boolean factors inside the `r=5` source-service
integral.  Its role-residue signal is nonzero in one certified split-field
image.  The apparent complete `F_7 x F_13` spectrum is, however, a rank-one
`delta_(ell=0)` lift: the endpoint owner factor forces cell zero, so this is
not genuine cell/residue mixing (MISTAKE-417).  It does not isolate an exact
address `C(a;X,m)`, identify arrival and source time, transplant `U_clock`,
exclude a row, or prove LRC(14).

## Inheritance and the changed object

The closest proved mechanism is THM-2471, equations (25), (28), and (31):

```text
U_u(y) = (P_(13^5) f)((y+u)/13),
V_q(y) = (P_(13^5) e)((y+q)/13),
w_u    = (y+u)/13.
```

THM-2594 already puts its cell, deep, and word factors inside this same
Boolean ancestry integral.  The corrected near miss was the source-aligned
seven-cell package: it has maximal `7 x 13` spectral capacity, independently
audited at `d2b4c10d0`, but its endpoint `AX/BY` values are scalars integrated
before the source contraction.  The least-used relevant sidecar is the
desheeting identity from the actual `U_full` endpoint engine,

```text
91t = 7a+r,
y   = r/7 = 13t-a,
t   = (y+a)/13 = w_a.
```

The new object uses that identity literally.  For each endpoint character
and guard shift it forms, before integration,

```text
Q(13y) exp(2*pi*i*57122*y)
 [sum_u U_u(y) E_u(y)] [sum_q V_q(y) E_q(y)] cell_ell(y).
```

Here `E_u(y)` is the actual shifted `U_full` Boolean endpoint factor at the
owner node `w_u`; it is not an atom weight or an `AX/BY` marginal.  The
source and endpoint factors therefore share one displayed base and are
multiplied before integration and Fourier inversion.  Information retained:
the collision-root indices, source-service weights, endpoint owner sheets,
endpoint character, guard shift, and seven-cell coordinate.  Information
still destroyed: the `a,e',b` ancestry sheets, exact address orbit, prescribed
clock comparison with `U_clock`, and source/arrival temporal copy.

## One exact joint coordinate and field

The inherited grids are incompatible:

```text
T_source   = 297836897838480,
T_endpoint = 483730250419703196.
```

Using either old field would silently project away rational breakpoints.  The
probe instead uses

```text
C = lcm(T_source,13*T_endpoint)
  = 9684279613402457983920,

L = lcm(T_source,13^2*T_endpoint)
  = 125895634974231953790960.
```

The prime

```text
p = 6L+1 = 755373809845391722745761
```

has the fully listed factorization of `p-1`.  The Lucas checks certify `23`
as a primitive generator, and `148035889=23^6` has exact order `L`.  Thus the
source endpoints, endpoint inverse branches, `F_7` and `F_13` characters, and
the oscillatory factor all live in one lawful split-field image.  The endpoint
frequency descends as

```text
742586 = 13^2*4394,
exp(2*pi*i*4394*(13y-j)) = exp(2*pi*i*57122*y),
```

because the discarded branch phase `4394j` is integral.

## Exact result and the cell-collapse hostile

The complete coupled character bank has digest

```text
b5246eb2a69f35e4dac7dabbf26b1703f21ed22bf803061399ebbf766b9a073d.
```

After inverse transformation to the thirteen `(1,0,t)` residue classes, the
formal seven-cell spectrum has support census

```text
(total,DC,F7-axis,F13-axis,mixed) = (91,1,6,12,72).
```

and double-centering reports all `72` mixed modes nonzero:

```text
(72,0,0,0,72).
```

Those census numbers do **not** certify a two-coordinate interaction.  The
endpoint set contains the owner factor `||13t||<1/14`.  On every desheeted
branch `13t=y+u`, hence

```text
OWNER(t)=1  =>  ||y||<1/14  =>  ell=0.
```

The exact character-bank support by cell is therefore

```text
(2197,0,0,0,0,0,0),
```

and both the inverse table and its ANOVA matrix have rank one.  If `R(t)` is
the nonzero thirteen-entry residue profile, the exact factorization is

```text
table(ell,t) = delta_0(ell) R(t),

ANOVA(ell,t)
 = (delta_0(ell)-1/7)(R(t)-mean(R)).
```

Thus the seven nonzero cell Fourier modes are repeated copies of the same
one-dimensional residue coefficient; output centering manufactures a full
mixed support census from a separable rank-one tensor.  At the fixed relation
class `(1,0,6)`, all seven `F_7` Fourier reductions are consequently the same
nonzero field value

```text
317699132065964946247468.
```

The six zero cell values are geometric characteristic-zero zeros, not merely
zeros modulo the split prime.  The nonzero residue and bridge reductions still
prove characteristic-zero nonvanishing of those one-dimensional quantities.

The two role values and their bridge are

```text
q_H  = 125385278409587426725290,
q_Q5 = 657486478079327229022863,
q_H-q_Q5 = 223272610175651920448188 != 0 mod p.
```

This differs from the source-erasure control, so the source service changes
the current rather than merely rescaling the endpoint-only answer.

## Controls and why they matter

Five literal guarded endpoint sets independently reproduce the rows obtained
by deleting and then restoring the guard through chamber safety.  Normal and
optimized transcripts are byte-identical.

The source profiles independently recover

```text
total numerator = 168908463464745122312762880,
total / DENC     = 169 I_5.
```

More sharply, every same-root product satisfies

```text
U_u(y)V_u(y)=0
```

on every exact source segment.  Hence the entire same-root sector vanishes
pointwise before endpoint restriction, oscillation, integration, or reduction
modulo `p`.  The positive signal is wholly off diagonal in the collision
roots.

Replacing both source profiles by `1` gives a source-erasure bank.  It has the
same delta-cell factorization and independently restores the literal endpoint
guards, so it is a positive engine control rather than an attempted causal or
mixing proof.  Its fixed-relation geometry matches the previous desheeted
endpoint calculation, now evaluated in the larger joint field.

## The tournament clue becomes an exact directed-cut statement

The raw thirteen-root support is not a tournament.  On all 33 source
segments it has exactly one of two forms:

```text
|U|  |V|  missing unordered pairs  one-way pairs  two-way pairs
 1    12             66                  12              0
 2    11             56                  22              0
```

In fact `supp(V)` is the complement of `supp(U)`, and every supported arc is
oriented from `supp(U)` to `supp(V)`.  Thus each local support graph is a
complete directed cut `K_(1,12)` or `K_(2,11)`, with all within-side edges
missing and no bidirectional edge.  Collapsing the two sides gives one
oriented edge; it does not give a four- or six-vertex tournament.  Any genuine
size-4/6 tournament carrier must therefore come from a further typed quotient
such as endpoint chambers, Walsh arms, cut modes, or proof obligations.  This
is a useful rejection criterion for cosmetic tournament encodings.

## Honest consequence and next obstruction

The computation advances one part of the previous stopping boundary: source
service and actual `U_full` endpoint factors now coexist as an atomwise
product on one owner-node base, and their thirteen-class role bridge is
nonzero.  It is therefore stronger than a contraction against preintegrated
endpoint atom weights.  It does **not** advance genuine seven-cell spectral
closure, because the owner factor collapses that coordinate before the source
weights can interact with it.

It is still only a candidate current for four precise reasons:

1. it uses the unsplit whole source packets, not an isolated exact address
   `C(a;X,m)`;
2. the source/arrival temporal distinction of THM-2471 (33)--(37) is not
   removed by evaluating both factors at an owner node;
3. `U_full` has not been intertwined with the required `U_clock` packet at a
   common prescribed time; and
4. the computation is for the canonical `r=5` host, not uniformly for all
   165 covering rows.

The cheapest decisive continuation has two independent gates.  First, replace
the owner-subordinate seven-cell label by a refiner not forced by `OWNER`, or
retain a deeper/word sheet whose support genuinely crosses cells, and demand
matrix rank at least two before reading a mixed spectrum.  Second, retain one
exact-address sidecar (or the inverse ancestry sheet that determines it) and
test whether the nonzero `(1,0,6)` residue signal survives.  Separately, test
whether `U_clock` descends to the same owner-node coordinate without changing
the source-time copy.

## Reproduction

```text
python -B 04-computation/lrc_r5_ufull_owner_node_common_base_current_probe_20260816.py
python -B -O 04-computation/lrc_r5_ufull_owner_node_common_base_current_probe_20260816.py
```

The corrected semantic digest is
`98a27d4540648377c544d8e1b86c3dd3df7bb16d3431f6a9471f4844ba2e6b9f`.
