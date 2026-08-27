# Sun parity, odd cycles, and the `W=0` planar-JC attachment frontier

Scratch-only research report, 2026-08-26.  Nothing here has been promoted to
maintained canon.  Dependencies are current THM-4230 and corrected THM-4241.

## Inheritance and concept board

- Closest proved mechanism: THM-4230 section 6.2 confines degree `34/42`
  collapse to finite marked-ratio sets; THM-4241 gives the full integral Hom
  lattice and its unique `F_4` glue line.
- Canonical hostile: THM-4230 `(U,W,Z)=(2,-2,1)` shows that attachment equality
  need not imply constancy.
- Corrected near miss: the `M=10,11,12` degeneration dual graphs have large
  cycle rank but no odd cycle.  Graph-cycle parity is not the obstruction.
- Least-used sidecar: the free involution `iota=tau^6` on the twelve labelled
  nodes and the order-three action `sigma=tau^8`.

Live concepts were: Sun's free/fixed involution parity, the four attachment
`C_3` orbits, the THM-4241 `F_4` glue line, the bivariate projection theta
series, and the finite marked-ratio sets.

## 1. Rigorous no-go: degeneration odd cycles are absent

After contracting rational chains, the exact dual-graph skeletons are:

- `M=10` (THM-4218): vertices `R,C,E`, ten parallel `R--C` paths and two
  parallel `C--E` paths; `b_1=10`.
- `M=11` (THM-4222/4232): vertices `R,C,T,V`, eleven parallel `R--C` paths,
  one `C--T` path and one `C--V` path; `b_1=10`.
- `M=12` (THM-4230): vertices `R,C`, twelve parallel `R--C` paths;
  `b_1=11`.

Every skeleton is bipartite, but that observation alone is not enough:
arbitrary edge subdivision need not preserve bipartiteness.  Here the exact
local charts supply the missing parity.  The repeated `R--C` paths have
lengths `30` at `M=10` (`A_29`), `330` in dense `M=11` (`A_329`), `132` on
the `U=0` `M=11` wall (`A_131`), and `12` at `M=12` (`A_11`).  Every repeated
path length is even, so a cycle formed by two such paths has even length.
The other `M=10` cycle is the pair of direct `C--E` attachment edges and has
length two; the `M=11` side paths are bridges.  Hence every graph cycle in
the full resolved dual graph is even.  Nevertheless the cycle ranks are
`10,10,11`.  Thus:

> **PROVED no-go.** Neither positive dual-graph genus nor the current
> `M=10,11,12` planar-JC obstruction can be detected by the presence of an
> odd cycle in the degeneration dual graph.

The hostile is `M=12`: it has the largest current graph contribution
`b_1=11` while remaining bipartite.  Vertex genus, labelled attachments,
component multiplicity, and Hom/monodromy are the indispensable sidecars.

## 2. The lawful odd cycles are attachment-action cycles

On `C_0:x^6+y^4=1`, the twelve contacts form one `tau` orbit.  The action

```text
sigma=tau^8
```

partitions them into four disjoint directed `3`-cycles, while
`iota=tau^6` pairs them freely into six antipodal pairs.  Every hidden map
`ell` satisfies

```text
ell o sigma=[omega]ell,       ell o iota=[-1]ell.
```

If a hidden map has one common value `P` on a full attachment orbit, the
`C_3` action gives `(1-omega)P=0`, while the free involution gives `2P=0`.
Since the two kernels have coprime orders, `P=O`.  This is the exact
odd-cycle/torsion mechanism already implicit in THM-4230.

Connection contract:

```text
source:      four directed C3 orbits plus six free C2 pairs on {Q_j}
target:      E_0 torsion
map:         Q_j -> ell(Q_j)
preserved:   sigma eigenvalue omega and iota sign -1
destroyed:   mixed visible-hidden half-coordinate and individual node values
sidecar:     M/(V+L)=F_4 and exact attachment denominators
hostile:     a C3 alone permits nonzero E_0[1-omega]; the C2 sidecar is essential
```

This is genuinely analogous to THM-4037's Sun parity ladder: an involution
acts freely on the relevant finite set, and the target contribution is forced
into a small fixed/torsion fibre.  It is not a transfer from tournament OCF;
the intrinsic graph here is the permutation graph of `sigma`, not a
manufactured tournament.

## 3. Exact bivariate projection theta census

For `m in M=Hom(J(C_0),E_0)`, put

```text
v=m+m o iota,          ell=m-m o iota.
```

THM-4241 gives `v in V`, `ell in L`, and orthogonality gives

```text
q(v)+q(ell)=4q(m).                                      (1)
```

Using the normalized THM-4241 basis `[u,f,g,h]`, where
`2h=v_0+(omega^2f+g)`, the exact coordinates are

```text
ell=(2b+omega^2d)f+(2c+d)g,
q(v)=16N(a)+4N(d).
```

The scratch certificate independently reproduces the full theta counts
`36,288` at degree `34` and `16,992` at degree `42`, and refines them as:

```text
d=34, q(ell):count
12:1536  24:2304  36:5952  60:5760  72:8064
84:5376  108:4992  120:1728  132:576

d=42, q(ell):count
12:672  24:288  36:2304  60:288  72:3456  84:3072
108:3744  120:1728  132:576  156:672  168:192
```

Every hidden projection degree is a positive multiple of `12`.  This is the
exact analogue of Sun's role-labelled multiplicity refinement: total degree
alone is replaced by the indexed fibre `(q(v),q(ell))`.  It preserves the
degree identity `(1)` but loses node evaluations, so it schedules resultants
rather than proving collapse impossible.

If `m` collapses all twelve attachments to `P`, choosing an `iota`-fixed base
point gives the necessary conditions

```text
ell(Q_j)=O for every j,       v(Q_j)=2P for every j.   (2)
```

## 4. New exact boundary exclusion: `q(ell)=12`

The hidden degree-`12` shell has exactly `24` vectors.  Postcomposition by
`mu_6` and precomposition by `T` split it into exactly two free orbits of
size `12`.

For a degree-12 hidden map, THM-4230's pole ledger gives

```text
6 deg(d_ell)-2 <= 24,
```

and oddness of the reduced `X/x` denominator sharpens this to
`deg(d_ell)<=3`.

For one orbit, direct characteristic-zero group law gives a nonzero scalar
times `t(t^2-1)`; its scalar has absolute resultant `65536`, and exact
numerator resultants at `t=0,+1,-1` are nonzero.  For the other orbit, the
independent symbolic probe gives a raw denominator containing
`t(t^2-1)` and nonzero numerator norms

```text
t=0:   4,477,456,
t=±1:  40,290,721,869,103,654,477,234,176.
```

Therefore the reduced denominator contains all three factors.  The analytic
upper bound forces equality:

```text
d_ell(t)=nonzero scalar * t(t^2-1)                    (3)
```

for both orbit representatives, hence for all `24` degree-12 maps up to the
attachment-preserving symmetries.  Four good reductions (`313,349,373,397`)
independently give reduced degree `3` and reciprocal gcd exactly `t^2-1` for
both orbits.

The simultaneous-pole/collapse test compares `d(t)` with
`t^deg(d)d(-1/t)`.  By `(3)` their only common roots are `t=±1`.  But

```text
Z/U=((t^2-1)/(2t))^2,
```

so `t=±1` is the excluded wall `Z=0`.  Hence:

> **PROVED relative to the THM-4230 attachment/pole interface.** No hidden
> degree-12 map collapses an admissible gate-interior attachment orbit.

Combining with `(2)` removes all `q(ell)=12` full vectors: `1,536` at degree
`34` and `672` at degree `42`.

There are two immediate exact fibre-degree deletions:

- at degree `34`, the `(q(v),q(ell))=(4,132)` row has `576` vectors; a
  nonconstant degree-four visible map cannot put twelve distinct points in
  one fibre;
- at degree `42`, the `(0,168)` row has `192` vectors and is precisely the
  pure-hidden degree-42 shell already excluded by THM-4230.

Thus the raw vector frontier decreases from

```text
degree 34: 36,288 -> 34,176,
degree 42: 16,992 -> 16,128.                           (4)
```

These are vector counts before source/target symmetry quotient.  Equation
`(4)` does not enumerate or empty `S_34,S_42`.

## 5. Sun multiplicity method: what transfers and what does not

THM-4036's useful move is to keep a role-labelled diagonal fibre and peel an
outer coordinate only after separating exact equality from off-diagonal
difference multiplicity.  The corresponding JC object is the bivariate
theta fibre `(q(v),q(ell))`, not the univariate degree theta series.

```text
Sun source: role-labelled additive fibre
JC target:  involution-labelled orthogonal degree fibre
map:        tuple -> (outer atom, tail sum)
            m -> (q(v),q(ell))
preserved:  exact total and indexed multiplicity
lost:       Sun height/sign; JC attachment image and marked ratio
sidecar:    square avoidance in Sun; exact node resultants in JC
hostile:    high multiplicity/local support can coexist with a Sun hole;
            a represented JC degree can coexist with failed attachment collapse
```

The next useful object is therefore a trivariate finite table

```text
(q(v),q(ell), attachment-resultant orbit)
```

rather than another scalar theta coefficient.

## 6. Precise next computations

1. Enumerate `mu_6 x <tau>` orbits separately in every surviving projection
   row of `(4)`, beginning with the sharp fibres `q(v)=12` and `q(ell)=24`.
2. Reuse the denominator-degree budget with map degree `q(ell)`; record the
   reciprocal-resultant polynomial in the marked variable
   `R=4t^2/(t^2-1)^2`, not only whether a finite-field gcd appears.
3. On the visible side, compute the analogous exact fibre polynomial for
   `v(Q_j)=v(Q_0)`.  Intersect its marked-ratio roots with the hidden roots.
4. Retain the `F_4` glue class when combining the two resultants.  Conditions
   on `2m=v+ell` alone can leave a nonzero two-torsion difference.
5. Use odd `C_3` attachment orbits only through the proved torsion map above;
   do not use odd cycles in the bipartite degeneration dual graph.

## Reproduction

```text
python -B .scratch/sun_jc_oddcycles/jc_w0_projection_histogram.py
python -B .scratch/sun_jc_oddcycles/jc_w0_hidden_degree12_attachment_audit.py
python -B -O .scratch/sun_jc_oddcycles/jc_w0_hidden_degree12_attachment_audit.py
python -B .scratch/sun_jc_oddcycles/symbolic_degree12_orbit0_probe.py
```

The normal and optimized degree-12 audit outputs byte-match.
