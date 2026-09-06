# Independent audit of the larger-unit-component LRC closure

**Status: PASS, independent proof and computational audit.** The complete
[producer report](overnight15_20260906_lrc_larger_unit.md),
including its native refinement (7), has been reviewed. No repair is
requested. Promotion to **PROVED RELATIVE TO THE STATED PROVED AND CITED
SUPPLIERS + INDEPENDENTLY AUDITED** is justified. The theorem's lower-
dimensional LRC supplier is CITED, not re-proved by these controls.

Frozen producer source SHA256:
`23f42086c09e367bf9ddc87958587901b338a8a27d712e34e6eaea64ba8ce935`.
Frozen producer output SHA256:
`3b8c091e3e387ce0cacc24e22690bd78ecc9c3d313fffd4fcc0f4f1127106b36`.
No producer or repository file was edited.

## 1. Types and cited inputs

The physical row is `tV union gU`, with thirteen distinct positive integers,
primitive overall and of sum at most Q^2, Q=91^6. The actual all-scale decoder
graph has two components of sizes a,b, where a+b=13 and a<=b. The primitive
shapes V and U have respective gcd one; U contains one and has maximum K.
The equality W=V_dec is required. It is not inferred from rank eleven alone
or from a proposed component split.

Because V and U are primitive and the whole row is primitive,
`gcd(t,g)=1`. The internal pair in the physical larger component is
`g,gK`, and its primitive height is K. The audited general decoder's
physical-box necessity argument therefore gives K<=Q. Since b>=7 and
the speeds are distinct, K>1, so these are separate coordinate columns.

The current
[joint-shadow theorem](lrc14_joint_shadow_empty_core_next_sep06.md)
was read directly. It is proved relative to cited lower-runner LRC and
independently audited. In a primitive thirteen-speed row with no weak
1/14-safe phase, every seven-speed subset has gcd at most 90. Since g
divides the gcd of any seven speeds in its larger component, hypothetical
failure forces g<=90. This implication does not discard common phase labels
or infer a phase-uniform lifting result at clock 90. A bound g<=Q would
already suffice for the arithmetic step; the inherited bound is stronger.

The primary preprint
[Sungkawichai--Trakulthongchai, arXiv:2604.23906v2](https://arxiv.org/pdf/2604.23906v2)
was opened and its Theorem 1.3 on printed page two checked directly. It
states LRC(k) for k<=12 nonzero speeds, with clearance 1/(k+1). Thus it
supplies both proper component phases used below, including the singleton
case and the b=12 case. This is the precise CITED input; the theorem does
not supply LRC for thirteen nonzero speeds.

## 2. A direct crossing construction verifies the scale inequality

Let v=min V and delta=gcd(g,v). Coprimality of t and g gives
`gcd(g,tv)=delta`. Set

```
c=g/delta,       x=tv/delta.
```

Under hypothetical failure, c<=g<=90<Q. Suppose that x<=Q(K+1). A direct
positive-box construction avoids needing the full minimal-coefficient lemma:

```
s=min(floor(x/K),Q),       r=x-sK.
```

If floor(x/K)<=Q then 0<=r<K<=Q. Otherwise s=Q and the assumed upper
bound on x gives 0<=r<=Q. Thus r,s are integers in [0,Q], and

```
c*(tv)-r*g-s*(gK)=0.
```

This is a literal support-at-most-three relation of height at most Q.
Its weighted partial sum on the smaller component is ctv>0, so it lies
outside V_dec, contradicting W=V_dec. No coordinate has been rescaled
without accounting for its coefficient. Therefore

```
tv/delta > Q(K+1),
t > delta Q(K+1)/v >= Q(K+1)/v.                     (1)
```

This strict inequality is stronger than the advertised t>QK/v and is
particularly useful at the endpoint of the proposed minimum-speed bound.

## 3. The actual phase grid and its endpoint margin

The cited b-speed theorem supplies a phase theta with U clearance at least
1/(b+1). Since u<=K for every u in U and distance to the nearest integer
is one-Lipschitz, the closed circular arc of radius

```
R=(1/(b+1)-1/14)/K
```

around theta is weakly 1/14-safe for U. Its full length is

```
ell=2R=a/[7(b+1)K].                                (2)
```

Its interior is strictly 1/14-safe. The factor two in (2) is essential;
the statement uses a full arc length, not its radius.

Fix a smaller-core phase eta with V clearance at least 1/(a+1)>1/14.
The physical phases

```
tau_j=(eta+j)/t,       j=0,...,t-1
```

preserve this smaller-core phase exactly. Their larger-core clocks are
`g eta/t + gj/t modulo one`. Since gcd(g,t)=1, these are a complete
translated t-grid. They are not merely t samples with potentially repeated
residues. Every circular point is within 1/(2t) of such a grid.

If

```
v <= floor(Qa/[7(b+1)]),                            (3)
```

then (1)--(2) give

```
t*ell > (K+1)/K > 1.
```

The closest grid point to theta lies strictly inside the safe arc. At its
physical phase tau_j the smaller component keeps its strict margin and
the larger component also has strict 1/14 clearance. This contradicts
hypothetical absence of a weak-safe full-row phase.

Thus the weak-safe conclusion is valid under (3). The conditional arithmetic
branch actually produces a strict witness. The overall theorem should retain
weak-safe scope unless the inherited g>90 exclusion is separately known to
give strict safety. No parity assumption or rounding of an odd grid size is
used. A closed threshold in (3) is harmless because the preceding scale
inequality is strict.

The producer's stronger native refinement is also valid. Without assuming
hypothetical failure or invoking the joint-gcd ceiling, the conditions

```
g/delta<=Q,          7(b+1)Kv<=a delta Q(K+1)
```

combine with the direct crossing inequality (1) to give t*ell>1. Therefore
they yield a strict full-row witness directly. Equality in the second
condition is allowed because (1) is strict. This retains the actual gcd
delta and the extra K+1, and can be stronger than the uniform minimum-speed
criterion in the remaining 6+7 branch.

## 4. The spanning-tree bound and singleton boundary

For a connected primitive a-coordinate decoder shape, choose a spanning tree
and root it. Write the ratio along each oriented edge as b_e/a_e in lowest
terms. Each numerator and denominator is at most 355. Multiply the root
coordinate by the product of all tree-edge denominators. At another vertex,
replace the denominator factor by its numerator on each edge of the root
path. This gives a positive integer realization of the entire tree ratios,
with exactly a-1 factors, each at most 355, in every coordinate.

The primitive shape divides this realization by its common gcd. Consequently

```
max V <= 355^(a-1),     hence min V <= 355^(a-1).      (4)
```

For a=1 the primitive singleton is {1}, and the empty product gives one.
The actual singleton still has physical speed t, which is retained in the
crossing and phase formulas. Its scale has not been inferred from a lossy
packet invariant.

The exact threshold comparison is:

| a+b | floor(Qa/[7(b+1)]) | 355^(a-1) |
|---:|---:|---:|
| 1+12 | 6,240,321,451 | 1 |
| 2+11 | 13,520,696,477 | 355 |
| 3+10 | 22,124,776,053 | 126,025 |
| 4+9 | 32,449,671,545 | 44,738,875 |
| 5+8 | 45,068,988,257 | 15,882,300,625 |
| 6+7 | 60,843,134,147 | 5,638,216,721,875 |

Thus (3) is automatic for a<=5. The 6+7 theorem retains the explicit
minimum-speed threshold floor(3Q/28). It is not automatically discharged by
(4). If both primitive shapes contain one, the smaller minimum is one and
the threshold holds for every split. The larger-shape unit assumption remains
part of the theorem.

## 5. Independent exact controls

The standalone companion imports no mathematical producer and passes 29,288
always-active gates. It checks 6,825 complete positive unit-pair boxes and
7,767 literal physical crossing constructions, retaining the cleared
coefficient and the coprime-scale gcd identity. Exact rational grid tests
give 1,620 translated-grid controls, including equality at the closed-spacing
boundary and strict interior with slack. A noncoprime t=4,g=2 hostile has
only two distinct phases, exposing why that hypothesis cannot be omitted.

Eleven actual thirteen-speed controls cover all six component splits and
include connected primitive smaller shapes without one. Set U to the first
b powers of three and t=2Q max(U)+1, g=1. The unit smaller shapes are the
first a powers of three. For a>=2 a second family uses powers of three
times two and times three, with total size a. The 2:3 atlas edge joins
these two strings, so their primitive smaller shape is connected and has
minimum two. Every row is checked positive, distinct, primitive and in the
physical box, with the actual two decoder components.

The independent dominance t>2Q max(U) proves absence of every bounded
crossing for these controls, hence W=V_dec. Their literal phase tau=1/4
has clearance exactly 1/4. It is also realized in the retained grid with
eta=3/4 and j=(t-3)/4. That chosen eta is an explicit strictly safe phase;
it need not attain the sharp lower-dimensional supplier value. The abstract
proof above separately uses the CITED supplier when an explicit phase is
not available.

In addition, all seven declared producer controls are independently
reconstructed, including both explicit lower-dimensional supplier phases,
their actual decoder graphs, physical box and entry-dominance proofs. The
audit uses the grid point nearest the arc center, rather than the producer's
first point beyond its left endpoint, and verifies every physical clearance.
The five-coordinate star is regenerated from 179,181,183,185. Its primitive
minimum is 1,013,861,907 and its physical sum is
90,168,220,083,627,413,785,737. The independent chosen phase is
`3785795013607/34072155122462`, and every full-row clearance is strictly
greater than 1/14. Thus the billion-minimum example has been checked as an
actual safe equality entry, not accepted merely from its size or a Boolean.

Run:

```
python -B 04-computation/overnight15_20260906_lrc_component_closure_audit.py
python -B -O 04-computation/overnight15_20260906_lrc_component_closure_audit.py
```

Both modes emit identical actual LF bytes. Independent source SHA256:
`72dc83b007422f079d5663686fa7312013669acf89dfabb348b7f7c326ee2bf2`.
Independent output SHA256:
`409ec49be081fa208382d37d2e065e36bc93620180c5966501e21fb5509f7979`.

The proof preserves the actual decoder equality, primitive component scales,
full arc length and labelled phase grid. The finite controls supplement
those arguments; they neither establish the cited lower-dimensional theorem
nor extrapolate safety to the remaining 6+7 shapes or nonunit larger cores.

**Filing:** root read the complete proof and independent audit. All source and
output bytes are retained; reproduction commands are relative to repository root.
