# Independent centered-division audit of the coprime-maximum entry criterion

**Status: PROVED INDEPENDENT REPROOF OF AN INCOMING SPECIALIZATION +
FINITE-EXACT + INDEPENDENTLY AUDITED.** The incoming
[signed-box criterion](overnight12_20260906_lrc_gcd_semigroup.md), read
directly from `origin/main` at the twelfth checkpoint before this note was
written, already proves a strictly stronger endpoint-gcd theorem. Its
Corollary 4 includes the entire closure reproved below. This note retains a
short independent centered-division mechanism, a distinct actual nonunit
entry control, and an independent audit of the stronger incoming lemmas.
It is not a new LRC closure, sharp coefficient-box radius, or theorem-ID
promotion. LRC(14) remains **OPEN**.

## 1. Exact inherited domain

Put `Q=91^6=567869252041`. Retain all of the actual-entry hypotheses in
[THM-3818, scaled inert cubeclass support-two pair packet](../../01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md),
Sections 6.4--6.5:

```text
n = tV disjoint-union g(p,q),
|V|=11, gcd(V)=1, t,g>0, gcd(t,g)=1,
n has thirteen distinct positive integer coordinates, sum(n)<=Q^2,
the actual decoder graph has exactly the components tV and g(p,q),
W_(Q,3)(n)=V_dec and rank(V_dec)=11.
```

Here `W_(Q,3)` is the span of **all** integer relations of support at most
three and coefficient height at most Q. The primitive pair satisfies
`1<=p<q`, `gcd(p,q)=1`, `p+q<=356`, and every prime in its sum is 2 modulo
three with exponent at most two. This is the actual 5,855-entry atlas;
in particular `p<=177`, `q<=355`. A chosen split, a connected selected
subgraph, or a rank computation on selected relations is not a substitute
for these hypotheses.

Let `K=max V`. The specialization considered here is

```text
there is u in V with u<K and gcd(u,K)=1.                  (1)
```

**Conclusion:** n has a phase of clearance at least `1/14`.
The literal normalized coordinate one is not required. Incoming
Corollary 4 already gives this whenever `gcd(u,K)<=76,388,115`.

The nearest preceding proof is the
[eleventh unit-component theorem](overnight11_20260906_lrc_unit_component.md).
The map below keeps two actual core columns and one actual outside column,
their integral scale, and the full coefficient budget. Its centered residue
choice discards other core labels only during the relation construction;
the complete core maximum K and a common phase grid remain in the gluing
step. The signed-box incoming result supplies the exact stronger radius.
The earlier primorial normalization hostile remains relevant algebraically,
but is already excluded as a strict failure by the large-subset gcd caps.

## 2. An elementary sufficient centered interval

Let `1<=u<K<=Q` with `gcd(u,K)=1`. Every integer x with
`2|x|<=QK` has a representation

```text
x=a u+b K,       |a|<=K/2,       |b|<=Q.                 (2)
```

Choose a centered representative of the solution
`a=u^(-1)x mod K`, so `|a|<=K/2`, and put `b=(x-au)/K`.
Then b is an integer and

```text
|b| <= |x|/K+u/2 <= Q/2+u/2 <= Q.
```

The endpoints are included. Signs and zero coefficients are allowed;
they do not increase support. This radius `QK/2` is only sufficient.
For example at `(Q,u,K)=(4,2,3)`, x=7 lies beyond it but equals `2u+K`.
Without coprimality, `(u,K,x)=(2,4,1)` is an elementary obstruction.

Incoming Lemma 1 gives the exact centered radius

```text
R=Q(u+K)-(u-1)(K-1) >= QK,
```

and its first missing integers are `±(R+1)`. Thus the centered construction
here provides an independent, simpler sufficient interval and does not
improve that result.

## 3. Actual crossing prohibition forces coherent gluing

The full equality branch of THM-3818, (15q), gives
`K/gcd(u,K)<=Q`; by (1), K<=Q. This implication uses the maximum column
and the actual equality branch, not just `gcd(V)=1`.

For any primitive positive pair p<q, choose

```text
z=j/(p+q),       p j=floor((p+q)/2) mod (p+q).
```

Both pair clearances are at least `1/3`. Hence its closed `1/14`-safe set
contains an arc of length `11/(21q)`. If `11t>=21q`, the symmetric cyclic
gluing in THM-3818 (15aa) closes the row: the t lifts of a core-safe phase
form a complete shifted t-grid in the pair clock because `gcd(t,g)=1`.
The weak endpoint is retained.

Otherwise `11t<21q` and therefore

```text
t <= floor((21q-1)/11) <= 677 < Q.
```

Put `delta=gcd(t,p)=gcd(t,gp)`, `c=t/delta`, and `x=gp/delta`.
If `2x<=QK`, equation (2) gives the literal physical relation

```text
c(gp)-a(tu)-b(tK)=0.                                    (3)
```

Its coefficients have height at most Q, its support is at most three,
and its partial sum on the pair component is the nonzero number cgp.
It is therefore in W but outside V_dec, contradicting the full equality
hypothesis. This includes cases a=0 or b=0. Thus equality entry in this
branch forces

```text
gp/delta>QK/2,
g>delta QK/(2p)>=QK/(2p)>42K,                            (4)
```

since `Q>84*177>=84p`.

Cited LRC for eleven nonzero speeds supplies a phase of V with clearance
at least `1/12`. The full core's clearance is K-Lipschitz, so its closed
`1/14`-safe set contains an arc of length `1/(42K)`. The g lifts of a
pair-safe phase form a complete shifted g-grid in the core clock. Equation
(4) makes its spacing strictly smaller than that arc length, yielding an
actual common physical phase. This proves the conclusion.

The lower-dimensional supplier is inherited from THM-3818 and the incoming
unit/signed-box notes; their cited source is arXiv:2604.23906v2, Theorem 1.3.
No lower-dimensional LRC proof is reproduced here.

The [joint-shadow gcd theorem](lrc14_joint_shadow_empty_core_next_sep06.md)
additionally implies t<=2 in any hypothetical strict counterexample,
because t is the gcd of an eleven-subset. The proof above does not need
that improvement, and does not assert t<=2 for every safe decoder entry.

## 4. Audit of the stronger incoming result and remaining support gap

The incoming signed-box note was read in full at `origin/main` before
finalization. Independent analytic audit: **PASS** for the following.

- Lemma 1: complementing x by `L=Q(a+b)` reduces the central interval to
  a nonnegative representation above `(a-1)(b-1)`. Its canonical residue
  representative has coefficients within the required box. The first gap
  `ab-a-b` proves exactness of R, and `2R>L` is retained.
- Lemma 2: every positive outside coefficient is a multiple of its minimal
  cleared coefficient c. If x is absent from the box, x>R>L/2, so every
  multiplier at least two leaves the entire box support. The conditions
  `b<=Q` and `c<=Q` are both necessary for this iff as stated.
- Corollary 4: for `D=gcd(u,K)<=H=floor(Q/(42*177))=76,388,115`, the small-t
  branch gives `c<=677H<Q`. Its exact radius satisfies `R>=QK/D`, so the
  crossing prohibition forces `g>42K`. It strictly subsumes Section 3.

The other incoming result,
[connected decoder gcd descent](overnight12_20260906_lrc_decoder_descent.md),
retains an actual boundary edge when extending a chosen subset: the gcd
loss divides its oriented coefficient, at most 355. Combined with our
seven-subset cap 90, this gives physical gcd caps `31,950`, `11,342,250`,
and `4,026,498,750` on core subsets of six, five, and four labels.
The five-label gcd cap is below H, but it is not a pair-gcd conclusion.
Compressing that identity into (3) would spend additional relation support.
No such compression is claimed here. The live question is whether actual
edge paths and full inherited subset profiles force an eligible pair or
another genuine support-three crossing.

## 5. Nonvacuity, hostiles, and exact controls

A fully typed entry without a unit label is

```text
V=(2,3,4,5,6,7,8,9,10,11,13), u=2, K=13,
t=1, g=2^45, (p,q)=(1,3).
```

The verifier rebuilds every decoder edge and obtains components 11+2.
The row is primitive, distinct and inside `sum(n)<=Q^2`. Also `g>2QK`.
A crossing relation with one pair label has a nonzero pair term of size
at least g and at most two core terms of total size at most 2QK. A crossing
with two pair labels has a pair sum that is either a nonzero multiple of g
or zero; neither can cancel one nonzero core term. Thus no bounded crossing
exists. Connected internal decoder rows span both component hyperplanes,
so W=V_dec. This is actual entry, not merely a selected-relation test.

Its literal safe phase and minimum clearance are

```text
x=5799621772885/70368744177664,
clearance=5026338869841/70368744177664 >=1/14.
```

The THM-4052 primorial shape is retained solely as an algebraic hostile:
it has gcd one, maximum `237907127334685115>Q`, all internal pair heights
at most 127, and no label coprime to that maximum. Its displayed subsets
of sizes seven, eight, nine, ten have gcds `392430,3810,30,5`. These violate
the current caps `90,30,9,4`, so it is not a surviving unsafe core.

The [standalone source](../../04-computation/lrc14_coprime_max_entry_empty_core_morning_sep06.py)
imports no mathematical producer. Its universe is 14,129 complete signed
centered-division tests for Qtoy=2,...,12; 2,304 cleared relation controls;
all 5,855 atlas pairs for the scale inequalities; 94 literal signed boxes
and 1,740 minimal-outside-multiplier controls for the incoming lemmas;
the actual nonunit entry, the normalization hostile, and 30 direct coherent
grid witnesses on two nonunit core shapes and both gluing branches. Generic
grid witnesses are not asserted to have the actual equality decoder graph.

Reproduction and hashes are recorded with the frozen output below. All
checks are explicit exceptions and remain active under optimized Python.

```text
python3 -B 04-computation/lrc14_coprime_max_entry_empty_core_morning_sep06.py
python3 -B -O 04-computation/lrc14_coprime_max_entry_empty_core_morning_sep06.py
```

Both runs pass **24,545 gates** and reproduce the
[saved output](lrc14_coprime_max_entry_empty_core_morning_sep06.out)
byte for byte. Semantic SHA-256:
`42e046c369798af9fb4a1fa1c6a0e371a90766c22033c6f861ceedf5cd1f6317`.

The root agent independently read the complete proof and accepted the
centered division, exact physical crossing relation, both gluing branches,
and typed nonunit dominance control: **PASS**. The root also independently
audited incoming Lemmas 1--2 and Corollary 4: **PASS**. The subsumption and
the distinction between an algebraic hostile and a surviving unsafe core
are retained explicitly.

Frozen raw LF SHA-256:

```text
source e65f9df0a7607e32a18b7e3aa2839c20bb054f00469e19d1b47a2a8663ef5fb0
output a104137ef7d6369b60d465fc1058300ffc41452becbd89c661e29e24fa49499e
```

All three owned files are frozen for root integration. No prior frozen
file, theorem ID, navigation file, or Git state was modified by this lane.
