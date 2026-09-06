# Quantitative retained pair packets and unbounded body suppliers

**Status: PROVED elementary + FINITE-EXACT controls; independent full
written-proof and replay audits PASS.** The source
is the old pack-clock count, now applied to a geometrically retained body
packet. It gives an actual safe-measure lower bound, a uniform `1/12` tail
margin, and an explicit family with unbounded body and tail heights that
defeats every clock of denominator at most fourteen. It does not supply
arbitrary chart entry or settle LRC(14).

## 1. Inheritance and the retained information

The closest proved engine is
[THM-737 — pack-clock sampling and measure dispatch](../../01-canon/theorems/THM-737-pack-clock-sampling-measure-dispatch.md):
an open danger arc meets a residue grid in a bounded number of points,
with gcd multiplicity retained. The exact body-phase input comes from
[THM-2072 — half-shift certificate](../../01-canon/theorems/THM-2072-fixed-owner-clock-bank-no-go-and-half-shift-certificate.md)
and the audited
[scale-three transport](lrc14_antipodal_consumer_transport_overnight_hexagon_sep05.md).
The canonical hostile is the actual ten-body Haar floor failure. The
corrected near miss is either to infer a scalar body floor, or to normalize
the tail without transporting the actual body phases. The least-used
sidecar is the **measure of a translated safe-set intersection**.

The concept board is: gcd-weighted branch counts; actual dyadic valuation;
geometric rather than arithmetic packs; translated safe-set intersections;
strict endpoint counts; residue-preserving body lifts. The map sends two
located body phases to six physical lifts. It preserves body safety and
each tail's residue multiset, but not an integer common-factor six-pack of
the body. The phase separation is the required sidecar.

Incoming
[THM-4442 — bounded ten-body parity-free completion](../../01-canon/theorems/THM-4442-lrc14-bounded-ten-body-parity-free-scale-three-completion.md)
already gives point existence for every ten-subset of `[13]` and every
ternary-unit tail triple. Its proof was read completely. No bounded body or
tail census is repeated here. The new consumer is quantitative measure,
and Section 4 has genuinely unbounded bodies. Exact searches for the packet
phase, modulus, and displayed constants found no prior identical family.

## 2. Actual-tail pair-packet measure theorem

For a finite positive body `C` and `0<lambda<=1/12`, write

```text
G_C^lambda={y in R/Z: ||cy||>=lambda for every c in C},
L_lambda(S)=mu{x: ||sx||>=lambda for every s in S}.
```

Let `T` consist of three distinct positive integers, none divisible by
three. Use the **actual**, possibly nonprimitive tails to define

```text
g=2^min_(w in T) v_2(w),     u_w=w/g,
e=#{w in T:u_w even},        0<=e<=2,
delta=1/(2g),
A_delta=G_C^lambda intersect (G_C^lambda+delta).       (1)
```

Then

```text
L_lambda(3C union T) >= (3-e) mu(A_delta)/6.           (2)
```

This is a lower bound for actual physical safe measure. If the translated
intersection consists only of isolated points its measure bound is zero,
but any one point still gives a qualitative safe witness at that threshold.

### Pointwise residue count

For `y in A_delta`, both `y` and `y-delta` are body-safe. Hence all six
physical points

```text
x_(epsilon,j)=(y-epsilon delta+j)/3 mod1,
epsilon=0,1,       j=0,1,2,                           (3)
```

are safe for `3C`. For actual speed `w=g u`, its phases on this packet are

```text
w y/3-epsilon u/6+j g u/3 mod1.                      (4)
```

Since `g` and `u` are both units modulo three, (4) is a full six-point grid
if `u` is odd, and a three-point grid with each point repeated twice if
`u` is even. The strict dangerous arc has length `2lambda<=1/6`. It meets
at most one distinct point in either grid; equality at an endpoint does not
increase the count because danger is open. Thus each odd `u` kills at most
one packet point and each even `u` at most two. At least

```text
6-[(3-e)+2e]=3-e                                      (5)
```

points survive all tails. This is precisely THM-737's retained gcd count,
not a new counting principle. When `g>1`, the six **physical** points need
not be uniformly spaced; the six-grid statement concerns each actual
tail's observed multiset (4).

### Integrating without losing a factor

Let `N(y)` count survivors in (3). For each fixed `epsilon`, the three maps
in `j` cover the physical circle once, with `dy=3dx`. Therefore

```text
integral_(A_delta) N(y)dy
 =3 integral 1_(tail safe)(x)
       [1_A(3x)+1_A(3x+delta)] dx
 <=6 L_lambda(3C union T).                            (6)
```

Every nonzero summand on the middle line is body-safe by the definition of
`A_delta`; all phase arguments are modulo one. Equations (5)--(6) prove
(2). This proves the general dyadic statement without assuming that
`A_delta` itself is `delta`-periodic.

For `g=1`, the actual physical packet is the clock-six grid
`y/3+k/6`; here `A_(1/2)` is half-periodic. Restricting `y` to a half-circle
gives the same factor `1/6`, not `1/3`.

### Qualitative stronger-margin form

If two actual body phases separated by `delta` have common body margin
`alpha`, apply the same count at

```text
lambda=min(alpha,1/12).
```

The full row has a physical point with margin at least that value. No
positive-measure body hypothesis is needed for this point statement.

## 3. Exact gains for the earlier hostile packets

The midpoint choices in the three positive antipodal intersections give:

| Body from the preceding sidecar | Phase `y`, with mate `y+1/2` | Common body margin | Uniform full-row margin |
| --- | --- | ---: | ---: |
| first height-thirteen Haar hostile | `167/2002` | `15/182` | `15/182` |
| height-thirteen minimum-mass hostile | `695/2912` | `33/364` | `1/12` |
| incoming clock-filtered height-eighteen hostile | `515/6188` | `39/476` | `39/476` |

Each uniform margin applies to every actual ternary-unit tail triple with
an odd coordinate. The first two point-existence families were already
subsumed by THM-4442; their stronger margins and uniform-in-tail-height
measure bounds are the useful quantitative survivor. The height-eighteen
body is outside THM-4442's bounded universe.

At `lambda=1/14`, their antipodal intersection masses are respectively
`8/1001`, `1/52`, and `1101/136136`, so (2) yields their exact explicit
floors after multiplication by `(3-e)/6`. These floors do not deteriorate
when the actual tail heights grow.

## 4. An explicit unbounded body-and-tail family

Put

```text
P=480480=165*2912,
B=(1,3,5,7,8,9,11,12,13,11650),
c_i=B_i+P n_i,             n_i arbitrary nonnegative integers,
C={c_1,...,c_10},
T={14,u,v},                                                (7)
```

where the three tails are distinct positive ternary units and at least one
of `u,v` is odd. Body and tail heights are independently unbounded. Then

```text
max_x min_(s in 3C union T)||sx|| >= 1/12.                (8)
```

Every resulting row has thirteen distinct speeds and is primitive, and
**every rational clock of denominator at most fourteen fails**.

### Proof of the packet and primitivity

Take `y=695/2912`. The common distances over `y,y+1/2`, multiplied by
2912, for the base body in its displayed order are exactly

```text
(695,629,563,497,264,431,365,396,299,1390).             (9)
```

Their minimum is `264/2912=33/364>1/12`. The last body speed is congruent
to two modulo 2912, so (9) is also the earlier minimum-Haar body's packet.
The increment `P` is even and divisible by 2912, hence all ten independently
chosen increments preserve both actual phase residues. Section 2 proves
(8), with `g=1` because at least one actual tail is odd.

The base residues are distinct modulo `P`, so the body speeds are distinct.
The two parts of the physical row are disjoint modulo three. Finally
`c_1=1+P n_1` is one modulo fourteen, so `gcd(3c_1,14)=1`; the complete row
is primitive without imposing a hidden primitive-body or primitive-tail
filter.

### Exact small-clock obstruction

The modulus preserves each of the following pins:

```text
q=2,3,4,6,8,12: 3c_i with B_i=8 is divisible by24;
q=5,10:         3c_i with B_i=11650 is divisible by10;
q=7:            3c_i with B_i=7;
q=9:            3c_i with B_i=3;
q=11:           3c_i with B_i=11;
q=13:           3c_i with B_i=13;
q=14:           the actual tail14.                     (10)
```

Consequently, for every `q=2,...,14` and every numerator `p`, some actual
speed has zero distance at `p/q`. Denominator one also fails trivially.
This is a constructed family of missed small-clock certificates, not an
unresolved family or a counterexample to LRC. The retained six-point packet
has denominator dividing 8736. It is deliberately a family-dependent
packet, compatible with THM-2072's no-go for one clock bank uniform over all
bodies.

### Explicit measure floors independent of tail height

Write `H=max(C)` and let `e` be the number of even speeds in `T`, so `e` is
one or two here. At the body packet (9), the slack above `1/14` is `1/52`.
The Lipschitz bound `| ||cy||-||cz|| |<=c|y-z|` retains both body phases in
two disjoint intervals of radius `1/(52H)`. Hence
`mu(A_(1/2))>=1/(13H)`. Applying (2) yields

```text
L_(1/14)(3C union T) >= (3-e)/(78H).                  (11)
```

At the stronger threshold `1/12`, the body slack is `2/273`, and the same
argument yields

```text
L_(1/12)(3C union T) >= 4(3-e)/(819H).                 (12)
```

These bounds depend on body height but not on either free tail's height.
The two intervals are disjoint because their centers are a half-circle
apart and both displayed radii are below `1/4`.

## 5. Sharp count, hostile, and stopping boundary

The survivor count `3-e` is sharp pointwise for the supplied information.
At `lambda=1/14`, the exact packets

```text
T=(1,5,7),   y=5/84:    three survivors;
T=(1,4,5),   y=83/420:  two survivors;
T=(1,4,10),  y=83/420:  one survivor
```

attain it. The transcript has corresponding controls at `lambda=1/12`
and after actual dyadic dilations `g=2,4,8`. No global optimality of the
Haar constants or of the family margin is claimed.

At the very same geometric body packet `y=695/2912`, the actual all-even
triple `(4,26,34)` kills all six half-shift points, already at `1/14`.
Its three dangerous pairs are `{1,4}`, `{0,3}`, and `{2,5}` in the uniform
clock-six order. The correct theorem instead uses its actual `g=2` and
asks for a body pair separated by `1/4`. Dropping the valuation coordinate
is therefore false, not merely inconvenient.

The theorem remains conditional on an actual translated body-safe
intersection. The preceding sidecar refutes a universal half-pair supplier
by an inclusion-minimal nine-core, and a ten-body extension. Neither the
constructed residue family nor the incoming bounded theorem supplies the
general body geometry or arbitrary chart entry.

## Reproduction and audit

```bash
python3 -B 04-computation/lrc14_quantitative_pair_packet_overnight_hexagon_sep05.py
python3 -B -O 04-computation/lrc14_quantitative_pair_packet_overnight_hexagon_sep05.py
```

The [script](../../04-computation/lrc14_quantitative_pair_packet_overnight_hexagon_sep05.py)
imports no carrier or network engine. Fifty-four literal physical interval
computations at both thresholds and actual dyadic scales are checked against
an independent integral of six-point survivor counts. Thirty-six controls
give positive certificate floors; zero translated intersections are retained
honestly. The unbounded arithmetic family is proved above, while named
literal residue controls reach body height `480480000011650` and multi-million
tails without an infeasible interval expansion. The
[output](lrc14_quantitative_pair_packet_overnight_hexagon_sep05.out) records
the universe, exact count controls and stopping boundary. Root and observer
message-level audits agree on the generalized dyadic proof and its factor
`1/6`. The observer's subsequent full written-package audit and independent
replay passed all 223 gates and 54 measure controls, including the modulus,
every divisor pin, distinctness, full-row primitivity, and both explicit
Lipschitz floors. Its semantic digest agrees with the frozen output.

Normal and optimized executions have identical output, with 223
optimization-live checks. Frozen source SHA256:
`47125ede363b024b70e36162307319722af48a1d7875ebdcd03195bb42778ca8`;
output SHA256:
`c530c5dd93478a34dfc2d8b6685d9319120b78417fb0dd5d93cf04f5411d4a4b`.
The root's separate full written-proof audit of Sections 1--5 passed,
including the residue family, denominator pins, primitivity and both
Lipschitz measure floors.
