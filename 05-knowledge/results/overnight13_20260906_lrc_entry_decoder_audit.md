# Independent audit: the complete actual-entry equality test

**Verdict: INDEPENDENT ANALYTIC PASS; FINITE-EXACT controls.** The theorem and
implementation in
[the thirteenth report](overnight13_20260906_lrc_entry_decoder.md)
correctly decide W=V_dec for a primitive positive distinct thirteen-speed row
in the Q² box whose actual decoder graph has components of sizes11 and2.
No repair is requested. This is an entry-equality theorem, not a safe-phase
classification or a proof of LRC(14).

The referee derived the all-orientation argument and independent arithmetic
path before reading producer code. After the independent engine passed, the
implementation was read through its entry routine. The domain rejection,
physical-coordinate sorting, witness labels, internal-height gate, and all121
support construction agree with the proof. In particular, an out-of-box row
does not receive a false negative equality verdict.

## 1. Complete proof and finite-box necessity

Under the diagonal rational map z_i↦n_i z_i, each weighted decoder edge row
becomes a nonzero incidence difference. Connectedness therefore identifies
the decoder span exactly with the two componentwise weighted kernels. Its
dimension is11, and it is contained in W because every actual edge has
primitive coefficient height at most355<Q.

Suppose an internal pair has primitive height greater thanQ. On that pair
and any opposite-component label, the (Q+1)³ nonnegative coefficient triples
have at most Q³+1 possible sums. The physical bound is used here explicitly.
Two triples coincide, producing a nonzero bounded relation. Its opposite-
component coefficient cannot vanish: that would yield a bounded internal
pair relation, contradicting the primitive height. Hence W contains a vector
outside the decoder kernel. The internal-height condition is necessary.

Once every internal pair has height at mostQ, the audited minimal-coefficient
lemma applies to every mixed triple. Global sign reversal makes the distinguished
coefficient positive; zero coefficients on the paired labels remain allowed.
The exact interval test is therefore equivalent to a bounded crossing on that
support. Any such crossing has a nonzero component partial sum and lies outside
V_dec. Conversely, every relation of support at most three meeting both components
is captured by one of these tests. A crossing with only two occupied labels is
captured by adjoining an unused label in the same component as one of them.

Thus all121 negative tests imply that every generator of W is internal, so
W⊆V_dec. Actual decoder rows supply the reverse containment. The count is
2·C(11,2)+11=121. The old110 supports omit precisely the orientation whose
distinguished coordinate lies in the eleven-component. No necessity/sufficiency
direction, sign, zero coefficient, or label orientation is lost.

The proof of this pointwise theorem needs no lower-dimensional LRC phase
supplier. The generic pigeonhole argument is reproved from the stated box;
the actual atlas and its weighted-kernel interpretation are inherited from
**THM-3818**, `01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`,
Section6.4, and the coefficient lemma from the audited twelfth report.

## 2. Independent arithmetic and literal rational spans

The referee does not call the producer's interval routine. For x≥0, it instead
sets M=Q(a+b)−x and seeks M=au+bv with 0≤u,v≤2Q. It chooses
v₀=b⁻¹M mod a, using v₀=0 when a=1, and u₀=(M−bv₀)/a. All solutions are
v=v₀+ah and u=u₀−bh. Intersecting these nonnegative upper and lower bounds
gives an independent membership decision and a signed witness r=Q−u,s=Q−v.
Negative x is handled by symmetry. This reverses the modular variable and
uses the nonnegative complement representation rather than the producer's
signed-coefficient base solution.

A separate literal engine enumerates every relation of support two or three,
up to overall sign: the first occupied coefficient is positive; a second
coefficient is enumerated and the last is solved by exact divisibility.
Fraction-free Gaussian elimination with integer content removal computes
the rational rank. The direct equality test checks both rank n−2 and zero
weighted partial sum on each proposed component. It never substitutes a
selected crossing family for the complete relation span.

The complete toy universes are primitive four-element subsets of1,…,13 at
Q=4,…,8 and primitive six-element subsets of1,…,10 at Q=5,…,9, restricted
to the physical Q² box. All2+2 partitions, once up to swapping components,
and all4+2 partitions are tested. This gives **3,105 rows, 18,939 partitions,
and335,957 literal relation rows**. Every arithmetic decision agrees with the
full rational span. These dense toy universes happen to contain zero equality
partitions; this is reported explicitly, and separate positive controls follow.
They are abstract partitions, not claims about the actual atlas components.

The actual decoder positive toy rows (1,3,319,957) at Q=53 and
(1,3,9,1819,5457) at Q=101 lie in their boxes. Their full literal relation
banks have34 and2,299 rows up to sign, respectively, and ranks2 and3. Every
relation is internal; every arithmetic orientation is negative. These give
literal positive equality controls independently of the large-Q certificates.

## 3. All large certificates and both hostiles

The referee replays every one of the **363 retained support certificates**
across the three in-domain thirteen-speed controls. It checks unique support
coverage, orientation, physical gcd clearing, canonical solution, both claimed
integer interval endpoints, the independent complement result, and every full
positive witness. An independent prime sieve and multiplicative admissible-sum
construction rebuild the actual graph; union-find and literal weighted-edge
rank recover its two components and dimension11.

The canonical unit and unitless rows have no positive support in either
orientation. Their independent dominance proof excludes every crossing, so
their actual equality claim does not rest on an isolated rank sample.

For t=3Q+1 and the scaled canonical core tU with pair(1,3), all110 original
tests are negative. The opposite orientation has exactly one positive test,
with relation t−1−3Q=0. The returned vector raises literal decoder-plus-witness
rank from11 to12. A separate small actual-graph example at Q=41 uses core
(124,496,744), pair(1,3), and the same coefficients(1,−1,−41). It lies inside
the box, passes every first-orientation test, and its complete literal relation
span has rank4 while its decoder span has rank3.

For the producer's powers-of355 core with g=2Q·355¹⁰+1, dominance proves
actual equality, yet the internal pair height exceedsQ. Its sum exceedsQ²,
so the algorithm correctly reports out-of-domain. The smaller actual graph
row (1,3,9,127,381) at Q=7 demonstrates the same failure with a complete
literal relation bank: equality holds at rank3, although internal height9>7.
This independently discharges the physical-box firewall.

An additional in-box thirteen-speed height-gate control, proposed separately
by the producer during the audit, was independently verified. Its core is
{355^j:0≤j≤5}∪{3,9,27,81,243}, K=355⁵, with outside pair g,3g and
g=1000K+1. It has the actual11+2 graph and lies in the box, but K>Q. The
literal relation g−1000K−1=0 has height1000≤Q and raises rank11 to12.
This tests the height-failure branch directly; no minimal-coefficient lemma
is applied to the inadmissible large primitive pair.

## 4. Frozen inputs and reproduction

The independent script checks these producer raw-byte hashes before use:

```
source 9b8fe6804a037a0f93a29396940e510691ee3b2765be1eeca1271ed28e5b7b6c
output 8932a324412f5c0b7ce8cfc243146af374454a2239aa6059030f00af261bd57d
JSON   af177ea9703390a8169d9466d04a5e4d8bb53f66833cfbfd223c94a0006cb847
```

[Independent source](../../04-computation/overnight13_20260906_lrc_entry_decoder_audit.py)
and [output](overnight13_20260906_lrc_entry_decoder_audit.out):

```
python -B D-computation/overnight13_20260906_lrc_entry_decoder_audit.py
python -B D-computation/overnight13_20260906_lrc_entry_decoder_audit.py
```

Both runs pass **31,754 always-active gates** with byte-identical LF output.
In addition to the universes above,3,708 signed integers independently validate
the complement membership routine against full coefficient sets. The source
supports adjacent outside-folder artifacts and repository filing with source
and JSON in `04-computation` and output in `05-knowledge/results`.

Referee source SHA256:
`a580b7865ff4691ed481f55772d70ce5426c0c40eb864f19ebd5476974ed78a1`.
Output and optimized-output SHA256:
`7861e9ddc8bd460025b01b77728430ac51472cd52e4c955733541a24ab4f7123`.
The finite computations support the all-height proof and its scope boundaries;
they are not an exhaustion of the inherited speed box. No Git, navigation,
producer code, or earlier frozen artifact was modified by this referee.

**Filing:** root integrated these independently audited artifacts in the thirteenth
checkpoint. Reproduction commands are relative to the repository root.
Transcript hashes refer to filed LF bytes; Windows CRLF captures were normalized.
