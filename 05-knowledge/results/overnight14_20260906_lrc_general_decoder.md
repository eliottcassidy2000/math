# A complete decoder-span equality test for thirteen-speed rows

**Status: PROVED ELEMENTARILY ON THE STATED DOMAIN + INDEPENDENTLY AUDITED; FINITE-EXACT controls.**
This is a general-component corollary of the
[audited thirteenth entry decoder](overnight13_20260906_lrc_entry_decoder.md),
with its [independent referee](overnight13_20260906_lrc_entry_decoder_audit.md).
It decides equality of two relation spaces for any row in the declared box.
It is not a safe-phase criterion, a classification of all ranks in every graph
branch, or a proof of LRC(14). No theorem identifier or priority claim is made.

## 1. Inheritance, target, and the recovered singleton coordinate

Fix Q=91^6=567,869,252,041. Let n be thirteen positive distinct integers with
gcd(n)=1 and sum n_i≤Q². Construct the actual all-scale graph of
**THM-3818**, `01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`,
Sections1 and6. For an unordered pair, reduce its coordinates to coprime p<q.
The edge exists exactly when p+q≤356 and every prime divisor of p+q is2 modulo3
with exponent at most two. This is the inherited5,855-ratio atlas. It supplies
primitive edge coefficient heights at most355<Q.

Let c be the number of actual graph components. Let V_dec be the rational
span of its weighted primitive edge relations, and W the rational span of
all integer relations of support at most three and coefficient height at
mostQ. The target is the literal equality W=V_dec. The diagonal coordinate
map z_i↦n_i z_i identifies V_dec with the direct sum of ordinary zero-sum
spaces on the graph components; thus dim V_dec=13−c and V_dec⊆W.

The closest mechanism is the thirteenth two-component crossing test. The
canonical hostile is its in-box row with110 negative supports and one
opposite-orientation crossing. The corrected near miss is dropping the
physical box before rejecting an internal pair height. The least-used sidecar
is a singleton component: packet data may lose its physical scale in
THM-3818 Section6.4, but this routine receives the full physical row and
therefore retains that scale in every distinguished-coordinate test. It does
not reconstruct missing speeds from packet data.

The five-concept board is: graph-component count; the physical height box;
bounded relation rank; both crossing orientations; and singleton scale.
Anchor: exact LRC entry data. Niche: two-generator signed boxes. Wildcard:
the distinction between decoder connectivity and connectivity of the full
bounded relation span. The new sharp hostile in Section5 has three decoder
components while the full bounded span has rank exactly eleven.

## 2. The general classification

**Theorem.** In the declared domain:

1. If c=1, then W=V_dec and both have dimension12.
2. If c≥3, then W≠V_dec. The general bound is11≤dim W≤12; component
   count alone does not choose between these two values.
3. If c=2, with component sizes a,b and a+b=13, equality holds if and only if
   every internal primitive pair has height at mostQ and every mixed triple
   has a negative signed-box test below. There are exactly11ab/2 such triples,
   at most231. Equality gives dim W=11; failure of either gate gives dim W=12.

For a mixed triple, select the two labels in one component with physical
coordinates A<B and the distinguished coordinate Y in the other. Put

    D=gcd(A,B),  a₀=A/D,  b₀=B/D,
    δ=gcd(D,Y),  c₀=D/δ,  x=Y/δ,
    r₀=a₀⁻¹x mod b₀,  0≤r₀<b₀,  s₀=(x−a₀r₀)/b₀,
    lower=max(ceil((-Q−r₀)/b₀),ceil((s₀−Q)/a₀)),
    upper=min(floor((Q−r₀)/b₀),floor((s₀+Q)/a₀)).             (1)

The test is positive exactly when c₀≤Q and lower≤upper. It is evaluated
only after the internal height gate ensures b₀≤Q. A positive test returns
a literal crossing by choosing h between the endpoints and using coefficients
−(r₀+hb₀),−(s₀−ha₀),c₀ on A,B,Y. The distinguished coordinate can belong
to either component, including a singleton.

## 3. Proof and the role of the finite box

The full lower rank bound is the generic pigeonhole argument in **THM-2052**,
`01-canon/theorems/THM-2052-finite-height-forces-high-rank-bounded-relation-code.md`,
Sections1–2. On any three labels, the (Q+1)³ coefficient vectors in {0,…,Q}³
give at most Q³+1 sums, so there is a nonzero bounded relation on those labels.
If W⊥ had dimension at least three, a matrix whose rows span W⊥ would have
three independent coordinate columns, contradicting that relation. Therefore
dim W≥11. Positivity of n gives dim W≤12. This argument is unconditional
for the declared box; no lower-dimensional safe-phase theorem is invoked.

For c=1, the decoder edges already span dimension12, proving equality.
For c≥3, dim V_dec≤10<11≤dim W, proving strict inclusion.

For c=2, an internal primitive pair above heightQ has no bounded pair relation.
Adjoin any label in the other component and apply pigeonhole. The resulting
bounded relation must have a nonzero coefficient at that distinguished label;
otherwise it would be a forbidden pair relation. Its component partial sum
is nonzero, so it lies outside V_dec. Consequently the internal height gate
is necessary, and its failure adds the twelfth dimension.

When all internal heights pass, the audited minimal-coefficient lemma makes
(1) equivalent to a bounded crossing on that mixed triple. Any crossing of
support at most three has two labels in one component and one in the other,
after adjoining an unused label if its support is only two. Such an unused
label always exists here, since the total size is thirteen. A global sign
change makes the distinguished coefficient positive. Thus no sign or
support-two relation is lost, even for a1+12 split.

All negative tests exclude every crossing generator of W. Every remaining
generator belongs to a componentwise weighted kernel, hence to V_dec.
The reverse inclusion comes from decoder edges, proving equality. A positive
test explicitly adds a vector to the rank-eleven decoder span, so dim W=12.

The number of mixed triples is

    C(a,2)b + C(b,2)a = ab(a+b−2)/2 = 11ab/2.                (2)

Since a+b=13, the maximum ab=42 occurs at6+7, giving231. The six unordered
size types have counts66,121,165,198,220,231. The1+12 case has66 tests;
there are no internal pairs in the singleton, and the other orientation
still captures every crossing. The argument classifies all possible actual
component counts, not just the original11+2 branch.

## 4. Actual positive controls for every two-component size type

For each a=1,…,6, set b=13−a and

    A=(1,3,…,3^(a−1)),  K=3^(a−1),  g=2QK+1,
    full row = A ∪ g(1,3,…,3^(b−1)).                         (3)

Every row is positive, distinct, primitive, and in the physical box. Its sum
is less than Q·3¹²+3¹³<Q². Each group is connected by actual ratio1:3 edges.
There are no edges between the groups: if a≥2 then g≡1 mod3, so every reduced
cross ratio retains a coordinate at leastg; if a=1 its only small coordinate
is1 and the same conclusion is immediate. In both cases its sum exceeds356.

The relation-space equality has an independent proof. Every nonzero partial
sum on the large group is an integer multiple ofg. It cannot be canceled by
at most two bounded terms from A, whose total absolute value is at most2QK<g.
A zero large-group partial sum cannot cancel a nonzero single small term.
Thus no bounded crossing exists, and connected internal edges give W=V_dec.

All six actual rows pass every mixed-support test, including all231 at6+7.
The verifier checks label rotation and reversal, including the singleton
orientation, and separately checks membership through the nonnegative
complement representation Q(a₀+b₀)−x. These controls show that the maximum
support count is attained on an actual equality row; no claim is made that
231 is the minimal possible algorithmic decision complexity.

The connected row(1,…,13) supplies the c=1 case. The thirteenth opposite-
orientation hostile and its additional in-box internal-height hostile supply
the two distinct rank-twelve rejection paths. An out-of-box powers-of355 row
is rejected as outside the domain, not assigned a false equality verdict.

## 5. Three components can give either bounded rank

The following control prevents a false inference from c≥3 to the rank-twelve
case. Let

    A=(1,3,…,3⁹,1000),  K=3⁹=19683,
    g=2QK+1,  n=A ∪ (g,3g).                                  (4)

The actual graph has components of sizes10,1,2. The powers of3 are connected.
The coordinate1000 is coprime to each power and has ratio sum greater than356
with it. It also has no edge to the large pair; this follows from the exact
ratios and is independently reconstructed by the verifier. The row lies
inside the physical box. Its decoder rank is10.

The bounded relation1000·1−1000=0 joins the isolated coordinate to the ten-
component and adds one dimension. Together with the decoder edges, these
literal rows have rank11. On the other hand, g>2QK forbids every bounded
crossing from the eleven coordinates A to the pair(g,3g). All of W lies in
their two componentwise weighted kernels, of total dimension11. Hence

    dim W=11,  dim V_dec=10,  W≠V_dec.                         (5)

For the other endpoint, take A=(1,3,…,3¹⁰), K=3¹⁰, g=QK+1, h=2g+1,
and n=A∪{g,h}. This is another primitive in-box row, with actual graph sizes
11,1,1 and decoder rank10. The two bounded relations

    g−QK−1=0,             h−2g−1=0                            (6)

add two independent dimensions, yielding dim W=12. Both records therefore
belong to the c≥3 non-equality branch, while their full bounded ranks differ.
The implementation correctly returns only the interval[11,12] in that branch.
Neither branch is labeled safe by this theorem.

## 6. Implementation and exact reproduction

[Standalone implementation and verifier](../../04-computation/overnight14_20260906_lrc_general_decoder.py)
and [output](overnight14_20260906_lrc_general_decoder.out):

```
python -B 04-computation/overnight14_20260906_lrc_general_decoder.py
python -B -O 04-computation/overnight14_20260906_lrc_general_decoder.py
```

The routine first checks the physical domain, then the actual graph. Connected
and at-least-three-component cases use the proved rank bounds. A two-component
row gets all internal height gates and at most231 exact mixed-support tests.
Each test uses a fixed number of integer gcd, inverse, division, and rounding
operations; bit cost depends on operand size. No Q-sized enumeration is used.

The verifier regenerates the actual atlas, checks literal weighted-edge ranks,
replays the1,001 mixed supports across all six positive size types by an
opposite-variable complement calculation, and checks all stated controls,
including explicit rank witnesses in both c≥3 cases. The universal theorem
comes from the proof above, not this finite bank. Both normal and optimized
runs pass **4,744 always-active gates**, with byte-identical LF output.

Source SHA256:
`2bf80cad3bef06ef5c2b4ccb5e524190bfd09079591ac3aa2a1c278fd627126e`.
Output and optimized-output SHA256:
`10ed81b1cc784a6a2ff941bcc730bc463540641d4ca0a3a296420c7b5c8498b6`.

No Git, navigation, scarce identifier, or prior frozen artifact was modified.

The [independent audit](overnight14_20260906_lrc_general_decoder_audit.md)
accepts the full theorem and implementation with11,832 additional gates.

**Filing:** root read these proofs and audits and integrated the fourteenth
checkpoint. Reproduction commands are relative to the repository root.
