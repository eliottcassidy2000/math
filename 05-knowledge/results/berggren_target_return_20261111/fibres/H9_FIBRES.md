# The complete decorated H9 fibre and its Berggren responses

**PROVED in this note; FINITE-EXACT independent controls.** This classifies
one-apex/transitive-base tournaments with a named apex and ordered base,
not all tournaments of count9. It uses the proved append-letter law from
[Session256's graph lift](../../session256_chan_pak_triples_20261111/spinor/PROOF.md).
The general ordinary one-apex classification is older work, recovered in
`../../session26_random_bridges_20261022/session26-tournament/APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md`;
the specialization and exact branch-padding responses below are the present
connection to107,73,17 and30000.

## 1. Exactly three decorated cores, at every order

Let T(w) have a totally ordered base with all arcs forward; the apex beats
exactly the O letters and loses to the I letters. Every word of count9 is
uniquely

```
I^a v O^b,    a,b>=0,    v in {OIII,OIOI,OOOI}.                 (1)
```

Proof. Leading I vertices are global sources and trailing O vertices are
global sinks. Deleting them preserves H. The remaining word starts O and
ends I; if it is empty H=1. For a nonempty core of length l, the Hamilton
paths with apex between the first and last base vertices contribute
`2^(l-2)`, and the two endpoint-apex paths contribute2. These are disjoint,
so `H>=2+2^(l-2)`. H9 forces l<=4. The short closed-form enumeration gives
OI at l2 with H3, OOI/OII at l3 with H5, and at l4 precisely OIII,OIOI,OOOI
with H9; OOII has H11. This proves (1) for every length, not only the
tested lengths.

Consequently the decorated generating function by graph order is
`3z^5/(1-z)^2`: there are `3(n-4)` decorated words at order n>=5.
OIII and OOOI give the same ordinary strong tournament with different
apices, whereas OIOI gives the other strong five-vertex H9 class.
Their exterior SCC positions are intrinsic, giving ordinary generating
function `2z^5/(1-z)^2` for this subfamily. The missing H9 source type
`C3 -> C3` has feedback vertex number2 and is outside this one-apex chart.
These are contributions to A000568, not a replacement for that sequence.

The word OIOIO in the user's example has one global sink. Thus107,73,17
come respectively from `(v,b)=(OIOI,1),(OOOI,0),(OIII,0)` under A.

## 2. Complete A response, including arbitrarily much padding

Write the state after the substituted core A(v) as (h,p,e), where p is
the binary integer O=1,I=0. In core order OIII,OIOI,OOOI these pairs are

```
(h,p)=(17,16), (33,36), (73,84).                            (2)
```

A fixes I and replaces O by OI. Leading I remains harmless, while every
trailing O becomes an active OI block. Appending OI sends

```
(H,P) -> (H+2P+2,4P+2).
```

Solving this recurrence gives the full all-order response

```
H(A(I^a v O^b)) = h + [2(3p+2)(4^b-1)+6b]/9.             (3)
```

Equivalently the three rays are `(c_i*4^b+6b+d_i)/9` with

```
(c_i,d_i)=(100,53),(220,77),(508,149).                    (4)
```

The b term is real: an affine dyadic ray alone loses the repeated-block
contribution. Padding preserved the source count, but the substitution
turns it into feedback and changes the branch response without bound.
This is why a scalar-H quotient cannot support this branch operation.

## 3. Complete B and C responses

For a substituted core u define p(u), h(u), its length L, and
`alpha(u)=sum_{j:u_j=I}2^(j-1)` using positions j starting at1. Since u
starts O, appending u to a state adds `alpha(u)*P+h(u)-1` to H and replaces
P by `2^L*P+p(u)`; the old last-bit coordinate disappears at that first O.

For C, the left padding becomes (OOI)^a and the right padding remains
harmless. The complete response is

```
H(C(I^a v O^b)) = h(Cv)
  + [(42*alpha(Cv)+24)(8^a-1)+28a]/49.                    (5)
```

For B, first set

```
h_* = h(Bv)+[(35*alpha(Bv)+10)(8^a-1)+28a]/49,
p_* = p(Bv)+2^L*5(8^a-1)/7.
```

Then

```
H(B(I^a v O^b)) = h_*+[2(3p_*+2)(4^b-1)+6b]/9.          (6)
```

Proof: OOI sends `(H,P)->(H+4P+4,8P+6)`; OIO sends it to
`(H+2P+2,8P+5)`. Sum the geometric series, apply the core response, and
then use the OI block law. These formulae retain both exterior positions;
there is no exchange of independent branch labels.

## 4. A complete fixed-translation obstruction

No two values in the A image of the entire decorated H9 fibre differ by
30000. This is an all-order statement about (3), not a global tournament
or prime obstruction.

For unequal pairs (i,b),(j,c), put M=max(b,c). The leading difference in
(4) has absolute value at least `27*4^M`: at equal exponents the minimum
coefficient difference is120; at adjacent exponents the least absolute
normalized difference is `|100-508/4|=27`; at a gap at least2 it is at
least `100-508/16>27`. The lower-order difference has absolute value at
most `6M+96`. Thus a difference30000 requires

```
27*4^M-6M-96 <=270000.
```

The left side is increasing and already442230 at M7, so M<=6.
The exact21-value residual contains no pair at difference30000.
`h9_fibres.json` gives those values and the empty pair list. This combines
an infinite cutoff proof with a tiny complete integer check.

## 5. Reproduction and bridge ledger

Run `python -B h9_fibres.py` here. Independent subset Hamilton DP checks
all511 cut words of length<=8 and76 substituted small graphs. All729
formula comparisons for a,b<=8 pass. The full H9-word census through
base length10 matches (1); explicit graph isomorphism distinguishes the
two ordinary cores and identifies the apex switch.

Source: decorated H9 fibre. Target: exact branch counts. Map: A/B/C
substitution with exterior padding retained. Preserved: literal graph
construction and count formula. Lost by H alone: core, apex and two padding
coordinates. Sidecar: `(v,a,b)` or full word. Cheap hostile:107/73/17.
The30000 no-difference theorem says this *specific* branch fibre cannot
supply the requested translation; it does not block the separate pure-H
rooted constructor in the arithmetic lane or a growing-defect return.

## 6. The numerical relation107=73+2·17 and its first hostile

The three values also satisfy107=73+2·17, and their prime indices are
28,21,7 with28=21+7. Thus107 is outside the first25 primes. These finite
identities are exact but do not define a preserved graph operation.

Append the same b>=0 extra O letters to OIOIO, OOOI and OIII before
applying branch A. All three source counts remain9. The formulas above
give the first response minus the second minus twice the third as

    [172(4^b−1)−12b]/9.

It is zero at b0 and positive at every b>=1, since4^b−1>=3b. At the
first padding the responses are401,243,51, giving discrepancy56.
So even neutral padding of the actual sources destroys the numerical
relation; the cut word remains necessary.

For arithmetic translation put y_p=M+p. The correct affine identity is
`y_107=y_73+2y_17−2M`. At M30000 the correction is60000. Differences
are preserved by translation, while a linear relation whose coefficients
sum to−2 acquires that nonzero constant.
