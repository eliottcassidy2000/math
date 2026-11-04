# Independent audit: coloured endpoint-forest transport

2026-11-03. **INDEPENDENT ANALYTIC / SOURCE / REPLAY AUDIT PASS.**
I read the final proof and the complete standalone source, then independently
ran both normal and optimized modes. All 216,070 always-active gates pass;
both 474-byte outputs agree with the frozen producer output. No correction
remains. Primary status promotion is left to the parent.

## 1. Exact state and all-order birth mechanism

The objects are actual unordered spanning covers by nonempty directed paths.
The anchored signature retains paired endpoints, not separate start/end
marginals. An active physical name is present at most once in the signature,
except for the two ends of its own singleton. Its absence records an interior
vertex. Fully anonymous paths retain their ordered pair of endpoint colours
as a monomial variable. This state correctly allows a named vertex to become
interior without deleting it or losing the physical sidecar.

The signature recurrence is correct: the least active name is absent, a
singleton, paired with one of q anonymous colours in either direction, or
paired with one of the other active names in either direction. These disjoint
possibilities give the displayed `(2+2q)` and `2(p-1)` coefficients. The stated
coefficient-count bound is an upper bound, not a reachability or rank claim.

Deleting the newborn from an actual cover gives a unique old cover and zero,
one or two selected old components. In the two-component case they must be
different. Colour-uniform incidence makes this selection depend only on the
retained endpoints. The source decrements the anonymous exponent after the
first selection, so selecting one type twice contributes `k(k-1)` rather than
`k^2`; selecting the same anchored component twice is separately rejected.
The resulting endpoints give precisely

```
sum_(b in I,c notin I) sum_(a,d) z_ad * partial_ab * partial_cd.
```

Differentiation precedes multiplication. A new vertex may become interior,
leaving its name absent from the resulting anchored signature. For one
anonymous colour the mixed term vanishes because that colour cannot lie on
both sides of the cut. Thus the old first-order grammar is recovered without
assuming its false extension to several differently oriented populations.

The inverse deletion argument proves the rule for every finite palette,
active-name boundary and component degree. Finite tests are not being used
to extrapolate this all-order statement.

## 2. Other operations and exact losses

Forgetting is a pushforward on paired endpoints; a newly fully anonymous pair
adds its correct variable. Whole-colour maps act on both ends and exponent
indices, including noninjective merges. They cannot split an occupied colour.
Converse reverses both anchored and anonymous endpoint pairs.

Named joins select the unique component ending at one active name and the
different component starting at the other. This pays the cycle exclusion.
The signed flip formula has the actual arc orientation as a sidecar. Its cut
argument, join commutativity and square/opposite-join zero identities remain
valid after retaining endpoint colours. In particular, mixed flip terms are
not dropped, and a reversed second flip uses the reversed sign. The proof
does not declare arbitrary signed intermediate states realizable.

The one-colour alias has actual outputs H=1 and H=3 after a forbidden split;
two colours distinguish its sources before that operation. The same-component
hostile and the paired-colour hostile identify two different information
losses: falling multiplicity and endpoint pairing. Both are correctly kept.
The fact that the same two paired signatures occur as covers of a transitive
four-vertex graph does not identify the whole graph's state with either one
of those single-cover signatures.

The literal graph/word remains necessary for physical adjacency, old-name
identity and native Hamilton-complement legality. No operation resurrects a
forgotten physical name or permits nonuniform future incidence within one
colour.

## 3. Named defects, strongness and the long control

For nonempty graphs, the one-component readout is H. At one or two components
there is a unique assignment to the ordered slots when the first component
must end at a fixed active name r; this gives F_r. The joint observer uses
different components. Consequently the actual changing defect is F_r-H,
and the inherited bad-join identity gives its positive interpretation.

The proof correctly distinguishes two facts: F_r-H remains readable even if
historical r-incidences are mixed inside a colour, while the incoming-colour
sum interpretation requires whole-colour incidence. A forgotten root loses
its named observer. No scalar reset is inferred.

As an additional independent small control, literal permutation/cut counting
on the two strong tournaments with rows `(6,4,8,3)` and `(6,12,8,1)`, retaining
root 3, gives H=5 for both but F_3=5 and 7 respectively. They differ by one
active arc reversal and have defects zero and two. Thus retaining H alone
already fails at four vertices, even among strong actual sources.

The SCC birth guard is exact: the newborn must reach the first old SCC and
be reachable from the last. Its proof covers arbitrary old SCC count, the
one-SCC case, and the empty initialization separately. It is properly not
called a native-complement guard.

The displayed infinite control is literally C3({0},{1},TT_(k+1)). Alternating
the new colour preserves the transitive order inside the third block. The
Hamilton paths consist of the `2^(k+1)` cuts around the adjacent pair 0,1 and
the additional order `1,chain,0`. The endpoint deletion calculation gives
F_0=2^(k+2)-1. Hence the same retained physical root has the stated changing
defect for every length; no two-digit closure property is claimed.

## 4. Complexity contract

The fixed-p,q corollary is valid with its stated input contract: initialization
is empty, a fixed explicit seed, or a supplied exact initial profile. It does
not promise to compute an arbitrary tournament's starting profile cheaply.
Every operation has a bounded number of local choices per monomial for fixed
p,q, and polynomially many coefficients in the final order. The coefficient
bound `N! 2^(N-1)` overcounts all unordered covers, giving polynomial bit size.
Each named join selects at most one pair from each physical cover, so the
three-term signed flip has absolute coefficient sum at most three times the
old total. Positive birth partial sums are bounded by the resulting cover
count. The graph/name/colour sidecar and SCC test also cost polynomial time.
Fixed p,q bounds include the active name introduced during a birth, before
any subsequent forgetting. This is exact word evaluation, not word discovery,
graph recognition, global spectrum decision or return certification.

## 5. Source, independent replays and pins

The direct source evaluator is independent of the transport operators: it
chooses the component containing the least remaining vertex and enumerates
its directed paths by literal permutations. This counts every unordered
cover once, retaining internal path multiplicity. A separate endpoint subset
DP checks H, ordered two-slot covers, all retained F and joint observers.
The complete profile, not only its small-degree readout, is compared after
every finite operation in the declared universe.

The raw universe sizes are independently consistent with the declared loops:
76 labelled graphs of orders 0–4; 5,422 active/red/blue assignments; 67,732
typed births, including duplicate incidence descriptions when a colour is
empty; 3,530 actual active flips; 14,282 forget operations; and 5,422 merges
and converses. No isomorphism filter is hidden in those counts. Full-profile
checks extend to order five through births. The long word's full profile is
independently recounted through order six, its readouts by DP through order
twelve, and its explicit formula is proved for every length.

Independent commands, from the worktree root:

```
python3 04-computation/multiport_grammar_20261103/transport_colours.py > /tmp/transport_colours_audit_normal.out
python3 -O 04-computation/multiport_grammar_20261103/transport_colours.py > /tmp/transport_colours_audit_optimized.out
```

Both outputs are byte-identical to `verification.out`, with **216,070 gates**.

| Artifact | Bytes | SHA256 |
|---|---:|---|
| Primary proof, pre-promotion | 20,855 | `fafb04a1ac43a8dbd15b6d6f3818c07357e2cc2ec289176e4488c5a5b07020e2` |
| `transport_colours.py` | 13,609 | `860928f457f1227d017d677f1412b2e9da017e2604f8ee7d4a6d513f618ee77a` |
| Producer and both independent outputs | 474 | `d415cc1fc86727fb986a7a094c5701f99629372d750deaf8c325d36782b34d7e` |

The semantic receipt is
`c9601f087a73eb440bf18e940e31b0f9ad2e25dfdb255344a412f6388b8afc2a`.
All producer files were preserved during this audit. **PASS.** The accepted
result is an unbounded typed evaluator with its changing inner defects; actual
reusable two-branch closure and global spectrum coverage remain separate open
obligations.
