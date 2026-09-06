# Connected decoder edges force bounded gcd loss

**Status: PROVED, independently audited; FINITE-EXACT sharp safe controls.**
The [independent referee](overnight12_20260906_lrc_decoder_descent_audit.md)
accepts the full proof without repair and supplies61,138 additional live gates
with byte-identical normal/optimized output. The elementary graph lemma below
is unconditional. Its LRC
consequences inherit the actual THM-3818 decoder graph and the proved, audited
joint-shadow gcd restrictions. They do not close LRC(14), construct a failure,
or turn a five-label gcd identity into a three-label bounded relation.

## 1. Inheritance and the retained coordinate

The inherited graph is **THM-3818**,
`01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`,
Sections 1 and 6. Its primitive edge ratios a<b satisfy a+b≤356, gcd(a,b)=1,
and every prime in a+b is 2 modulo3 with exponent at most two. There are
5,855 such ratios; the largest oriented coefficient is355. Its Section6.1
already uses products along a connected tree to bound primitive component
heights. The present contribution starts at an arbitrary selected subset and
retains its exact gcd loss when one actual boundary edge is crossed.

The closest incoming supplier is the proved and independently audited
`05-knowledge/results/lrc14_joint_shadow_empty_core_next_sep06.md`, together
with its exact JSON profiles. It improves the earlier recursive-sieve caps
96,32 to90,30. Under a primitive thirteen-speed strict failure, every subset
of size12,11,10,9,8,7 has gcd at most1,2,4,9,30,90, respectively. Its set of
seven-body gcds has42 entries and maximum90. This note uses that completed
proof, not its namespace reservation.

The canonical hostile is the primorial component from **THM-4052**,
`01-canon/theorems/THM-4052-lrc14-affine-component-width-escape-cones.md`:
gcd one and small primitive pair heights do not imply a small component maximum
or a coprime pair. It is already excluded by the incoming large-subset gcd cuts.
The corrected near miss is the virtual pair split whose purported components
were crossed by a bounded three-label row. Thus the actual all-scale graph and
W=V_dec must be retained. The least-used coordinate is the oriented coefficient
on the edge that leaves a selected subset, rather than a global pair-height cap.

The six-concept board is: connected decoder edges; prime-adic gcd valuations;
hereditary shadow profiles; physical versus primitive scale; exact coefficient
boxes; and relation support. Anchor: genuine LRC entry. Niche: the
[audited signed-box criterion](overnight12_20260906_lrc_gcd_semigroup_audit.md).
Wildcard: reuse a graph boundary to control valuation descent.

The map sends S to S∪{v} along an actual edge u–v with u∈S. It preserves
divisibility of the enlarged gcd into the old gcd, while losing the quotient.
The retained sidecar is the edge orientation and the exact loss multiplier.
The decisive hostile is to ask whether merely knowing gcd(all coordinates)=1
supplies a suitably small edge; it does not. The actual graph supplies it.

## 2. Exact boundary-edge lemma

Let positive integers label a connected graph. Suppose every oriented edge
u–v has

    u/gcd(u,v) ≤ B.

For a nonempty proper vertex subset S, put d=gcd{w:w∈S}. Choose any edge
u–v leaving S, with u∈S and v∉S. Such an edge exists even if S itself is
disconnected. Then

    d' = gcd(d,v),       ℓ = d/d',
    ℓ divides u/gcd(u,v),       d=ℓd'.                 (1)

**Proof.** At a prime, let α=v_p(d), β=v_p(u), γ=v_p(v). Since d|u,
α≤β. The valuation of ℓ is max(α−γ,0), at most max(β−γ,0), the valuation
of u/gcd(u,v). This holds at every prime, proving divisibility. In particular
ℓ≤B. There is no claim that an arbitrary newly adjoined vertex obeys this bound.

Iterate actual outgoing edges until a k-element subset has grown to r elements.
If every r-element subset has gcd at most M, then

    gcd(S) ≤ M B^(r−k).                               (2)

More precisely, along any chosen growth path,

    gcd(S)=gcd(T) ∏ℓ_i,   ℓ_i divides a_i,             (3)

where T has r elements and a_i is the oriented primitive coefficient on that
step's actual edge. Formula(3) is stronger than the scalar envelope(2).
It retains divisibility, the endpoint subset, and an explicitly realizable
growth path. No bounded-coefficient Bézout representation is asserted.

## 3. Actual eleven-core consequence

In a hypothetical primitive strict failure with the actual THM-3818 11+2
decoder partition, apply (2) to the eleven **physical** core coordinates with
B=355, r=7, M=90. Every core subset has the following bounds:

| Subset size | Physical gcd upper bound |
|---|---:|
| 6 | 31,950 |
| 5 | 11,342,250 |
| 4 | 4,026,498,750 |

For general k≤7 the bound is90·355^(7−k). If the physical core is tV with
gcd(V)=1, its primitive subset gcd bounds are divided by t. The incoming
eleven-subset cap additionally gives t≤2 in a strict failure.

The new six-subset bound restores a finite clock universe in this actual
connected branch. The unrestricted six-body/seven-tail shadow problem has no
such scalar cutoff. Let C₇ be the incoming42-element seven-body gcd set.
From (1), a physical six-subset gcd must lie in

    C₆ = {ℓc : c∈C₇, ℓ divides an oriented atlas coefficient}.   (4)

The exact atlas contains every integer1,…,355 as an oriented coefficient.
Thus (4) equals {ℓc:1≤ℓ≤355,c∈C₇}, with **6,121 distinct values** and
maximum31,950. This is a necessary arithmetic universe, not a sufficient
entry condition or a classification of actual component phases.

At Q=91^6, the four-subset bound is belowQ. The five-subset bound is below
H=floor(Q/(42·177))=76,388,115 from the signed-box endpoint criterion. This
does **not** trigger that criterion: it requires a gcd on a selected pair,
and an identity using five labels does not supply a three-label relation.
Both support size and coefficient height remain separate debts.

## 4. Sharpness within the retained arithmetic and entry conditions

The following three finite safe examples attain the bounds. Set

    U₄=(1,3,4,9),  U₅=(1,3,4,9,10),  U₆=(1,3,4,9,10,16),
    D_k=90·355^(7−k),
    W_k = D_k U_k ∪ {90·355^j : 0≤j≤6−k} ∪ {2,4,6,12},
    K=max W_k,   g=2QK+1,   full row = W_k ∪ {g,3g}.             (5)

Each W_k has eleven distinct coordinates and gcd2. Thus the physical core
has t=2 and primitive shape V=W_k/2, while gcd(t,g)=1. The top k coordinates
have gcd exactlyD_k. The full thirteen-coordinate row is primitive, distinct,
and inside the actual sum≤Q² box. These are exact finite statements checked
below; the displayed largest physical maximum is36,238,488,750.

The graph is connected within W_k: its top factors connect to1 through ratios
1:3,1:4,1:9,1:10,1:16 as needed, whose sums4,5,10,11,17 are allowed.
Successive scale levels connect through1:355, since356=2²·89 is allowed.
The lowest scale90 connects to12 through primitive ratio2:15 with sum17.
The remaining small vertices2,4,6,12 connect through ratios1:3 and2:3.
The outside pair has ratio1:3. All these are actual decoder edges.

There are no crossing decoder edges or any other bounded crossing of support
at most three. A nonzero outside partial sum is an integer multiple ofg, of
absolute value at leastg>2QK; at most two core terms cannot cancel it. For
two outside labels and one core label, a zero outside sum cannot cancel a
nonzero core term either. Internal connected decoder rows already span rank
eleven, so the full bounded relation space is exactly W=V_dec. This checks
actual entry, not just its necessary numerical filters.

The controls also pass **every retained joint-shadow profile**, for every
subset of the full row of sizes7,…,12: compute its body gcd c and the sorted
word gcd(c,outside coordinate), then check the exact inherited profile table.
In particular this includes the joint exclusions and the restricted c=90
words, not only the scalar caps. The maxima for sizes7,…,12 are respectively
90,6,6,3,2,1. Passing these profiles is only passing the stated arithmetic
relaxation; no implication of unsafe behavior is made.

Indeed x=1/7 is literally strictly safe for each full row. Every core integer
in (5) is nonzero modulo7, Q is divisible by7, and g≡1 modulo7. Every speed
therefore has clearance at least1/7>1/14. The normalized core contains a unit
and is already covered by earlier closure results. The examples prove sharpness
of (2) under the retained graph/profile conditions; they are not remaining
unsafe candidates and do not prove a sharp boundary for LRC.

The first draft used base96 and later base90 with small vertices1,3,4,9. The
new incoming profile table rejected the latter at a c=90 complement word
(1,1,2,3,3,9), despite all scalar caps passing. Replacing these vertices by
2,4,6,12 repairs the exact profile condition and exposes why the word sidecar
must be kept. No rejected control is treated as surviving current constraints.

## 5. Reproduction, provenance, and remaining question

[Standalone verifier](../../04-computation/overnight12_20260906_lrc_decoder_descent.py),
[output](overnight12_20260906_lrc_decoder_descent.out), and
[frozen incoming profiles](../../04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json):

```
python3 -B 04-computation/overnight12_20260906_lrc_decoder_descent.py
python3 -B -O 04-computation/overnight12_20260906_lrc_decoder_descent.py
```

Both runs pass **221,166 always-active gates**, with identical LF outputs.
The verifier regenerates the atlas, checks the elementary divisibility lemma
for u,v≤180 and every divisor d|u, constructs all three actual-entry controls,
checks every full profile at sizes7,…,12, and retains actual outgoing growth
paths from every core subset of sizes4,5,6 to size7. It does not rerun the
incoming shadow classifier and does not use the finite bank as a universal proof.

The profile bytes were retrieved read-only from
`origin/main:05-knowledge/results/lrc14_joint_shadow_empty_core_next_sep06.json`
after the parent fetched384e5eca4. They are pinned by SHA256
`935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f`.
Verifier source SHA256:
`a2e9d6052a77522b2cb9c09afaec4866ba4b92c503b72dffa0c402bb6f64b80e`.
Output and optimized-output SHA256:
`13e462e4b569f3929ff05ade52646aa95bb5577ce30a4ebd56043a9622cda466`.

The remaining question is whether the exact edge-loss paths and large-subset
profiles force a suitable **pair** or a genuine support-three crossing, with
all coefficient budgets retained. This note provides a finite six-subset
clock target and a sharp obstruction to improving its scalar envelope using
only the present graph and profile assumptions. No Git, navigation, theorem
namespace, or earlier frozen producer artifact was modified.

**Filing:** root integrated these independently audited artifacts in the twelfth
checkpoint. Reproduction commands are relative to the repository root.
