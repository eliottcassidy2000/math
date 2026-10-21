# Recovered proof: grounding the actual optimal frontier-credit owners

Recovery date: 2026-10-21.

## Recovery status and evidence boundary

**This document is reconstructed from retained conversation and worker reports. It is not a byte-identical recovery of the lost packet, not a replacement for its missing source data, and not a newly executed verification.**

A workspace reset removed the local repository and the Session21 staging directories during an approval wait. The recovery inspection found neither `math/` nor `session21-owned-credit-children/`. No commit of this packet had been made according to the parent's recovery notice. This reconstruction does not imply that any historical receipt, input file, audit report, or staged deliverable is currently available.

The mathematical arguments below are preserved derivations. Every numerical result described as historical must be reproduced after the pinned inputs are recovered. The old hashes are identification targets, not proof that a file with that hash presently exists. Only this recovery document has been created by this recovery task; no remote write, route discovery, annulus enumeration, or mathematical replay has been attempted.

The lost packet's scope was:

- **FINITE-EXACT, historically verified:** all 233 scale-specific owner records in the frozen optimal N=4096 and N=8192 selections, representing 212 distinct source/pair addresses.
- **PROVED application of an inherited compiler:** each selected owner and its paid child has a strict fully odd ROOT word, extending along its guarded j+9k ray.
- **OPEN and unchanged:** uniform native 153N coverage, uniform payment of larger sibling banks, a uniform bound on the remaining R', and arbitrary-source Collatz entry.

No new general theorem, priority claim, negative interaction, or additional extracted credit was claimed.

## 1. Inputs, attribution, and what the bridge certified

### 1.1 Historical input locations

The local repository was `math/`, at revision

    5ce22d4ecb752c1a38df8bc739aa7751932755fd.

The frozen local source inputs were:

1. `math/05-knowledge/results/session18_native_annulus/frontier.json`.
2. `math/05-knowledge/results/native_cancellation_105_20261018g.json`, the frontier's original producer receipt.
3. `math/05-knowledge/results/session18_native_annulus/README.md`.
4. `math/05-knowledge/results/session18_native_phase_lifts/README.md`.
5. `math/05-knowledge/results/source_gap_bursts_20261020k/README.md`.
6. `session21-cross-scale-credits/receipt.json` and its `README.md`.
7. `session21-native-remainder/receipt.json` and its `PROOF.md`.

The incoming compiler and comparison receipt were read from `eliottcassidy2000/math` at immutable revision

    881a9903995884c02a03ffb45c68e9174f8b6c58,

from these paths:

- `05-knowledge/results/credit_children_20261021a/README.md`.
- `05-knowledge/results/credit_children_20261021a/check.py`.
- `05-knowledge/results/credit_children_20261021a/receipt.json`.

The exact paid-prefix/ROOT compiler, including authentication of the parent's first letter before replacing that letter, is inherited from this incoming packet. Its verified historical source URL was:

https://github.com/eliottcassidy2000/math/blob/881a9903995884c02a03ffb45c68e9174f8b6c58/05-knowledge/results/credit_children_20261021a/README.md

The primitive inverse, source-gap debt, first-arrival convention, and affine scale transport are earlier inherited mechanisms. The bridge changed the explicit certification target: it grounded the owners actually selected by our two optimal banks.

### 1.2 Historical overlap boundary

The incoming packet explicitly grounded a different 380-source bank. The lost checker compared exact pairs `(parent, source)`, not just method descriptions:

| Selection | Owner records | Already in incoming 380 | Not in incoming 380 |
|---|---:|---:|---:|
|N=4096|79|16|63|
|N=8192|154|42|112|

The two intersections have union 42. The two selections have union 212 distinct source/pair addresses. Therefore 170 addresses were newly certified **relative to that explicit incoming receipt**. Every overlapping ROOT word, paid child, and child ROOT suffix agreed exactly in the historical replay.

This is not a claim that those 170 integers had never previously been shown to converge. The separate Session21 grounded-child packet concerned a fixed 67-source Mersenne bank and a cap language, not this selected-owner universe. A historical check found zero overlap between its 67 original source integers and these 212 selected source integers. That comparison is background, not an input needed for the bridge proof.

No original bank receipt was edited. The bridge authenticated selected addresses, weights, both source rows, terminal, first-arrival depth, full native arrival clock, and strip count against both frozen banks.

## 2. Definitions and native source identity

Let

    B(x)=4x+1,
    U(x)=(3x+1)/2^{v2(3x+1)}.

A primitive source is a positive odd integer not congruent to 5 modulo 8. Literal Python booleans are not integers for this API.

The native primitive map P starts at primitive s and chooses

    a=2 if s=1 modulo 8,
    a=1 otherwise.

It forms z=(3s+1)/2^a, then repeatedly replaces z by (z-1)/4 while z=5 modulo 8. Let j be the number of these removals and u the resulting primitive endpoint. The native P-edge has endpoint u, edge depth 1, and integer division/braid clock a+2j. It is not the fully odd U-step count.

Since z=B^j(u),

    U(s)=B^j(u),
    s=[2^a(3u+1)4^j-(2^a+3)]/9.                    (2.1)

The inverse API must authenticate all of the following:

- u is a literal primitive positive integer.
- a is a literal integer in {1,2}.
- j is a literal nonnegative integer.
- The numerator in (2.1) is divisible by 9.
- The resulting s is primitive.
- Recomputing its actual native edge returns exactly `(u,a,j)`.

This is an exact source identity, not a residue-based convergence argument.

## 3. Compile all 627 parent ROOT words from the existing rank DAG

### 3.1 Frozen frontier obligations

The inherited frontier V contains 627 primitive sources, with maximum 63665. Its supplied rank map has domain exactly V, has rank(1)=0, and assigns positive integer rank to every other vertex. The lost checker required exact equality between the frontier and rank arrays in `frontier.json` and their original producer receipt.

For each v in V other than 1, independently recompute its native edge `(u,a,j)` and verify

    u in V,
    rank(u)<rank(v).

This is the only traversal used to construct the parent words. It is finite induction over an already certified DAG, not a new forward search for Collatz routes and not enumeration of an annulus.

### 3.2 Literal inherited splice

A strict ROOT word is an ordered list of the exact valuations on the fully odd orbit, stopping at the first 1. Set w(1) to the empty list.

Process vertices v other than 1 in increasing `(rank(v),v)` order. The rank check guarantees that w(u) is already available. Authenticate the entire parent word by replay **before** replacing its first letter.

If j=0, use

    w(v)=(a,w(u)).

If j>0 and u>1, write w(u)=(b,w_tail), where b=v2(3u+1), and use

    w(v)=(a,b+2j,w_tail).                            (3.1)

If j>0 and u=1, use

    w(v)=(a,2+2j).                                 (3.2)

The proof is the identity

    3B^j(u)+1=4^j(3u+1),

which gives exact next valuation b+2j and endpoint U(u). The u=1 case uses its actual valuation 2 but does not append a fictitious parent step after ROOT. For nontrivial v, the case j=0,u=1 cannot occur: branch2 gives v=1 and branch1 gives a nonintegral inverse.

### 3.3 Exact word replay and backward reconstruction

For `replay(n,word)`:

1. Require n to be a literal positive odd integer and word to be a literal list.
2. Initialize x=n and debt D=0.
3. For each letter a:
   - Reject if x is already 1.
   - Require a to be a literal positive integer.
   - Require a=v2(3x+1).
   - Replace x by `(3x+1)>>a`.
   - Update debt by the exact rational recurrence in Section 5.
   - Retain the endpoint and debt.
4. Require the final endpoint to be 1.
5. Reconstruct backwards from z=1. For each reversed letter a, require `(z<<a)-1` to be positive and divisible by 3, set z to that quotient, and require z odd.
6. Require the reconstructed z to equal n.

This checks the source, every valuation, strict stopping, and the full word. It prevents a forged parent first letter from being silently discarded by (3.1).

### 3.4 Historical parent results

The missing receipt historically contained one record for each of the 627 parents, including source, rank, source registers, ROOT word, ROOT step count, ROOT valuation sum, and actual native-to-ROOT edge/clock totals. Non-ROOT records additionally retained their native parent, branch and strip count.

Historical maxima:

- Frontier rank: 64.
- Actual native-to-ROOT edge count: 64.
- Strict fully odd ROOT word length: 65.
- Fully odd ROOT valuation sum: 113.

These numbers need rerunning from the recovered pinned frontier.

## 4. Explicit clock completion corollary

This is a derived identity of the inherited compiler, not a new general clock-loss theorem.

For a primitive v>1 whose strict P-chain reaches 1, let d_P(v) be the actual number of P edges and t_P(v) the sum of their native edge clocks. Then

    length(w(v))=d_P(v)+1,
    sum(w(v))=t_P(v)+2.                            (4.1)

Proof: the final native edge to 1 has j>0 and literal word `(a,2+2j)`, so its odd length is 2 rather than native depth 1, and its valuation sum exceeds native clock a+2j by 2. At each preceding edge:

- j=0 prepends a, adding one odd letter and native clock a.
- j>0 replaces the parent's b by `(a,b+2j)`, adding one odd letter and native clock a+2j.

Induction proves (4.1). ROOT v=1 is the explicit exception: both words and both clocks are zero.

The historical checker verified that every supplied frontier rank equals actual native depth, not merely an upper bound. Thus for a selected source first arriving in frontier terminal u at stopped depth d_s=1 and native clock t_s=a+2j,

    ROOT length = d_s+rank(u)+1,
    ROOT valuation sum = t_s+t_P(u)+2.              (4.2)

The selected parents were all greater than 1, although the compiler also tested the u=1 boundary. The certificate retained both the stopped frontier-arrival clock and the completed native-to-ROOT clock. Equation (4.2) does not make those clocks interchangeable.

## 5. Strict paid children of the actual selected owners

### 5.1 Ownership and first-arrival guards

Each selected record has old endpoint u and new endpoint s satisfying

    u in V,
    u<4N<=s<256N,
    s not in V,
    P(s)=u with clock a+2j,
    beta<0.                                        (5.1)

The old source has stopped route `(u,clock0,depth0)`; the new source has `(u,clock a+2j,depth1)`. Because s>u, the incoming compiler's payment guard holds.

The quadratic owner is `(u,s)`. Its newly certified paid child does not replace either source, its terminal, or its charge.

### 5.2 Exact debt and payment proof

The original-source debt starts at zero and evolves along actual odd valuation letters by

    D'=3D/2^a+max(0,1-4/2^a).                       (5.2)

All debt arithmetic must use exact rational numbers. Strict payment at endpoint x means

    x<s and D<s-x.

For branch a=2, the first endpoint z=(3s+1)/4 is strictly smaller than s for s>1, and its debt is zero. Thus step 1 is strictly paid.

For branch a=1, j=0 would give s=(2u-1)/3<u, contradicting (5.1). Therefore j>=1. The second endpoint is U(u), with

    b=v2(3u+1),
    D2=1-4/2^{b+2j}<1.

Already j>=1 gives

    s>=(8u+1)/3>(3u+1)/2>=U(u).

Both s and U(u) are odd, so their positive gap is at least 2. Hence D2<s-U(u), and step 2 is strictly paid.

After choosing the paid endpoint, independently replay its strict ROOT suffix from the paid child and invert that suffix backwards to the same child. This is distinct from carrying the original-source debt through the whole word; the suffix's ROOT validity does not require reassigning the original credit.

### 5.3 Historical exact counts

| N | Owner records | Distinct parent integers | Step1 paid | Step2 paid | Unchanged extracted mass |
|---|---:|---:|---:|---:|---:|
|4096|79|74|41|38|273|
|8192|154|131|74|80|474|

Across both scales:

- 233 scale-specific records.
- 212 distinct `(parent,source)` pair addresses.
- 212 distinct new source integers.
- 184 distinct parent integers.
- 21 pair addresses shared across the two scales.
- Among unique sources: 108 first-step and 104 second-step payments.
- Maximum selected-source strict ROOT length: 63.
- Maximum paid-child strict ROOT suffix length: 61.
- Selected base strip counts ranged from 2 through 8.
- 212 distinct `(u,a,j modulo9)` ray classes.

Counts 233 and 212 have different meanings. The shared addresses cannot supply extra independent credits by being mentioned at another scale. Their historical weight comparison was:

- 17 shared pairs had weights (-3,-3).
- 2 had (-3,-2).
- 1 had (-6,-2).
- 1 had (-3,-1).

The two entries in each weight pair refer to N4096 and N8192. Therefore even the weight of a shared physical address can depend on the scale.

### 5.4 Preserved example

The selected N4096 owner `(u,s)=(95,130161)` had branch2, j=5, native first-arrival terminal95, depth1, and clock12. Its first paid odd child was97621, not95. Its strict odd ROOT word had39 steps and valuation sum79, with paid-child suffix length38 and valuation sum77. The interaction remained -6.

Historical source rows were:

    old: [95,5,25,95,0,0,1,-1,3]
    new: [130161,6,31,95,12,1,-1,1,3]

The row schema is given in Section 7. Its payment debt was0 and payment margin32540. These example values are retained context, not freshly replayed values.

## 6. Exactly the guarded j+9k rays

For a selected base record, keep u,a and its supplied parent word fixed. Replace j by j+9k for literal integer k>=0. Set

    F=64^{3k},
    N_k=NF.

Then

    s_k=F s+(2^a+3)(F-1)/9.                        (6.1)

This follows directly from the full-denominator inverse (2.1). The same primitive parent is reached by an actual edge. The annular guard persists:

- Lower bound: s>=4N gives s_k>=4NF.
- Upper bound: s is odd, so s<=256N-1; since `(2^a+3)/9<1`, formula (6.1) gives s_k<256NF.
- Frontier exclusion: k=0 is a checked base guard; for k>=1 the source exceeds the frozen frontier maximum63665.

The native first-arrival depth remains1. The full native edge clock is

    a+2j+18k.                                     (6.2)

The strict fully odd word is obtained by the same splice. Its second valuation and total valuation sum increase by18k, while its odd step count is unchanged. Its prefix remains paid in step1 for branch2 or step2 for branch1. Clock completion (4.1) continues to hold.

The inherited one-scale source transports are

    s ->64s+35 for branch1,
    s ->64s+49 for branch2.

They preserve the native annular band. Three such lifts restore the new-source charge phase; the fixed old-source charge is unchanged. Thus each particular selected interaction weight is unchanged on (6.1), preserving the chosen finite masses273 and474 along their respective three-scale rays.

**This does not imply that those finite masses pay the enlarged sibling banks at all lifted scales.** It also gives no all-scale bound on R'. The rays of a fixed finite parent bank do not provide a number of owners proportional to N within one annulus.

### 6.1 Source and child registers

Retain `(mod15,mod35,mod105)` for the source, parent and paid child. Require

    r15 modulo5 = r35 modulo5,
    (r15+15*((r35-r15) modulo7)) modulo105 = n modulo105.

This is a compatible mod105 state, not an independent 525-state Cartesian product.

For a lifted source, compute the numerator

    2^a(3u+1)4^{j+9k}-(2^a+3)

modulo `9m`, then divide the resulting divisible-by-9 representative by9, for m=15,35,105. Reducing modulo m before this division is invalid.

Three-scale charge restoration does not necessarily restore these registers. For T_c(s)=64s+c, the inherited exact joint-register period is3 on its special modulo5 stratum and6 otherwise. Under three-scale lifting those periods become1 and2 respectively. Actual integer clocks still increase even when registers return.

The historical checker tested k=1,2,10 for every selected record, giving699 lift records with source sizes up to201 bits. At N4096 the largest probe had200 bits; at N8192 it had201 bits. The all-k conclusion follows from the algebra above, not extrapolation from those probes.

## 7. Exact cross-authentication and charge algorithm

### 7.1 Charge and interaction definitions

Use the native palettes

    L=((1,2),(-2,-1),(1,-1)),
    U=((0,1),(-1,-1),(1,0)).

For primitive s, compute A_N(s) by summing over its B-chain points m<4N:

- Add `L[m modulo3]` if m<2N.
- Otherwise add `U[m modulo3]`.
- Replace m by4m+1 and continue until m>=4N.

Use

    Q(x,y)=x^2-xy+y^2,
    beta((x,y),(a,b))=2xa-xb-ya+2yb.

For each owner, the negative interaction must be recomputed as

    beta(A_N(u),(-1)^a A_(64N)(s)).                  (7.1)

The old charge must also satisfy A_(64N)(u)=A_N(u).

### 7.2 Both source-row receipts

The source-row schema is

    [source,mod15,mod35,terminal,full_native_clock,P_depth,
     charge_coordinate_1,charge_coordinate_2,Q(charge)].

The cross-scale receipt contained `cases`, with each case's `N`,
`selected_count`, `extracted_credit`, `selected_credits`, and other prior
bank data. Each selected credit retained `old`, `new`, `branch`,
`strippings`, `beta`, `old_row`, and `new_row`.

The remainder receipt contained `cases`, each with `N`,
`selected_owners`, `retained_owner_source_rows`, and the original ledger.
Each selected owner retained `old`, `new`, `beta`, `terminal`,
`old_clock`, `old_depth`, `new_clock`, `new_depth`, and `strippings`.

For N=4096 and8192, in that order, the reconstructed checker must:

1. Require exactly79 and154 selected records, respectively.
2. Require that all selected `(old,new)` addresses are distinct within that case.
3. Compare sorted `(old,new,beta)` triples between both receipts.
4. Require negative mass273 and474, respectively, equal to the cross-scale case's recorded extraction.
5. Recompute each old row as

       [u,u%15,u%35,u,0,0,*A_N(u),Q(A_N(u))].

6. Recompute each new row as

       [s,s%15,s%35,u,a+2j,1,*A_(64N)(s),Q(A_(64N)(s))].

7. Require each recomputed row to equal both the cross-scale row and the corresponding retained row in the remainder receipt.
8. Require the remainder-owner terminal, depths, clocks and strip count to equal the literal edge data.
9. Check (7.1), strict negativity, source bounds, frontier membership/exclusion, and source identity from the inverse compiler.
10. Compile the ROOT word and paid child; retain the result under the original scale-specific owner address.

The bridge does not need to rebuild the full sibling bank or enumerate all sources. Its input is the already selected frozen bank. Prior optimality and remainder claims are inherited from the pinned source packets, not reproved by this small checker.

### 7.3 Lift records and historical certificate structure

Each base owner certificate retained at least:

- N, source, parent, branch, strips, beta.
- Source, parent and paid-child registers modulo15,35,105.
- Native edge clock, native first-arrival terminal/depth/clock.
- Completed native-to-ROOT edge count and clock.
- Payment step count, paid child, exact rational debt, exact payment margin, payment valuation sum.
- Full ROOT word, ROOT step count, ROOT valuation sum.
- Paid-child ROOT suffix, its step count and valuation sum.
- Both authenticated source rows.
- Whether the exact pair was already explicitly grounded in the incoming380.
- Lift probes for k=1,2,10, retaining the corresponding derived fields; full lifted words were reconstructed during checking rather than duplicated in every stored probe.

For each lift, independently verify the inverse, strict ROOT replay, paid suffix, affine source identity, annular/frontier guard, native and odd clocks, unchanged odd length, unchanged charges and beta, and full-denominator residues.

The historical top-level receipt included:

    summary
    input_sha256_lf
    parent_certificates
    owner_certificates
    cross_scale_shared_pairs
    incoming_grounded_overlap_pairs

The precise historical JSON byte formatting is not guaranteed by this reconstruction. The lost producer used sorted keys, indentation2, and a final LF for its generated receipt.

## 8. Hostile controls and required API discipline

All acceptance checks used explicit exceptions, not Python assertions, so optimized mode did not disable verification. Literal integer tests used `type(x) is int`, rejecting bool aliases.

Historical negative controls rejected:

1. Boolean branch: `compile_paid(3,True,8,w(3))`.
2. Boolean strip: `compile_paid(3,1,True,w(3))`.
3. Nonintegral strip type: `compile_paid(3,1,8.0,w(3))`.
4. Forged discarded parent first letter: `compile_paid(3,1,8,[999,4])`.
5. ROOT-padded parent: `compile_paid(1,1,1,[2])`.
6. A legal inverse that violates the payment guard: `compile_paid(11,1,0,w(11))`; its source7 is below parent11.
7. Empty non-ROOT word: `compile_paid(3,1,8,[])`.
8. Nonprimitive parent: `compile_paid(5,1,0,[4])`.
9. Negative lift index in the full-denominator register reader.

Positive boundary controls included:

- Empty ROOT parent: `(u,a,j)=(1,1,1)` gives source3 and strict word `[1,4]`.
- Legal zero-strip branch2: `(u,a,j)=(7,2,0)` gives source9 and paid child7.

The historical independent reviewer reported ten hostile mutations in its own distinct audit. That is not a contradiction with the producer's nine named rejection controls; the two verifiers had different test programs.

## 9. What remains unchanged in the energy ledger

The prior decomposition remains

    I=R'+(J+C).

For these fixed selections, C is0 at depth0 and equals -273 or -474 after depth1. Extracting a negative interaction from the remainder adds its positive magnitude to R'. Grounding its source or child does not create a second negative term and does not improve the physical energy by itself.

Thus the bridge changes no original pair weight, no source diagonal, no arrival cut, no extraction mass, and no prior remainder bound. It proves neither all-scale availability of a growing parent supply nor universal entry.

The later independent audit historically compared514 ledger cuts, corresponding to257 cuts at each of the two scales, and found the source ledgers unchanged. Those files and that audit output are currently missing and require recovery or rerunning.

## 10. Historical verification and known lost-output hashes

### 10.1 Producer verification history

The final producer checker passed normally and under Python `-O`, with identical summary output, at854236 production checks. Verification against the saved receipt adds its separate immutable-equality gate. This count is historical; no present claim is made that the program or receipt is still available.

The final stable producer hashes reported before the reset were:

| Lost file | Historical SHA256 |
|---|---|
|`check.py`|`dfbd2729e2004def7817088557897bfed578a804aa34bb883613320b725c2c4e`|
|`receipt.json`|`9199c63a9668f2ed5332f1d5cf8fba4541f8244890e927bff1c1dbadb33a0dbc`|
|`PROOF.md`|`6caa72e00a09638f1d38134310151fe6c7ea25df53dc39a494480e7b458eda71`|
|`pins.json`|`2c3b1ca77f325e309d94ed0b1a23ae43267beb1800d47d68812ac7765a7e9303`|
|`MANIFEST.json`|`59ddd2b356ab4184aa50f42eaab83bf3654f2f719fb2e664e0b376b0ef6cf948`|

The manifest remained labeled `READY FOR INDEPENDENT AUDIT`; no optional prose/status edit was made after the auditor requested freezing. The separate PASS report superseded that label.

Earlier, superseded receipt hashes appeared during construction:

- `49a847de83ea9436bd76f0cbe9c6a035833ccc79b1a433ab2ca5caec8f19c114`.
- `3744d77931055794d2c57b2e0f89bf9f56d0628d6bafc066f918877bb5c40cc5`.

They must not be mistaken for the final clock-corollary receipt beginning9199c63a.

The lost deliverable set was:

    PROOF.md
    README.md
    check.py
    pins.json
    receipt.json
    normal.out
    optimized.out
    MANIFEST.json
    sources/incoming-README.md
    sources/incoming-check.py
    sources/incoming-receipt.json
    sources/incoming-reviewed-pins.json

Hashes for files not listed above are not recoverable from this worker's retained context and must not be invented.

### 10.2 Historical independent audit report

The independent reviewer reported **PASS**, with805885 checks in matching normal and optimized runs. Its stated coverage was:

- 627 parent records.
- 233 owner records.
- 699 lifted records.
- 514 unchanged ledger cuts.
- 12 source pins.
- 10 hostile mutations.

It explicitly accepted the optional clock induction and requested no mathematical correction. It then froze the producer packet unchanged. No hash for that separate review packet was retained in this worker's conversation, and its existence cannot currently be asserted.

## 11. Historical source hash targets

All SHA256 entries below refer to LF-normalized bytes. These are exact retained identifiers for reacquisition and validation, not assurances of present availability.

| Role | Historical SHA256(LF) |
|---|---|
|Frontier JSON|`391f6d3380bcae4afb9a1bbd9522b74b1e9006012fa10a92be5bf0bbf938754f`|
|Original frontier producer JSON|`d2d4d152593a03eb653d1fd6fecfbc13f573a6eb74e1310f59b7a6cb78b01721`|
|Native-annulus proof|`f5edff0fe9009171f00b23e9bb38f5732eac72b682e58f9ae38f05e8fbff16fb`|
|Native phase-lift proof|`20d2cf37ed149345e259b1107509b7a9480de709461bad326df57b396f1b61ff`|
|Source-gap burst/debt proof|`473e64947de2f06badd7bbe86cb8d5f43f7f8fba8103bd4639b88fd7f09f9c2b`|
|Cross-scale credit receipt|`d2f996a98d6a12b8cc203223cc62f8bc18ec97b105cf9246f21ef612ebde0dd2`|
|Cross-scale credit proof|`33006feb994515cc1592e9d4b01566172e6d9654e87039b4931239ae78a2bc45`|
|Native remainder receipt|`0d4d647606269d7a89e2c1d2c12a56073934ff0a8b9f6317b04e8e997366ad69`|
|Native remainder proof|`e13e3f9a22947ead26bc3e2a25473041589a52fe48448cb0b6a6dc02117e3611`|
|Incoming compiler proof|`bcfc629a5f2d7083045aa254e0af42c20ddd04f3b1e4c85f4152c8f7ba6764b3`|
|Incoming compiler Python|`a806a25e5f619d4bb61c7670c61f3836f7338e15767af75d4946a6e7793751ab`|
|Incoming380-source receipt|`8a04dbcb81a5d5f2755d2549df956a7836b9a8e8495618ae58734a8d3b93d082`|

The incoming files' Git blob identifiers at881a990 were:

| Incoming file | Git blob SHA1 |
|---|---|
|README.md|`8c12fe9a6dc0725d9d96a860c95b54d8e4714c2f`|
|check.py|`4f14fc7e0f1222e6c82d8cdaa2b64acc4018f395`|
|receipt.json|`6d856336ca474b7609f2b74e3e2ee376f04913d8`|

The receipt's historically fetched content length was489372 bytes, ending in a single LF. The initial staging write temporarily added an extra LF; that was corrected before the above pins were computed. Git blob verification used SHA1 of `b'blob '+decimal_length+b'\0'+raw_bytes`. Preserve exact content when validating the blob, and separately apply the declared LF-normalization rule for SHA256.

## 12. Recovery and rerun plan

The immediate blocker is missing input bytes, not a failed mathematical argument. No current ownership certificate can be certified solely by remembering its count or hash.

After the parent recovers the repository or uploaded tree:

1. Recover the original source inputs and producer files where possible. Authenticate the repository revision, all12 input hashes, and incoming immutable blobs.
2. Prefer a byte-identical recovery of the final producer packet. Its final receipt must match9199c63a..., not either superseded draft hash.
3. If only source inputs survive, reconstruct the checker from Sections2–8. Treat it as a new implementation unless it exactly matches the historical check.py hash.
4. Rerun the bounded rank-DAG compilation and authenticate every selected owner against both frozen banks. This requires only627 parents and233 owner records, not a full annulus census.
5. Replay every base ROOT word, paid prefix and child suffix; independently invert each word. Recompute source registers, charges, interactions, native clocks and clock-completion identities.
6. Recheck all699 lift probes and prove their all-k continuation using (6.1), rather than relying on probes alone.
7. Recompute all overlaps and distinct counts, including shared addresses and incoming380 comparison. Preserve their different meanings.
8. Rerun all hostile and positive boundary controls with normal Python and `-O`. Compare deterministic outputs.
9. If recovered original producer/receipt bytes match their historical hashes, use immutable-receipt verification. If reconstruction changes the program or byte serialization, issue a fresh clearly versioned receipt and new audit rather than claiming the old hash.
10. Rebuild the independent audit, including the514 unchanged ledger cuts, and obtain a fresh PASS before any publication claim.

The historical default commands were:

    python -B check.py
    python -B -O check.py

The lost CLI supported `--repo`, `--cross-receipt`, `--remainder-receipt`, and `--write`. Default input locations were adjacent `math/`, `session21-cross-scale-credits/`, and `session21-native-remainder/` directories. `--write` wrote only this packet's derived receipt. It did not mutate repository files or original receipts.

No step in this recovery document authorizes remote publication, data substitution, fresh unbounded route search, or a new full-source enumeration. Those actions remain with the parent's current task and authorization decisions.
