# Endpoint-589 q50 residual response cover: independent 19/95 audit

Status: **FINITE-EXACT, INDEPENDENT AUDIT PASS** for the fixed response
representation described below.  This is not an optimum computation, a
carrier exchange, a physical-entry theorem, or an LRC(14) proof.

## Result

Let `P` be the frozen thirty-label pool.  The inherited endpoint-589 failure
ledger contains 20,025 rank-nine bodies at `(50,589)` and 11 at `(96,589)`.
Apply the already frozen q96 peel masks

```text
0000b3a5  0220932c.
```

At `(50,589)` the first mask is inactive; the second is active and hits 891
failures.  The exact residual `U` therefore has 19,134 bodies and ordered
FNV `0e67c635fbc71d3b` (the hash adds original ordinal, then body).

For a rank-eight or rank-nine mask `R`, define a response edge

```text
E_R = {B in U : R is active at (50,589) and R intersect B is empty}.
```

Activity is the literal exact inequality

```text
63 * mu(G_((P\R) union {50,589})) - 4 >= 0.
```

The independent audit proves

```text
19 <= tau({E_R : |R| in {8,9}}) <= 95.                 (1)
```

Here `tau` is for the **19,134-body post-peel residual**.  For the unpeeled
20,025-body q50 failure family, adjoining the active peel mask `0220932c`
gives the immediate corollary

```text
19 <= tau(full q50 response family) <= 96.              (2)
```

The 95-mask artifact does not by itself cover the 891 peeled bodies, and
`0220932c` is not one of its masks.  Thus “95 covers all q50 failures” would
be a false restatement of the result.

### Lower certificate

The nineteen residual bodies are

```text
013c6401 071c6400 052c6401 07186401 031c5600
27046401 23165400 27146008 171a5000 15107401
35105401 07485401 31106601 0518f008 236c4001
055c4408 2504f400 23943001 14386401
```

The audit streams all

```text
C(30,8) + C(30,9) = 5,852,925 + 14,307,150
```

masks.  It reconstructs the wall geometry directly, finds 480,409 active
rank-eight and 4,112,383 active rank-nine masks (ordered active-stream FNV
`b79e52255c5b3522`), and verifies that every active response hits at most one
of the nineteen bodies.  Exactly 172 rank-eight and 17,049 rank-nine masks hit
one.  The nonzero hit-stream FNV is `430fdb51e2ee1fa1`.  Unit weight on each
body is therefore an integral dual of value nineteen.

### Upper certificate

The ordered 95-mask artifact is `input_cover95.csv`; all masks are distinct
and rank nine.  The verifier recomputes each literal activity margin and then
replays each declared greedy gain directly against the 19,134 bodies.  It
covers 19,134/19,134, with mask-stream FNV `fa44f9bfad76cfe7`.  Its smallest
activity surplus is 1,550,209,054,968 ticks at mask `0a624049`.

The complete ordered cover is:

```text
24811216 08a21059 1a128888 10a410d2 28422886 0031132c 04e12242 006095d0
0040a559 280a0a92 080a2296 1806aa80 04299124 0a820e81 0090855a 00a0b183
0070003f 04098750 29420a82 08d20249 0401023f 0260a889 0681a320 08070a89
00f03250 0c038889 04708924 26611050 18028a89 12140836 2000303f 08321825
06618308 30c4210a 10c83016 28229c40 0ea22082 108c0586 0c530888 02500d58
0c21122c 0ab28042 0a624049 0441a324 08322a81 30540416 00aaa442 02d09908
21808127 3e060880 190a0a82 22802616 300420d9 00b8212c 04510925 0a46c880
09822ac0 048b2089 0c8306c0 38860608 08623212 2c4b0880 00d8812c 034089a4
1282b082 24b0012c 0802403f 0a606890 08c88cc0 0ac11c08 1084a125 26210216
38241068 0288c132 10d48128 04414750 2014b016 00ec1092 03a01164 086260c2
12452a80 06019324 12982a02 0c8e2880 1382008b 24048329 2a421c04 2110832c
0220641e 04923083 2008329c 001ac213 01809552 030a0296 0421b124
```

### Exact local obstruction at the incumbent

The 95-mask cover is inclusion-minimal and, more strongly, admits no
two-for-one response exchange.  For each of its `C(95,2)=4,465` mask pairs,
remove the pair and collect the bodies left uncovered.  Any single replacement
must be disjoint from the union of those bodies.  The exact audit enumerates
all 40,468 rank-eight/rank-nine subsets of the resulting complements that can
possibly qualify and evaluates activity from the wall classes; none is active.
This proves only a local minimum under `2 -> 1` exchange, not global optimality.

## Independence and provenance

This audit imports neither the capacity agent's 82 MB signature atlas nor its
greedy implementation.  `q50_response_common.hpp` rebuilds all literal walls
and the 2,383 rank-at-most-nine failure classes.  The audit then evaluates
activity from those classes and response incidence from the raw disjointness
predicate.

The three packet inputs are linewise identical to the frozen capacity inputs.
Their original capacity-packet raw SHA-256 values were:

```text
endpoint589_failures.csv                 9cf04fe2bd9c361e68321d0ae8bad3f165fd5f5a8c8ed093589fa795868ce973
q50_core512_find19_solution_O3.csv       7b2405f3fb6584c1d62a9146157ee3bea440693db00064262069ca3fc4d0e3a1
q50_seeded_greedy_cover_O3.csv           0893d5ce51e73642ddb661909c37a9c8db092a318c83b0276e86456264192221
```

This packet normalizes text to LF and therefore freezes its own raw hashes in
`SHA256SUMS`.  O2/O3 exhaustive transcripts and normal/optimized local-exchange
transcripts are respectively byte-identical.

## Map and lost data

Source: the full 20,025-body q50 carrier-failure ledger.  Map: delete the 891
bodies hit by the active fixed peel mask, then associate to every active
rank-eight/rank-nine mask its disjoint residual bodies.  Preserved predicate:
response sufficiency.  Lost data: q96 obligations, inactive carrier masks,
simultaneous deletion/exchange behavior, higher-rank responses, full literal
mass assembled across different response masks, and physical-entry data.

That loss matters.  The separate direct literal-mass theorem already closes
all endpoint-589 carrier failures without responses.  Consequently `(1)` is
a theorem about the complexity of this sufficient-certificate representation,
not a measure of physical difficulty and not the endpoint-589 closure proof.
