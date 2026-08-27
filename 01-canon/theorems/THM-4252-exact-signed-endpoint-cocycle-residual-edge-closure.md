---
id: THM-4252
title: "Exact signed endpoint-cocycle residual-edge closure"
status: >
  PROVED RELATIVE TO THM-4228/4233/4242 + FINITE-EXACT + INDEPENDENTLY
  AUDITED FIXED-THREE-EDGE CLOSURE. For every nine-body in the displayed
  thirty-label pool, adjoining either pair (466,699), (616,769), or (721,769)
  leaves Haar-safe mass at least 4/63. All three pairs fail the THM-4233 scalar
  oscillation gate; closure uses exact signed endpoint phases. Removing these
  three edges from the post-THM-4242 proof residual leaves 181,123 edges with
  maximum endpoint 768. This is not arbitrary pair entry, and LRC(14) remains
  OPEN.
source: root/cross-frontier-bridge/2026-08-26
depends_on:
  - THM-4228-common-gcd-two-outsider-periodic-observable-haar-ray
  - THM-4233-pair-specific-primitive-observable-oscillation-haar-charts
  - THM-4242-fixed-fifty-direct-r590-tail-and-twenty-three-label-chart
related:
  - THM-4245-primitive-observable-component-floor-and-cofinal-gate-redundancy
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
  - THM-2162-signed-endpoint-cocycle-and-bv-component-split
  - THM-4207-two-newcomer-sharp-depth-transition-base-surplus-composition-and-variable-pool-chart-number
primary_script: 04-computation/lrc14_exact_signed_endpoint_cocycle_residual_edge_closure_primary_thm4252.cpp
primary_output: 05-knowledge/results/lrc14_exact_signed_endpoint_cocycle_residual_edge_closure_primary_thm4252.out
independent_audit_script: 04-computation/lrc14_exact_signed_endpoint_cocycle_residual_edge_closure_independent_audit_thm4252.cpp
independent_audit_output: 05-knowledge/results/lrc14_exact_signed_endpoint_cocycle_residual_edge_closure_independent_audit_thm4252.out
residual_postprocess_script: 04-computation/lrc14_exact_signed_endpoint_cocycle_residual_postprocess_thm4252.py
residual_postprocess_output: 05-knowledge/results/lrc14_exact_signed_endpoint_cocycle_residual_postprocess_thm4252.out
primary_script_sha256: ead75230a5abdfb177220c377470a0dbf1b1a2bd0832e464fc0da9e770be89c9
primary_output_sha256: 96536059edffcfa58710907aea6d485c130a1aa8feb1cd5f03ccf40d4397fb35
independent_audit_script_sha256: 5bf28a44af861c29e97729dbade7ab7e81e2e242a60980193db62ad4f7d82a77
independent_audit_output_sha256: c351ceb60dfcdb60a262a10bc40db452758377647fa47ece9ff53b5f45d8b72f
residual_postprocess_script_sha256: b174282cff9d3bd0a6e8ce3f54fcbdcacdde308eebae14d9054d46e33a66c74c
residual_postprocess_output_sha256: c881ad7a8e98fcbcc237bb15b0fa8e2bec26847669cc016b94bdcdd87d41250a
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. The primary computes the exact signed endpoint
  cocycle on every fixed-pool cell, applies a labelled superset transform to
  all 5,852,925 repairs, and exhausts all 14,307,150 bodies for each named
  pair. The clean-room audit instead builds fresh full joint walls, directly
  integrates every repair in the emitted prefix, and recursively enumerates
  every body. Both -O0/-O3 transcript pairs are byte-identical. ASan+UBSan
  replays reproduce both transcripts with empty diagnostics; leak detection
  is unavailable on the audited Apple toolchain. No load-bearing assert,
  floating point, or sampling is used.
---

# THM-4252 -- exact signed endpoint-cocycle residual-edge closure

**PROVED RELATIVE TO THM-4228/4233/4242 + FINITE-EXACT + INDEPENDENTLY
AUDITED FIXED-THREE-EDGE CLOSURE; LRC(14) REMAINS OPEN.**

## 1. Statement and inheritance

For a finite positive set `A`, write

```text
G_A={x in R/Z:min_(a in A)||ax||>=1/14},
alpha=4/63.
```

Retain the fixed pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,264,286,290}.  (1)
```

Then, for every

```text
(q,r) in {(466,699),(616,769),(721,769)}
```

and every `B in binom(P,9)`, one has

```text
mu(G_(B union {q,r}))>=4/63.                            (2)
```

The closest proved mechanism is THM-4228's exact component transfer, refined
by THM-4233's pair-specific primitive observable. The canonical hostile is
THM-4207: two lawful newcomer marginals need not compose. The corrected near
miss is MISTAKE-524's prohibition on consuming a reserved theorem; THM-4252
is promoted here only after its complete proof artifacts. MISTAKE-464 fixes
all component counts below as **positive-length circular components**, with
isolated safe points excluded. The least-used sidecar is the labelled signed
endpoint chain of the unprojected repair set.

THM-4245 is now proved and shows that THM-4233's scalar gate is wholly
cofinal-redundant. It is related explanatory context, not a dependency of the
exact endpoint calculation below.

The live concept board is

```text
primitive H | signed endpoint addresses | exact repair hypergraph
transversal number | residual top layer.                (3)
```

## 2. Lossless integer endpoint cocycle

Fix a named pair and write

```text
(q,r)=g(u,v),       gcd(u,v)=1,       N=14uv,
D=18,241,159,416,480.                                  (4)
```

Put `A=G_u intersect G_v`. Sort the unique joint-comb wall ticks

```text
v(14k+-1),       u(14l+-1)       modulo N
```

together with `0,N` as

```text
0=e_0<e_1<...<e_s=N.                                   (5)
```

Let `A_j in {0,1}` be the indicator of `A` on
`(e_j/N,e_(j+1)/N)`, and set

```text
d_j=e_(j+1)-e_j,
S=sum_j A_j d_j,              beta=S/N,
C_0=0,
C_(j+1)=C_j+d_j(N A_j-S).                              (6)
```

For THM-4233's centered primitive

```text
H(t)=integral_0^t(1_A(s)-beta) ds,
```

equation `(6)` gives

```text
H(e_j/N)=C_j/N^2.                                      (7)
```

For a repair `R in binom(P,8)`, put

```text
U_R=G_(P\R).
```

Every `U_R` is empty in a neighborhood of the chosen origin, so its
positive-length circular components have nonwrapping oriented lifts

```text
[a_i/D,b_i/D],       0<a_i<b_i<D,
m_R=sum_i(b_i-a_i).                                    (8)
```

For a fixed-pool endpoint tick `t`, let `p_t` be the representative of
`gt mod D` in `[0,D)`, and choose `j` with

```text
D e_j <= N p_t <= D e_(j+1).
```

Define

```text
Z_g(t)=D C_j+(N p_t-D e_j)(N A_j-S).                   (9)
```

Continuity makes `(9)` independent of the adjacent-cell choice at equality.
Direct affine integration on the primitive cell gives

```text
H(gt/D)=Z_g(t)/(D N^2).                                (10)
```

Therefore each oriented component has the exact centered error

```text
epsilon_i=(Z_g(b_i)-Z_g(a_i))/(g D N^2).               (11)
```

Summing THM-4228/4233's exact component identity yields

```text
K_R=sum_i(Z_g(b_i)-Z_g(a_i)),

mu(U_R intersect G_q intersect G_r)
  =(g N S m_R+K_R)/(g D N^2).                          (12)
```

Thus exact activation is the single non-strict integer predicate

```text
63(g N S m_R+K_R)>=4g D N^2.                           (13)
```

No endpoint ambiguity affects Haar mass. The programs nevertheless retain
the endpoint addresses and signs until after the exact sum.

The older scalar inequality is obtained only after the lossy replacement

```text
K_R>=-c_R D(max_j C_j-min_j C_j),                       (14)
```

which is equivalent to discarding the signed phase data in favor of
`-c_R omega/g`. Equation `(14)` destroys endpoint addresses, signs, phase
correlation, and cancellation; equation `(12)` does not.

## 3. Exact repair hypergraph and the body implication

For a named pair define the eight-uniform labelled hypergraph

```text
E(q,r)={R in binom(P,8):
        mu(U_R intersect G_q intersect G_r)>=4/63}.     (15)
```

The computation ranges over all

```text
binom(30,8)=5,852,925                                  (16)
```

repairs, including repairs rejected by THM-4228's scalar-surplus filter.

Suppose that every `B in binom(P,9)` is disjoint from some `R in E(q,r)`.
Then

```text
B subset P\R,
B union {q,r} subset (P\R) union {q,r},
G_((P\R) union {q,r}) subset G_(B union {q,r}).         (17)
```

Equations `(15)` and `(17)` give `(2)`. Equivalently, the transversal number
of `E(q,r)` is greater than nine. Notice that the inclusion in `(17)` points
from the more constrained active repair set into the target body set; reversing
it would invalidate the proof.

## 4. Exact three-edge certificate

Every named edge satisfies `beta>=66/91` but fails THM-4233's scalar condition
`omega/g<=1650/28710227`.

The primary exact census gives:

| pair | primitive data | active repairs | active-deck FNV |
|---|---|---:|---|
| `(466,699)` | `233(2,3)`, `beta=16/21`, `omega/g=67/205506` | 5,114,702 | `33a4fecc55e46462` |
| `(616,769)` | primitive, `beta=36/49`, `omega=11629/11605748` | 4,632,061 | `e7fb7af64ae2edcc` |
| `(721,769)` | primitive, `beta=36/49`, `omega=17862/27168001` | 4,762,470 | `059c1899134692d7` |

Here and below an FNV ledger means standard FNV-1a-64 over each unsigned mask
serialized as one little-endian eight-byte word, in the displayed order.

Order active repairs by

```text
(SplitMix64(mask xor 0x4245422842334245), mask).        (18)
```

The first `K` repairs in this fixed order already have transversal number
greater than nine:

| pair | `K` | prefix FNV | order-minimality body | final repair |
|---|---:|---|---|---|
| `(466,699)` | 481 | `a5461a33cd4e1e8c` | `{88,132,168,170,176,193,240,252,264}` | `{10,15,42,60,80,126,145,190}` |
| `(616,769)` | 709 | `f07e55a1f5c2b1b6` | `{80,85,88,95,170,176,190,252,264}` | `{8,16,42,63,132,143,145,290}` |
| `(721,769)` | 672 | `88a16a55ce12cb5a` | `{8,85,88,143,145,176,190,240,252}` | `{16,20,60,63,132,168,170,193}` |

All `14,307,150=binom(30,9)` bodies are covered. In each row the displayed
body intersects the first `K-1` repairs and is disjoint from repair `K`.
Thus the prefix is minimal **within the frozen order**. No globally
minimum-cardinality certificate is claimed.

For the final repairs, the exact positive margins are

```text
(466,699): 63mu-4=205871087047/16865833904920,
(616,769): 63mu-4=69858259059/112909716920,
(721,769): 63mu-4=228482347209/3713369504845.            (19)
```

The signed-component controls are:

| pair | components | positive / negative / zero errors | exact total centered error | scalar absolute-value budget |
|---|---:|---:|---|---|
| `(466,699)` | 196 | `94 / 102 / 0` | `-68912209443349/22313498256209160` | `134/2097` |
| `(616,769)` | 212 | `100 / 112 / 0` | `-26703447540083/5925389034244680` | `616337/2901437` |
| `(721,769)` | 186 | `98 / 88 / 0` | `-1504924197793/674304215379795` | `3322332/27168001` |

Their absolute ratios are respectively about `4.8331%`, `2.1215%`, and
`1.8250%`. The scalar failure and exact success are therefore explained by
signed endpoint cancellation, not by an improved oscillation ceiling.

Equations `(12)--(19)` prove `(2)`.

## 5. Primary and clean-room algorithms

The primary implementation:

1. reconstructs the `7,133` fixed-pool open cells on denominator `D`;
2. computes the primitive cell word `(e_j,A_j,C_j)` exactly;
3. evaluates `H(g right)-H(g left)` on every fixed-pool cell and checks it
   against a separately evaluated direct pair-safe prefix integral;
4. groups the exact cell masses by their labelled pool-failure masks;
5. applies an ordinary-colex superset transform to all repairs;
6. thresholds with `(13)`, freezes the full active-deck ledger, and exhausts
   every labelled body; and
7. reconstructs the positive-length components of the final repair and checks
   their signed sum against the transformed mass.

The clean-room implementation reads only the emitted prefix masks. For each
named pair it:

1. constructs a new common-denominator wall arrangement for all thirty pool
   speeds together with `q,r`;
2. classifies every open joint cell directly at its exact midpoint;
3. integrates every prefix repair literally, without `H`, atom grouping,
   colex ranks, or a zeta transform;
4. cross-multiplies the final literal mass against the primary component mass;
5. enumerates bodies recursively rather than by Gosper masks; and
6. verifies universal coverage and the order-minimality body.

The referee does not independently reconstruct the complete multi-million
edge active hypergraph or the positive/negative component census. Those are
primary exact sidecars; the smaller closing prefixes, their literal masses,
and the full theorem quantifier are independently checked.

All C++ gates use explicit throwing `require` functions that remain active
under `NDEBUG`. There are no load-bearing `assert` statements.

## 6. Proof-graph consequence

THM-4242 leaves an exact residual `E_4242` with

```text
|E_4242|=181,126,
FNV(E_4242)=bdf59726990a6c92,
max endpoint=769,
top layer={(616,769),(721,769)}.                        (20)
```

All three edges proved here belong to `E_4242`. Put

```text
E_new=E_4242\{(466,699),(616,769),(721,769)}.           (21)
```

Exact postprocessing gives

```text
|E_new|=181,123,
FNV(E_new)=6ec03ed4c4dc841b,
SHA256(E_new)=9a9b6fbe14db00e9d7f8f08ecddaa1e3d263fd063c6b3c003e18c210b3334ef8,
max endpoint=768,                                      (22)

endpoint-768 layer={(616,768),(721,768),(744,768),
                    (750,768),(765,768),(766,768)}.     (23)
```

Consequently the combined proof-graph cutoff becomes `769`: every outsider
pair with maximum endpoint at least `769` is now proved. This cutoff is exact
for the combined certificate graph, not claimed literal-minimal. No surviving
edge is asserted unsafe.

## 7. Connection contract and related controls

```text
source:       primitive joint-comb word, gcd g, and labelled signed boundary
              chain of every repair U_R
target:       exact active-repair hypergraph E(q,r), then the proof residual
map:          pair H_(u,v) o m_g with the signed boundary, sum exact mass,
              threshold, then test nine-body transversals
preserved:    primitive ratio, scale, endpoint phase/sign, cancellation,
              repair/body labels, exact mass, and body consequence
lost by sum:  null endpoint membership; phase grouping loses component pairing
              but preserves exact mass
lost by gate: surplus magnitude after activation
lost scalar:  all endpoint addresses, signs, correlations, and cancellation
sidecar:      signed phase multiplicity nu_R(p)=#right(p)-#left(p)
hostile:      (466,699)=233(2,3), the closest scalar hostile from THM-4245
decisive:     exact tau(E(q,r))>9 for each named pair.   (24)
```

The same discovery engine and a fresh-wall replay also closed the residual
controls `(319,638)=319(1,2)`, `(200,600)=200(1,3)`,
`(540,720)=180(3,4)`, and `(41,533)=41(1,13)`. These are **AUDITED RELATED
EXAMPLES**, not statements of this theorem, not dependencies of `(2)`, and
not removed in `(21)`. A non-load-bearing audit summary is retained in the
promotion packet rather than expanding the theorem universe.

## 8. Reproduction and raw hashes

From the repository root:

```bash
clang++ -std=c++20 -O0 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_exact_signed_endpoint_cocycle_residual_edge_closure_primary_thm4252.cpp \
  -o /tmp/thm4252-primary-O0
/tmp/thm4252-primary-O0 > /tmp/thm4252-primary-O0.out

clang++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_exact_signed_endpoint_cocycle_residual_edge_closure_primary_thm4252.cpp \
  -o /tmp/thm4252-primary-O3
/tmp/thm4252-primary-O3 > /tmp/thm4252-primary-O3.out
cmp /tmp/thm4252-primary-O0.out /tmp/thm4252-primary-O3.out

clang++ -std=c++20 -O0 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_exact_signed_endpoint_cocycle_residual_edge_closure_independent_audit_thm4252.cpp \
  -o /tmp/thm4252-referee-O0
/tmp/thm4252-referee-O0 /tmp/thm4252-primary-O0.out \
  > /tmp/thm4252-referee-O0.out

clang++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_exact_signed_endpoint_cocycle_residual_edge_closure_independent_audit_thm4252.cpp \
  -o /tmp/thm4252-referee-O3
/tmp/thm4252-referee-O3 /tmp/thm4252-primary-O3.out \
  > /tmp/thm4252-referee-O3.out
cmp /tmp/thm4252-referee-O0.out /tmp/thm4252-referee-O3.out

python3 -B \
  04-computation/lrc14_exact_signed_endpoint_cocycle_residual_postprocess_thm4252.py \
  > /tmp/thm4252-residual.out
```

The raw-LF hashes are:

| artifact | SHA-256 |
|---|---|
| primary source | `ead75230a5abdfb177220c377470a0dbf1b1a2bd0832e464fc0da9e770be89c9` |
| primary output | `96536059edffcfa58710907aea6d485c130a1aa8feb1cd5f03ccf40d4397fb35` |
| independent source | `5bf28a44af861c29e97729dbade7ab7e81e2e242a60980193db62ad4f7d82a77` |
| independent output | `c351ceb60dfcdb60a262a10bc40db452758377647fa47ece9ff53b5f45d8b72f` |
| residual postprocessor | `b174282cff9d3bd0a6e8ce3f54fcbdcacdde308eebae14d9054d46e33a66c74c` |
| residual output | `c881ad7a8e98fcbcc237bb15b0fa8e2bec26847669cc016b94bdcdd87d41250a` |

The stored `-O0` and `-O3` copies of each C++ transcript are byte-identical.
An Apple Clang ASan+UBSan build of each source also reproduces its frozen
transcript with empty diagnostic stderr. Apple AddressSanitizer does not
support leak detection, so leak checking is not claimed.

## 9. Scope firewall

This theorem does **not**:

1. close any pair other than the three in `(2)`;
2. promote the four related controls in Section 7;
3. prove a dilation ray or monotonicity in `g`;
4. claim globally minimum certificate cardinality;
5. assert that any surviving residual edge or failed scalar gate is unsafe;
6. replace the fixed pool, prove arbitrary-pair entry, or prove physical entry;
7. prove that the new cutoff `769` is literally minimal; or
8. prove LRC(14). **QED.**
