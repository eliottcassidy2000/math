# Independent audit of the sharp full-profile sheet-edge bound

**ACCEPTED: PROVED actual-entry reduction and sharpness; FINITE-EXACT
exhaustive quotient certificate.** This independently audits
[the selected-edge note](third_20260906_decoder_profile_graph.md).
The audit imports no repository producer. It reconstructs a larger
scalar-only quotient universe, then consults the complete words only at
seven states. No mathematical correction to the theorem or actual control
was needed. LRC(14) remains OPEN.

## 1. Scope and faithful quotient reduction

Keep the producer's actual balanced hypotheses and its inherited-profile
supplier. The physical row is `tV union gU`, with six and seven labels,
`gcd(V)=gcd(U)=gcd(t,g)=1`. Its two actual strict-atlas components are
connected. The hypothesis that every physical subset passes its full
inherited gcd/complement profile is retained; this is a necessary
condition for strict failure, not a characterization of strict failure.

Put `d_i=gcd(t,u_i)`. For a nonempty index set I of size k<=6,

```text
gcd(tV union {g*u_i:i in I}) = gcd(d_i:i in I) = c_I.
```

Since `c_I` divides t and is coprime to g, each remaining complement entry
is `gcd(c_I,g*u_j)=gcd(c_I,d_j)`. Thus the entire required profile is
exactly the stated pair `(c_I,sort(gcd(c_I,d_j):j outside I))` at level
`7-k`. Also `gcd(d_1,...,d_7)=1`. These identities hold for every subset
I, rather than for a selected spanning set or only the largest gcd.

For an actual atlas edge ij,

```text
gcd(t,gcd(u_i,u_j)) = gcd(d_i,d_j).
```

If all actual edges had sheet gcd at least seven, the state graph
containing every pair of distinct positions with gcd at least seven
would contain that connected actual graph. It would therefore be
connected. This is the only direction of connectivity used. Primitive
ratios and additional state edges are lost by the quotient, and no
converse realization assertion is needed. Zero quotient survivors
therefore imply an actual atlas edge with sheet gcd at most six.

Repeated state values are allowed throughout. Sorting forgets only
label order; subset enumeration still counts distinct positions and
retains repeated entries in their words.

## 2. Independent exhaustive reconstruction

The sole mathematical data dependency is
[the inherited profile JSON](../../04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json),
pinned by raw SHA256

```text
935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f.
```

The script derives the exact allowed scalar gcd sets from the profile
pairs and verifies they equal the stored scalar projections. It does
not use the producer's partial-word tables or partial full-word filter.

For threshold s, begin with every allowed level-six state at least s.
Build all connected sorted multisets by adjoining every state joined to
an existing state at gcd threshold s. During growth, impose only that
every k-subset gcd lies in the exact scalar projection of level `7-k`,
and that the full seven-state gcd is one. The full complement words
are not consulted in this enumeration.

This is exhaustive. Every connected graph has a vertex ordering with
connected prefixes, obtained from a spanning tree. Every subset of a
valid final row satisfies the inherited scalar restrictions. Therefore
some growth ordering of every possible final row survives all scalar
filters. All partial multisets and all extensions are explicitly
included, modulo permutations only.

The independent scalar-only layer counts are:

| Threshold | size 1 | size 2 | size 3 | size 4 | size 5 | size 6 | size 7 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 7 | 36 | 104 | 377 | 1,031 | 1,923 | 1,700 | **104** |
| 6 | 37 | 144 | 788 | 2,929 | 6,407 | 5,985 | **493** |

For each complete state row, the audit independently verifies
connectivity, then applies literal full profile membership for every
subset size 1 through 6. Threshold seven has no survivor. Threshold six
has 26, and the complete independently reconstructed set agrees exactly
with the producer's frozen list.

Every one of the 104 threshold-seven rejections already occurs at a
singleton state's full six-entry complement word. The same is true of
the 467 rejected threshold-six rows. Thus, after the exact scalar
projections have been retained, these finite certificates require only
the singleton full words; the other full words are still checked on
each survivor. This is a finite compression observation, not an asserted
classification at arbitrary graph size.

A stronger hostile than the producer's coarse-cap example is

```text
D=(8,8,9,9,10,60,72).
```

It passes every exact scalar allowed-gcd set at all subset sizes and
has a connected threshold-seven graph. At its state 60 the full
complement word is `(3,3,4,4,10,12)`, which is excluded. This shows
that complete words contribute information beyond even the exact
scalar projections, not merely beyond their numerical maximum caps.
The hostile is a quotient row, with no actual realization claimed.

## 3. Independent actual sharpness verification

The audited control is

```text
V=(1,2,3,4,5,6),
U=(12584,14872,117,9999,98890,132990,10296),
g=1, Q=91^6,
t=360*(1000Q+1)=204432930734760360.
```

The physical row has 13 positive distinct labels, gcd one, and sum
`4293091545430247308<Q^2`. The states are exactly
`(8,8,9,9,10,30,72)`, one of the 26 threshold-six survivors.

The atlas is rebuilt by a prime sieve and multiplicative generation:
each allowed primitive sum is a product of primes congruent to two
modulo three, each with exponent at most two, and is at most 356.
All coprime ordered sizes `p<q` at those sums give exactly 5,855 atlas
ratios. This differs from the producer's per-sum factorization test.
It accepts 25 and rejects the exponent-three hostile 125, in agreement
with **THM-3818**,
[scaled inert cubeclass pair packet](../../01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md),
Section 1.

All 78 physical pairs are inspected. The exact graph has precisely the
six- and seven-label components claimed. Its weighted edge relation
matrix has rational rank eleven under a separate exact matrix engine.
The seven-component has exactly these six edges, in the displayed U
index order:

```text
(0,6): ratio11:9,   sheet gcd8
(1,6): ratio13:9,   sheet gcd8
(2,6): ratio1:88,   sheet gcd9
(3,6): ratio101:104,sheet gcd9
(4,5): ratio29:39,  sheet gcd10
(5,6): ratio155:12, sheet gcd6.
```

Consequently the selected-edge ceiling six is attained by the minimum
over all actual seven-component edges.

The audit establishes absence of every mixed bounded support relation
without importing the producer's signed-box algorithm. For two tV
labels and one outside label u, a relation with outside coefficient c
would imply `t | c*u`. Therefore `t/gcd(t,u)` divides c. Each of these
seven divisors exceeds Q, ruling out nonzero `|c|<=Q`, independently
of the gcd of the selected V pair. This modulo-t argument retains the
required divisibility even when that pair gcd has additional primes.

In the reverse orientation, any nonzero tV term has magnitude at least
t times its primitive label, whereas the sum of the two bounded U terms is at most
`Q*(u_i+u_j)`. The literal inequalities
`t*v>Q*(u_i+u_j)` hold for every pair and every v. Thus all
`binom(6,2)*7+binom(7,2)*6=231` mixed supports are excluded, including
relations whose other coefficient is zero.

Every internal relation lies in the weighted kernel of its component.
The connected decoder edges span those two kernels and have coefficient
height at most 355, below Q. Therefore the mixed-relation exclusions
and rank calculation prove both inclusions in
`W_(Q,3)=V_dec`, of rank eleven. Internal primitive pair heights are also
checked directly and are at most Q.

Finally, all 4,095 physical subset profiles of sizes seven through
twelve are recomputed, including those discarded by the seven-state
quotient. Every complete word passes. Their maximum gcds are
`(72,10,9,3,2,1)`. The exact phase `1/7` has clearance `1/7`.
This verifies sharpness on the stated profile-surviving actual-entry
domain. The example is safe; it does not establish sharpness after
adding the distinct hypothesis of strict failure.

## 4. Reproduction and limits

The [independent source](../../04-computation/third_20260906_profile_edge_audit.py)
and [frozen output](third_20260906_profile_edge_audit.out) reconstruct the
quotient, compare all 26 survivors, and verify the complete actual
control. The computation passes 5,599 always-active exact gates.

```bash
python3 -B 04-computation/third_20260906_profile_edge_audit.py
python3 -B -O 04-computation/third_20260906_profile_edge_audit.py
```

Normal and optimized outputs match byte for byte. Raw LF hashes are:

```text
source 010a08b6396db09eca2350accf7e6480c5a59106d6e23b5c876c5f8b7e98ab00
output 6ed810cd3aaea51cdcfab2d754e4b4e7bdd905c8870decc6d8e23fb0f8728ba9
semantic 64779d70a6b19516a1b3317c034e56a6e217bfe8bb6598d27e1a896e57d65986
```

The first round's universal parameter proofs and the translated-grid
consumer are separate objects. This audit accepts the edge supplier
and its actual sharpness; it does not infer a pair ratio, a universal
phase obstruction, or closure of LRC(14).
