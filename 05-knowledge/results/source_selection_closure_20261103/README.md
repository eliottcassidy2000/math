# Constructive source selection: reusable resets and the remaining branching core

**PROVED scoped construction and obstruction theorems / FINITE-EXACT small
certificates. Full prime and 7p realization, and full native metagraph
classification, remain OPEN.** This continues the
[two-port reset session](../two_port_return_induction_20261103/README.md).
The target remains odd primes p!=7 and 7p for odd primes p!=3, under
[THM-4094, two-prime-lane completeness](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md).
The exclusion of 7 and 21 is inherited; incoming independent recertifications
are reconciled in [INCOMING](INCOMING.md).

The main new positive result is a **closed branching family of actual
primitive sources with repeatable native resets**. The main limitation is
equally precise: even all block lengths in that balanced architecture are
too sparse for the complete target set. The other results identify where
an induction must retain source information and make genuinely new choices.

## 1. An actual reset can now be repeated

For an actual strong balanced rooted profile (h,s,q,F)=(b,t,t,b), a cycle
of adaptive births, selection of the oldest physical leaf, and two further
births returns to another actual balanced profile:

    (b,t) -> ([3*2^k+4*3^(k-1)+1]b-2t, (2^k+1)b), k>=1.

For k>=2 the provisional newest-leaf pair has continuation divisor b;
moving that same physical leaf to the oldest position changes it to
gcd(b,t). This erases arbitrarily large introduced factors. An adjacent
swap compiler reduces the native cost to 2(k-1) moves, with k-1 temporary
births/deletions and complete named reentry. The endpoint pair is ready for
the next cycle. Temporary intermediate graphs have their own counts and
are not silently treated as neutral.

For k=3 and k=5 the balanced maps are

    C3(b,t)=(61b-2t,9b),     C5(b,t)=(421b-2t,33b).

After the first step from rooted C3, keep gcd(b,t)=1, 3|t and 3∤b.
C3 always preserves this guard. C5 does so exactly when t/b!=7 in
P1(F11), and sends that ratio to zero. Eight projective control states
suffice; the unbounded integers, actual source, and physical pair are also
retained. This is a proved branching source family, not a finite bank of
fixed integer seed states.

The short branch C3 then C5 reaches (76147,5973). Its next birth has

    H=146321=7*20903.

The 18-vertex graph is independently counted. A subsequent C5 oldest exit
gives 17203249=7*2457607 by the proved source transfer; its 24-vertex graph
and programs are recorded without claiming a 24-vertex DP recount.
Both quotients are proved prime by trial division.

[NATIVE_SELECTION](NATIVE_SELECTION.md) proves these results and the general
position law: for k children at a rooted source (h,s,q,F), the pair formed
by the old root and a chosen leaf has content gcd(h,s,q,F) at the
oldest position, and gcd(q,F,h-s) at every other position. The two values
can coincide; no inherited common factor is erased by this selector.

## 2. Source changes repair the H31 branch obstruction

The old 37 H31 marked states come from exactly three adaptive ordinary
classes. Even allowing every future old vertex as a root cannot give any
one of these three sources both positive branches 95 and 97. Merely keeping
more names therefore does not fix this particular source obstruction.

Two concrete operations do improve it:

* A newborn losing to both retained ports, already supported by the same
  universal seventeen-coordinate observer, gives one named H31 source all
  three successors 93,95,97. Each child has a valid marked C3 reentry when
  the discard choice is correct. Complete continuation through H295 shows
  that none of its available children repeats this branching property in
  the nine proper-cut command language.
* One H-neutral native move carries an actual H31 source supporting 95 to
  another supporting 97, with the same two physical port names. The second
  source has adaptive core order five, so its ordinary source changes even
  though the universal observer still applies. The exact fixed-name neutral
  component is a six-cycle over two rigid classes. Its class-return action
  is C3, although the class-staying native groups are trivial.

[BRANCH](BRANCH.md) gives literal rows, paths, terminal programs, complete
finite universes, and the distinction between a one-level repair and a
closed recursive family.

## 3. The four-coordinate grammar has a canonical inverse boundary

Every strong rooted source satisfies h-s>=q. In the stronger physical
coordinates

    x=(h-s-q,q-s,s,F-h-q+s)>=0, with s>0,

all 28 inherited maps remain nonnegative. Fifty-four exact dual certificates
show that each full-rank map, 041 and 403, has an image disjoint from every
other last-map image on this cone. A profile therefore has a unique maximal
full-rank suffix. Those maps preserve primitivity; their suffixes can be
peeled without making a source choice.

At the remaining singular boundary, each predecessor fibre has at most two
free integer coordinates and a necessary finite first-birth bound. This
gives a terminating membership decision procedure for the full 21-entry/
28-map profile grammar, returning literal programs. It also gives a finite
target-directed scalar decision procedure. No larger coverage scan was run.

A fixed budget K of singular steps yields only
O_K(X^(log_6 2)(log X)^K) states below X. Hence a successful prime/7p
construction **within this grammar must use unboundedly many singular
choices**. This does not exclude a finite-state controller with unbounded
execution. [GRAMMAR](GRAMMAR.md) contains the decoder, proof and controls.

## 4. Closure and sufficient target supply remain separate requirements

Even allowing every block length k>=1 in the balanced reset architecture,
its declared adaptive and return outputs satisfy

    # {outputs <=X} < 30 (X/b0)^(2/3)

from an initial count b0. The proof sums reciprocal growth weights over all
words; it is not a finite-value test. Thus this architecture misses infinitely
many primes and 7p values. The construction remains a reusable reset
mechanism, but the target selector must leave this architecture or enlarge
its declared output language. Temporary native intermediate counts are
explicitly outside this bound. [BALANCED_BOUNDARY](BALANCED_BOUNDARY.md).

The wildcard test explains the cost of one possible enlargement. Admitting
both polarities of births cannot have a universal finite-dimensional linear
observer, even on strong one-port sources. An exact constructive replacement
uses eighteen endpoint-cover polynomials and first derivatives. The degree
records the forgotten path components that the old observer could discard.
This carrier preserves all 24 two-polarity responses, while actual sources
and names remain sidecars. [GRADED_OBSERVER](GRADED_OBSERVER.md).

## 5. Procedural board and next proof obligations

Anchor = productive actual-source selection. Niche = native transport of a
physical port through a homogeneous block. Wildcard = graded attachments
when the operation alphabet changes. The live board was source ancestry,
terminal names and port age, continuation content, inverse-image cones,
return clocks, and output supply. The closest proved mechanisms were the
common seventeen-query observer, four-coordinate sealing, and traced arc
compiler. Hostiles were H31 branch incompatibility, the weak H5 predecessor,
and the failed latest-two-leaf reset. The previous false extrapolation from
three neutral examples was retained as a warning to derive an all-order law.

Each round changed a specific obligation: extra names failed at H31; mask3
repaired one level; the native bridge changed the relevant source core;
port position supplied an exact gcd selector; balanced renewal made it
repeatable; growth exposed its insufficient target supply; and cone
separation localized the remaining inverse choices. The only broader
computations were the stated small complete closures, not a new scalar record.

The next productive construction must combine actual source-changing moves
with unbounded singular choices or a richer output language. It needs to
preserve a provable entry program and ordered pair, choose a primitive next
interface, and supply the requested target with a legal reusable child.
Neither a primitive vector, a complete modular orbit, nor a source family
closed under its own operations proves that last condition.

[AUDIT](AUDIT.md) states reproduction and proof-review scopes. Five frozen
jobs replay normally and under -O; source pins and a packet manifest retain
provenance. New results are relative to the recovered workspace lineage,
without an external literature-priority claim.
