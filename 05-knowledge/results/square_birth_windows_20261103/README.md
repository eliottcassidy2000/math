# Square-sum birth windows: incidence, unused edges and future memory

Status: **PROVED scoped / FINITE-EXACT; same-agent audit with independent
partition refinement.** This packet extends the complete 25-to-26 interface,
not the unrestricted tournament carrier or A000568 enumeration.

The new exact minima are **761 states for 26 -> 27** and **903 states at 25
for 25 -> 26 -> 27**. The operation language retains current named square-sum
Hamilton paths, every legal actual path complement, and the specified
orientation of each new square-sum edge at each birth. Global converse also
descends. Arbitrary re-marking outside this path bank is not included.

| Source | Last permitted birth | Smallest full-source carrier | Status here |
|---:|---:|---:|---|
| 25 | none | 33 | inherited and replayed |
| 25 | 26 | 97 | inherited and replayed |
| 25 | 27 | **903** | new finite-exact |
| 26 | none | 71 | inherited and replayed |
| 26 | 27 | **761** | new finite-exact |
| 27 | none | 637 | inherited and replayed |

These include the entire orientation cube, including sources with no current
path. At 27, 56,808 of the 62,952 active orientations restrict to sources with
no Hamilton path at 26. Recursion through only old Hamilton paths misses them.
One old guard bit, on edge 2--7, separates two sources identical under the
entire fixed-26 observer once the same ear 22 -> 27 -> 9 is added.

The user's edge-arrival idea recovers an already proved arithmetic clock:
`b(n)=floor(sqrt(2n-1))-floor(sqrt(n))`. Its count can stay constant while its
actual ports change. At 9 and 25 a square exits exactly when another enters.
Thus 8/10 and 24/28 should be compared using root names and parents, not only
equal edge counts. [GROWTH.md](GROWTH.md) supplies the exact comparison.

Two further elementary mechanisms clarify the orientation proposal:

- Every acyclic orientation admitting a Hamilton path is the unique forward
  orientation of its path order. If the graph has any extra edge, its actual
  path complement creates a directed cycle. This small model is not closed.
- At 25 a path leaves 8 of 32 edges unused; only edge 1--15 is unused by
  every path. **56** is the first order with more edges unused than used:
  111 edges versus 55 path edges. The inequality persists at all larger
  orders; a literal Hamilton path at 56 certifies nonvacuity.

The [first-25 ledger](LEDGER.md) distinguishes integer Q_k, integer Q_(p_k),
and the square-sum graph on actual prime values. The 25 primes below 100
have 24 square-sum edges, but three isolated vertices. The 15 Ogg primes
have residue-bipartition sizes 6 and 9, so cannot have a square-sum Hamilton
path. The user's three sets are instead the three components of Q18 with
13,14,15 removed. Equal cardinality does not identify these objects.

Read [BIRTH_WINDOW.md](BIRTH_WINDOW.md) for the general finite-horizon reduction
and exact certificates, [RECOVERY.md](RECOVERY.md) for inheritance and typed
connections, and [AUDIT.md](AUDIT.md) for reproduction and boundaries.
The all-size smallest complement/re-marking carrier and the tournament
Hamilton-count gap conjecture remain **OPEN**.
