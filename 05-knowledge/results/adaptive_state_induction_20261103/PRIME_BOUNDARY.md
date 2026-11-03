# The five prime groups and the boundary of future named contacts

**FINITE-EXACT for the25-prime ledger; PROVED elementary boundary and
square-scaling transport. No spectrum or minimum-pathwidth theorem.**

The [recovered five groups](../prime_lane_ports_five_groups_20261102/prime_groups.md)
partition all25 primes up to100, with sizes6,5,4,5,5. Their square-sum graph
has24 edges. The earlier packet proved that no group is an equitable class:
the group name does not determine the individual prime's neighbors. Here we
connect that loss to the new live-port representation through a precise
sequential boundary, not an imposed orientation of square-sum edges.

For any ordered finite graph, the future boundary after a prefix consists
of old vertices having a neighbor in the remaining suffix. In a generator
where every later contact must explicitly name its old endpoint, and which
permanently forgets a vertex name when it discards that port, these
are exactly the names that cannot yet be discarded. This is the same set
construction as L_t in the adaptive root word, applied to the undirected
parent-incidence tree instead of the full tournament. It forgets tournament
arcs and Hamilton counts, so neither realization theorem transfers merely
from a shared boundary size. It is also not a lower bound for algorithms
allowed to reconstruct discarded names later.

For primes in ascending numerical order through100, the maximum future
boundary is **8**, attained after47 and after53. At either cut the live
primes are

    2,3,5,11,17,29,41,47.

The receipt supplies a future square-sum partner for each at each cut,
and checks every cut. Removing the even prime2 gives an ascending boundary
maximum **7** on the remaining24 vertices. These are prescribed-order
widths, not minima over all vertex orders. Five colors do not identify the
eight physical contacts at the full-prime boundary.

There is also a literal realization bridge. The earlier square-carry
construction supplies distinct tournament classes labelled H=49p for every
odd prime p in this list, using the H343 carry atom when p=7. The relation
between these classes is the declared secondary relation that their H
values sum to a square. It is not an arc inside a tournament. Since

    p+q=s^2 iff 49p+49q=(7s)^2,

the increasing map p->49p preserves every edge, every nonedge, and every
prefix boundary. Thus the width7 ledger transports exactly to these
realized odd Hamilton-count labels. For p!=7 the construction can use an
ordered join with H49; the p=7 vertex uses its separate actual H343 witness,
since H7 is forbidden. The prior23–2–47 obstruction still prevents simply
adjoining an odd H label for2 under the same odd-square dilation.

The connection contract is explicit: ordered prime vertices map to actual
Hamilton-count-labelled classes; square adjacency, numerical order and
live boundaries survive; internal tournament structure and construction
port slots do not. Named primes, square roots and chosen actual witnesses
remain necessary sidecars. No group-theoretic or prime-infinitude conclusion
is extracted from the numbers7 or8.

`prime_boundary.py` independently enumerates edges by pairs and by square
roots, reconstructs all prefix boundaries, and checks every scaled odd-prime
boundary. Normal and optimized runs must match prime_boundary.json.
