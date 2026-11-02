# Clocked attachment carriers, Fermat neighbours and modular-genus certificates

**Status:** PROVED scoped results and CITED inherited inputs, with same-agent
proof review and exact normal/optimized checks. No independent-agent or
formal audit, external novelty, unrestricted enumeration or universal
Collatz result is claimed.

- [Corner carrier](CORNER_CARRIER.md): with h endpoint edits remaining,
  endpoint-chord queries need exactly (h+1)(h+2)/2 orientation bits,
  conditional on the shared ordered path and clock. Only h bits of a new
  attachment row matter. Arbitrary re-marking destroys this compression.
- [Fermat neighbours](FERMAT_NEIGHBOURS.md): orders 15 and 17 both have 16
  multiplier orbits of circulant presentations. More generally the pairs
  (3,5),(15,17),(255,257),(65535,65537) have matching orbit counts. This
  follows from Akin's free-action theorem and the classical Fermat product
  identity. A binary-carry proof classifies these as the ONLY collisions
  within the full free-action family, even allowing unknown Fermat primes.
  Different unit groups and composite block symmetries obstruct
  the naive symmetry identification.
- [Genus bridge](GENUS_BRIDGE.md): for a Fermat prime p=2^(2^m)+1, m>=2,
  g(X0(p))=(p-5)/12 and 3g+1=2^(2^m-2). Thus 1,21,5461 have exact odd-step
  Collatz certificates. At primes beyond 17 these genera are nonconstructible
  polygon orders. Four residue channels give the larger exact family.
  The pasted Ogg/Fermat bijection is repaired to correctly typed intersections.
- [Collatz geometry](COLLATZ_GEOMETRY.md): an explicit binary inverse grammar
  already inside the basin of 1 has 2^r distinct labels within 7r native steps.
  This rules out bounded-dilation injective placement in a fixed honeycomb.
  A dyadic upper-half-plane layout keeps a constant hyperbolic macro length
  and retains each macro edge's native time.
- [Trees and ports](TREE_PORTS.md): the already proved equivariant tree
  trijection connects marked objects, not A000055 to A000568 class counts.
  Exact bridge responses use one endpoint count; two ports require joint
  endpoint data, with a four-vertex cycle hostile to marginal compression.

[Concept board](BOARD.md) · [Audit and reproduction](AUDIT.md) ·
[Next tests](HANDOFF.md) · [Source pins](SOURCES.json).

Incoming concurrent work now proves sharp weighted stability for n=3..9
and n>=42; orders 10..41 and full raw-H spectrum completeness remain OPEN.
The new cutoff belongs to the [capacity-potential packet](../tensor_primorial_ports_20261102/STABILITY_42.md),
not this session. Our contribution is to the representation and composition tools.
