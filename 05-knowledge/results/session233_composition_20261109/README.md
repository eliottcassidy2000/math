# Session233: constructive composition and endpoint universality

**Checkpoint2: three independently audited results. Research is continuing;
Chan–Pak completeness remains OPEN.**

- [Coprime interfaces](graft/PROOF.md) combine the h37 and h45 moving-core
  constructions to interlace their dilated progressions into consecutive odd
  windows. The explicit asymptotic width/vertex ratio is37 times the inherited
  single twelve-bit core's ratio, with a large finite overhead.
- [Complete endpoint spectrum](endpoint/PROOF.md): every positive integer N is
  an endpoint count and endpoint-stabilizer quotient of a strong tournament
  on at most max(5,ceil(log2(N+1))+3) vertices. Binary terminal weights provide
  a real selector; the exact response sidecar explains why this does not
  automatically select full Hamilton counts.
- [Three-state bridge transfer](transducer/RESULT.md) gives the full packing
  polynomial for complementary transitive blocks. Its complete bounded-value
  census has1374 values through10000. Both canonical run-window endpoints,
  even closed under ordinal sums, miss11.

Canon: THM-5000, THM-5001 and THM-5002, with full filenames linked in the
respective proofs. Independent audits: [composition](graft_audit/AUDIT.md),
[endpoint](endpoint_audit/REPORT.md), [bridge transfer](transducer_audit/AUDIT.md).
The next experiment tests whether shorter nonunit moving-core loops can recover
unit spacing after composition. No result from that experiment is claimed here.
