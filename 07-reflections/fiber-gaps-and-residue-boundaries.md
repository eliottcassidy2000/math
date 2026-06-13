# Fiber Gaps and Residue Boundaries

Session: codex-2026-05-30

This session wandered through four threads that initially looked unrelated:

- good-cut height and the missing bucket `1`;
- quotient transport checksums;
- endpoint-transfer parity boundaries;
- projection defects and old-coordinate residue.

The common object is a quotient map.

```text
q : total space -> bucket space
```

Everything else is what the fibers do.

## Four Faces

**Good-cut height.** The quotient is `goodCutBucket`.  The fiber over `1` is
empty, and after THM-355 the finite image is exactly `{0} union {2,...,n-1}`.
The missing bucket is not a poetic absence; it is a zero row and zero column in
every transport matrix that uses this quotient.

**Transport.** The checksum says every nonempty row has a conserved half-line
budget.  THM-355 adds the complementary fact: empty fibers are silent.  They
cannot emit, and nothing can land in them.

**Endpoint transfer.** THM-347 says endpoint insertion sees child fibers by
column sums, and modulo 2 it sees the odd-fiber boundary.  This is not the same
as good-cut support, but the representation is parallel: first determine
fiber support/parity, then ask whether the transfer rows remember enough
information.

**Projection defect.** Deletion, support shadow, even-graph projection, and
path-homology old projection all ask what survives after a quotient forgets
coordinates.  A projection kill is a gap after deletion; a near-kill is a tiny
residue fiber.

## Hypotheses Generated

1. **Gap columns explain many spectrum absences.** Fixed-n H-value gaps should
   behave like good-cut bucket gaps when viewed as fibers of `T -> H(T)`.
   Their transport columns should be exactly zero; neighboring columns should
   reveal the obstruction geometry.

2. **Parity boundaries are support derivatives.** Endpoint-transfer mod-2
   boundaries, self-complementary merged nodes, and odd automorphism 2-adic
   graph orbits may be different parity shadows of the same fiber-support
   calculus.

3. **Projection defects are residue masses.** Instead of treating
   `delta_proj` as a standalone feature, treat it as the mass left in fibers
   after a projection deletes or identifies coordinates.

4. **TDA features should be quotient-native.** A robust feature block should
   include `fiber_support`, `gap_count`, `boundary_parity`, `transport_escape`,
   and `residue_mass`, rather than only scalar invariants such as `H`, score,
   or `g`.

## Feedback Loop

The loop that feels productive is:

1. Choose a quotient because it has semantic content.
2. Prove or compute its fiber support.
3. Use transport to see which fibers communicate.
4. Use projection defects to ask what residue survives coarser views.
5. Feed the residue back as a new quotient or feature.

This is not a replacement for the older stories about parity, staircase
geometry, or the two-sheeted cover.  It is a lower-level grammar in which those
stories can be compared.  The two-sheeted cover says score fibers have
monodromy.  Endpoint transfer says quotient towers expose fiber boundaries.
Good-cut height says interval coverage has missing states.  Projection defect
says forgetting coordinates leaves residue.

All four are about fibers refusing to be uniform.

References: THM-347, THM-349, THM-352, THM-353, THM-355, HYP-1783,
`04-computation/lean/TournamentH7/TournamentH7/BucketBalance.lean`.
