# Independent audit of the complete tournament response

2026-10-26. **PASS**, after one prepublication wording clarification. This
reviews [information/RESULT.md](../information/RESULT.md), its
[checker](../information/check.py), and the parent's
[independent audit](../audits/INFORMATION_AUDIT.md). It does not promote
unrelated historical statements or turn finite controls into all-order proofs.

## All-order proof review

The coefficient calculation is correct with the displayed orientation
convention. For two incident edges, a directed vertex word using both must
contain their three vertices as a consecutive block. Its two traversals
give the same product of the two reference signs, and there are
`2(n-2)!` words. Dividing the path-indicator expansion by `2^(n-1)` gives
exactly `epsilon*2(n-2)!/2^(n-1)`. For disjoint edges, reversing one
two-vertex block is a fixed-point-free involution of the contributing
words and reverses exactly one selected sign. Their coefficient is zero.
There is no odd-order hypothesis in either argument.

The shift identity for `R_T(F)=H(T xor F)` follows by changing variables
in the Boolean character sum. Its sign is `x_e*x_f`, not its negative.
The graph on arc coordinates with adjacency for a shared endpoint is
connected for `n>=3`. Consequently equality of all those products forces
all coordinate ratios between two tournaments to be equal. The converse
direction follows by reversing the entire Hamiltonian vertex word.
This proves the full-response fibre is precisely a global converse pair.

The `m-1` code is a sharp **fixed-length binary encoding**: products on
a spanning tree have all `2^(m-1)` possible sign patterns and each has two
oriented lifts. This does not prove an `m-1` bound on the number of H
queries, nor an optimal adaptive query bound. Computing the displayed
correlations directly uses the full response table. The packet expressly
retains this distinction.

Full signed curl separates the two converse lifts because each triple curl
is a nonzero odd integer. The aligned pair therefore determines the
tournament for `n>=3`; the order-two boundary correctly fails. Passing to
simultaneous vertex relabeling gives ordinary tournament isomorphism,
whereas taking separate quotients discards the relative alignment. The
`C3 -> v` versus `v -> C3` example has the stated equal curls, converse
isomorphism, and distinct score multisets, and order four is minimal.

The path-start selector is also sound given the recovered curl-fibre
classification: any Hamiltonian path visits the ordered strongly connected
components consecutively and in that order. Its starting component fixes
the initial cut of their circular order. The promise that the marked vertex
starts a valid Hamiltonian path is essential; a free vertex mark has no
such information. Retaining the full valid word gives the genuine tiling
operation, with only non-path arcs complemented.

## Clarification made during this audit

An earlier sentence suggested the repair for the curl intervention problem
was "an I2 orientation bit." Read without its full I2 context, this could
mean that one oriented arc added to curl suffices. The producer replaced
it with the complete response (or its `m-1` correlation code) **together
with** one oriented arc, or the aligned curl/response pair.

The distinction has a minimal exact hostile. With edge order
`01,02,12` and mask bits denoting reversals, order-three masks `0` and `6`
both have full curl `K_012=1` and the same oriented arc `0->1`. They are
different tournaments. Thus the deleted ambiguity would have supported
a false inference; the current theorem and repaired sentence do not.
This was a draft clarification, not a defect in the coefficient proof.

## Independent computations and replay

[check_information_audit.py](check_information_audit.py) imports neither
the producer nor the parent implementation. It uses three independent
paths:

* Direct expansion of degree-two terms over **all directed vertex words
  through order eight**, with their exact reference signs. This avoids
  both the producer's subset-DP and the parent's path-subcube/FWHT method.
  It checks incident and disjoint edge pairs. The magnitudes at orders
  `3,4,5,6,7,8` are `1/2,1/2,3/4,3/2,15/4,45/4`.
* Literal Hamiltonian-word enumeration and **actual complete response
  arrays** for every labelled tournament through order five. Their fibres
  are exactly converse pairs. It computes all triple curls, rather than
  only the triangles anchored at vertex zero. Aligned curl/response pairs
  distinguish all 1,024 labelled tournaments at order five.
* Independent directed reachability determines every SCC, reconstructs
  all cyclic rotations, and compares these with every full-curl fibre
  through order five. Every valid path-start key identifies one tournament.
  Direct isomorphism canonicalization of actual response arrays through
  order four gives four aligned types but only three separately quotiented
  pairs at order four. It does not use the response/converse theorem as a
  proxy for these response-array comparisons.

The independent run has 1,877 always-active checks. Normal and optimized
outputs are byte-identical:
[normal receipt](information_audit_receipt.json),
[optimized receipt](information_audit_optimized.json), SHA256
`128392e7e32a513123fbd02c59283751c178b9ea4501fb726b77109bfd2c6e80`.

```text
python check_information_audit.py > information_audit_receipt.json
python -O check_information_audit.py > information_audit_optimized.json
```

The refreshed producer checker was also replayed successfully with 667,959
checks, including the path-start and six-cycle topology additions. The
parent checker independently replays the order-five historical fractions
`165/32`, `1935/112`, `135/8`, `129/133`, and `18/19`. Its exact coefficient
probe checks magnitudes; this audit additionally checks every coefficient's
orientation sign through order eight. The historical Fisher calculation,
the different projection spaces, and the cube-versus-one-skeleton Euler
distinction agree with the elementary formulas in the current packet.

Reviewed producer checker SHA256:
`d6809002bdbc1a3cab3a73e673cf8abd818c80f5a6958bd0e0e8dee7dbeb79d7`.
The proof remains all-order because of the arguments above; these bounded
computations are independent controls, not a replacement for those arguments.
