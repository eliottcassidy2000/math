# Independent audit of the ordered-end obstruction

2026-11-03. **INDEPENDENT ANALYTIC / SOURCE AUDIT PASS.**
The theorem is unbounded in both old and new tournament sizes. Its scope is
exactly the two physical operations in [PROOF.md](PROOF.md), with old port
memberships preserved, one common nonempty transitive prefix, and strong full
attachments. It is not a no-go theorem for all native or internal constructions.

## Analytic review

I read the full proof and source. The S182 signature and attachment conventions
are retained, especially the `eS` and `gF` terms in the joint update and the
actual inner defect `j`. The terminal-SCC strongness condition is used on the
correct block: appending requires a terminal port in the new block, while
prepending leaves the old terminal SCC unchanged.

The actual cut injection gives `f>=|W|h`. Rédei's odd Hamilton count and the
fixed-point-free slot swap give odd `h`, even `p`, and `p>=2h`. Excluding an
empty first slot gives `f<=p-h<p`. These are constraints on actual tournaments;
they cannot be imposed only after an abstract signature calculation.

The append equation `(H-3)f+pF=0`, with `F>=H>=1` and odd `H`, forces `H=1`.
Then `pF=2f` and `0<f<p` force `F=1,p=2f`. The proof that `h=1` implies a
transitive tournament is valid: insert the remaining vertices into each of
the three distinct directed-triangle Hamilton orders, retaining the old
relative order. A transitive terminal-port block has `F>=2^(n-1)`, so `H=F=1`
forces the marked singleton. The converse and displayed full updated state
are exact.

The prepend equation `(P-3)f+Fh=0`, with positive even `P`, forces `P=2`, hence
a singleton. Its unmarked case cannot satisfy the equation. The marked case
is exactly `f=h`; the cut injection then gives one old port and consequently
`j=0`. The displayed negative digit and full updated state follow.

Both first steps have `(h',p',f')=(h,2p,3f)`. A second prepend would require
`3f=h`, impossible. A second append would require `p=3f`: this contradicts
the first append's `p=2f`, or the first prepend's `f=h` together with even
`p` and odd `h`. This exhausts both adaptive choices without a block-order
cutoff and without a restriction on the additive digit. The theorem does
not assume a fixed or zero inner defect.

The actual marked singleton gives a first step of digit -4, so the theorem is
not vacuous. The inherited +4 example changes an old port. The inherited +2
example retains old ports but inserts a cyclic block through an old vertex,
creating a new-to-old arc forbidden by an ordered-end extension. Both are
correctly kept as scope hostiles. The formal odd-p example shows precisely
where realizability is used. No global sign assertion about append digits
is inferred from the finite census.

## Independent finite reconstruction

I inspected all of `closure.py` and its frozen producer replay record. The
source uses always-active exceptions and distinct literal-permutation and
subset-DP implementations. I did not rerun its 4,066,300 raw extension
comparisons, because the all-size proof and a different exact reconstruction
already pay that obligation.

The independent
[`audit_component.py`](../../../../04-computation/multiport_grammar_20261103/audit_component.py)
now reconstructs the finite first-step census using its own literal
permutation/cut endpoint counts and separate reachability. It recovers:

| Quantity | Independently recovered value |
|---|---:|
| Actual tournament/port states, orders 1–5 | 33,866 |
| Terminal-port guarded states | 27,475 |
| New-block states through order 3 | 74 |
| Guarded new-block states | 43 |
| First append states | 775 |
| First prepend states | 311 |
| Raw extension universe size | 4,066,300 |
| Guarded second-extension universe size | 127,062 |

The independently recovered append offset multiplicities are
`{-42:120,-20:120,-16:144,-12:88,-10:150,-4:153}`; prepend multiplicities are
`{-48:40,-30:120,-28:120,-16:24,-10:6,-4:1}`. They match the producer output.
The reconstruction applies the independently audited exact first-step iff;
it does not purport to be a second raw extension scan. Separately, its 484
literal two-atom products verify all seven transfer coordinates on actual
graphs, including both joint cross terms.

The independent verifier's normal and optimized outputs agree, with 18,770
gates across its full component and closure audits. The producer reports
213,973 always-active gates in each mode; its two frozen 844-byte outputs
have identical bytes and the declared hash. Those producer mode runs are
credited to the producer, not misreported as a fresh independent full replay.

## Frozen pins and conclusion

| Artifact | SHA256 |
|---|---|
| Primary proof, pre-promotion | `c1f4478b145bba4d3f4095a14dee843765bac876493c8f559cf1bf6df38809b5` |
| `closure.py` | `3878a30e9e20ec1ae29caa5a523f0bd1df62a312c458c6e2014d11223cd6c962` |
| Producer normal/optimized output | `d26f2e4f531305d54134a5762c736b134d3ddd0e94a3cff695a1706e451b53c3` |
| Independent verifier | `d96549dccd67aba2c98f788f8ba3a9d0dfc81fc28b278172cc923666f04ed353` |
| Independent receipt | `b1719d7a7e787ef2947dfe8cddc4efb34fb2eb1eb5add80b1a78e2f9d3581d93` |

**PASS; no correction remains.** Any reusable two-branch grammar must leave
this exact ordered-end/old-port-preserving frame. The theorem supplies that
necessary design constraint, not a return construction or a global spectrum
classification. Status promotion and Git integration remain with the parent.
