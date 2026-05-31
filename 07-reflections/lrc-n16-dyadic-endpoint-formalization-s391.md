# LRC n=16 Dyadic Endpoint Formalization

codex-2026-05-31 S391

The useful move this session was to stop treating "endpoint debt" as a mood
and count it exactly.

For a pure dyadic owner `u=2^k` at denominator `16`, the endpoints are

```text
(16m +/- 1)/(16u).
```

A protector `p` sees them only through

```text
p*(16m +/- 1) mod 16u.
```

That makes the local law almost embarrassingly 2-adic.  If
`p=2^j q mod 16u`, then the count is decided by `k-j`, with only the two
shallow cases remembering `q mod 16`:

```text
super-gate: all endpoints
j>=k, non-super: none
k-j>=3: fixed deep count
k-j=2: only q=+/-1,+/-3 mod 16
k-j=1: only q=+/-1 mod 16.
```

This became THM-367.

The surprising part was not the formula; it was what the formula says about
the proof shape.  The owner `u=16` is the first pure dyadic layer whose lower
protectors can cover every endpoint.  That local closure costs exactly nine
residues:

```text
1,3,5,7,8,9,11,13,15.
```

At `u=16`, those nine are not just sufficient.  They are forced by private
endpoints.  This is the clean private-leaf picture S390 was hoping for.

But one layer higher, the picture changes.  Owners `u=32` and `u=64` still
have exact lower-cover number `9` in the bounded solver, and there is a
constructive nine-cover for all pure dyadic owners `u>=16`.  The private
endpoint certificate is a first-layer phenomenon, not the whole proof.

So the global proof target sharpened:

```text
not "every dyadic owner has a private leaf",
but "every attempt to close a dyadic owner exports a debt flow that cannot
balance across only fifteen speeds".
```

This is the Cayley-Dickson analogy in its more honest form.  The sedenion row
does not merely break by exposing one obvious zero-divisor.  It creates a
self-similar closure move whose norm has to be tracked globally.  The local
private leaf at `16` is the first crack.  The higher dyadic covers are the
reason the proof still needs a conservation law.

The odd payload scans point in the same direction.  For `u=16w`, the maximum
capacity in each protector `v2` layer scales by `w`, but the actual counts
inside shallow layers vary by odd residue.  That smells like a CRT-side
constraint rather than a plain tensor product.  In pure `n=16`, there is no
odd torsion in the denominator, but speed owners can still carry odd payloads,
so the future debt ledger has to record them.

Current best proof object:

```text
dyadic endpoint-flow network
vertices: owner layers and labelled endpoint classes
edges: strict protection arrows with THM-367 capacities
potential: demand - incoming capacity, weighted by maximality and half-turn
target: every labelled endpoint cycle has positive divergence.
```

THM-366 gives the entry gate.  THM-365 says a counterexample must contain a
labelled endpoint cycle.  THM-367 gives exact dyadic capacities.  HYP-1859 is
now the missing conservation law between them.
