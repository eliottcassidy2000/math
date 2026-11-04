# Additional independent bit-width and gate audit

The reviewer independently checked NOTE.md's distinction between arithmetic gates, bilinear rank, and exact state bits. No mathematical correction is required.

## Bit-width proof and explicit implementation

After i inputs, P is the represented i-bit integer, so P<=2^i-1. Before step i, P+L<=2^(i-1). Starting at h=1 and summing these nonnegative possible increments gives h<=2^i. Therefore registers of width n+1 suffice for every intermediate value in a fixed-length n-input circuit. These are O(n) bits, not constant-sized Boolean signals.

An explicit, conservative implementation uses two width-(n+1) ripple adders, one complemented input bit, and n+1 AND gates per step. Each ripple bit uses five B2 gates: p=x XOR y, g=x AND y, t=p AND carry, carry'=g OR t, sum=p XOR carry. P'=2P+c is a wire shift and insertion, and L'=c is a wire. This gives n*(11*(n+1)+1) gates with constant zero/one wires available. If the exact gate convention requires constructing constants, for n>=1 a shared b_1 XOR b_1 and its negation add at most two gates. This is a multioutput binary encoding of H, not a new single-output NP characteristic-family bound.

The implementation is intentionally conservative, with unoptimized gates and zero-padding. It confirms the claimed O(n^2) upper bound; no optimality is asserted. Product by 1-c is bitwise gating, not a generic multiplication of two n-bit integers.

## Independent finite controls

`independent_bit_audit.py` does not import the producer's checker or recurrence. It implements literal ripple carries and compares their integer output with the apex-deletion combinatorial formula for all 2,047 binary words through length10. It verifies the exact gate inventory and absence of every carry/shift overflow. Separately, it uses the combinatorial formula to verify both suffix probes and their difference 2P on all 8,191 prefixes through length12.

Normal and optimized outputs match in `independent_bit_audit.json` and `independent_bit_audit.optimized.json`. Explicit exceptions keep all controls active under optimization.

The complete derivative graph is not a gate count. Exact cut rank two concerns sums of products of linear forms over Q only. The one-way k-bit lower bound concerns deterministic exact recovery with two suffix inputs, and is consistent with the above quadratic Boolean upper bound. No transfer to unrestricted B2 lower bounds is inferred.
