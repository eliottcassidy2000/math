# Pair midpoints explain the shifted radix digits

**Status: elementary proved size algebra, independently reviewed. This is not a new existence proof for graph pairs.**

Here odd/even order means the number of vertices, not the number of edges. A Shorter/Longer pair has orders n and n+1 and midpoint x=n+1/2. Its integer-order parity and its per-label positional phase are different data.

Let Q=2h+1 be odd. The certified square-sum lifts use an output order

    N=Q*n+d,    h<=d<=h+Q-1.

Setting r=d-h gives exactly

    N+1/2 = Q*(n+1/2)+r,    r in {0,...,Q-1}.

Thus the apparently shifted digit window is the ordinary digit range in **half-integer midpoint coordinates**. The half-integer lattice is preserved by odd Q and integer digits. An even Q would send such midpoints to integers under this particular formula; this is not an impossibility theorem for every conceivable even-radix graph construction.

For Q=49 the raw digits are24,...,72; for Q=25 they are12,...,36. Carrying from the last digit to the next parent is a statement about these size coordinates, not a guarantee that chosen path realizers fit together cheaply.

## Uniform threshold formula

A zero digit maps a base index a to Q*a+(Q-1)/2. Iterating l times gives

    T_l = ((2a+1)*Q^l-1)/2.

For mixed radices replace Q^l by their product. Size products commute, although the actual certified prefix operators may depend on their order.

Examples:

- Original nice-pair base a=41, Q=49: first recursive index2033.
- Weak-B/ninja base a=42, Q=49:2082.
- Weak-B/ninja base a=42, Q=25:1062.
- A merely formal base a=32, Q=49:1592.

The last numerical identity supplies no constructor: both actual pair classes have empty fibers at32, despite G32 having a Hamiltonian cycle. The necessary initial graph certificate cannot be supplied by arithmetic relabeling.

The constants83,85,107 and225 appearing in the inherited thresholds are simply2a+1 for base indices41,42,53 and112. This explains their common algebraic form; it does not independently verify the relevant graph interfaces or their base data.

## Parity bookkeeping

If Q is an odd square, Q is1 modulo8, so h=(Q-1)/2 is divisible by4. Therefore

    N modulo2 = (n+r) modulo2.

The output pair still has opposite order parities. But compatibility additionally requires the positional parity of every common vertex label to agree between its two paths. The latter cannot be reconstructed from this single size bit.

Chung–Graham's rounded shift ceil(phi^2*n) is a separate numeration operator with a separate legality condition. Both systems can address the same integers; their digit operations are not interchangeable graph-building instructions. Ordinary successor, numeral carry, and graph-witness repair remain distinct operations.
