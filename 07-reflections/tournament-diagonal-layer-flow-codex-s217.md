# Tournament diagonal-layer flow - S217

The useful outcome is that the user's `k^2+k` inter-layer lines have a precise
rank law.

Between diagonal layers of sizes `k` and `k+1`, the lines form
`K_{k,k+1}`.  If a line records an XOR/equality observation between tile bits,
the system is the `GF(2)` coboundary `delta(x)(uv)=x_u+x_v`.  The block has
`k(k+1)` lines but rank only `2k`.  A spanning-tree basis of `2k` lines
determines every other line by

```text
L(a,b)=L(a,b0)+L(a0,b)+L(a0,b0).
```

The surplus `k(k-1)` dimensions are exactly rectangle laws:

```text
L(a,b) + L(a,b') + L(a',b) + L(a',b') = 0.
```

That is the algebraic efficiency.  Do not count the duplicated lines; classify
the cycle-space generator.  Globally, adjacent-layer flow has rank `#tiles-1`,
so it reconstructs all tile bits up to global complement.  For a full
`n`-tournament this is `2*C(n,3)` line observations compressing to
`C(n,2)-1` independent coordinates.  The redundancy decomposes as
`2*C(n-1,3)+C(n-2,2)`: local rectangle cycles plus extra hourglass cycles
linking adjacent layer bridges.

The A000568 connection is now cleaner.  A000568 is the orbit count of all
binary upper-triangular half-tilings under the full `S_n` relabelling action.
A fixed Hamiltonian path gives a smaller cube with `C(n-1,2)` free nonpath
bits, and it surjects onto A000568 by Redei's Hamiltonian-path theorem.  The
fiber over a class `T` is `H(T)/|Aut(T)|` for ordered Hamiltonian paths,
verified through `n=6` by the S217 script.  So fixed-path half-tilings are a
presentation atlas with uneven fibers, not the quotient.

The diagonal-symmetry guess also sharpens: path reversal plus converse is a
real `Z2` quotient on fixed-path bits, but it is far too weak to be A000568.
The quotient sizes already exceed A000568 at `n=4` (`6` versus `4`), `n=5`
(`40` versus `12`), and `n=6` (`544` versus `56`).  It is useful as a
controlled-forgetting sidecar, not as the final orbit law.

For LRC the translation is direct.  Inter-layer lines should carry endpoint
owners, barcode-bar support, active bottleneck owners, route/status labels, or
proof obligations as sidecars.  If every rectangle law stays zero, the sidecar
descends to layer potentials.  If a rectangle law fails after adding LRC
payload, the failed rectangle or hourglass residue is the hidden coordinate,
exactly like the pair-good decoy generator and Haar/fixed-margin `zeta`
lessons.  Caveat: the rank law is only for the XOR/equality shadow; once a
line carries owner labels, direction, or active bottleneck data, the graph is
only the skeleton and the payload needs sidecar equations.

Assumption challenge: I considered runners, arcs, diagonal layers, line
carriers, rectangles, fixed paths, path-reflection orbits, A000568 classes,
edge/cycle/clique perspectives, barcode bars, active bottleneck owners,
endpoint-owner packets, and proof obligations.  The vertices that served the
route were proof carriers.  The preserved predicate is parity/fiber flow plus
the LRC controlled-forgetting rule; the destroyed data is labelled identity,
specific Hamiltonian-path presentation, and owner/route/status coordinates
unless attached as sidecars.
