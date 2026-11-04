# Independent audit of one-exit transport and persistent old names

**PASS — 2026-11-03.** The coordinator independently checked the proof and
standalone source of
[ONE_EXIT_AND_PERSISTENCE.md](ONE_EXIT_AND_PERSISTENCE.md). The accepted
results are analytic; the finite controls support, rather than replace,
their unbounded quantifiers.

The expansion at a sole-exit auxiliary partitions paths according to
whether they end there or continue immediately to its unique out-neighbour.
Coefficient extraction gives sigma=e_b(T−a) and
tau=[W_w E_b]_(T−a), with two nonempty ordered ending slots. The repeated
endpoint case vanishes. For distinct endpoints, cutting an ending-b
Hamilton path immediately after w is injective and proves tau>=sigma.
No strongness assumption was imported into this identity.

The inherited endpoint-monotonicity and permanent-certificate statements
were checked in their actual source, NOTE §§I1–I2 of
[`balanced_ternary_returns_20261103/interfaces`](../../balanced_ternary_returns_20261103/interfaces/NOTE.md).
The zero alternative of endpoint monotonicity is retained. Each of the
three certificate types covers its stated future-degree cases, and all
named triangles are disjoint from the deleted names. For the actual family,
the witness endings 28*2^k+15 and224*2^k−44 are correct, as are H_(j−1)
and2H_(j−1) for subsequent pairs. The saved-vertex obstruction therefore
holds at every clock and every zero-channel depth under every induced
extension. It does not exclude new auxiliaries or changed old arcs.

The clone block formula is inherited and correctly attributed. Its
application gives H_a=h+2^a(h+2), F_a=(1+2^a)f; the physical step at
a=1 to3 is exact. The exponent equation has only that positive solution.
All clone vertices have the listed permanent certificates. This is one
continuation, not reusable closure.

The source
[`algebra_one_exit.py`](../../../../04-computation/nonzero_preparation_20261103/algebra_one_exit.py)
was read and independently rerun in normal and optimized modes. Both pass
**44,726 always-active gates**, and both outputs match the frozen
[`algebra_one_exit.json`](../../../../04-computation/nonzero_preparation_20261103/algebra_one_exit.json)
byte for byte. The source/receipt SHA256 pins are respectively
`a2fa146f99ed10102d4a41cc699967575d8d2d2bdf726d820b72df4a253218c0`
and `9e6a4dceca0bdaaf72697966313440331cc27764f68146ac49fc26cbb70c6570`.
The finite universe includes all5434 one-exit marked labelled cases through
order5, actual clocks3–5, zero depths0–2, clone lengths1–5 and every one-birth
cut on the first actual eleven-name output. The last cut universe checks
surviving certificates, not all Hamilton responses; that distinction is
stated correctly.

The independently authored no-import
[`root_audit.py`](../../../../04-computation/nonzero_preparation_20261103/root_audit.py)
additionally tests314 marked one-exit cases through order4 using literal
permutations and the two actual C2 births to recover the joint response.
Normal and optimized modes agree with its frozen receipt; the combined
induced/algebra audit has154382 gates. Reproduction commands and pins are in
the adjacent [induced audit](../induced/AUDIT.md).

No complete arbitrary-auxiliary classification is attributed to this note;
that extension has a separate proof and audit. Status/link-only promotion
does not change this audit's source or receipt pins.
