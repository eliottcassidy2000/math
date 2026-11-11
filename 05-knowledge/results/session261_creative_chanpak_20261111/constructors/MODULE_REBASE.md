# Profile-preserving module rebases, hidden port menus, and content barriers

**PROVED:** the all-host marked replacement theorem, content obstruction, and
stable-port obstruction below. **FINITE-EXACT:** the explicit profile collisions
and bounded host searches. **OPEN:** embedding the trade into a positive unit
seed and universal same-q rearming. No Chan–Pak completion is asserted.

## 1. Interface and inherited input

Use Session258's selected seed statistics H,W,U, with D=H−W and multiplier
q=U. All paths in W have selected first endpoint and second start; U allows
an empty first path and requires a selected second start. The target is D=1.
The inherited unmarked substitution kernel is THM-3121,
`01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`.
The marked extension and the consequences for rearming are given here.

For a nonempty tournament M, let C_r(M) count **ordered nonempty** spanning
r-path covers. Thus C_1=H and C_r=0 for r>|M|. This convention differs from
Q_2, which allows empty slots: Q_2=C_2+2H.

**Marked module replacement theorem.** Suppose M is an unselected module of
(T,O): O is disjoint from M and every exterior vertex has a uniform orientation
toward M. Replacing M by any M' with C_r(M')=C_r(M) for every r preserves
(H,W,U) exactly. Hence it preserves both the unit-seed property and exact q.

Proof: compress each maximal nonempty run of module vertices in a Hamilton
path, or in either of the two paths defining W or U, to a module symbol.
The exterior vertices, path break, and selected endpoint conditions remain
fixed. No endpoint condition specifies an interior module vertex, because all
marks are exterior. If there are r run slots, their order in the first path
then second path labels them canonically. Expanding them is precisely an
ordered spanning r-path cover of M. Thus each statistic F has the form

    F(T,O)=sum_{r>=1} a_r(F) C_r(M),   a_r(F) nonnegative integers.

The sum is finite. Empty first paths in U introduce no exception. The same
coefficients apply to M', proving the assertion. This is a physical graph
replacement, not a relabelling of scalar response.

**Content obstruction.** Put d(M)=gcd(C_1(M),...,C_|M|(M)). Then d(M) divides
H,W,U and D whenever M is an unselected module. In particular every unselected
module in a unit seed has d(M)=1. More sharply, only the reachable run numbers
need be retained in this gcd.

**Stable-port obstruction.** Suppose replacing an unselected vertex v by TT_m
keeps D=1 for every m>=1. Replacing that same vertex by an arbitrary tournament
M instead gives D=H(M). In particular a nontransitive module cannot be inserted
there while retaining a unit seed.

Proof: the signed skeleton coefficients d_r=a_r(H)−a_r(W) give
D(TT_m)=sum d_r r! S(m,r). For m=1,2,... the triangular surjection basis
forces d_1=1 and d_r=0 for all r>=2 (the skeleton has bounded run number).
The same coefficients evaluated on M give C_1(M). The deletion boundary m=0
is separate and is not used. This explains why a new module trade cannot be
installed at an old stable port merely because it preserves a scalar profile.

## 2. A genuine profile trade that changes the future port menu

Outgoing-bit rows, zero-based vertices:

    A=(2,4,9,19,39,15),
    B=(2,4,57,3,11,27).

Both have the exact full cover vector

    (C_1,...,C_6)=(17,190,828,1752,1800,720).

Its gcd is1, so it survives the content test. Yet the rooted clone responses
differ. Replacing vertex1 of A by a transitive m-chain gives

    H_A(m)=9*2^m−1,  m>=1,

whereas replacing vertex0 of B gives

    H_B(m)=8*2^m+1,  m>=1.

At m=1 the modules agree at17. The first physical clone births give35=2*17+1
and33=2*17−1 respectively. Thus a profile-preserving rebase really can change
which arithmetic operation a subsequently chosen port implements. This
escapes the full-reversal symmetry: reversing a whole module transports all
its rooted transitive-inflation profiles unchanged, up to vertex relabelling.

The all-length formulas are certified in the finite run basis: for A's port,
c_1=17,c_2=9,c_r=0 for r>=3; for B's port c_1=17,c_2=8,c_r=0 for r>=3.
Their deletion values are5 and9 respectively, separate from m>=1. The checked
certificate is exact, not an extrapolation from the displayed samples.

This is not yet a unit-seed rearm: the replacement theorem needs a unit host
with marks exterior to the traded module. Exhausting the82 rooted quotient
presentations on3,4,5 vertices (all unlabelled quotients, every root, every
exterior selected subset) found no such host for A. This is only a bounded
negative result. Moreover B itself has a unit marking O={1,2}, q=48, whereas
A has no unit marking. Marks inside a module cannot be discarded from the
replacement interface.

## 3. First hostile and wider positive evidence

The first collision found by the differing-port-menu test at order5 is

    (2,4,25,3,11) and (2,4,9,19,7),
    C=(9,66,186,240,120).

Its gcd is3. The content theorem excludes its embedding as an unselected
module in **any** unit seed, at every host size. This separates the failed
first candidate from the gcd1 six-vertex candidate, which remains open.

The order6 census has56 unlabelled classes and24 full cover profiles. Among
gcd1 profiles, three collision fibres have different rooted inflation menus:
H=17,29,37. Explicit members and rooted values are in
`six_module_collisions.json`. This is a small source bank of exact physical
replacements, not a claim about arithmetic coverage.

## 4. Provenance, transfer ledger, and next gate

Anchor: exact multiplier return outside the old L/R/G carrier. Niche: modules
with the same ordered cover vector. Wildcard: boundary-count kernels from a
random historical source. Live board: unit defect / exact multiplier / marked
boundary / cover content / hidden root menu / physical substitution.

Random seed261771 sampled four files from the sorted2312 historical result
Markdown paths, excluding sessions260/261; `random_sources.json` records the
actual pool and paths. The useful unrelated source was
`unit_equality_quantized_guards_20261011e.md`: its boundary-count matrix and
integer quanta suggested retaining the full extension-count interface and
checking divisibility before seeking a count-preserving switch. The random
simplex/cyclotomic note was used as an idea source only, not as proof of its
historical broad claims. The parent additionally suggested module-flip fibres
from the corrected intransitivity work.

Transfer ledger: source = ordered path-cover module profile; target = marked
seed host; map = literal module replacement; preserved predicate = exact
(H,W,U), provided all marks are exterior; lost state = root identity and future
rooted response; sidecar = vertex-labelled clone profiles; cheapest test =
profile equality, then content gcd, then a unit host. The first two tests pass
for A/B; the actual unit host is still missing.

Two concrete next gates remain: find a gcd1 trade with an admissible unit host
beyond quotient order5, or retain selected-endpoint tensors inside traded
modules so a mark can cross the interface without losing D=1. A joint trade
may cancel marked defects even when either individual trade cannot.
