# Residue Feedback Loop

**Session:** kind-pasteur-2026-05-30-S3
**Status:** synthesis pass after random thread walk
**Related:** HYP-1779, HYP-1780, THM-351, THM-352, THM-354

The random walk landed on five threads that did not look like the same object:

1. Boolean-cube bucket checksums.
2. Good-cut height and SCC condensation.
3. OCF-derived Redei parity.
4. H=63 / THM-025 odd-cycle residues.
5. Paley versus Interval support shadows.

The common shape is a projection followed by a surviving residue.

## Bookkeeping Versus Residue

The bucket-balance story is the cleanest bookkeeping case. Boolean xor by
nonzero masks is reversible and fixed-point-free, so the row checksum is forced:

```text
2*internal + escape = bucket_size * mask_count.
```

Once that checksum is removed, all remaining structure is in the residual
transport distribution: which buckets receive escape mass, which channels are
spine/ribs/sea, and which lines are invisible to one quotient but visible to
another.

This is a useful separation. Conservation is not the interesting data. The
residue after conservation is.

## Good Cuts: Coordinate Becomes Defect

The good-cut count used to look like a coordinate artifact of a chosen base
path. THM-354 changes the meaning:

```text
g_P(T) = n - #SCC(T).
```

So good cuts are not a tiling decoration. They are the residue of SCC
condensation. Bad cuts are exactly the component boundaries. This turns the
old "merged classes are pure in g" experiments into a theorem: SCC count is an
isomorphism/complement invariant.

It also changes future transport language. A nonzero `Delta g` is really a
nonzero `Delta #SCC`, hence it cannot stay inside a merged tournament class.

## OCF: Parity as the First Residue

The OCF thread now has the same flavor. Lean already proved `H` is odd from
OCF. This session added the explicit residue:

```text
H(T) % 2 = 1.
```

That is the first projection of the OCF partition function. Higher questions
then ask what survives after stronger projections: which alpha-vectors are
possible, which residue profiles would force H=7 or H=21, and where the first
nontrivial independent packets appear.

The phase-transition theorem `alpha_k > 0 -> H >= 3^k` says the first
nonempty k-packet residue costs exactly `3^k`. That makes the H-spectrum look
less like a list of forbidden numbers and more like a rank filtration.

## Exact Kill, Near Kill

H=63 and THM-025 are the sharpest contrast.

For H=63, deleting the core vertex kills all odd cycles. The residue is empty,
so `Omega=K31` and `H=1+2*31`.

For THM-025, deleting a near-core vertex kills almost everything, but a tiny
two-cycle residue remains. That small survivor has enough disjointness to keep
the unique independent triple and break real-rootedness.

The difference between exact kill and near kill feels like a useful diagnostic:
empty residue gives rigid complete conflict graphs; small dangerous residue
gives failures of naive positivity or real-rootedness.

## Paley/Interval: Same Shadow, Different Fiber

Paley and Interval on 7 vertices share a support shadow but differ in
multipity and disjointness. The visible projection is the same; the residue in
the fiber is not. This is the same pattern in a softer form:

```text
same shadow, different residue geometry.
```

That suggests a feature hierarchy for TDA work:

- first record compressed invariants (`H`, score, support shadow);
- then record residue features (support excess, max multiplicity, deletion
  loss, SCC count, transport escape rank);
- then ask which failures or gaps appear only at the residue layer.

## Feedback Rule

When a new invariant appears:

1. Ask what support family it is counting.
2. Identify the projection or quotient that makes it cheaper.
3. Measure what survives after the projection.
4. Use the residue to generate the next theorem or feature.

This loop produced THM-354 from good cuts, and it suggests HYP-1780: the next
obstructions should be stratified by residue rank, not by raw support size.

## S355 Loop: Rank Becomes a Feature

The next pass made the slogan operational in `tournament_tda.py`.

The extractor now records both the size and the shape of the max-loss deletion
residue:

```text
keep_v*,
alpha_1(R_v*),
alpha_2(R_v*),
rank_res(v*),
I(R_v*,2).
```

The first probe gives a clean three-way split:

- H=63 single-core classes are exact kills: `rank_res=0`.
- THM-025 is a small dangerous near-kill: `alpha=(2,1)`, `rank_res=2`.
- Paley/Interval are broad fiber residues: rank 2 appears, but not at small
  residue size.

That creates HYP-1785: first Omega obstructions should be enriched among
small nonempty max-loss residues with `rank_res>=2`, while complete-Omega
unlocks remain in the exact-kill stratum.
