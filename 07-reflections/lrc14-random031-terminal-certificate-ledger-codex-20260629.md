# LRC14 Random031 Terminal Certificate Ledger

The useful new angle was to stop treating the latest random031 objects as
parallel descriptions and force them into one terminal ledger.

HYP-3486 says `242` cells hit endpoint-rank-2 gates and `40` are free holes.
HYP-3511 says the `40` free-hole cells are locally bracketed.  HYP-3510 says
the whole seam complement is connected after branch order and survivor ports
are retained.  HYP-3490 says projection-current deletion is blocked.  HYP-3521
joins those facts and makes the sharper terminal partition:

```text
282 = 230 ordinary rank-2 route cells
    +  40 bracketed free-hole cells
    +  12 pure bypass owner-boundary cells.
```

This clears up a small but important ambiguity: the older `242` routed-cell
count is correct, but it is `230 ordinary + 12 bypass`.  The bypass should not
be swallowed by the ordinary endpoint-rank lemma, because it is the only
gate-routed hard-component owner-boundary certificate.

The free-hole packet also became smaller as a formal target.  There are `14`
free-hole components, but after HYP-3511 doublet collapse there are only `12`
free-hole certificates:

```text
10 ordinary-bracketed singles
2 same-branch doublets, with cell sizes 4 and 10.
```

The fiber PGF sidecar got sharper too.  The `24` zero-exit fibers are pure
free-hole fibers, while the mixed `free_hole+ordinary` fibers have one escaped
ordinary sheet and one bracketed no-gate sheet.  That suggests a clean lemma:
an ordinary escaped sheet can bracket its missing same-fiber companion without
using the forbidden seam.

The proof target now feels like:

1. Prove the ordinary-route lemma for the `64` ordinary components.
2. Prove the free-hole bracket lemma once for singles and once for doublets.
3. Prove the pure bypass owner-boundary lemma for the single `12`-cell packet.
4. Use HYP-3490 to forbid falling back to projection-current deletion.
5. Use HYP-3486 to forbid vertical half-turn quotienting unless class sidecars
   are retained.

That would turn random031 from a named hard packet into a compact terminal
certificate theorem: `64 + 10 + 2 + 1 = 77` terminal certificates.
