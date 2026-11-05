# Independent review: fixed-five-port binary–ternary bridge

**PASS. No correction required.** Reviewed 2026-11-05. This review binds the frozen files in `session197-hybrid-constructions/binary-ternary-bridge/`:

- `PROOF.md`: `65c57c4f232dc77713b98048765ac901e4039ee761d518066817c84395ef3fc7`
- `bridge.py`: `48a18bea9eee3d09cf839b3181de6aa43584aa4b30e9959179c8eff65048ca7d`
- `graphs.py`: `3204eeb1e9d2b05cac9bcc21c98f25697d10f925d6ff743393e945ce4e73615c`
- `verify.py`: `62e36a8b1edb08031210d8faa86beac00813e62c8e54848697096fe0b9449072`
- `VALIDATION.json`: `e6229bb9d0fe2555166a19ce831890c187202d66a0831d4e450dba611810379d`

The original files were not changed. Source replays and independent controls are in `binary-ternary-audit/` beside this report.

## Mathematical review

1. The port order `[0,1,4,5,6]` makes right mask7 mean that the newborn beats exactly0,1,4. Left mask24 makes it lose exactly to5,6. Drop5 forgets only the newborn as a port. Every old port persists, and no operation silently re-exposes a forgotten vertex.
2. The right newborns, together with original7, form an increasing transitive module. Left newborns form a reverse-birth-order transitive module and beat every right-module vertex. Thus the explicit relabeling is correct for every interleaving. In particular an L newborn beats an R newborn in either chronological order: this is physical graph commutation, not merely equality of scalar formulas.
3. The root crossing graph has precisely the old edges1->2,1->3,4->5,4->6, covered by singleton vertices1,4 uniformly in the bridge-module size. The bridge crossing graph is covered by the old singleton4 uniformly in root-module size. Every internal fragment gives a disjoint crossing edge, so these covers bound positive run indices by3 and2 respectively. They do not count an entire variable module as a singleton cover.
4. Surjections onto ordered nonempty transitive runs give j!S(s,j). For positive s, the spans are therefore bases1,2,3 and1,2. Recovering coefficients in one variable by an invertible two-point matrix establishes the joint tensor span; the second Vandermonde has determinant12. Thus the six physical counts identify an all-order theorem, not an extrapolation from a grid.
5. The boundary t=0 is different: the deleted-module contribution is not recovered by evaluating the positive-variable exponential expression at zero. The negative run vector gives60*3^m-1 there. The first L operation is exactly2H+3. Only after t>=1 are the returning operations2H-1 and3H-2. Their affine compositions commute and give2^l*3^r*(H-1)+1.
6. The source has the displayed Hamilton cycle. Adding the first bridge vertex has both incidence signs with the strong source. Nonempty transitive substitutions preserve strongness, including m=1 and t=1. No source/sink padding is used.
7. Unique factorization gives the exact interior criterion N-1=2^a*3^b*5 with a>=3,b>=2, and parameters m=b-1,t=a-2. Graph order7+m+t is correct. This is a complete inverse only for the displayed family. Neither a reverse sign gate nor unrestricted spectrum universality follows.

## Independent computation

The independent checker does not import the producer's graph builder, birth routine, endpoint DP or run DP.

- Rebuilt every graph as a quotient blow-up of the literal nine-vertex source.
- Counted19 graphs, including the complete tensor grid, boundary cases and off-grid pairs through order14, by signed inclusion–exclusion over all walks of length n-1. This is distinct from the source subset endpoint recurrence. All counts agree with60*2^t*3^m+1 for t>0 and60*3^m-1 for t=0.
- Independently solved the six-by-six rational tensor system, obtaining coefficients `(1,0,0,0,0,60)` in the ordered bases `(1,2^m,3^m) tensor (1,2^t)`.
- Enumerated permutations of the outside vertices and all permissible/forced root slots combinatorially. Recovered complete run vectors `(59,179,180,60)`, `(121,361,360,120)`, and `(179,361,180)`, including the separate deletion entries and all higher zeros.
- Checked canonical graph identities and strongness for all255 R/L words through length7, with the original port labels fixed.
- Checked820 crossing-cover instances with module sizes through20 and the exact inverse criterion for every integer through200000.
- Copied-source normal and optimized Python replays match each other and the frozen validation receipt byte-for-byte.

These finite controls are checks of the stated all-order structural proof. They are not the basis for extending the formulas beyond the grid.
