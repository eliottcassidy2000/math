# Internal adversarial review and reproduction contract

2026-10-26. This is a proof review and separate computational formulations by
the producing session, **not an independent-agent audit**. No new numbered
canon theorem is reserved or promoted by this packet.

## Status and exact scope

- R1: **PROVED**, application of the classical two-weight discrete Hardy
  criterion, with a self-contained constant4 proof. Supremum bounds all input
  vectors; one observed current alone need not saturate it.
- R2: **PROVED**, Hilbert-space refinement identities and explicit hostile.
  Conditional expectation acts on observables; mass pushforward sums.
- R3: **PROVED**, a rational identity on every finite positive odd Collatz
  prefix, with an explicit n-dependent input cost. No uniform cheap input
  estimate, stopping result, or native ROOT observation is obtained.
- R4: **PROVED conditional extraction theorem**, assuming an actual energy
  budget and using actual source masses. The polynomial alias tail is inherited
  for the specific R6 prior. No replacement of native ROOT weights is permitted.
- T1: **PROVED**, all-order closure of rank/forest frames under marked
  relocation, the two binary gluings, and full converse. The stated compactness
  counts combinatorial bits only. Free-edge complement has an exact exception
  update; T2 proves it does not preserve the forest subclass even after reframing.
  The explicit seven-vertex witness needs exactly one exception cycle.
  Ordered joins prove exactly r necessary cycles on7r vertices after one flip;
  the certificate is checked for r=1,...,12 independently of the all-r proof.
- M1: **PROVED corollary** of audited THM-4808: order-sqrt(H) modulus scale for
  uniform positive endpoint-vector decoding of branch character. Scope does
  not include projective slope-only observers or unbounded height without a
  sidecar. No new spectral theorem is claimed.

## Exhaustive and finite controls

`verify.py` is dependency-free Python with exact integer/Fraction arithmetic.
It raises explicit exceptions on failure; optimization cannot remove checks.
Normal and `-O` runs must produce byte-identical JSON receipts.

| Lane | Declared universe and checks |
|---|---|
| Forest frames | All forests on the fixed rank order 0,...,n−1, n=1,...,6: counts 1,2,7,38,291,2932. Other labelled rank orders are transported by relabelling. Every permutation is tested for Hamilton legality; 57,495 legal paths. |
| Marked relocations | All cut/insert positions for those legal paths: 2,044,019 comparisons of at most three guards with a full path check. Exact arc-table reindexing and named source-weight/port pullback are checked for each legal path. |
| Composition | 200 ordered J/star joins of all fixed-rank forest frames of orders1–3, using deterministic insertion paths. Forest, cross arcs and concatenated-path legality checked. This finite test supplements the disjoint-union-plus-bridge proof. |
| Complement hostile | At every tested legal path, exact frame XOR update checked against literal edge reversal. The first fixed-frame forest failure is n=3, D={02}, P=012, updated D=K3. This does not test nonexistence of an alternative rank frame. |
| All-rank forest obstruction | Explicit forest predecessor of Paley7 under the legal path0123456; all5,040 rank orders independently checked for absence of a forest frame. Minimum feedback size7; explicit natural-order unicyclic repair. Fourteen triangles/two incidences per edge supply a separate all-ordering proof. Joined exact-r-cycle witnesses checked for r=1,...,12 (up to84 vertices). |
| Weighted renewal | All u,v entries in {1/2,1,2}, lengths1–4, plus 240 deterministic rational signed-gain examples of lengths1–12: 7,620 cases total. Exact positive-definite elimination of 4B*diag(v)−R*diag(u)*R and explicit lower-bound test vectors. The gain kernel is also built independently by interval products. |
| Refinement | 27 coarse and729 fine vectors, one declared split into2,3,1 children with nonuniform rational masses. E J=I and both weighted isometries checked. Equal-mass false microclock refinement checked by exact Rayleigh ratios. |
| Collatz expansion | Explicit source 2^(L+1)−1, L=1,...,32. All first L valuations equal1, endpoint2*3^L−1, clock2L. These are finite prefixes, not divergent orbits. |
| Passive energy | 32,768 rational single-branch identities across every odd n≤4095 and a16-element error/forcing grid;256 actual positive odd sources≤511, each followed for exactly32 odd steps, with a full rational telescoping identity and physical clock. The finite runs do not assume or test universal convergence. |
| Adaptive source | 192 finite alias controls: r∈{1,3,6}, source indices0–15, q∈{16,32,64,128}, fixed signed responses supported on0–127. Three nested schedules for G_t=9^t, epsilon_t=2^(−t), t=0,...,16, marked source m=7. |
| Berggren | Every fixed modulus M=1,...,256 in the inherited hostile word; all9,841 words through depth8 for pole determinant, its mod3 image, and dyadic positive-coordinate lifts. These do not redo the spectral census. |

## Hostile implications and repaired forms

1. **Mass conservation implies energy stability:** false already at two equal
   children, energy ratio5/2. First failure: treating refinement children as
   sequential transitions. Survivor: J R E at the same clock is isometric on
   copied inputs and norm preserving globally; actual fine defects need payment.
2. **One boundary scalar gives bounded gluing:** false, cut norm squared
   a(L−a) for unit-gain/unit-weight integration. Survivor: retain the scalar
   and its cut product. New question: does the actual signed current pay it?
3. **Relative contraction makes unweighted renewal uniformly bounded:** false
   on arbitrarily long all-one Collatz prefixes. First failure: contractions
   approach1 with the source height. Survivor: the paid passive identity, with
   cost (3n+1)^2/(6n+1), and the supplied-weight Hardy criterion.
4. **A forest frame remains a forest after the prescribed complement update:**
   false at the triangle. The stronger any-rank repair proposal fails at the
   explicit forest-to-Paley7 transition. Survivor: exact D symmetric-difference
   path update with a cycle budget; one unit is necessary and sufficient at P7.
5. **A finite residue can decide parity without a height/history coordinate:**
   refuted for every M by THM-4808. Survivor: endpoint precision of order sqrt(H)
   or the already proved fixed-modulus pole-frame detector.

## Proof review details

The Hardy bound uses strictly positive u,v and nonzero gains within segments;
zero gains split segments before taking the supremum. Real negative gains are
removed by an isometric sign gauge. The Schur sum uses right-endpoint integral
bounds in the correct direction. Lower tests certify necessity, not merely
sufficiency. The value4 is not asserted optimal for each finite weight choice.

Source refinement is not a Markov restart. Weighted fine-only modes can be
arbitrarily large even if every copied input is preserved. Changing q changes
observational precision only. Its budget must additionally preserve the native
valuation word and desired endpoint bits; native time is h+S, not log2(q).

The passive identity retains a feedthrough-containing output r. Replacing r
by a nonnegative "dissipation" score would drop its interference term. Positive
definiteness requires n>0; negative cycles are outside that positivity claim.

Forest frames retain full labelled tournaments, not only path counts, a minor,
or a clone quotient. An active path is distinct from the frozen rank order.
Repeated composition may use interior vertices as hinges, so bounded local
two-port width is not asserted. Guard evaluation counts are not native move
costs. Shared port identities persist; no resource-capacity theorem is inferred.

## Reproduction

From this directory:

```sh
python3 verify.py > receipt.json
python3 -O verify.py > receipt.optimized.json
cmp receipt.json receipt.optimized.json
LC_ALL=C shasum -a 256 -c SHA256SUMS
```

Source provenance: inherited local claims were read at `31c118ba7e12`;
the incoming Session72 synthesis at `31002464ae` was read before integration.
The source files and exact slugs are linked in README; their claims, not their
chronological summaries, supply the dependencies. No literature novelty is
claimed for Hardy inequalities, isometric refinement, or contraction dilation.

## Remaining decisive obligations

1. Insert actual fixed-source signed ROOT currents into (1), (6), or (7),
   retaining their original masses and proving the complete forcing budget.
   A positive cell lower bound must be independent of the desired conclusion.
2. Bound actual fine-mode defects, not just copied coarse responses.
3. Find a uniformly economical legal reframe after changing-path complement,
   or prove a necessary growth law for its minimal exception budget. The local
   triangle witness only tests a frame; T2 rules out the entire forest class,
   while permitting a compact carrier with a nonzero cycle budget.

These are not discharged by the finite checks or by choosing a larger modulus.
