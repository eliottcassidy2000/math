# Internal proof review, hostile controls, and reproduction

2026-10-26. This packet has **internal adversarial review with distinct
computational formulations**, not an independent-agent audit. No numbered
canon theorem is reserved or promoted.

## Proof and status gates

1. **Native coordinate, PROVED.** Starting from the literal definition, rather
   than an inherited displayed identity, gives 2^t(3y+1)=3(3s+1)+2^a. It is
   checked on 37,501 primitive seeds through 100001. The old actual-clock note
   had an additive error3, witnessed at 3→1; its support inequality and inverse
   compiler survive. The error is repaired at its source and in MISTAKES.
2. **Paid forcing, PROVED scoped.** Eisenstein polarization proves the two-port
   identity. The exact costs lambda and 1/lambda distinguish the retained sign
   gauge. The global cost inequality is proportional to the **actual integer
   clock**, not to a compressed odd-step count. The constant 1/9139 is sharp
   locally; no sharp aggregate constant for all prescribed source families is
   claimed. The forcing payment is not a bound on the terminal synthesis.
3. **Dyadic pair ledger, PROVED.** Each pair belongs to exactly its first
   differing bit, in the odd-source index coordinate. Terminal labels and
   the integer sign remain fixed. Positive part is taken after the complete
   level sum. Uniform control is an OPEN sufficient target, not an equivalence
   to the observed final energy and not proof of universal entry.
4. **Module formula, PROVED all sizes.** The lower bound applies to every
   order, including interleaved clones and reversed internal orders. The
   equality case in triangle counting is essential for the extra 1 in the
   cycle-rank bound. Sorting clone labels is a graph relabelling followed by
   deletion of internal backward edges, so finite interleaving tests do not
   exclude a potentially better order. Separate edge/cycle minima are retained.
5. **Rerouted grammar, PROVED one step / OPEN arbitrary words.** A cyclic
   module baseline is an exact alternative representation. Its one-complement
   exception update is algebraic. Transitive cloning is not proved to be in
   the old native four-map alphabet. Therefore the quadratic family does not
   refute the narrower native-only linear-budget target.
6. **Norm-seven twist, PROVED.** Multiplication is injective of index 7 over Z,
   with the image congruence and the changed ROOT readout explicitly retained.
   The normalized inequalities are unchanged. The source-shell and inverse-
   lift connections are exact; the aliquot/tournament gap comparison has no
   predicate-preserving number map. Only the complete small fibres are tested.

## Finite universes and independent paths

All computations use Python integers and Fractions. Decimal forcing summaries
are not inputs: bounds are certified by flooring/ceiling each edge price after
multiplication by 10^12. Aggregation uses the original nonnegative source norm,
so the final enclosing intervals are valid. No checks rely on assert.

- Native single edges: all positive odd s≤100001 with s not 5 mod 8,37,501
  seeds; full affine coordinate, two-port identities, and exterior clock bound.
- Frontier: exact byte copy of the inherited627-source file, maximum63665;
  every native image lies in the file, every nonroot supplied rank decreases,
  rank(1)=0, and every primitive seed below1009 is present.
- Word carries:49 finite prefixes from seven declared starting seeds (skipping
  a nonprimitive seed if any), each of length at most16, stopping on V. Exact
  affine carry, product multiplier and rational forcing bound are checked.
- Source reconstruction: N=4096,5997,11002,16384,65536, every primitive source
  s<4N. The original three real rows and direct two-coordinate B-chain sums
  are constructed separately and compared. Zero charges remain named sources.
- Routes retain terminal, full clock, native depth and exterior entry port;
  every table hashes these source records. Repeated states or a10,000-edge
  unresolved route raise an error rather than assigning ROOT weight zero.
  No such cap is reached in the declared finite universe.
- Source refinement: all dyadic odd-index moduli through source separation.
  Direct bin energies and an independent pair-polarization evaluation agree
  at the first six cuts. Their telescoping identity retains the source norm.
- Prefix scopes: all distinct arrival-depth cuts and, separately, all distinct
  integer-clock cuts at N=4096,5997,11002,16384. Between arrival events the
  result is unchanged. No every-N prefix or joint rectangular-cut claim.
- Complete scale scan: every integer 4096≤N≤65536,61,441 scales. A separate
  three-threshold-event engine starts from N=0 and updates original charges;
  it agrees with the full-chain construction at all five shared scales.
  All terminal residues and the finest source energy are checked at every
  scale. The sufficient threshold1/2 was chosen before this complete scan.
- Cyclic modules: every one of 59,083 colour interleavings across all 20 sorted
  triples1≤a≤b≤c≤4; all 720 unrestricted vertex orders in the 2+2+2 case;
  all 119 candidate nontrivial proper Paley7 modules. Edge and cycle minima
  are computed independently as distinct objectives.
- Norm-seven twist:225 integer vectors in[-7,7]^2; each tested through six
  repeated twists, with norm, image congruence and covariant readout.304
  native inverse lifts retain divisibility guards and the exact extra 6 clocks.
- Aliquot controls: additive proper-divisor sieve through 400, the complete
  inverse cutoff for targets2,5,7,21. The table is finite and does not claim
  exclusion of other prime untouchables.

## Hostiles and strongest survivors

- **Native/ordinary-map confusion:** the previous ordinary relative gain
  3n/(3n+1) cannot be assigned to the compressed seed map. Repair with X=3s+1
  and the literal positive carry2^a. Old braid-clock and ordinary-step clocks
  remain distinct.
- **Cheap forcing implies cheap terminal readout:** false as an abstract
  implication. Two unit inputs summed at one port have energy4 from source
  energy2. The actual source-specific coherent pairs need (8), and the existing
  N=4096 all-frequency hostile already excludes a free unitary readout.
- **Refining sources deletes coherence:** false. At full source resolution,
  diagonal energy is S_arrived; restoring the requested terminal sum adds
  exactly the telescoping dyadic pair debt. Changing modulus identifies owners
  but pays no positive term by itself.
- **Minimum edges minimizes cycles:** false at C3[T3,T3,T4]. Minimum edges9
  forces cycle rank 4; a10-edge frame attains the true cycle minimum2.
- **Paley's linear copies suggest every clone-closed rank budget is linear:**
  false. Balanced cyclic clones force(t−1)(t−2) cycles. The strongest survivor
  is the smaller native-only target and a cyclic-module baseline with at most
  n−1 exceptions after one complement.
- **Norm shells transfer spectrum gaps:** false. The twist changes the
  marked readout unless transformed covariantly.7 and 21 have nonempty aliquot
  fibres; tc exclusions remain bounded-order statements. H=63 also refutes a
  blanket ban on multiples of the missing H-value 7.

## Reproduction and provenance

From this packet directory:

```sh
python3 verify.py > receipt.json
python3 -O verify.py > receipt.optimized.json
cmp receipt.json receipt.optimized.json
LC_ALL=C shasum -a 256 -c SHA256SUMS
```

Normal and optimized outputs must agree byte for byte. `frontier.txt` is an
unaltered source copy, so the checker does not depend on sparse-checkout
hydration elsewhere. `SOURCE_PINS.json` records inherited commits/blobs,
including the pre-correction clock note for lineage and incoming material
read through 64a7985e70. The corrected note is committed with this packet.

Research moves used: recover the precise source object before choosing a
norm; retain the selected side/readout under a basis change; account for
every refinement pair once; replace a poor baseline before enlarging its
exception budget. These are applications of existing META-PATTERNS, not a
new universal meta-theorem.

## Remaining obligations

Uniform paid forcing at unbounded clocks, the all-N/all-depth dyadic debt
bound, universal source entry, and a compact all-word reparsing rule remain
OPEN. Finite test success, a larger modulus, or a coefficient norm twist
does not discharge any of them.
