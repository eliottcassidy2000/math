# Session 133 independent adversarial audit

**Verdict: PASS in the scopes below.** These audits do not establish that every
odd Hamiltonian-path count other than 7 and 21 occurs. That completeness
question remains **OPEN**. The all-order exclusions of 7 and 21 are inherited,
not a new result or novelty claim. In particular cite the full path
`01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`, since the theorem
number also occurs on an unrelated file.

**Closest pointed-fibre predecessor recovered before final publication:**
`feedback_word_hamilton_spectrum_20261020/BINARY_PRODUCT_DENSITY.md` already
proves `W(X)=X ln(2)+O(sqrt X)` and, in Section 6, uses the same
complement-reversal involution with `J(X)=O(sqrt X)`. Those mechanisms and
the density law are inherited. The audited refinement is the exact G(q),
fixed-H injectivity/parity, and
`J(X)=k-1+1[G(k)<=X]`, where `k=floor(sqrt X)` and `X>=1`.
That last formula follows immediately from the proved disjoint square shells.
The final provenance update changed no executable or computational receipt.

This directory was written by a separate auditor. Producer directories were
read only. Checks use explicit exceptions and are repeated under `python -O`.
Run from this directory, with the five sibling lanes present:

```sh
python3 replay.py
sha256sum -c SHA256SUMS
```

`REPLAY.json` records receipt hashes. `INPUT_PINS.json` identifies the exact
reviewed inputs; `SHA256SUMS` covers this audit directory. Every normal/optimized
receipt pair is byte-identical. No network, third-party package, absolute
workspace path, or repository checkout is required for these scripts.

## Executed independent checks

| Area | Independent evidence | Deliberate boundary |
|---|---|---|
| Refined compiler | 120 natural-semantics AST cases against baseline/refined IR and full integer first-eligible fraction scans; six self-loop guard probes | Finite stress tests supplement the primitive/optimization invariants, not a proof by sampling |
| Universal H evaluator | All 1,100 labelled tournaments of orders 0–5 against independent subset DP; all 12 of orders 0–3 through actual emitted fractions; five whole-integer runs, including invalid input | Large native tests are not represented as completed large FRACTRAN runs |
| Seven orbit counters | Direct full-permutation Burnside, signed unordered-pair orbits and complement-twisted fixed sets, orders 0–7; 28 IR, 21 emitted-fraction and seven literal-integer cases | Pair layers and vertex colours remain distinct objects |
| One-apex 32-row program | Every 511 binary words of lengths 0–8, full integer ordered scan, independent Hamiltonian DP and explicit odd-cycle enumeration | One-apex is a restricted family; output is not a general spectrum characterization |
| Core 14/20-row programs | Whole integers through order 21; exact terminal shape and exact step count; independent component-product DP through 150 | Tournament interpretation imports the audited strong-floor theorem |
| Linear streams | 32 small matrices, full integer scans; 397-row recurrence matrix and inventory; 67/43-row tails, exact seeds and six actual integer rounds each | All-order tournament meaning uses the complete bounded strong inventory |
| Compact tail outputs | All aliases through order 1,000; unused prime 197 absent from seed and every fraction numerator/denominator | Equal counting coordinates do not identify isomorphism classes |
| 56 forbidden obstructions | Independent disjoint orbit cover of all 32,768 order-six tournaments; 56×64 order-seven augmentation; exact same 29+13 strong classes; all 56 witnesses and deletion lists | Completeness above the finite atlas uses the proved Moon/Busch cutoff |
| Pointed fibres and G | Cycle-based f checks through 4,096; all odd fibres through 2,047; square-shell bounds through q=4,096; independent word-state reconstruction of the full delivered bounded-fibre census | The apex mark is essential; this is not the unpointed tournament census |
| Extremizer/G counters | Extremizer source and IR through N=30 and nine literal runs; G source through q=4,096, IR through 31, actual emitted rows through 10, literal integers through 5 | Bounds on tested executions are reported separately from totality proofs |

The seven replay programs are `audit_compiler.py`, `audit_path_counter.py`,
`audit_linear_stream.py`, `audit_direct_programs.py`, `audit_forbidden_basis.py`,
`audit_seven_counters.py`, and `audit_pointed_fibres.py`.

## Structural proof audit

* **First-eligible semantics.** Every applicable compiled row retains its
  unique control-prime denominator. A positive DECJZ row precedes its zero
  fallback. A redirected self-target needs the explicit bridge, or reduction
  would cancel the control guard. Arithmetic shortcuts preserve named source
  states; only jump threading is full-register trace equivalence up to silent
  steps. Division tests positivity, and reachable exact divisions have their
  separate combinatorial integrality proofs.
* **Linear update.** The first pass consumes all old x registers into y=A x
  before any x is repopulated; the second pass moves y to x and clears scratch.
  Hence there is no accidental in-place update. A round takes exactly
  `1+2d+2 sum(x)+2 sum(Ax)` fractions. For every finite natural input it returns
  to the checkpoint, and repeated checkpoints implement the infinite stream.
* **Direct programs.** The one-apex read phase consumes the input exponent in
  pairs, preserving multiplicities such as the denominator factor 5². The
  final sentinel 1 is not processed as a word bit. Core loops transfer every
  output unit into three scratch units, then back; finite loop counters prove
  termination and the stated exact runtime. Repeated prime powers must never
  be collapsed to Boolean presence.
* **Exact core bounds.** Given the cited Busch minimum, replacing each SCC by
  C3/S4 blocks weakly increases size and weakly decreases H. This proves the
  exact budget K(M), not an exact-H-fibre maximum. In particular K(21)=7 is
  witnessed at H=15, while H=21 is impossible. The eventual g(N) floor and
  equality classification pass: n≥6 components incur a strict cost decrease;
  n=5 components strictly increase replacement size, ruling them out at the
  monotone block floor. Orders 0,1,2,5 retain their explicit exceptions.
* **Obstruction bounds.** For B≥1 a minimal bad tournament has no singleton
  SCC. Unless it is a triangle chain, Moon deletion leaves a singleton-free
  core with H≤B, giving order≤K(B)+1. A chain of r triangles is minimal exactly
  when 3^(r−1)≤B<3^r and proves the sharp universal maximum 3r. The complete
  forbidden-basis plateau criterion compares sets, not their cardinalities.
* **Threshold recognition.** For a supplied tournament, the K(B) rejection
  leaves total SCC size O(log B); Busch bounds each `2^s` by
  `2 B^(3 log_5 2)`. Positive capped sums/products retain exactly
  `min(H,B+1)`. Summing subset-DP costs yields the stated
  `O(B^(3 log_5 2) log³(B+2))` bit bound, separate from adjacency/index work;
  memory is `O(B^(3 log_5 2) log²(B+2))`. This does not solve existence at H=h.
* **Pointed square shells.** The corrected fibre domain is h≥3 odd, or the
  inverse range must explicitly exclude n=0. For q≥1,
  `q²<G(q)≤(q+1)²<G(q+1)` proves strict increase; upper equality is exactly
  q=2^a with a≥1. Converse fixes at most one pointed class in an H fibre, so
  fibre parity detects that fixed class. This does not survive forgetting x.

## Corrections and adversarial boundaries

1. The initial inverse-fibre wording omitted its h=1 endpoint: n=0 would
   falsely count a nontrivial strong object. The auditor identified this;
   `fibres.inverse` now returns empty for h<3 and explicitly uses n≥1.
2. **Marked/unmarked hostile:** at H=5 the pointed words `100` and `110` are
   interchanged by converse, so R_5=2 and no pointed class is fixed. Their
   underlying unpointed tournament is the unique strong order-four class,
   which is self-converse. Erasing the apex invalidates the parity claim.
   The optional inherited root-forgetting repair also passes a direct check:
   in a reduced word, deleting an interior chain vertex retains the endpoint
   triangle; deleting an endpoint removes every triangle exactly for
   `100...0` or `11...10`. These are the same double-apex family, with
   H=`2^(n-2)+1`. Apart from C3, its two apices cannot be exchanged by an
   automorphism because their arc would reverse. This justifies the precise
   unpointed correction, rather than transferring pointed parity unchanged.
3. The spectrum producer self-review corrected checkpoint 5811bb4's false
   H=15 packing-profile hostile: `[1,7]` and `[1,7,0]` represent the same
   polynomial. The repaired H=9 control compares `1+4x` and `1+2x+x²` after
   removing trailing zero coefficients. The core count, strong inventory,
   forbidden basis and structural conclusions are unchanged. See
   `../spectrum/CORRECTION.md` for lineage.
4. Source/sink padding preserves H but destroys induced minimality; arbitrary
   transitive-module contraction need not preserve H. These are different
   operations from the canonical singleton-SCC deletion kernel.

All-order claims here rely on their stated proofs and cited inputs. Exact
finite receipts certify the stated execution universes only. There is no
remaining mathematical blocker within the audited scope.
