# Independent audit: the original character bank obstructs both selected rectangles

23 October 2026. **VERIFIED:** the exact charge-order certificates, signed
rectangle algebra, disjoint-target obstruction and distinct-port bound below.
**FINITE-EXACT:** the stated original order-eight data. **OPEN:** legal moves
outside these fixed rectangles, arbitrary multi-source blocks, and realization
of any altered aggregate network by full native wiring.

The audit imports no producer code and runs no producer linear solver.
`audit.py` reads the authenticated original edge list, builds its weighted
adjacency and diagonal, and independently evaluates `D z - A z`. It also
reconstructs the two canonical source tournaments and enumerates all `8!`
permutations at each source. The four original byte pins are retained in
`SOURCE_PINS.json`. The producer's full image-witness JSON is pinned separately
by its hash in `receipt.json`.

## 1. Exact charge orders, not a guessed complete bank

Write `C8 = D^0/(L Z^V + 8 D^0)` for the original connected graph. In the
standard identification with zero-sum vectors modulo eight, this is the
zero-sum module modulo `im L`. The inherited, already audited annihilator
identity makes character evaluation equivalent to image membership. Here
the independently checked certificates give more directly:

| Difference charge | Checked character | Character difference mod 8 | Image witness |
|---|---|---:|---|
| `e1238-e2246` | `2Y0` | 2 | `L z = 4(e1238-e2246) mod 8` |
| `e881-e2762` | `2Y0` | 2 | `L z = 4(e881-e2762) mod 8` |
| `e2650-e2843` | `2Y1` | 2 | `L z = 4(e2650-e2843) mod 8` |

`Yi` denotes zero-based entry `i` of the actual original electrical receipt's
`radical_lifts_mod4`. Each vector is checked against the full integer
Laplacian modulo four, and `2Yi` is checked modulo eight. The supplied
3528-entry image potentials have root coordinate 3527 equal to zero. The
audit reproduces their entire integer correction `(Lz-4d)/8`, not just a
solver success flag or a parity reduction. All 438,216 original pair rows
are consumed; loops have zero Laplacian contribution and all other weights
are exactly `blue + black`.

For each charge `d`, the image witness proves `4[d]=0`. The displayed
harmonic character evaluates `2d` to four modulo eight, so `2[d]` is not
zero. Consequently each charge has **exact order four** in `C8`. This proof
does not need to assume the selected bank generates the whole harmonic
module. The original single frozen primitive-lift character `X` has
difference zero modulo eight on all three pairs, so its blindness is real.

## 2. Signed convention and the general source-gap obstruction

Let `S != R`, `u=eS-eR`, and let `v` be an integer zero-sum target vector
with `vS=vR=0`. For a nonnegative integer `b`, define

    Delta L = 2^b (u v^t + v u^t).

This is a symmetric zero-row-sum matrix. For `v=a(eU-eV)` on four distinct
vertices, it means removing edge weight `2^b a` from `(S,U),(R,V)` and
adding it to `(S,V),(R,U)`. Negative amplitudes reverse the operation.
These are signed algebraic changes; their positivity and native path
realization require separate checks.

Suppose an original harmonic character `x mod 2^k` has source difference
`u.x` of valuation `s<k`. If this character remains harmonic after the
perturbation, then at each off-source coordinate `j`,

    0 = (Delta L x)_j = 2^b v_j (u.x) modulo 2^k.

Therefore every `v_j` is divisible by `2^max(0,k-b-s)`. If the source gap
vanishes, use `s=k`, giving no restriction. This is exact scalar
divisibility, not cancellation of a nonunit in a dyadic ring. The audit
exhausts 18,582 scalar cases through modulus 32, including zero gaps,
negative coefficients and factors already zero modulo the modulus.

Disjoint source/target support is essential to this formulation: with
`u=v=(1,-1)`, `x=(0,4)` and coefficient one modulo eight, the two tensor
terms cancel to zero, although `v` is odd and the source gap has valuation
two. The audit retains this hostile. More generally, source rows still
require `2^b(v.x)=0` for every original harmonic character; the target
coordinate divisibility alone is not sufficient.

The supplement's exact general criterion also passes: writing
`g=gcd_j |v_j|` with `g=0` for the zero vector, all original characters
remain harmonic if and only if both `2^b v` and `2^b g u` lie in the
original image. The source rows impose the first condition. The target
rows impose membership of every `2^b v_j u`; these are equivalent to the
gcd condition by integer Bezout and closure of the image under integer
scalars. The audit independently exhausts 200 cases on full finite
image/kernel sets of two connected five-vertex weighted graphs modulo
four, including zero and nonprimitive target vectors.

An actual-original-bank hostile makes that latter boundary concrete. Keep
`u=e1238-e2246`, set `v=2(e0-e5)`, and use the original character
`Z=eighth_radical_lifts_mod8[0]`. Although every target coefficient is even,
`Z0-Z5=7 mod 8` and the black-factor-two tensor has nonzero current four
at each of sources 1238 and 2246. This is a counterexample to a proposed
general sufficiency claim, not an asserted legal port switch.

## 3. Why the two named rectangles have an exact parity criterion

Take the shared source charge `u=e1238-e2246`, and

    v1=e881-e2762,       v2=e2650-e2843.

The original paired-black convention gives

    Delta L = 2(u (a v1+b v2)^t + (a v1+b v2) u^t).

The four target vertices are distinct from one another and from both
sources. A target row is thus individually `+/-2a u` or `+/-2b u`. Since
`[u]` has order four, preserving all original harmonic characters requires
both `a` and `b` even. Equivalently, the single authenticated character
`2Y0` already forces this at the disjoint target coordinates because its
source gap equals two modulo eight. No cancellation between these target
rows is available.

Conversely, if `a,b` are even, every target row is a multiple of `4u`
and each source row is a sum of multiples of `4v1,4v2`. All three are
explicitly certified image vectors. Hence every row annihilates every
original harmonic character, proving sufficiency. Thus, for this signed
two-rectangle family, retaining all original modulo-eight guards is
**equivalent to both amplitudes being even**. The audit checks the finite
window `-5 <= a,b <= 5` against the authenticated characters as a control;
the proof and image certificates give the unrestricted integer statement.

“Retaining the bank” here means `ker(L mod 8)` is contained in
`ker((L+Delta L) mod 8)`, with the same grounding. It does not claim that
the new harmonic kernel is equal to the old one, or that the new critical
quotient or native wiring is unchanged.

## 4. Native paths, shared pins, converse and the amplitude cap

Independent complete permutation enumeration recovers exactly 409 paths
at source 1238 and 523 at source 2246. Every path's original source and
target owners are authenticated using the original forward-word owner
map. Every reversed path is checked directly in the converse tournament,
including both merged owners. Complete automorphism enumeration gives the
identity only at each source, so no omitted nontrivial automorphism/kernel
lift supplies an extra copy of a scarce source path.

All eight selected paths are checked, including each rectangle's opposite
retargeting orientation, fixed endpoint owners, indexed internal arc,
actual frozen phase and clock costs `2/3` or `3/3`. Direct original-edge
lookup reproduces pure-black weights `4,2,4,4` for the first rectangle and
`4,2,4,6` for the second. Adding the reversed converse partners has the
declared black weight two, not weight one and not another factor of two.

The exhaustive source-1238 census has exactly one path to target 2762:

    (7,1,0,5,2,6,4,3),

and exactly one path to target 2843:

    (7,0,6,4,1,5,3,2).

This is a global count at that source, before additional pins. Either
orientation of amplitude two would need two distinct paths on that scarce
side, in the old or new selection. Hence the declared distinct-observer
model has `|a|<=1` and `|b|<=1`. The literal amplitude-one selections give
the matching lower feasibility bound, but fail the full original bank.
Combining the exact parity restriction with these capacity bounds leaves
only `a=b=0`.

The conclusion closes these two fixed selected-port rectangles and all
their integer superpositions under the stated distinctness rule. It does
not rule out repeated-observer models, other target pairs, additional
sources, overlapping supports, target-multiset-preserving labelled-slot
trades, or different legal move menus. An altered aggregate adjacency is
still not a full native ownership bijection.

## 5. Reproduction

From the workspace root, with both packets present:

    python3 -B session44-audit-rectangles/audit.py --repo session35-publication/clean-repo --producer session44-rectangle-bank
    python3 -O -B session44-audit-rectangles/audit.py --repo session35-publication/clean-repo --producer session44-rectangle-bank

For a published layout, pass the repository and producer packet paths
explicitly. Both modes use explicit exceptions, not validation-critical
assertions. The independently generated receipts are byte-identical.
The final checked code, input pins, output and producer snapshot hashes
are recorded in this packet's manifest.

Final results: **898,404 independent checks** in each mode, with identical
output. The producer's final checker was separately rerun in both modes;
both 8,820-check outputs exactly match its frozen receipt. Final proof
review found no mathematical correction necessary. A minor solver-label
clarification was adopted: there is one wholly unpivoted column, while
lower-modulus pivot variables retain higher-bit freedoms. No vector or
equation changed. The final image-witness SHA256 is
`5bd986838a40035fa617f567ff9e09845c407f9aa94fc0d0ae7fb7fee3505c54`.
