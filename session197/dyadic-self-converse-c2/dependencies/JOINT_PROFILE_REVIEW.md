# Independent audit: dyadic physical profile universality

**PASS, after the local L>=2 scope clarification already incorporated in the frozen proof.** This is an all-order constructive congruence theorem, not an exact-target coverage theorem. The proof and numerical evidence are separated below. No source files were changed by this reviewer.

## 1. Bound snapshot

The reviewed packet is `session197-prime-atoms/dyadic-profile-universality/`:

- `PROOF.md`: `49bfc941828bf680dcaca1d5e4ec448f3b689d0f72c71f164bcb8b6bf0088c64`.
- `verify.py`: `e797180e25f059a7924de13683f593605509febdcf1f63dc4a42cc8c45bc0479`.
- Both source receipts: `8cc7db3752c35169996a74ebb1651ad2fb945e768881980ff3e20e615bdfe3ae`.

The proof's named ROOTED and modular-gadget dependencies were already reviewed. The new argument rederives the relevant loop law and explicitly maintains actual minimal binary width. The earlier proof hash is superseded: its strong-source paragraph did not explicitly exclude L=1. The corrected version does so. The construction in the joint theorem always had L>=2, so this was a boundary statement clarification, not a change to the construction or computed outputs.

## 2. All-order proof review

The formula Q=1+revc(bin(L-1)) follows from a literal first-selected-base-vertex decomposition of the reserved root-starting path. It uses the actual minimal width; padding the integer with a leading zero would change Q. The prescribed leading k-bit prefix starts with one precisely because the requested Q residue is odd. Appending later bits preserves that prefix, so it preserves the low k bits of Q independently of the low bits of L.

Starting at n=0 modulo 2^k, appending 1^a0^k changes f by 2-2^a modulo 2^k. The first zero has an odd predecessor and supplies the extra one; the later predecessors are even. This calculation remains valid from any integer with the reset residue, not only the literal integer zero. Each bit of the required half-correction is supplied by the stated one or two loops. The k=1 reduced modulus is one and requires no correction loop.

The terminal k-bit suffix changes f by f(r)-1 modulo 2^k because the recurrence uses only the previous residue and its parity. This is valid because the modulus is even. Initial prefix, reset, correcting loops and suffix therefore prescribe L, g and Q independently in the stated parity classes, with B divisible by 2^k. The stated quadratic length bound is conservative and correct. Inserting 10^k after the reset preserves all three residues without recomputing any correction, and makes the actual integer strictly larger.

For the strong corollary, the exceptional vertex a has sole incoming arc b->a. A Hamilton cycle must traverse b,a, then the transitive base in order, and return to b through the last base vertex. This gives strongness and exactly one cycle when the base is nonempty. Deleting b is transitive. A path starting at b must immediately use a; a path ending at b traverses the transitive deletion before b. Thus both endpoint counts and the deletion count are exactly one. The inherited one-seam readout yields the displayed H and R. Taking L=1 modulo 2^k makes H congruent to g, while B=0 makes R congruent to Q. The witnesses have unbounded H, since H>=L.

Appending sinks preserves H and the root-start count, forces root-ending count zero, multiplies R by two per sink, and leaves deletion count one. The valuation statement correctly restricts k>t and does not claim independently arbitrary four-coordinate control.

## 3. Independent finite controls

`check.py` imports only the submitted constructor as the system under test. It does not import its recurrence, graph builder or path counter. Independent evaluations use a direct sum over final binary zero positions. A separate physical evaluator assigns every transitive-base vertex to the portion before or after the apex. Physical graph profiles are counted by enumerating every simple directed word on every vertex subset, rather than either source subset-DP direction. The exceptional vertex is labelled last rather than first in the independent graph construction.

Both normal and optimized runs independently check:

- 4,182 joint residue triples, including exact integers through modulus bit length 128.
- 8,340 neutral-loop insertions into already constructed words, with the original correction left untouched.
- 680 loop controls from arbitrary reset integers and 16,368 arbitrary-reset suffix controls.
- All 4,096 literal two-slot Hamilton counts for n from 0 through 4,095.
- 512 literal gadget and strong-source profiles for every L from 2 through 257, including both cover responses, endpoints, deletion count and cycle count.
- 128 literal retained-root sink profiles.
- The explicitly excluded L=1 transitive-edge profile, whose ending count and cycle count are zero.

These controls support the formulas and implementation. They are not the proof of the all-order congruence conclusion.

## 4. Reproduction and audit hashes

Run from this directory:

```
python3 source/verify.py > source-normal.json
python3 -O source/verify.py > source-optimized.json
python3 check.py > independent-normal.json
python3 -O check.py > independent-optimized.json
cmp source-normal.json source-optimized.json
cmp source-normal.json source/receipt-normal.json
cmp independent-normal.json independent-optimized.json
```

The scripts do not modify source inputs. Python may create its ordinary import cache. The `source/` directory contains exact copies of the four named source files.

- `check.py`: `9f57fab5547228430f6cca4367430132391b547f2125932aec18e8848d0e643c`.
- Both independent receipts: `3a6167a3447183a6476f83f74fbd30d34c474d3d942181857ecb8389f03294f2`.
- Both source replay receipts: `8cc7db3752c35169996a74ebb1651ad2fb945e768881980ff3e20e615bdfe3ae`.

The exact-target spectrum, primality, arbitrary odd moduli and arbitrary four-coordinate realization are outside this theorem and were not inferred from these controls.
