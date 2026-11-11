# Recover the joint-residue engine, then retain the height coordinate

Session251, 2026-11-11. **PROVED-IN-NOTE local image classification and height obstruction; FINITE-EXACT new automaton controls. Chan–Pak OPEN.** No canon promotion or literature-firstness claim.

## Useful outcome first

A much stronger modular tool than a numerical residue census was already hidden in `session197/modular-gadget-universality/PROOF.md`: the binary one-apex compiler realizes **every joint residue** of its cover count and its odd Hamilton count. This closes the modular part of essentially every fixed-prefix one-apex-join approach. The unresolved issue is exact height, not finding another residue class.

The new usable refinement is a complete local-image formula for **every** rooted prefix, including nonprimitive ones. A new small executable extracts shortest binary witnesses. An elementary counting obstruction explains why upgrading the joint construction to linear height cannot possibly work: there are quadratically many joint targets. Any spectrum proof using this tool must choose the unused cover-count residue freely rather than prescribe it first.

## Inheritance and concept board

Closest proved mechanism: the Session197 joint modular loop, combined with THM-4875's rooted one-arc join. Canonical hostile: singleton-prefix image misses the exact count13 while realizing its residue modulo every integer. Corrected near miss: THM-4883's signed responses fill residue classes but do not prescribe the positive parent count. Least-used sidecar: the integer size of the binary address, rather than merely its residue or bit length.

Live concepts: (1) positive rooted joins; (2) binary finite-state transducers; (3) local congruence images; (4) exact target height; (5) signed-response clocks. The output below connects all five without identifying positive counts with signed offsets.

Read first for current context:
- `05-knowledge/results/decomposition16_odd_spectrum_20261110/README.md` (threshold95 extension criterion remains open).
- `05-knowledge/results/session250_sixteen_20261110/SYNTHESIS.md` (positive languages and named ports, not moments alone).

Recovered source addresses:
1. `session197/modular-gadget-universality/PROOF.md`, Theorems1–2 and section5. All joint residues modulo even M, O(M²)-bit general constructor, O(log² M)-bit dyadic constructor.
2. `01-canon/theorems/THM-4875-old-mark-digit-lower-bound-lorentzian-gadget-adjoint-and-one-apex-joins-realise-every-odd-number-to-2-32-except-7-and-21.md`, partD. Positive prefix formula and inherited finite coverage; its addendum reports coverage to2^35.
3. `session197/path-polynomials/ONE_CORE_BINARY_INVERSE.md`, sections2–5. Physical binary words, distinct prefix/suffix conventions, and finite target inverse. This precludes calling a generic bounded inverse new here.
4. `01-canon/theorems/THM-4883-signed-one-apex-responses-exact-2-adic-structure.md`, T1–T5. Signed functional universality, parity sidecar and unique source-independent clock. This is useful context, not a positive-count proof.
5. `01-canon/theorems/THM-4887-one-apex-energy-at-scale-two-to-the-n-over-three.md`, two-block identity. The same binary address controls multiplicity/energy. Its averaging bound does not bound the least address of a prescribed count.

No inherited source is audited again. Files absent from sparse checkout were read selectively with `git show HEAD:path`.

## 1. Conventions and physical map

Put f(0)=f(1)=1 and, for q>=1,

    f(2q)=f(q)+q+(q mod2),   f(2q+1)=f(q).

Let G(L)=f(L-1), L>=1. A one-apex rooted gadget has the inherited triple (G(L), L mod2, L). For a physical rooted prefix with triple (h,e,R), the one-arc join count is

    F(L)=h G(L)+R L-e(L mod2).                         (1)

Here h is odd, 0<=e<=h, R>=h, and R=e mod2. We need only h odd and the parity relation for the local statement. In particular the abstract theorem below does not assert that arbitrary triples satisfying these inequalities are physical.

**Recovered joint universality.** For every even M, every u mod M and every odd z mod M occur as

    (L mod M, G(L) mod M)=(u,z)

for infinitely many actual positive L. The source gives literal words and hence tournaments, not just formal matrix states.

Source -> target map: feed a joint residue witness into (1). Preserved predicate: the output is a Hamilton-path count of a literal tournament and has the specified residue. Destroyed information: exact output value and minimum order. Retained sidecar: full binary L, prefix triple, and integer F(L). Cheapest decisive test: evaluate (1) over the reachable finite state graph and compare its residue image with the formula below.

## 2. Complete local image of every prefix

**Lemma 1 (PROVED-IN-NOTE).** Fix an even M and let

    d=gcd(2h,2R,M).

Then the set of F(L) modulo M is exactly

    { x mod M : x=h mod d or x=h+R-e mod d }.          (2)

Every residue in (2) is attained for infinitely many L. This is a necessary and sufficient local criterion, not a sufficient exact-integer criterion.

**Proof.** Write u=2v+p and z=2w+1, p in{0,1}. Formula(1) becomes

    h+2hw+2Rv+p(R-e) mod M.

For each fixed p, w and v range independently modulo M/2 by the recovered joint universality. Their additive image is the subgroup generated by2h and2R, namely the multiples of d modulo M. This gives both cosets, and joint universality supplies infinitely many actual L for each admissible pair. QED.

**Odd-modulus corollary.** For odd m, put D=gcd(h,R,m). Then the image is exactly

    { x mod m : x=0 mod D or x=-e mod D }.            (3)

Proof: apply(2) at2m and project; h and R vanish modulo D and2 is invertible modulo m. The parity lift is necessary: L mod an odd modulus does not determine L mod2.

**Primitive corollary.** If gcd(h,R)=1, F realizes every parity-compatible residue modulo every positive integer, infinitely often. In(2) d=2 and both offsets are odd; in(3) D=1.

**Physical hostile.** The directed triangle rooted anywhere has (h,e,R)=(3,1,3). Its final-gadget image modulo3 is exactly{0,2}; counts1 mod3 cannot come from this fixed prefix. This is stronger than the weaker common-content test gcd(h,e,R): that gcd is1 here, yet a genuine local obstruction survives. Retaining e separately is essential.

For a finite prefix bank, take the union of(2) for its members. A single gcd(h,R)=1 member eliminates every fixed-modulus obstruction to the whole bank, without implying completeness.

## 3. A precise local/global hostile

The singleton prefix has (h,e,R)=(1,1,1), so

    F(L)=G(L)+L-(L mod2).

It has no congruence obstruction by the primitive corollary. Yet it misses13 exactly. It is enough to check1<=L<=13; for L>=14, G(L)>=1 gives F(L)>=15. The executable performs this terminating check, and lists its exact gaps through101. The count13 is realizable by other tournaments; it is a gap of this constructor, not of the global spectrum.

Thus even the statement “for each modulus there are infinitely many physical witnesses congruent to13” supplies no approach to an exact13 witness in this family. The recovered signed-clock theorem has the same missing height coordinate when used as a positive-count proposal.

## 4. New shortest-word compiler and its sharp information barrier

For even M, the state graph is

    (n,z) -> (2n+1,z),
    (n,z) -> (2n,z+n+(n mod2)),                      (4)

with coordinates modulo M, starting at(0,1). Reading MSB-first produces (n,f(n)). There are exactly M²/2 reachable states by joint universality. Breadth-first search therefore returns a shortest binary witness for every joint residue. The predecessor chain itself is the certificate; integer reevaluation checks it independently of residue transitions.

**Lemma 2 (PROVED-IN-NOTE, address-count barrier).** Any constructor that supplies every joint pair modulo even M must use some address n>=M²/2-1. Equivalently its worst-case binary length is at least ceil(log2(M²/2)), for M>=4.

**Proof.** Each nonnegative integer n supplies only one pair. The interval0<=n<M²/2-1 contains fewer than M²/2 integers. For a maximal bit length K there are only2^K distinct nonnegative integers of length at most K (including0). QED.

Consequently there cannot be a uniform O(M)-address constructor for *all joint pairs*. This does not rule out O(M) addresses for a single scalar equation F(L)=a mod M, because that equation allows many choices of the unused cover residue u.

This is the actionable source-target bridge:
- Source: complete joint modular reachability, with word certificates.
- Target: exact positive count N via(1).
- Map: choose a state on h z+R u-e(u mod2)=N mod M.
- Preserved: literal gadget realizability and the congruence.
- Lost: integer height.
- Required sidecar: keep the least address, and allow u to vary over the target line.
- Cheapest next test: measure **minimum address on each scalar target line**, not the diameter over every joint state. Requiring all joint states to have linear height is provably the wrong strengthening.

For an exact certificate one can retain integer F(L) and choose M larger than its permitted range, so a congruence forces equality. This needs an actual upper bound placing the constructed F(L) in that range. The current modular proof supplies no such bound for a requested N. No claim of exact-target progress is made merely from small modular diameters.

## 5. New finite controls and limits

Run from repository root:

    python3 -B 05-knowledge/results/session251_creative_completion_20261111/gems/local_image.py

`receipt.json` records all64 even moduli2..128, four prefixes per modulus, 256 exact local-image comparisons, and one reconstructed farthest integer word at each modulus. It also checks the primitive exact13 hostile. This is a new finite automaton experiment, not a rerun of the tournament census.

At M=2,4,8,16,32,64,128, the maximum shortest bit lengths are respectively

    1,4,8,10,12,15,18.

At M=128 all8192 joint states are reached; the information lower bound is13 bits and the measured diameter is18. These observations suggest studying logarithmic shortest-word bounds, but no all-modulus logarithmic theorem is asserted. Even such a theorem with a constant larger than1 would only give polynomial address height, insufficient by itself for an exact-target spectrum induction.

Normal and optimized replays are byte-identical. The full requested Chan–Pak construction remains open. The strongest deliverable here is the completed local classification, its useful primitive/nonprimitive distinction, and the reason the next positive construction must optimize a scalar fibre rather than the whole joint-state atlas.
