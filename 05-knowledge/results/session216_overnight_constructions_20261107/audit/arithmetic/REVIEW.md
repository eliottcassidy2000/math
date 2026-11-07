# Independent audit: eleven common-coefficient contexts

Status: APPROVED. The packet proves the stated arbitrary-nonempty-module identity and the eleven-value dyadic interval family. No eventual-spectrum theorem, unrestricted coverage extension, single-edge staircase, or complete classification of cores is approved.

## Mathematical proof

For a fixed marked core, deleting all module vertices from a Hamilton word leaves an order of the outside vertices and a subset of occupied gaps. Every occupied gap carries a nonempty directed module path. Assign those paths in gap order. This gives a bijection with the product of a permitted skeleton and an ordered spanning path cover of the module. Thus each coefficient cj counts actual skeletons and is nonnegative. No division by j! occurs.

For nonempty U, F_U(1)=H(U) and F_U(2)=Q_U(2)-2H(U): the two omitted contributions are precisely the bipartitions with one empty bin. Every listed core has c1=a and c2=36 with higher coefficients zero. Therefore its substitution count is

    a H(U)+36 F_U(2)=36 Q_U(2)+(a-72) H(U).

Negative coefficients d=a-72 in this alternate basis do not indicate cancellation of physical paths; the path-cover basis has positive coefficients. This proves the arbitrary-input identity for all eleven contexts.

For a transitive k-vertex module, k>=1, H=1 and Q(2)=2^k. The resulting eleven counts are exactly36*2^k+d for d=-7,-5,...,13, all consecutive odd numbers within that window. Each substituted graph has order n-1+k<=k+7.

Every core has a checked spanning directed cycle, hence is strong. Replacing its marked vertex by nonempty U preserves strongness: quotient paths connect outside vertices and the module; any module vertex can leave via the quotient cycle through the mark and return to any other module vertex. This does not require U to be strong or to contain a Hamilton cycle. Nonemptiness is essential for this argument.

The empty-module extension genuinely fails in general: at d=3, deletion leaves count33 whereas the displayed nonempty formula, if misapplied with H(empty)=Q_empty(2)=1, gives39. Some contexts coincide numerically at the empty boundary, which does not license uniform extension.

## Independent coefficient method

The producer sums gap polynomials over all external permutations. independent.py instead contracts a possible module visit into each transition between external vertices, creating polynomial transition weights

    direct-edge indicator + z*(edge-to-hole indicator)*(edge-from-hole indicator).

Optional initial and final visits contribute analogous degree-one factors. Subset DP over external vertices computes the entire polynomial without enumerating permutation/gap pairs. All eleven complete coefficient vectors agree with VERIFY.json, including c0 and every zero higher coefficient. Hence this independently verifies the finite certificates establishing the all-order identity.

All adjacency rows, vertex ranges, self-loop exclusions, tournament orientations, root indices, and displayed spanning cycles were checked. The offsets are exactly the eleven specified odd values.

## Literal graph verification

The independent generator places module vertices FIRST, unlike the producer's outside-first labeling. Seed H and Q(2) values are obtained by literal permutation enumeration on induced subsets, not the producer's dynamic program.

It builds all891 specified composites: eleven contexts times all75 labeled nonempty tournament seeds of orders1..4 and six transitive seeds of orders1..6. literal.cpp independently counts their Hamilton paths using a pull recurrence from predecessor-ending paths. All891 targets pass in optimized and debug builds. All891 actual graphs also pass bitset strong-reachability checks. Maximum order is13; all counts are bounded by13!, safely below unsigned64-bit range.

Normal and optimized independent coefficient/generator receipts agree byte-for-byte. Producer verify.py was replayed normally and with -O; both reproduce the published receipts byte-for-byte. The proof comes from the complete fixed-core coefficient calculation and the bijection, not from finite module sampling.

## Provenance and scope

THM-3121-path-cover-walk-content-substitution-kernel.md, the two S203 module-rebuilding proofs, and S215 transfer/cover methods are correctly treated as inherited techniques. This packet contributes a literal bank of eleven common-coefficient contexts, without claiming scalar-count novelty.

Every supplied core row list occurs exactly in the inherited tournament bank. PROVENANCE_CHECK.json records the canonical current path, source hashes, and row locations. The bank's completeness is irrelevant to existence: all required graphs and all certificates are supplied. The currently available bank is nested under session203/induced25-51-rebuilding/session194-weak-parent-certificate; the proof's session194 reference is historical provenance.

The dyadic windows leave gaps. Choosing a different context changes d; it does not yield a common-graph edge flip or a recurrence for arbitrary represented inputs. Substituting arbitrary U produces an arithmetic progression with step2H(U), not necessarily consecutive odd values. Only the transitive special case sets H(U)=1. No uniform all-target selector follows.

## Frozen artifacts and reproduction

frozen/ contains exactly the five reviewed producer files. FROZEN_SHA256SUMS pins them. Independent artifacts include INDEPENDENT.json, GRAPHS.txt, LITERAL.json, LITERAL_DEBUG.json, and PROVENANCE_CHECK.json.

    python independent.py > independent-replay.json
    cmp independent-replay.json INDEPENDENT.json
    g++ -O3 -std=c++17 literal.cpp -o /tmp/session216-arithmetic-audit
    /tmp/session216-arithmetic-audit < GRAPHS.txt
    cd frozen && sha256sum -c ../FROZEN_SHA256SUMS

No source files were edited and no publishing was performed by this auditor. Concurrent parent changes to session216-stage were observed and left untouched. Parent remains the sole publisher.
