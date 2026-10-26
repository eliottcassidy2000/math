# Classifying composition is not yet closing the metagraph

Status: inherited PROVED flux law, derived all-order family nonclosure, and independent FINITE-EXACT endpoint checks. This note uses ordinary classes and marked Hamiltonian paths. It does not assert black-core connectivity, a closed grammar for all prime tournaments, or reconstruction of the complete line ledger from its flux distribution.

## 1. Recover the exact native operation

Pinned at fc3768ac5fcd87022fde9e6d7e3f06c0e5bc22e8, THM-785 (`THM-785-cyclic-triangle-flow-on-blue-black-metagraph.md`) and THM-646 (`THM-646-line-metagraph-score-complement-skeleton.md`) already give the law used below. Exact source copies are in `inputs/`. It is recovered, not claimed new.

For a tournament T of order n>=2 and a directed Hamiltonian word P, its literal line partner U reverses every nonpath edge and retains every consecutive path edge. Equivalently, start with T's converse and toggle P's consecutive arcs back. Write f and l for the distinct first and last vertices, and s for outdegree. Then

    s_U(v)=n-1-s_T(v)+[v=f]-[v=l],
    c3(U)-c3(T)=s_T(f)-s_T(l)-1.

Proof: each internal path vertex has one incoming and one outgoing toggle, cancelling its degree change. The endpoints retain one change each. Global converse leaves the number of cyclic triples unchanged. In c3(T)=binom(n,3)-sum_v binom(s_T(v),2), changing the two endpoint degrees of T's converse changes the latter sum by s_T(l)-s_T(f)+1. This proves the second identity. At n=1 the endpoints coincide, so that separate-endpoint derivation does not apply. The degree identity and the literal complement involution still hold.

These edges are not the auxiliary legal relocations used to identify different presentations of one ordinary tournament. Relocation stays within a class; literal complement usually does not. Converse merging is a third distinct quotient.

## 2. The classified chirality families are not line-closed

Let A be a singleton dominating C3; its internal scores are 3,1,1,1. Its converse has scores 0,2,2,2. Let R be any regular tournament of odd order m and substitute one A-type block of order four at each quotient vertex. Every vertex has exactly 2(m-1) outgoing exterior neighbours.

In the uniform all-A tournament the scores are 2m+1 and 2m-1, all odd. In the uniform all-A-converse tournament they are 2m-2 and 2m, all even. Therefore the endpoint score difference is always even, and EVERY literal complement partner changes c3 by an odd number. In fact its change lies in {-3,-1,1}.

Consequently neither uniform node has a literal line partner within its fixed-C3-count family, let alone its exact cyclic-triple family. This applies both to the carousel prime families and to the iterated-C3 composition family: their undeformed quotients are regular. A complete classification of those families, even with exact automorphism groups and every path-cover coefficient, is not a closed classification under metagraph lines.

The cyclic-triple count in every order-4m chirality decoration of a regular R is

    m + 4^3 c3(R) = m + 64*m*(m^2-1)/24 = (8m^3-5m)/3.

There are m internal block triangles. A triple with two vertices in one module and one outside is transitive. The remaining cyclic triples are exactly transversal choices over a quotient cyclic triple. This count does not determine the quotient: the prime carousel R9 and C3[C3,C3,C3] are both regular but have different modular trees and automorphism groups.

## 3. Exact port evidence for the 81/27 witness

The script independently checks the addressed-degree law, triangle flux and partner involution on every directed Hamiltonian word of every labelled tournament through order five: 7,887 words. Separately it computes the full endpoint Hamiltonian-path matrices for the two order-12 witnesses by an exact subset dynamic program. No isomorphism oracle or floating arithmetic is used in this endpoint computation.

Both have H=462915. Their flux histograms differ substantially even though their labelled C3 fields and all path-cover counts coincide:

| Delta c3 | all-A labelled paths | one reversed A labelled paths |
|---:|---:|---:|
| -4 | 0 | 4293 |
| -3 | 73710 | 44064 |
| -2 | 0 | 86535 |
| -1 | 267624 | 132084 |
| 0 | 0 | 120258 |
| 1 | 121581 | 64017 |
| 2 | 0 | 11664 |

Dividing the respective columns by 81 and 27 gives ordinary normalized tiling-port counts, because automorphisms act freely on complete vertex words and preserve each flux category. This gives respectively 5715 and 17145 total ordinary presentations. Do not call these values the degrees of converse-merged nodes without separately accounting for their two ordinary fibres.

The example path [8,11,9,10,0,3,1,2,4,7,5,6] in the all-A witness produces a partner with c3=68, rather than67. The mixed witness's zero-flux paths only establish an unchanged TRIANGLE COUNT; they do not prove its partners share the whole C3 hypergraph or belong to the decorated family. That additional membership test has not been claimed.

## 4. Classification data and deliberate losses

- A canonical decorated modular tree determines the unmarked tournament. Child size alone or a prime-factor list is not a substitute for child isomorphism type and its placement.
- A coherent choice of child models gives a semidirect-product description of its automorphisms. The section is not canonical without those identifications.
- To classify marked presentations and literal lines, add a Hamiltonian word, or an equivalent exact run/allocation/port encoding. The unmarked classification does not choose this datum.
- Endpoint scores give this one exact line-flux observer, not the target's full binary relation. Keeping the full relation and path permits literal reconstruction; reducing to the flux discards it.

This is the precise role of the simple-group analogy: classify compositional atoms and retain extension/decorating data. It is not a claim that modular-prime tournaments correspond to simple groups. The prime carousel R9 already has cyclic automorphism group C9, which is not simple; conversely a factor's symmetry can disappear when its children cease to be isomorphic, as the 81/27 witness demonstrates.

## Reproduction

    python3 check_line_boundary.py > /tmp/line-boundary.json
    python3 -O check_line_boundary.py > /tmp/line-boundary-O.json
    cmp /tmp/line-boundary.json /tmp/line-boundary-O.json
    cmp /tmp/line-boundary.json line_boundary.json

The finite enumeration is supporting evidence for the stated law, not a proof at untested orders. Source script checks use explicit exceptions and remain active under -O.
