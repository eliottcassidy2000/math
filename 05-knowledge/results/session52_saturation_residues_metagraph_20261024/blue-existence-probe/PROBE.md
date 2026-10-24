# Internal blue-port existence: exact reduction and recovered literature lead

24 October 2026. **PROVED:** the reduction and low-order cases below.
**CONDITIONAL:** all-order existence, pending inspection of the original
Salvi-Zagaglia theorem. This is not a new self-contained all-order proof.

## Result of the bounded probe

The missing lemma appears to be a published theorem, not a new conjecture.
Farrugia's 1999 research survey, Section 2.30, printed page 59 (PDF page 75),
attributes the following stronger statement to Salvi-Zagaglia: in a regular
self-converse tournament of order 2p+1 with a specified involutory
antimorphism, each invariant arc belongs to an invariant directed cycle of
each available length, except possibly 4 and one of 3 or 5. For order at
least seven this includes Hamilton cycles without exception.

Bibliographic lead: N. Salvi-Zagaglia, *Sull'esistenza di particolari circuiti
di ogni lunghezza di un torneo regolare auto-opposto*, Istituto Lombardo
Accademia di Scienze e Lettere, Rendiconti A 114 (1980), 169–178;
MR 84h:05059.

Discovery reference (secondary, NOT a verified primary import):
https://alastairfarrugia.net/sc-graph/sc-graph-survey.pdf

The survey adds the exact sentence (12 words): “In particular, T must contain at least k directed Hamiltonian circuits.”
Its last sentence does not explicitly repeat symmetry. The argument below
therefore uses the main quantified clause only at N>=7, and independent
small-order controls at N=3,5, rather than extrapolating the last sentence.

The primary paper itself was not obtained. The publisher's archive at
https://istitutolombardo.it/attivita/pubblicazioni/elenco/5 lists Volume 114 A
(1980), but its row has no downloadable link. The current journal archive
https://www.ilasl.org/Scienze/issue/archive begins with 2011. Searches found
no primary full text. No account, paid retrieval, or external communication
was attempted. Do not promote the all-order statement to a primary-verified
CITED theorem until the paper is inspected or an independent proof replaces it.

## Exact reduction, proved without literature

Let T be a tournament of order 2p, with p vertices of score p and p vertices
of score p-1. Call them high and low. Fix an anti-involution q of T.
It has no fixed vertex: an anti-automorphism has at most one fixed vertex,
and the number of fixed vertices of an involution has the parity of the order.
Moreover d(q(v))=2p-1-d(v), so q exchanges high and low vertices.

Adjoin z, orient z toward all high vertices, and all low vertices toward z.
Every vertex now has outdegree p, and this is the unique regular extension
on the prescribed new vertex. Extend q by q(z)=z. The new incident arcs
are reversed under q because high and low are interchanged; therefore the
extension is still anti-involutory.

Suppose this regular extension contains a q-invariant directed Hamilton
cycle. Read it beginning at z as

    z, v1, ..., v(2p), z.

Since q fixes z and reverses the cyclic direction, uniqueness of predecessor
and successor on the cycle gives q(vi)=v(2p+1-i) at every position. Deleting z
therefore yields a reflection-fixed Hamilton path in T. Its first vertex is
high because z points to it; its last is low because it points to z. Thus
its fixed-path word is blue and internal by the balanced endpoint criterion.
Conversely, every q-symmetric high-to-low Hamilton path closes with z to
such a cycle. These operations are inverse, including the specific q.

On an odd q-symmetric cycle there is exactly one q-invariant arc: on the
rooted reading it is vp -> v(p+1)=q(vp). Hence the stronger cited lemma
would realize EVERY paired arc of q as the central arc of an internal blue
presentation. Different paired arcs give distinct labelled paths.

An arbitrary self-converse tournament always has some anti-involution by an
elementary group argument. A tournament has no nonidentity involutory
automorphism, since a swapped pair would require its arc to be preserved
with reversed endpoints. Thus its automorphism group has odd order by
Cauchy's theorem. If g is an anti-automorphism, g^2 has odd order m, and g^m
is an involutory anti-automorphism. This establishes that the specified-q
statement covers every self-converse class, rather than just a subclass.

## Consequences, conditional only on the primary literature lemma

For every p>=1 and every self-converse balanced T of order 2p, every
anti-involution q admits a high-to-low symmetric Hamilton path. More strongly,
each of its p paired arcs occurs centrally in such a path. Orders 2 and 4
are established independently below, so the exceptional small cycle lengths
in the reported literature statement cause no gap.

There are at least p labelled q-symmetric downhill paths. This bound is NOT
asserted for the quotient number of blue words of a class: automorphisms
can identify paths. The immediate quotient conclusion is positivity only.
No connectivity, bridge-freeness, minimum degree 17, or all-order count follows.
The order-ten minimum 17 remains a finite observation.

## Independent finite controls

Run `python3 check.py` and `python3 -O check.py`. Both reproduce `normal.json`
byte-for-byte. The checker enumerates all matrices anti-invariant under the
fixed reversal q, filters the balanced scores, then checks all 2^p p!
half-path choices. This is complete for pairs (T,q) up to relabeling because
any fixed-point-free involution is conjugate to reversal. It uses no theorem
about existence to filter its universe, and uses explicit exceptions rather
than assertions.

At n=2,4,6 the full anti-invariant universes contain 2,16,512 matrices;
2,8,112 are balanced. Their minimum labelled symmetric downhill-path counts
are 1,2,5, respectively. Every paired arc is realized centrally in every
case. In particular, n=2 and n=4 close the two small cases needed by the
conditional argument. This computation does not establish the general lemma.

## Inheritance and scope

The exact extension and cycle/path bridge are inherited from
`sources/balanced-metagraph-theorem.md`, Section 5 (pinned from
`session52-balanced-metagraph/THEOREM.md`), itself crediting
Session51 regular completion. Section 2 there supplies the internal HL iff.
The present note makes the anti-involution-preserving restriction explicit;
it claims no literature priority for regular completion or Hamiltonicity.
The closest canonical anti-involution reference located is
`sources/anti-involution-context.md` (original canonical file
`THM-4743-converse-extended-gluck-for-tournaments-primitive-classification-paley-tournaments-and-c3-towers.md`),
which also states anti-involution existence. The elementary proof above avoids
importing its unrelated classification claims. A targeted search found no
Salvi/Zagaglia reference in the current theorem or reference directories.

Stopped after the precise reduction, published-theorem lead, and small-case
certificate. No large enumeration or existing packet mutation was performed.
