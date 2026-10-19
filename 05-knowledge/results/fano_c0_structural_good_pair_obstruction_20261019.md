# A c=0 critical pair need not force any good pair

Date: 2026-10-19. Scope: a single explicit family of small counterexamples to a stronger structural criterion. These are not counterexamples to balance or Gold and do not settle the seven-thin c=0 frontier. In particular, root degree bounds must not be mistaken for ambient thinness.

## Construction

Use the cyclic Fano incidence relation on indices modulo7: point i is incident with line j exactly when j−i belongs to {0,1,3}. Take minima x,L0,...,L6 and maxima R0,...,R6,y. Set

    x<Rj for every j;
    Li<y for every i;
    Li<Rj exactly at a Fano incidence.

There are no other strict relations. Every point is incident with three lines, every line with three points, and distinct points (respectively lines) have one common incident line (respectively point).

The pair (x,y) is critical: D(x)=∅⊆D(y) and U(y)=∅⊆U(x). No third vertex is incomparable with both: every Li is below y and every Rj is above x. Thus c=0 and there are seven exclusive neighbors on each side.

## No good or very good pair

Use Zaguia's definition: an incomparable ordered pair (a,b) is good if, in P or its dual, D(a)⊆D(b), U(b)\U(a) is a chain, and p(a<b)≤1/2. His very good condition is equal downsets and two chain upset differences. [Definition1 and Definition5, author-deposited paper](https://arxiv.org/pdf/1610.00809v3).

Because this poset has height two, any subset of a single level is a chain only when it has at most one element.

Distinct Fano minima have both upset differences of size2; the corresponding dual statement holds for distinct Fano maxima. For x and a Fano minimum Li, U(Li)\U(x)={y}, while U(x)\U(Li) has size4. Thus the only structurally eligible orientation on that level is (x,Li). The dual possibilities are exactly its reflected upper-level counterparts. No distinct upper vertices have nested downsets, so there are no other same-level structural candidates in P; the dual argument handles the other direction.

Every remaining structural candidate is a nonincident minimum–maximum pair. Its forward probability is strictly greater than one half by the critical-pair swap injection, with strictness witnessed by placing all minima before all maxima: swapping the indicated minimum and maximum violates another incident relation. Thus none is good.

The exact enumerator verifies all candidates, with the following probabilities:

- p(x<Li)=p(Rj<y)=216/389>1/2;
- p(Li<Rj)=3359/3501>1/2 when Li and Rj are nonincident;
- p(x<y)=3476/3501>1/2.

There are exactly72 structurally eligible orientations across P and its dual:14 of the first type,56 of the second,2 of the third. All fail the good-pair probability condition. There are no very good pairs even before applying a probability test.

The total extension count is3,557,240,064. The numerator for x<L0 is1,975,228,416. The checker computes all pair counts by prefix/suffix ideal recurrences and cross-checks that numerator by separately adding x<L0 and recounting.

## Precise degree and Gold scope

The root incomparability degrees are8. Each Fano point and line has degree11, comprising seven same-level vertices and four nonincident opposite-level vertices. Hence the poset is **eleven-thin**, not eight-thin or seven-thin.

Simultaneously rotating point and line indices is a nontrivial automorphism. Consequently Gold already follows from Peczarski's automorphism theorem. [Section3, publisher PDF](https://cdm.ucalgary.ca/article/download/62349/46785/178291). This construction refutes only the assertion that a c=0 critical pair forces existence of a good or very good pair. It neither refutes nor proves a general c=0 balance/Gold theorem. It cannot satisfy the hypothetical no-balanced-pair assumptions used in maximum-majority-span arguments.

## Requested fourteen-point truncations

Delete L0 and one line. Both representative cases retain a c=0 critical root pair and six exclusive neighbors on each side:

- Incident deletion R0:51,528,960 extensions, no good pairs in either orientation or dual, and no very good pairs.
- Nonincident deletion R2:55,064,448 extensions, no good pairs in either orientation or dual, and no very good pairs.

Both have root degree7 but ambient maximum degree10. Surviving Fano vertices have degree9 or10. Thus these do not supply seven-thin counterexamples to the structural criterion. No girth/Moore-bound or seven-versus-eight ambient threshold conclusion follows from these tests.

## Reproduction

Run `verify_fano_c0_no_good_pair_20261019.py` normally and with `python -O`; the receipts match. `verify_fano_truncations_20261019.py` checks just the two explicit induced constructions. This is not a census and does not claim minimality in vertex count.
