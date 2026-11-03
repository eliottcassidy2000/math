# Scope and verification audit

**Status: PASS after same-agent proof review and independent computational
paths.** This is not an independent researcher audit. The normal and -O
receipts must be identical; failure checks do not use Python assertions.

## Mathematical obligations

1. D1 quantifies over EVERY integer k. The proof is the six-candidate
   restriction D6--D8, with all small cases listed. A large parameter
   scan is a hostile check, not the reason the theorem is uniform.
2. D2 classifies rational PREPERIODIC points with fixed parameter
   denominator sixteen. The denominator argument is p-adic. No assertion
   equates real boundedness with preperiodicity for arbitrary rational
   starting points. It does equate them on the invariant odd-quarter lattice.
3. The eight-point example, the period bound in this denominator class,
   and the 8->4->3 fibre obstruction have prior sources. No novelty claim
   is made for those facts. D1's explicit formula and proof are self-contained.
4. The hourglass is the user's clarified diagonal-layer relation. At full
   order four it uses four of six POTENTIAL vertices and four of eight
   LINE edges. At fixed-path order four there are three free potentials
   and no hourglass. These source types are not interchanged.
5. The integral basis proof uses an explicit spanning tree and a unit
   triangular change of basis. It proves an iff over every abelian group,
   including all composite moduli, rather than inferring this from mod2 rank.
6. The cube chart is explicit and marking-dependent. Its fold is a kernel
   relation, not a commuting complement automorphism. The two-bit minimum
   includes the stated next-value observer; the current source needs three
   bits. The chosen quadratic update fails on both tournament quotients.
7. The finite FRACTRAN source predicate has one state prime. The clock23
   and source tag29 never occur in a denominator. Multi-control integers
   are deliberate off-domain hostiles. The CA retains only time modulo3.
8. Wheel formulas have a genuine exception at pure powers of two: a cycle,
   not a union of paths. Joint-offset masks preserve labels and local
   forbidden residues. Neither local unit masks nor repeated residues test
   primality or bound integer dynamics.
9. Prime-stream saturation of finite sum rows is positive inherited
   information. It does not imply every row is nonempty, and modular
   complementary residue pairs may be aliases of larger integer sums.

## Exact universes and reproduction

From this directory:

```powershell
python -X utf8 verify.py > normal.json
python -X utf8 -O verify.py > optimized.json
python -X utf8 audit.py > audit.normal.json
python -X utf8 -O audit.py > audit.optimized.json
python -X utf8 retained_control.py > retained_control.json
python -X utf8 seal.py
```

* Producer: exact finite-graph pruning for all12,514 c=5 mod8 from-99
  through100005;12,502 generating-function/support-formula checks.
* Independent path: walks G_k directly in a deliberately loose interval
  for every k=-20..1500 (1,521 parameters), with memoized finite-orbit
  decisions;3,150 denominator controls for other odd numerator classes.
* Integral layer basis: triangular integer unit pivots for layer sizes
  1,...,m, m=2..12. Independent determinants for all117 layer-size tuples
  of lengths2..4 with entries1..3, including singleton boundary layers.
* All256 binary and6,561 ternary line words on the order-four full
  layer graph are checked against actual potential reconstruction.
* All eight quadratic cube states, all24 vertex permutations per tournament,
  H/Aut fibres, both quotient failures,72 source/clock machine steps,
  and a multi-control hostile are checked directly. An independent exact
  determinant verifies det(I-zA)=1-z^3 and twelve return traces.
* Wheel producer: every even m=4..2000, plus1160,1890,2310,5510,30030;
  independent circular-block scanner:496 moduli with odd-prime factors.
  Joint-port check:762 nonempty masks, including repeated residues and
  pure powers of two. Finite primality uses exact trial division.
* The inherited actual81/27 tournaments are regenerated, their automorphism
  groups recomputed, and their67 labelled cyclic triples compared. H is
  462915 for both; the named port8 probe after a sink gives2390013 and1754919.
  This replay imports the prior verifier and is explicitly inherited, not
  an independent implementation of that control.

The source pins and manifest use UTF-8 bytes with CRLF normalized to LF,
so Git's Windows checkout conventions do not change the mathematical seal.
The seal also checks local Markdown links and matching normal/-O receipts.

## Repairs made before final filing

The first exploratory assumption interpreted the hourglass as the strong
four-vertex tournament; the user corrected it to the diagonal-layer relation
before theorem filing. The free-bit order was then aligned with the actual
diagonal layers: central bit(0,3), leaves(0,2),(1,3). Both quotient hostiles
and all receipts were recomputed. A command initially ran the verifier from
the repository root; it failed without executing, created two empty redirected
files, and was rerun from this directory. Only those confirmed empty files
were removed. None of these failed commands supplied a certificate.

Incoming order-ten carrier data have two separate certificates: the earlier
lower-bound packet alone does not prove equality. The later
`carrier10_upper/PROOF.md` supplies the56-column integer kernel in each of45
endpoint blocks and hence the2520-dimensional upper correction. We cite that
later proof explicitly; no new replay of its large certificate is claimed.

No known false canon statement is changed here. The failures tested here are
research candidates and explicit scope controls, not established results
being retracted. No new theorem ID was reserved or promoted.
