# Lawful reset transport: an infinite physical bridge

**PROVED:** a general rank-one transport lemma, an exact O/Z/J reset family,
and its prime-output avoidance consequence. **FINITE-EXACT:** shortest resets
and literal graph checks. **OPEN:** Chan–Pak coverage; classification of the
entire reset ideal. No new literature-priority claim is made.

## 1. Random recovery, corrections, and scope

The reproducible selection uses commit `0e140250d07c`, Python Random seed
24620261110, and sorted canon theorem paths containing one of renewal,
continuant, confluence, controllability, semigroup, invariant-quotient,
automaton. The pool has 18 files; sampled indices are 0,10,3,13,16.
`RESET_CHECK.json` records exact paths and SHA256 hashes. This is a
mechanism-filtered random search, not a uniform sample of the whole repo.

Recovered sources, in selection order:

- THM-1770, `01-canon/theorems/THM-1770-the-localisation-lemma-as-first-return-renewal-and-the-pair-only-closure.md`: only primitive least return survives; scalar cancellation defeats atomwise inference. Its retracted conclusions are NOT used.
- THM-4250, `01-canon/theorems/THM-4250-odd-core-suffix-lift-automaton-and-density-one-compiler.md`: a reset acts on retained carries, while a suffix certificate is sufficient, not necessary.
- THM-3186, `01-canon/theorems/THM-3186-full-exterior-continuant-path-convolution-and-cancellation-wall.md`: support alone loses signed cancellation. Following its proved dependency leads to THM-3182, `01-canon/theorems/THM-3182-factorial-gauss-manin-rank-one-reset-and-two-transverse-smith-bands.md`, section 2: an actual rank-one modular transfer plus matching future transfer matrices proves all-tail descent.
- THM-4852, `01-canon/theorems/THM-4852-primegame-phase-law-both-published-lists-agree-selector-cellular-automaton-in-base-29-primorial.md`: corrected denominator-sensitive precision warns that an observer quotient must respect the next operation.
- THM-4878, `01-canon/theorems/THM-4878-graded-spine-grammars-cut-state-automaton-hankel-and-seam-ranks-cutwidth-censuses.md`: lawful Hamilton-count transfer is available, but low rank alone proves no coverage.

A targeted rank-one-reset search then recovered THM-4204,
`01-canon/theorems/THM-4204-rule30-debruijn-reset-and-dyadic-prefix-saturation.md`:
Rule 30 inverse-boundary words 0010 and 1010 reset to constant maps.
Its spatial predecessor theorem is not a temporal-center theorem.

Closest tournament inheritance is Session245 `descent/PROOF.md` and its literal
O/Z/J compiler. The present reset identities are direct consequences of those
matrices, not a claim to invent the matrices. Session193 response rigidity is
an earlier, different two-apex continuation theorem. No claim is made that
finite observability itself is new.

Concept board: legal transition; full-response reset; retained scalar weight;
matching future action; prime avoidance. The canonical hostile is equal
Hamilton count without equal response: the seed and its O extension both
have count 1, but adjoining J gives counts 3 and 5.

## 2. The common principle (elementary, exact)

Let physical objects carry vectors x in a module over a commutative ring.
Suppose each legal letter a has the exact linear update M_a and the desired
observation is l(x). Fix a legal word r with matrix

    M_r = v lambda,

where v is a fixed column and lambda a row. Assume v itself has a physical
representative, and every suffix w under discussion is legal both from that
representative and after r, with the same transfer matrices. Then

    l(M_w M_r x) = lambda(x) l(M_w v)                         (1)

for every finite suffix w. Proof: associate the three matrix factors.
The assertion covers every finite prefix of any infinite continuation; it
does NOT assert existence of an infinite Hamilton path count or a limit.

Call this **lawful reset transport**. The name describes an elementary
mechanism, not a new abstract theorem. Rank one without a legal reset word,
a physical representative, or matching future transitions is insufficient.
Over integers the retained multiplier can obstruct arithmetic universality.
Over a finite field it can be a unit or zero, which is a different outcome.

Exact source-to-target maps:

- Factorial source: x=(M,X,D), ring F_p[v], r is the singular step at p−1,
  v=(1,1,1), lambda=(0,0,s). THM-3182 matches subsequent indexed transfers
  after shifting p and replacing p+s by s. Preserved predicate is the full
  state congruence and therefore scalar-moment congruence. Lost information
  is the two transverse incoming coordinates modulo p. Sidecars are prime,
  parameter shift, unit condition p∤s, and index; this is NOT integer equality.
- Rule 30 source: x is inverse-boundary data and r is 0010 or 1010. A constant
  map E_i=e_i 1^T is a rank-one reset. Nonnegative column-stochastic transfer
  retains unit total mass, so the reset has trace one, even after further
  spatial letters (possibly changing its image boundary). The preserved
  predicate is unique periodic spatial predecessor. Lost information is
  incoming boundary choice. A spatial ring word, boundary interpretation,
  and trace observer are essential sidecars; a named temporal site is absent.
- Tournament target: physical rooted source state x=(D,E,U,d), reset below,
  ring Z, observer D. Preserved predicate is exact Hamilton-path count
  factorization for every legal future word. Lost information is E,d and
  the independent D/U direction. Sidecars are the rooted physical graph,
  legal O/Z/J letters, and the multiplier, which may not be discarded.

These are exact instances of (1), not a map from Rule 30 configurations to
tournaments or a transfer of either problem's unproved claims.

## 3. Exact physical reset family

Use Session245's physical source seed t=(1,0,2,1) and chronological letters

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D−E+U,d,2U−d,d),
    J(D,E,U,d)=(D−E+U,D−E,D+U,D).

Literal construction starts with arc 0→1 and root 0. Each letter adds a
vertex beaten by every old vertex, except Z/J make it beat the old root;
J designates the new vertex as root. Thus all words here are physical.

First, OJ and JZ are the SAME full linear map:

    OJ(x)=JZ(x)=(D+2U,D,D+2U,D).                         (2)

For m≥0 and k≥1 consider

    r=(OJ or JZ) Z^m O^k J.

Then the exact rank-one identity is

    M_r x = (D+2^(m+1)U) v_k,
    v_k=(1+2^k,1,1+2^k,1)=M_(O^(k−1)J) t.             (3)

Proof: put a=D+2U. After the initial pair the state is (a,D,a,D).
Each Z replaces the equal first/third coordinate b by 2b−D, retaining
second/fourth coordinate D. After m copies it is
(D+2^(m+1)U,D,D+2^(m+1)U,D). The positive O block resets the second coordinate
and multiplies the third by 2^k. Final J gives (3). The explicit formula
for v_k proves its physical realization, including k=1.

Consequently for every physical prefix u, suffix w, and these m,k,

    H(u r w)=[D(u)+2^(m+1)U(u)] H(O^(k−1)J w).          (4)

This is a genuine infinite family of physical reset words and future continuations.
In particular OJOJ, OJJZ, JZOJ, JZJZ all have identical matrix

    (3,1,3,1)^T (1,0,2,0).

The four are precisely the length-four rank-one words, and no shorter
word has rank one. The exact rank histograms at lengths 1,2,3,4 are
{3:3}, {2:2,3:7}, {2:11,3:16}, {1:4,2:40,3:37}.
This finite exhaustion proves minimal reset length, not a classification
of all longer reset words. The rank-at-most-one words form a two-sided
semigroup ideal since rank(AMB)≤rank(M); nonzero physical count prevents
the matrix from being zero. The infinite displayed sublanguage and its
OJ↔JZ rewrites are contained in that ideal. Completeness is OPEN here.

## 4. Prime-output avoidance and the physical positive witness

Reachable states satisfy D≥1, U≥2, d≥1, D≥E≥0, U≥E, and U≥d.
Check at t; O preserves them; Z uses U≥d and U≥E; J uses D≥E.
Each operation leaves D unchanged or increases it by U−E≥0; U never
decreases (for Z use U≥d). Thus the first factor in (4) is at least 5,
and the second is at least 1+2^k≥3, even after arbitrary suffixes.

**Prime pruning theorem.** No O/Z/J word with prime Hamilton-path count
contains a contiguous factor in the displayed reset language, or any
word obtained from one of those factors by OJ↔JZ replacements.
Proof: split at the occurrence and use (4); each factor exceeds one.
This gives a target-aware exclusion at ALL lengths, not a finite prime
census. It applies to this grammar, not to arbitrary descriptions of
all tournaments. It neither realizes primes nor excludes their realization.

The literal word OJOJ constructs the six-vertex tournament with row masks
[54, 60, 56, 17, 32, 8]. It is strongly connected and has H=15=5·3.
The last J provides a directed Hamilton cycle: take the unique
root-starting Hamilton path before that J, then close through the new
vertex. Alternatively the verifier checks reachability directly.
Therefore the arithmetic factorization is NOT a decomposition into
nontrivial strongly connected components. It is loss of response memory
under legal rooted attachments. For all suffixes w, its continued count
is 5 H(Jw); this positive example is not an abstract integer state.

The reset simplifies infinite continuation but charges an unavoidable
multiplicative factor. This explains why a successful physical renewal
mechanism may help certify families while being unsuitable for prime
coverage. Prime coverage within this grammar needs a route that avoids the displayed
reset sublanguage (and its rewrites). No exclusion theorem for the entire
rank-one ideal is claimed.

## 5. Reproduction

Run `python3 -I -B verify_reset.py`, then `python3 -O -I -B verify_reset.py`.
Outputs agree byte for byte with `RESET_CHECK.json`. `--fail` must exit
nonzero. The mathematics uses only the Python standard library; Git and this repository
(with the pinned base commit available) are required only for sampled-source
provenance. Git commands are anchored to the script’s resolved repository root
and use full-tree paths, so the output is independent of the working directory.
Standard-library rational elimination certifies ranks; literal
subset DP counts separately constructed adjacency lists. Checks include
all suffix matrices through length five, 18 literal graph counts, exact
reset matrices for a finite m,k grid, and the hostile equal-count pair.
The all-m,k/all-suffix proofs are (1)–(4), not the bounded tests.
