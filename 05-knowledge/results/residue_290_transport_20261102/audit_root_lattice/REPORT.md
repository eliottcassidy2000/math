# Independent audit of the ROOT/lattice reconstruction

**Verdict: PASS U1–U3, P1 and E1–E7 in their stated scopes, after three
precision repairs.** The audit was performed by a different agent, with
an independent implementation that imports none of the producer. It does
not use the disappeared October26 receipts or their reported gate totals.

Reviewed candidate: `root_lattice/PROOF.md`, SHA256
`fd8bc3e5ed62977ed959dc3742b426e972461a859e48069b395e5763e8510138`.
The source revision and every exact source pin are recorded in the manifest.
All ten producer source pins were verified against their Git blobs.

## Proof audit

**U1.** The inherited Beta(1,2) mixture is the original W, not a freely chosen
probability on trajectories. The exact counters L=tau−1 and
K=sum floor((a_i−1)/2) are valid only on completed strict ROOT words. The
adaptive-mixture proof gives the first two base generations and the geometric
contraction of later generations; the current source still states the
original13/3 nonroot bound, with later improvements explicitly distinguished.
Clearing the B2 denominator gives numerator r²+2r³+2r⁴−r⁶. The difference
from r²(1+r+r²)² is r⁴+2r⁵+2r⁶, coefficientwise nonnegative. The geometric
concavity estimate then integrates r^(−2alpha), exactly for alpha<1/2.
The constants imply M_(1/3)<239/6<40. Occupation identifies the lifetime
moment, with no replacement of odd-step time by the literal clock.

**U2–U3.** Telescoping the actual b factors proves the signed mode equation.
Strict first hit implies distinct nonroot odd vertices, so the product is
bounded by the factors at3,5,...,2tau+1. The logarithmic comparison yields
|h|^6<=tau+1, and U1 supplies genuine infinite-support ell2(W) membership.
Every rooted coordinate is determined by its ROOT value, giving precisely
one kernel direction. The maximal-domain decomposition and the failure of
forcing-only ROOT/source coercivity follow. Neither a coercivity estimate
that retains initial storage nor an estimate on the zero-trace complement
is disproved. The passive energy is paid by that original source storage.
For the finite-support obstruction, z_k=(64^(k+1)−1)/3 is an actual omitted
ROOT predecessor divisible by3, with the stated polynomial W weight and
exponentially growing1/d. ROOT evaluation is continuous because W(1)=1,
which proves the graph-norm non-density statement. No self-adjointness or
universal ROOT-entry conclusion is needed.

**P1.** Subtraction of the corrected affine identity removes each common
carry. Its multiplier is3^d/2^T, so the exact integer braid clock T consumes
T source-address bits. All common-word sources share those low T bits.
The bin relabel statement therefore requires **k>=T**; that restriction
was added during this review. For lower k there is no within-class pair
split. Original charges, source occurrences, terminal labels and arrival
filters are retained. Iterated-map controls in the checker are explicitly
not a census of first-V arrivals. The actual17/35 collision at3 has different
words and clocks4/5, so it does not contradict same-word injectivity. No
cross-word or unequal-arrival debt estimate follows. The complete-level
positive-part convention and the uniform debt target's OPEN status survive.

**E1–E5.** Simple-root lifting, the unit determinant beta−alpha and evaluation
prove the paired ring map and norm. Equality of ideal and evaluation kernels
uses inclusion plus the common index7^k, not a hidden unique-factorization
assumption. The one-way exchange at7 is not a two-port involution; its exact
failure locus at higher precision follows from v7(delta−2)=1. Two consecutive
color observations recover m,n moduloq. The constructive CRT lift proves
that the positive coprime opposite-parity Euclid source realizes precisely
primitive vectors. Thus the stated minimal observer cardinality and unequal
color masses have both lower and upper bounds. The real norm and the word
determinant losses are separate. During audit the local U and A branch
formulas were added to make the root-first word calculation self-contained.
Identity versus the shear is a transport-spectrum counterexample, not a
statement about tournament adjacency spectra. Lattice rotations preserve
norm but fail the positive-source guard as stated.

**E6.** At2,5,29 the Eisenstein polynomial is irreducible, so each local
quotient is a field. Conjugation is its Frobenius, and each nonzero norm
value has p+1 preimages. The character-sum proof for odd p is valid because
p!=3 and its completed-square discriminant is nonsquare. CRT gives112
scalar unit values with540 preimages each, totaling60480 ring units. A
unital evaluation to Z/290 is impossible already modulo5. The inert-prime
valuation argument excludes145 and290 as integral Eisenstein norms; the
Gaussian products instead verify both Gaussian norms. No prime clock or
quadratic-dynamics map is inferred from the repeated numerals.

**E7.** Modulo145,12²=−1 and24 is a unit, so the paired Gaussian colors are
invertible and multiply to the norm. The four roots are exactly
12,17,128,133. At2 the Gaussian factor is the dual-number ring, whose two
units both have norm1; the odd split factors have p−1 preimages per nonzero
norm. Therefore the global unit norm fibres have224 members, totaling25088.
The nonzero nilpotent145(1+i) rules out an isomorphism with(Z/290)^2; that
latter ring is reduced because290 is squarefree. Inertness and ramification
are two different exact obstructions, correctly separated in the note.
The final wording limits the obstruction to a **faithful split ring
presentation**. Nonfaithful observations can still support selected
operations: evaluation at±17 modulo290 has kernel{0,145(1+i)}, and that
kernel is stable under the D shear. The independent checker includes this
hostile so the failed ring-isomorphism claim cannot silently become a
no-observer theorem.

## Independent finite controls

`check.py` uses explicit exceptions under normal and optimized Python. It
implements the U mode by backward edge propagation, verifies every native
passive prefix for1023 odd sources below2048, clears the generation bound as
an exact polynomial, and checks original-W Pascal/tail payments. It compares
all source pairs in common P map-word classes through depth3 below1024.

For split lattices it enumerates the literal image of multiplication by
pi^k and compares that set with the evaluation kernel at7,49,343; inverts all
paired observations; checks exact exchange loci, norm defects and symmetries;
and realizes every primitive residue at49 by a genuine positive Euclid lift.
Full branch matrices, rather than endpoints alone, verify the determinant
hostile. For inert fields it constructs explicit inverses and norm fibres.
For Gaussian145 it checks all145² port states, and at290 it enumerates all
290² vectors, finding exactly two square-zero elements and the norm fibres. It also
checks the nonfaithful±17 observer and its well-defined D transition on
all290² source vectors.

The two independent receipts are byte-identical. Producer normal and `-O`
replays also match their frozen receipts; those replays are not counted as
additional independent gates. Exact totals and SHA256 bindings are in
`MANIFEST.json`.

## Correction lineage and limits

The initial candidate was
`a517258dc64efa5f7b1ddb9150d97ecc2d3676dbc0f4b6de244ffac68e38b3dc`.
The producer's `REVISION.json` records the k>=T/arrival-scope repair, local
Euclid convention repair, and root-authorized additive E7 section. The subsequent `WORDING_REVISION.json` records the restriction from a
failed general two-color observer to a failed faithful split ring
presentation. Its repaired candidate above was reread in full; the tests
bind it explicitly. Status/link-only promotion is separately inverted in
the audit manifest, recovering exactly these reviewed candidate bytes.

The all-order conclusions use the audited arguments and scoped inherited
measure facts. Finite controls alone prove neither a moment bound over the
infinite tree nor universal Collatz convergence, a uniform primitive ROOT
debt estimate, abc, RH, or a field-with-one-element theory. None is claimed.
