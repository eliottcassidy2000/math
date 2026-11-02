# Historical recovery: RH-adjacent ideas, five vertices, and 5/11

Source commit: `291b4d9757842c02aa451beb65acfade2848609b`.
Dates below are source-claimed session dates, not git-author dates. All source
paths are recorded with git blob pins in `source_pins.tsv`. Historical prose
is not inherited without the stated scope/correction. No novelty claim.

## Inheritance board

- Anchor: Session138 exact weighted stability is proved at orders <=9 and
  >=1500; the finite interval 10..1499 is not closed by this lane.
- Closest proved mechanism: weighted odd-cycle/determinant expansion; the
  old zeta-versus-incidence distinction (June 15).
- Hostile: the exact weighted cospectral pair in WEIGHTED_ZETA_BOUNDARY.md.
- Corrected near miss: treating OCF as an ordinary Euler product; treating
  common Bernoulli/Paley/5/11 labels as a transported inequality.
- Least-used sidecar recovered: edge/complement triangle incidence M23,
  the weighted digon analogue of the old disjoint-triangle packing defect.
- Other live concepts: triangle curl, determinant p-adic tangent, transfer
  discriminant, Paley character positivity, marked metagraph continuations.

## Source ledger, oldest first

### 1. March 6: low-order polynomial real-rootedness

`01-canon/theorems/THM-021-discriminant-real-roots.md`, opus S15:
independence polynomial I(Omega(T),x) has real nonpositive roots for binary
n<=8. At n<=5 its degree is at most one; at n6,7 Turan's triangle-free
bound on the disjoint-triangle graph gives the quadratic discriminant.
This is finite combinatorial real-rootedness, not the classical RH. The
historical sentence "for n>=9 degree>=3" must read "degree can be >=3":
transitive tournaments have constant polynomial at every n.

Connection to weighted n5: both are governed by support budgets. But the
weighted model introduces 2-cycles in determinants; the first hidden
2+3 incidence appears at five, whereas two odd cycles need six. The map
is through cycle collections, not an identification of their zero sets.

### 2. March 13: Paley-11 spectral flatness

`01-canon/theorems/THM-162-paley-flat-spectrum-H-maximization.md`, kps S60:
a VERIFIED finite circulant census at p7,p11, not a universal extremal
result. Later normalization repair in THM-1926 must supersede its rough
spectral text: nonprincipal adjacency eigenvalues (-1±i sqrt(p))/2
have modulus sqrt(p+1)/2; skew eigenvalues have modulus sqrt(p).
The phrase "2-transitive property" is not inherited: an orientation-
preserving tournament automorphism group cannot be transitive on ALL
ordered distinct pairs, because it preserves which ordered pairs are arcs.
Flat common-neighbor counts are the relevant design datum.

No transfer to the sharp weighted Delta/V constant follows from spectral
flatness alone. The weighted hostile retains even the entire spectrum.

### 3. March 16: genuine 5/11 discriminant, then analogy overreach

`01-canon/theorems/THM-224-golden-exceptional-points.md`, kps S115:
for M(x)=[[1,2x,0],[0,0,1],[1,x,0]], its characteristic polynomial
f=lambda^3-lambda^2-x lambda-x has

    disc_lambda f=4x(x^2-11x-1).

The two nonzero exceptional parameters are (11±5 sqrt(5))/2;
the double eigenvalues solve d^2+d-1=0. This is an exact map from a
specified transfer matrix to a discriminant. Here 11 and 5 have a common
carrier: disc_x(x^2-11x-1)=125=5^3. They are NOT the vertex orders five
and eleven, and the source gives no Delta-stability-preserving map.
The source's search-absence/priority claim is not inherited.

`04-computation/two_five_eleven_s90u.py`, opus S90u, March16:
explores Cunningham chain 2,5,11,23,47 and k-nacci discriminants.
Its exact survivor is recurrence/discriminant algebra, not a tournament
forbidden-value theorem. It hardcodes a wrong tetranacci value -283;
for lambda^4-lambda^3-lambda^2-lambda-1 the exact discriminant is -563.
For the pentanacci polynomial it is 9584=16*599; parametrically,

    disc(lambda^5-lambda^4-x(lambda^3+lambda^2+lambda+1))
    =-16x^3(4x^4-20x^3-14x^2-553x-16).

Its array entry disc_samples[99] is not x=1 (the grid is 0.01..5 with
200 points), so that earlier printed line must not be inherited as a
point evaluation. The later direct x=1 calculation is distinct.

### 4. March 17: explicitly RH/zeta-adjacent reflection

`07-reflections/riemann-zeta-tournament.md`, opus S74, with source/output
`04-computation/riemann_zeta_tournament_connections.py` and
`05-knowledge/results/riemann_zeta_connections.out`.
Recovered ideas: Bernoulli denominators; carry valuations; odd-cycle
collections; Paley Gauss sums; k-nacci traces; complement/functional-
equation analogy. None is an RH proof or an equivalent RH criterion.

**Current correction banner, August25 / MISTAKE-506:** denom(B6)=42 is
exact, but a common Bernoulli/tournament/Hurwitz obstruction is unsupported;
Hamilton's quaternion algebra does not ramify at 3 or7, and B114 repeats
denominator42. Shared prime labels transport no tournament invariant.

Other explicit non-inheritances:
- OCF is a sum over vertex-disjoint odd-cycle collections, not the ordinary
  product over all odd cycles; the independence/overlap constraint matters.
  The true dynamical Euler product uses primitive closed walks instead.
- k-nacci traces hitting 7,21 do not cause their forbiddenness (the original
  reflection itself records this correction). The same sequence hits
  attainable H values. Binary forbidden values are combinatorial theorems.
- The output's equation "(p-1)|2k iff p splits completely in Q(zeta_2k)"
  is not correct; e.g. p3,k3 is ramified. It is not used here.
- Reversal invariance is an involutive symmetry, not a map to the zeta
  functional equation's analytic variable or to its nontrivial zeros.

### 5. June 15: the useful spectral/correlation mechanism

`07-reflections/the-zeta-function-and-the-ocf-read-complementary-halves.md`
(monad-explorer), THM-499 (June13), THM-502/505 (June15):
Z_T=1/det(I-uA) records primitive closed-orbit totals; H's OCF additionally
records disjointness/overlap. Binary n5 is spectral, n6 has its first hidden
triangle-pair incidence. The reflection's same-day correction states that
simple cycle counts alone no longer suffice at n9: pair overlaps and triple
packings detach. The later THM-505 frontmatter further repairs its dimension
claims: the trace-basis count six overcounts the n9 packing-carrier rank
five, and packing-vector dimension is not the dimension of the single H
readout (which factors through cycle-count level sums). No unqualified
non-spectral dimension formula is inherited. Therefore no all-n simple-
cycle-only sidecar is inherited.

**Actual transfer in this lane:** weighted digons lower the first boundary
from 3+3=6 to 2+3=5; exact enhanced map (P,M23)->H proves it, with a
continuous rational orthogonal-conjugation hostile.

### 6. July 14–16: metagraph preservation and spectral symmetrization

`00-navigation/METAGRAPH-PRESERVATION-AVENUES-2026-07-14.md` requires
preserved predicate, continuation, and repair sidecars. MPA-04,06,12 already
warn that weighted node shadows can fail Markov composition and current
state equality can fail future equivalence. MPA-11 recovers mixed Haar
signs missed by margins. This is methodological precedent for M23, not an
unproved equality between it and a particular metagraph cocycle.

`00-navigation/METAGRAPH-ATLAS.md` routes to
`01-canon/theorems/THM-924-walk-reciprocity.md`:

    cpK(y)=2^(n-1)[cpA((y-1)/2)+(-1)^n cpA((-y-1)/2)].

It follows directly from A+A^T=J-I and hence extends to weighted A.
Symmetrization about Re(z)=-1/2 is an exact finite-matrix identity; it
neither identifies this line with Re(s)=1/2 for Riemann zeta nor restores
lost M23. The new pair is already cospectral before symmetrization.

### 7. July20–21: eleven spectra and the dynamical zeta

`THM-1455-why-eleven-spectra.md` is distinct from the colliding
`THM-1455-skew-pfaffian-parity-and-a-mod-16-law.md`. Its eleven is a count
of skew spectra at n7, not arithmetic modulo11. The parent lane is
recovering its F2 coboundary, first support gap, and erroneous "every
4-subset odd" Paley wording; do not duplicate an unreviewed promotion here.

`THM-1926-tournament-zeta-euler-product-over-strong-core.md`:
P factors over SCCs and is 1 on acyclic blocks; Z is the primitive-closed-
walk dynamical zeta. Oct19 correction keeps c3,c4,c5 spectral in binary
models but c6 nonspectral at every n>=6; Oct23 repairs adjacency/skew
normalization. The common rational pair is the weighted analogue.

### 8. August25: p-adic tangent and external margin boundary

`THM-4093-rational-edge-diagonal-gauge-and-padic-tournament-zeta-tangent.md`
proves W=DAD^-1 leaves closed-walk zeta unchanged; open-path weights
retain endpoint ratios. Its binary P=1-c3 u3+O(u4) gives
v_p(Z(x)-1)>=3v_p(x), equality iff p does not divide c3.
The weighted repair is P=1-Vu2-tau u3+O(u4), so degree2 is first unless
V vanishes. See the exact all-ring and p-adic statement in the theorem note.
This is not a theorem about Kubota--Leopoldt zeta or special L-values.

`THM-4089-hybrid-padic-zeta-margin-optimization-and-next-case-obstruction.md`
only optimizes the external manuscript's displayed real margin and gives
parameter-retuning obstructions. It marks the claimed irrationality theorem
AUTHOR-CLAIMED / UNREFEREED. No theorem about RH, or geometric/adelic input,
is inherited here. Its p11 cutoff boundary is not a tournament vertex count.

### 9. October17: a real arithmetic map, with sharp scope

`THM-4666-class-numbers-are-ranking-excesses-of-kronecker-tournaments-and-local-optimality-of-the-natural-paley-order-excludes-a-real-zero.md`
(audited, citation repairs MISTAKE-626): in the natural order of Paley P_p,
single-vertex insertion distance t changes ranking excess by -2S(t), where
S is the Legendre partial sum. If every S(t)>=0, summation by parts gives
Fekete polynomial positivity on (0,1), and the Mellin integral gives
L(s,chi_-p)>0 for real s>0. This rules out real positive zeros for the
stated condition, not all complex off-critical-line zeros or RH/GRH.
At p11 natural order has33 forward arcs while the optimal order has35:
local optimality is not global optimality.

Primary corroboration consulted Nov2: Alexander Kalmynin,
[Quadratic characters with positive partial sums, arXiv:2105.06910v2](https://arxiv.org/abs/2105.06910v2).
The primary abstract supplies the definition of L+ and its density upper
bound; the insertion/stability comparison here is not claimed by it.

## Connection cards (source / target / map / preserve / loss / test)

1. **Closed walks -> sharp n5 weighted margin.** Map A->P then enhance with
   M23; preserves V,tau,u,v, and enhanced map preserves H,Delta. Loses
   disjoint support incidence in bare P. Cheapest hostile: explicit A,B.
2. **Binary p-adic tangent -> weighted formal residues.** Map determinant
   expansion in the same parameter. Preserves coefficient identities in
   every residue ring; first nonzero coordinate changes degree3 to degree2.
   Loses real positivity under modular reduction. Cheapest test: fair pair
   P=1-u2/4, which is not tangent at degree3 at odd primes.
3. **March transfer matrix -> 5/11 arithmetic.** Map characteristic
   polynomial->discriminant. Preserves repeated-root locus; loses eigenvector
   incidence and does not encode Delta. Cheapest test: evaluate at x1:
   -44; compare the unrelated n5 weighted spectral pair.
4. **Paley ranking -> Dirichlet real-zero exclusion.** Map rank move->
   character partial sum->Fekete polynomial->Mellin integral. Preserves the
   implication nonnegative all partial sums => no positive real zero.
   Loses all complex zero locations and the full Hamiltonian path objective.
   Cheapest exact p11 check: residues{1,3,4,5,9}; all partial sums nonnegative,
   but natural ranking not globally maximum. No Delta map supplied.
5. **Metagraph quotient -> weighted stability.** The reusable operation is
   controlled forgetting with incidence/continuation sidecars, not a new
   claimed numerical equivalence. Cheapest test: keep P fixed while M23
   changes; the explicit pair proves the repair is needed for exact H.

## Stopping boundary

The recovered sources support an exact first-blindness theorem, a useful
negative test for spectral-only reconstruction, and a corrected all-modulus
coefficient/tangent map. They do not settle 10..1499, RH/GRH, global weighted
stability beyond the existing frontier, or universal prime-label roles.
The strongest next legitimate question is whether a local complement-
incidence hierarchy gives an inductive lower bound, rather than whether
an unmarked spectrum determines the answer.
