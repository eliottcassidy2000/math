# Session72 cross-domain synthesis: coherent frames, affine carries, and realizable observations

**Result first.** There is a precise positive unification beyond “preserve information”: **joint cyclic objects are classified by compatible relative phases, and the minimum sufficient collection of pair tests is governed by prime-power-layer connectivity.** This applies simultaneously to necklaces, chosen-normal-basis Frobenius orbits, tensor factors of finite fields, and an explicitly realized all-arity carousel-tournament family. The proof and finite controls are in [THEOREM-CYCLIC-ALIGNMENT.md](THEOREM-CYCLIC-ALIGNMENT.md).

A second exact subcarrier joins source-preserving arithmetic and dissection chains: **affine displacement is a cocycle, while legality is an independently composed domain.** A linear part, a trace, a common fixed centre, or a formally closed seam need not preserve the target predicate. The third recovered mechanism is older and sharper than an analogy: **the factorial moment functional has no fixed finite-dimensional multiplicative/linear carrier in characteristic zero**, even though a p-dimensional carrier exists in characteristic p. These give both constructive theorems and limits for a grander theory.

No new solution of Collatz, root-27 survival, LRC(14), planar Jacobian conjectures, Penrose infinite extension or dissection optimality is asserted. No literature priority is claimed. Theorem CA1–CA3 and GA below are elementary proofs; finite controls are separately labelled. Native Penrose geometry and other inherited theorems retain their existing scoped audits, not an assertion that this lane reran all of their code.

## 1. Provenance, corrections and prior synthesis

Repository pin: `4950589d056a64f8807774cf7da893dcafcd63a0`. `SOURCE_PINS.json` records exact Git paths, blob identities, SHA256s and unchanged snapshots. All authored files are in the isolated `session72-cross-domain/` lane; no shared-repository edits or commits were made. The parent owns publication.

Startup read AGENTS, START-HERE, relevant CURRENT-FRONTIER, ACTIVE-GUARDRAILS, RESEARCH-PROTOCOL, CORE-PAPERS and META-PATTERNS. The inbox command was attempted but could not run because `.machine-id` is absent; no identity was fabricated. This is a coordination-tool limitation, not missing mathematical evidence.

**Credit before extension.** The incoming reflection `07-reflections/why-the-advances-work-observation-composition-and-non-backtracking-macmini-20261026.md` already gives `(X,M,O,rho)` and five questions: reconstruction, descent, composition, cancellation and level. Earlier presentation-first and missing-coordinate frameworks are credited there. This lane does not rename that schema as a new theorem. It adds concrete all-arity alignment proofs and an economical observation criterion, and recovers the underapplied factorial-annihilator boundary.

The reflection is provenance, not a proof authority. In particular, at this pin:

- `THM-4807-untraced-arcs-are-cycle-closing-chords-elementary-legal-path-reconstruction-endpoint-decoder-one-guard-separation.md` is **PROVED and independently AUDITED**, with Rédei/Camion inputs and concurrent Session70/BJR-based reconstructions credited. This lane neither re-proves it nor claims its discovery.
- `THM-4808-berggren-tree-is-non-backtracking-branch-character-pole-frame-every-modulus-obstruction-cusp-pairing.md` uses “PROVED” in its filing text **and explicitly says “Pre-audit at filing.”** Its non-backtracking/spectral package is therefore **PRE-AUDIT / NOT A DEPENDENCY HERE**, not promoted because a reflection calls it proved. We use the older audited THM-3333 Gaussian-square theorem instead.
- Session66's correction overlay overrides the historical claim that the root-27 rank is determined solely by the first disagreement. Full matches require the residual state valuation; only fixed-size frontier optimality is established.
- MISTAKE-658 forbids interpreting a numerical coincidence or modular syntax as a target-preserving bridge. MISTAKE-587 repairs the Collatz/frieze reading; THM-4604 is read in its corrected form. Session71 keeps actual tile placements distinct from infinite extension and formal length closure distinct from geometric legality.

Inheritance portfolio:

| Role | Recovered mechanism | Hostile / boundary |
|---|---|---|
| Anchor | Session68 normal-basis/Frobenius and two-necklace relative phase | Separately classified necklaces lose gcd-many alignments |
| Niche | THM-2810 factorial-Hankel faithfulness | `s-s^14` vanishes in C13 exponent quotient but has factorial readout `-87178291199` |
| Wildcard | Session66 source-fixed affine switch cocycle, Session71 partial-affine domains | Macro86 source7 fuel `5 -> 2 -> 9`; equal formal length identities have unequal domains |

The live concepts were: action intertwiner; relative frame; cocycle; semantic realization; coefficient-field rank; native guard. Their precise roles appear below rather than being conflated into one undifferentiated obstruction.

## 2. Positive all-arity theorem: coherent relative phase really can be enough

For cyclic factors of periods d_i, the joint quotient is

    (product_i Z/d_i) / simultaneous translation.

Its classes are exactly edge residues `r_ij=x_j-x_i mod gcd(d_i,d_j)` satisfying

    r_ij+r_jk=r_ik mod gcd(d_i,d_j,d_k).

The number of classes is `product_i d_i / lcm_i d_i`. Keeping only edges of a graph E is complete **iff every set of factors divisible by a given prime power induces a connected subgraph of E**. This both proves pairwise sufficiency in the correct category and specifies exactly which pair tests can be discarded.

Examples matter:

- Equal periods d: an ordinary spanning tree suffices.
- Periods `(3,9,9)`: a star through the period-3 factor loses the phase modulo9 between the others.
- Periods `(15,21,35)`: all three edges are necessary, because primes3,5,7 require different pair edges. There are105 joint classes. This stays inside the odd-period carousel category.

The native construction is `B_e=A^{e1}->...->A^{er}`, with `A=1->C3`, substituted into an odd carousel R_m. Strong-component sizes decode the r-bit colour; maximal modules recover the outer blocks. Thus r fixed necklaces of periods d_i give exactly `product_i d_i/L` tournament types on4rm vertices, each with automorphism order `3^(rm)m/L`. No arbitrary tournament classification follows.

**The general relative-frame carrier.** For homogeneous G-sets `G/H_i`, aligned joint states form `G\product_i G/H_i`, not `product_i (G\G/H_i)`. For two factors the joint quotient is the double-coset set `H_1\G/H_2`: move the first representative to its reference and retain the second relative frame, modulo the two stabilizers. Session71 makes this literal for two named matching involutions a,b: their joint conjugacy classes are `C_G(a)\G/C_G(b)`, with stabilizers of all retained local labels when those labels vary. Keeping only ab is weaker still: on two side names, `(id,id)` and `((01),(01))` both have product id, but their matching diagrams have two components and one component respectively. No polygonal realization is implied for an arbitrary involution pair. CA1 and CA2 solve an economical special case of this general alignment quotient, not every G-action.

Over finite fields, the normal-basis coordinate map conjugates rotation to Frobenius; the same phase quotient indexes factors of

    tensor_i F_(q^d_i) = product_(product_i d_i/L) F_(q^L).

This is a genuine action-and-stabilizer transport. It is substantially stronger than equality of cardinalities and substantially narrower than identifying tournament substitution with field multiplication.

**Why this is not a contradiction of Session70/71.** A deterministic phase state is not a probability law on hidden placements. Coherent named pair differences can recover the former, while arbitrary pair marginals do not recover the latter. Parity-law hostiles keep that distinction explicit.

## 3. Guarded affine transport: a common exact subcarrier

Let K be a field, V a K-vector space, and let `(A,b,D)` mean the partially defined affine map

    f(x)=Ax+b for x in D subset V.

Take two such maps f and g=(B,c,E), in chronological order f then g.

**Theorem GA (PROVED).** Their exact composition is

    g o f = (BA, B b+c, D intersect f^(-1)(E)).           (GA1)

Composition is associative as partial maps. The homogeneous-coordinate carrier is

    [[A,b],[0,1]],

but it must still retain its domain. Forgetting b yields the linear action; forgetting D yields only formal affine transport. Neither forgotten coordinate is generally recoverable from the other.

**Proof.** Substitute `g(f(x))=B(Ax+b)+c` and require both `x in D` and `f(x) in E`. Associativity follows from composition of actual functions with their exact domains; expanding either parenthesization gives the same affine formula and the same three domain conditions. QED.

When `A` depends on a word, GA1 says b is a twisted additive cocycle. Changing the origin by `x=y+r` replaces b by `b+(A-I)r`. The formal total affine extensions therefore have a common fixed origin iff every displacement is simultaneously of the form `(I-A)r`. For an actual common fixed point of the partial maps, one must also require `r in intersection D`. In dimension one, two total affine maps with slopes unequal to1 commute iff they have the same centre; to infer commutation of restrictions from that formal identity, equality of the two composite domains is additionally required. Empty or lower-dimensional domains can also make distinct formal maps agree as partial maps; no unrestricted converse is intended. These are elementary affine facts, not new priority claims. The independent auditor supplied the sharp guard hostile: `f(x)=2x` on `(0,1)` and `g(x)=3x` on `(0,2)` share formal centre0, but their composite domains are `(0,1)` and `(0,1/3)`, and0 is not an actual fixed point in either domain.

**Scope of minimality.** General domains can be lower-dimensional, so distinct triples can define the same partial map there. We identify actual partial functions, not all triples. On nonempty open real intervals, distinct affine coefficients are distinguishable; all empty-domain triples represent the same empty partial function.

### Arithmetic specialization: retain the original source and the carry

THM-4604's odd Collatz valuation word has

    G_w=[[3^|w|,B_w],[0,2^(sum w)]],
    B_(uv)=3^|v| B_u+2^(sum u) B_v.

The diagonal/trace sees only length and exponent sum. It loses B_w, which decides merges at integer sources. The actual words `(1,7)` from483 and `(7,1)` from469 merge at17; their diagonal data agree but their carries differ. No convergence statement follows.

Session66 applies GA to the different native guarded process `T(z)=(3z+(z mod2))/2` with credit q. A word is executable only on its exact parity cylinder and while its guards hold. Balanced blocks P have centre `r_P=-C_P/(3^L-2^L)`. Pulling that centre through an executed prefix w of length n gives

    rho_(w,P)=(2^n r_P-C_w)/3^n,
    h_P(T^n(s))+n=v2(s-rho_(w,P)).                        (GA2)

The source s remains fixed; a change of shape changes the comparison centre, not s. For P=100,Q=101000, the centre defect is684. On the actual source7 orbit at macro86, P consumes fuel5 to2 and changing the reference to Q raises fuel to9. It is an equal-valuation cancellation, not a new source. The exact carrier GA2 is not a globally decreasing rank.

The corrected full-match hostile is especially small: word1000 matches both sources5 and21, with `Q_5=432` and `Q_21=1728`; ranks `4-v2(Q_z)` are0 and-2. The truncated longest-common-prefix data agree. Retaining residual valuation repairs the formula but proves no nonhalting.

### Dissection specialization: retain order and domain

A positive split of total L sends `x -> L-x` on `(0,L)`. GA1 produces precisely the sign, offset and admissible interval of Session71. Words `(2,2)` and `(3,3)` are both the formal identity, but have domains `(0,2)` and `(0,3)`. Input5/2 separates their legality. Thus even zero affine defect need not mean executable closure.

For no-flip Euclidean seams, the affine group is SE(2), and connected flat seam data develop placements up to common rigid motion. This is a frame-development theorem, not a nonoverlap theorem. Session71's five unit squares with successive origins `(0,0),(1,0),(1,1),(0,1),(0,0)` have a tree seam graph and legal local opposite-side seams, yet piece1 and piece5 coincide. Sum of areas5, union area4. Geometric acceptance remains a separate predicate after cocycle consistency.

This gives a definite proof-transfer rule: GA1 transfers along an actual affine intertwiner **and** a map preserving/pulling back the allowed domains. A matrix identity without domain preservation transfers formal motion only. The arithmetic parity cylinder, length positivity interval and planar nonoverlap condition are different native domains, not names for the same constraint.

## 4. Gaussian/Pythagorean norm: the representation is the sidecar

The audited THM-3333 uses

    Phi(m,n)=(m^2-n^2,2mn,m^2+n^2),
    <Phi(u),Phi(v)>_Lorentz=2 det(u,v)^2.                 (GN1)

On **canonical raw images of primitive projective integer rays**, Farey adjacency `|det|=1` is exactly Lorentz pairing2. Gaussian squaring is an actual quadratic representation map, and compatible changes of linear frame induce the symmetric-square action. It is not additive in u and it is not an equivalence between arbitrary null vectors and integer rays.

The norm-only quotient loses representation: `(8,1)` and `(7,4)` both have norm65, while their Lorentz pairings with Phi(1,0) are2 and32. Norm65 therefore does not decide the same adjacency predicate. Primitive normalization can change the raw pairing by content factors, and forgetting ordered legs can identify distinct adjacent vertices. The sidecar is ordered represented norm plus raw scale, not merely another scalar invariant.

This is the geometric counterpart of relative phase: a scalar invariant can correctly classify a single coarse object and still fail to classify its placement relative to another object. The equation GN1, not that general slogan, is the actual transport. No implication to minimum dissection pieces, golden inflation or LRC arrival is made.

## 5. Recovered older mechanism: finite radial carriers have a sharp obstruction

THM-2810, `01-canon/theorems/THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go.md`, is the deliberate older recovery. It concerns the characteristic-zero functional

    L:K[s]->K, L(s^n)=n!.

Its Hankel matrix H_d=((i+j)!) factors as `D P P^T D`, where D has factorial diagonal and P is the lower Pascal matrix. Thus

    det H_d=product_(j=0)^d(j!)^2 !=0,
    Ann(L)={q: L(qg)=0 for all g}=0.                     (FM1)

Every nonzero polynomial is eventually observable after multiplication. Consequently no proper polynomial-algebra quotient, no fixed finite-dimensional multiplicative carrier and no fixed-dimensional time-homogeneous **linear** state represents the full factorial tower. Matching n! through n=2d already requires at least d+1 linear-state dimensions.

The exact boundary is characteristic p:

    Ann(L_p)=(s^p),

so K[s]/(s^p) is the sharp p-dimensional multiplicative carrier. This follows from the same Hankel factorization on degrees0,...,p-1 and the vanishing of n! for n>=p. It does **not** say that a characteristic-zero factorial observer can be replaced by a fixed finite-field quotient, and does not substitute for the whole-face arithmetic argument of any Gaussian-moment theorem.

Two hostile controls make the distinction operational:

- In the cyclic height quotient s^13=1, `s-s^14` is zero but `L(s-s^14)=-87178291199`.
- `q=6s-s^3` has L(q)=0 while L(q^2)=504. A linear kernel is not automatically a multiplication-stable ideal.

A clock-dependent recurrence `(n+1)!=(n+1)n!`, a state whose dimension grows with truncation, and nonlinear state encodings are outside this no-go. For example `(n,v)->(n+1,(n+1)v)` is a two-coordinate nonlinear clocked realization; this theorem is not a universal ban on all finite-dimensional descriptions.

**Connection to the tournament/even-graph anchor.** A finite carrier can support exact Fourier linearization over a suitable coefficient field; that is a different question from whether an unbounded moment tower has a bounded carrier. Walsh linearization of an F2-torsor requires2 invertible, while FM1 has its own factorial characteristic-p boundary. These are distinct typed examples of coefficient-field dependence. They provide no universal special-prime explanation and no identification of their obstruction ideals.

## 6. Portfolio transfer ledger

| Source -> target | Explicit map and operations | Preserved predicate | Destroyed coordinate / required repair | Decisive test |
|---|---|---|---|---|
| r cyclic words -> joint phase class | Simultaneous rotation; r_ij modulo gcd | Joint isomorphism/stabilizer | Independent word quotients lose relative phase | CA1 triangle compatibility; CA2 prime-layer graph |
| Binary normal-basis coordinates -> finite field | Sum epsilon_i alpha^(2^i); rotation -> Frobenius | Orbit/degree/stabilizer | Minimal polynomial forgets chosen root; bad basis loses action | F8 normal basis(3,5,7); tensor rank27 |
| r-bit block colours -> carousel tournament | Ordered joins, then prime R_m substitution | Actual tournament type and Aut order | Separate row necklaces lose synchronized alignment | SCC decoder and module recovery |
| Guarded arithmetic words -> partial affine maps | Homogeneous matrix plus exact source cylinder/credit | Executed endpoint and clock | Trace loses carry; free rational map loses integer guard |483/469 merge; source7 macro86 refill |
| Split chain -> partial affine interval | Ordered composition `(A,b,D)` | Exact positive-length response | Formal identity loses admissible interval |(2,2) vs(3,3) |
| Local Euclidean frames -> developed placement | Products along seam paths, flatness | Consistent relative placement | Global nonoverlap/target union | Five-square path returns onto first |
| Gaussian rays -> Pythagorean raw null image | Phi, compatible symmetric-square action | Squared determinant/Farey adjacency | Norm or unordered primitive triple loses representation/scale | Norm65 gives pairings2 and32 |
| Radial polynomial -> multiplicative observation class | q~q' iff L((q-q')g)=0 for every g | All future moment readouts | No nonzero quotient ideal in char0 | Hankel full rank; cyclic13 alias |
| Hidden Penrose/tournament law -> low-order guard moments | Actual named subset moves/re-markings | Selected local expectations | Higher correlations | Parity laws; full native conjunction distinguishes |

None of these rows asserts that all their sources are equivalent. They share proof obligations with explicit successful maps or exact failures.

## 7. What is proved, what remains proposed

**PROVED in this lane:** CA1 all-arity cyclic compatibility/counts; CA2 sparse prime-layer completeness; CA3 native all-arity tournament realization; GA partial-affine composition. **INHERITED/AUDITED:** the source-fixed switch and rank correction, Gaussian/Lorentz map, Session71 geometry/guard protocols, factorial-Hankel obstruction. **CITED:** normal-basis existence and elementary finite-field Galois facts; external source links are pinned in EXTERNAL_SOURCES.json. **FINITE-EXACT:** the declared independent controls in check.py, not extrapolation of a census.

**Falsifiable next proposal, not a theorem:** for a concrete geometric assembly family, an economical faithful state might consist of local metric templates, a relative-frame orbit, affine domain guards and a finite global contact certificate. To test it, freeze the allowed motions and terminal target first, then search two realizable assemblies with identical proposed states but different next-seam legality or final target-union predicate. A flat five-square overlap already refutes the version lacking global contact data. No universal polynomial-size state or bounded locality claim is advanced.

**Source-fixed rank next obligation, still OPEN:** bound accumulated switch gains on the one original source7 ray while proving coverage by the selected legal blocks. Arbitrarily good local shape ranks cannot be summed without their actual cancellation defect. The full infinite centre bank merely re-encodes lifetime (Session66), so its existence is not independent progress toward nonhalting.

The actionable grander programme is therefore: **first intertwine the operation; then solve relative alignment; then certify the native realizability locus; finally ask the observer's algebra/field whether finite compression exists.** Each step has a proved positive model and a hostile control here.

## 8. Reproduction and exact boundaries

Run in this lane directory:

    python3 -B check.py > check.normal.json
    python3 -O -B check.py > check.optimized.json
    cmp check.normal.json check.optimized.json
    python3 -B check.py --inject-failure
    python3 -O -B check.py --inject-failure

The last two commands must exit nonzero. No repository research module is imported; only Python standard-library integer, Fraction, finite-field-bit and exact Bareiss arithmetic is used. Source bytes are checked offline before mathematics. All finite universes and counts are included in the paired JSON receipts and VERIFICATION.md. The tensor matrix test checks multiplication, not merely equal dimensions. The native graph test checks adjacency-defined modules, not merely a construction's own labels. General proofs are above and in the theorem note.

Meta-patterns used: controlled forgetting needs a sidecar; same representation is not same carrier; type every analogy; retain incidence; keep observer and recurrence class; compute the actual repair quotient; re-evaluate after fibre-changing operations. CA2 is concrete evidence for a refined card: **optimize aligned tests prime layer by prime layer, not on an unweighted connectivity graph.** Its counterindications are noncyclic actions, law-valued marginals, and geometric realization.
