# Session35 independent Lie/arithmetic audit

Status: PASS for initial Lie and arithmetic notes; no major errors found. Independent scripts and receipts are adjacent. No Git mutation performed. Metric-completion addendum, if issued later, is audited separately below.

## Lie checks

The invariant-form proof is valid. The algebra is perfect, hence its central z pairs to zero with every element. E7-equivariance kills B(V56,E7), since these irreducible modules are nonisomorphic. Invariance and [E7,V56]=V56 then kill B(V56,V56). Consequently invariant symmetric forms are precisely pullbacks of scalar Killing forms and their nonzero radicals are exactly H57. Zero forms of course have the whole algebra as radical.

Independently generated rational-coordinate E8 roots give grades (-2,-1,0,1,2)=(1,56,126,56,1). Complements in the grade-one set pair to theta and the only grade-two root sum is theta. Cartan accounting gives190 versus191, as stated. Counts alone do not prove bracket nonvanishing; the source note correctly separates that structural input.

Directly recovered corrected THM-4671 through Git object read. N57 enumeration and both expressions for coefficient A agree exactly for (2,115),(3,58),(4,39),(7,20). First-defect weight117 and zero/nonzero split are valid inherited applications. The note correctly refuses to distinguish all three nonboundary singular supports merely by the nonzero scalar. Boundary vanishing applies to the adjusted lift, not the unadjusted vector.

Primary source checks:
- Landsberg–Manivel, https://people.tamu.edu/~jml/LMsexpub.pdf , intermediate-algebra construction and printed p158 convention: E7 acts on the symplectic56, and its Heisenberg extension is the radical. Their H56 notation counts the symplectic space, unlike H57 total-dimension notation here.
- Dong–Mason, https://arxiv.org/pdf/math/0201318 , Theorem1 and §3: strong rationality requires strong CFT type, C2-cofiniteness and rationality. §2 explains strong CFT type by nonnegative integer grading, vacuum-only degree0, and L(1)V1=0. The cited full-V1 reductivity obstruction applies.
- Useful stronger consequence from the same source, §2 pp3–4: a simple strong-CFT-type VOA has a nondegenerate symmetric invariant form on V1. Thus the proved invariant-form theorem already excludes this full V1 without rationality or C2-cofiniteness. Simplicity and full-weight-one scope must remain explicit.
- Lee–Sun–Wang https://arxiv.org/html/2306.09230v1 , Conjecture1.1 is indeed conjectural, with positive levels1–5 and formal c=190k/(24+k). Kawasetsu https://arxiv.org/abs/1305.6463 verifies the intermediate vertex-subalgebra and modular-differential-equation statements. Neither source licenses silently replacing these objects by a strongly rational VOA with the specified full weight-one Lie algebra.

## Arithmetic checks

All elementary claims pass. For coprime opposite-parity e,o, L=o^e+e^o=a²+eb² with a=o^(e/2), b=e^((o−1)/2). Every prime divisor is coprime to ab, so −e is a nonzero quadratic residue. Square e is incompatible with a Blum prime; e=2s² forces Blum prime divisors to3 mod8. The displayed squarefree-d norm coordinates expand exactly.

For odd y>=3 with3∤y, N=2^y+y² is divisible3 and N/3 is3 mod8. The smallest permitted y is5, so the quotient is at least19 and cannot equal3. Thus Blum iff quotient prime is correct. Period18 residues mod9 and the y19 hostile factorization agree.

The ring map Z[sqrt(-2)]→Z/57 taking sqrt(-2) to13 is a well-defined surjection, kills5+4sqrt(-2), and the principal ideal has index57 by determinant. The quotient isomorphism follows. Both signed norm products and CRT idempotents are correct.

CRT square-root fibres have effective regular F2^s actions, where s counts nonzero local square targets. For a Blum prime, the quadratic-residue subgroup has odd order, so squaring is invertible there, providing the unique principal root. The explicit exponent in the note is correct. Primary source https://cacr.uwaterloo.ca/hac/about/chap2.pdf confirms Facts2.157–2.158 for principal-root uniqueness, with Fact2.160 the equivalent permutation/inverse assertion. Citation refinement was communicated to author.

Recovered the inherited literal decoder and physical central deletion from zero_heavy_ports_20261020g/microcosm.md §§2,4 and metagraph_fibre_network_20261020k/README.md §3. Independently instantiated all states for k2..5:3364 compatible parent/target comparisons verify actual paths, actual vertex deletion, and bijection of newest-bit lifts to every matching square-root fibre modulo57. No whole-metagraph map is inferred. Boundary stabilizers, ordered coordinates, fixed parent and target are essential; the note retains them.

Independent trial division reproduced the complete111-pair census below10^12 and all12 Blum rows. Completeness follows y<=39 since x>=2, and x<=10 since y>=x+1 and11^12>10^12. Independent modulo57 buckets yield9 four-root,10 two-root,1 singleton and37 empty targets. Jacobi symmetry and absence of a field of57 elements are correct.

## Reproduction

Run python3 session35-lie-arithmetic-audit/independent_checks.py from the workspace root. It writes independent_receipt.json, verifies rational-coordinate roots, coefficient identities, direct trial-divided census, and literal deletion/CRT comparisons. The finite checks support but do not replace the cited Lie/VOA inputs.

## Metric-completion addendum: independent PASS

Audited COMPLETION_ADDENDUM.md including generalization §7. Its finite-dimensional complex Z-grading, injective graded original Heisenberg bracket, and invariant nondegenerate symmetric form homogeneous of degree0 are all essential retained hypotheses. No unrestricted minimum or uniqueness is inferred.

At the lower pairing count247 the complete module profile is forced: L0=s, L±1=V,V*, L±2=C. Homogeneous nondegeneracy makes B|s nondegenerate. The solvable radical is stable under the grading derivation and thus graded; its degree0 part is a solvable ideal in simple s, hence zero. The nonzero s-valued bracket [V+,V−], forced by invariance and nontrivial s-action, excludes both irreducible degree±1 modules from the radical. Pairing the original [V+,V+] bracket with z− forces [V+,z−] nonzero, excluding R−2. These steps correctly leave only R=0 or Cz+.

For R=0, the innerness of derivations of complex semisimple Lie algebras gives D=ad h. Since [D,D]=0, ad(Dh)=0 and centerlessness gives Dh=0, putting h in L0=s. Since D|s=0, h is central in s and zero, contradiction. For R=Cz+, idealness and the degree profile make R central. Levi decomposition is then a direct sum with the central line, so no derived element can lie in that line, contrary to the original nonzero Heisenberg bracket. Thus247 is impossible.

E8 supplies the248 completion. The proof depends only on simple s and nontrivial irreducible symplectic V, so the general bound dim(s)+2dim(V)+3 is valid; attainment is not automatic. The parabolic's trace58 excludes any invariant nondegenerate form independently of grading. The coadjoint cotangent construction has dimension380 and the claimed homogeneous metric. Numerical checks are only dimension controls, not a substitute for Jacobi/Levi arguments.

Optional exposition refinement communicated: grading alone makes the radical graded; homogeneous B provides the tight opposite-grade dimension bounds. Their roles should be distinguished, though the theorem as stated is correct.
