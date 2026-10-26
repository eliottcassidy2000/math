# Independent audit: equinumerosity, symmetry, and operations

**Status:** PROVED elementary observations / FINITE-EXACT independent verification. Not a literature-priority claim. Prepared independently of the producer's script, against repository `4950589d056a64f8807774cf7da893dcafcd63a0`.

## 1. Inheritance and exact typing

Anchor: tournament/sign-even-graph correspondence. Niche: operation compatibility and cocycle trivialization. Wildcard: which generator families become indistinguishable only after taking isomorphism classes?

Concept board: (1) affine tournament torsor; (2) sign-even dual Fourier lines; (3) ordinary Eulerian cycle space; (4) stabilizer-containment Hall graph; (5) weighted operators versus Boolean support; (6) gauge/coarsening.

Closest proved mechanisms:

- `THM-4531-twisted-mallows-sloane-principle-and-no-natural-bijection-for-switching-classes.md`: every invariant edge subspace U; affine quotient T/U; signed monomial Fourier basis indexed by graphs in U-perp; stabilizer matching criterion.
- `THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md`: restricted odd-group Burnside, Sylow refinement, integral-localized Walsh basis. **The counting/Fourier basis is inherited.**
- `THM-4722-no-equivariant-bijection-tournaments-even-graphs-for-every-n-sylow-hall-violator.md`: no S_n-equivariant map on labelled tournaments inducing the sign-even graph type bijection for n>=5.
- `THM-4724-affine-fixed-point-criterion-restricted-twisted-burnside-for-every-submodule-even-n-switching-half-and-the-dc-obstruction-to-uniform-circulant-rules.md`: already gives the F2 affine fixed-point criterion and H^1 odd-group averaging for every U.
- `THM-4725-no-uniform-dilation-rule-for-prime-circulants-signed-set-burnside-and-the-p-mod-4-dichotomy-with-the-paley-rule.md`: already gives signed-set Burnside for arbitrary paired sets.
- `THM-587-per-level-signed-cycle-index-antipodal-euler-number.md`, and `THM-588-metagraph-laplacian-algebraic-connectivity-four.md`: ordinary tournament single-arc weighted quotient spectrum and gap4 already known.
- `THM-4069-even-graph-basis-dependence-and-canonical-cycle-envelope.md`, `THM-4078-even-graph-triangle-quotient-spectrum-and-boolean-noncommutation.md`: **different graph class**, namely ordinary even-degree/Eulerian graphs. The latter's native operation is cycle toggling, not arc reversal on tournament switching classes.

Canonical hostile: n=5 Hall deficit, plus n=4 nonclosure under XOR. Corrected near miss: unrestricted odd-H refinement fails at C15 on ten vertices (32 versus31); the all-odd-prime-labelled circulant correspondence is a surjection inducing a type bijection, not a labelled bijection. Least-used sidecar: operator multiplicity, affine cocycle, and the origin's stabilizer.

External fact checked directly against the primary abstract on 2026-10-26: Royle--Praeger--Glasby--Freedman--Devillers prove equality of **isomorphism-type counts** of tournaments and sign-even graphs; sign-even means every graph automorphism reverses an even number of edges. [arXiv:2204.01947](https://arxiv.org/abs/2204.01947), journal J. Algebraic Combin.57 (2023),515--524. No broader external novelty claim is made.

## 2. Three different graph state spaces

For n vertices and m=binomial(n,2):

1. All labelled graphs: F2^m, size2^m, a vector space.
2. Eulerian graphs Z=ker(boundary), size2^(m-n+1) for n>=1, a vector space.
3. Sign-even graphs E: automorphism-sign predicate, in general **not** a vector space.

The exact identities are

- tournament types = sign-even graph types;
- tournament switching-class types = sign-even Eulerian graph types;
- ordinary Eulerian graph types = ordinary signed-complete-graph switching-class types, from **untwisted** Fourier duality.

No pair in that list should be silently replaced by another.

Independent enumeration:

|n|labelled tournaments|labelled sign-even graphs|labelled Eulerian graphs|tournament/sign-even types|Eulerian types|sign-even Eulerian types|
|--:|--:|--:|--:|--:|--:|--:|
|0|1|1|1|1|1|1|
|1|1|1|1|1|1|1|
|2|2|1|1|1|1|1|
|3|8|4|2|2|2|1|
|4|64|20|8|4|3|2|
|5|1024|396|64|12|7|2|

At n=3 the two sign-even types are empty and P3; the Eulerian types are empty and C3. Their equal type counts are accidental, and the lists already differ.

At n=4, take X={01,02} and Y={01,02,03}. Both are sign-even: they are P3+isolated and K1,3. But X XOR Y={03} is odd. This is the smallest XOR-nonclosure witness (n<=3 checked exhaustively). It rules out treating all sign-even labelled graphs as a binary linear space.

## 3. What the general affine-Fourier envelope really gives

Let M be a finite abelian group, A an M-torsor, and G act compatibly on A and M. Choose a0 and write g(a0+v)=a0+c_g+rho_g(v), where c_(gh)=c_g+rho_g(c_h). Over C, Fourier characters of M form a monomial basis of functions on A. The stabilizer of a dual character chi acts on its line by the character chi(c_g) (or its inverse under the opposite pullback convention). The orbit contributes one invariant line exactly when this stabilizer character is trivial; otherwise it contributes none.

Thus dim C[A]^G is the number of admissible dual orbits. This is a **linear invariant-space correspondence**, not a point map A->dual types and not a state-count identity.

For a source-independent kernel k:M->C invariant under the linear G action,

    (C_k f)(a) = sum_m k(m) f(a+m)

preserves C[A]^G. The invariant Fourier lines have eigenvalues

    lambda_k(chi) = sum_m k(m) chi(m).

All such kernels form, on the invariant space, the full diagonal algebra on the admissible dual orbits. Proof: finite Fourier transform is invertible and commutes with the linear G action, so an arbitrary invariant function on the dual (in particular any specified values on admissible orbits and zero elsewhere) is the Fourier transform of an invariant k.

Consequently C_k=C_l on the quotient invariant space **iff** their Fourier transforms agree on all admissible dual orbits. This is the exact kernel-visibility criterion. Positivity or indicator-valued constraints on k are not preserved by arbitrary inverse Fourier interpolation; the full-algebra claim is over C, not a realization assertion for nonnegative walks.

For N<=M a G-submodule, A/N remains a torsor for M/N. Pullback retains chi in N-perp. Normalized fibre averaging kills all other characters. Both intertwine convolution with the **summed** quotient kernel

    k_bar(m+N) = sum_(n in N) k(m+n).

No average belongs in this last formula if the original operator is the unnormalized sum C_k. A normalized probability kernel pushes forward by the same sum.

These facts organize the producer's proposed theorem. They are standard Fourier/torsor mechanisms at abstract level; the useful research content is the exact combined operation/coarsening contract and its correctly typed applications, not a new name for Burnside.

## 4. Exact criterion for equivariant operation-compatible point maps

Let A,B be torsors for finite abelian G-modules M,M', with chosen origins giving cocycles c,c'. Let L:M->M' be a group isomorphism. A bijection F satisfying

    F(a+v) = F(a)+L(v)

exists and is G-equivariant iff L intertwines the linear G actions and

    L(c_g)-c'_g = rho'_g(t)-t             for every g

for some t in M'. In that case F(a0+v)=b0+t+L(v). Both directions follow by substituting this expression in F(ga)=gF(a). Thus the precise obstruction is equality of cohomology classes after transport by L, not equality of dimensions/cardinalities.

**Special case B=M with its linear action.** The following are equivalent:

1. an arbitrary G-equivariant set bijection A->M exists;
2. A has a G-fixed point;
3. c is a coboundary;
4. a G-equivariant M-translation-compatible bijection A->M exists.

Indeed 1 pulls the fixed point0 back to a fixed point, and a fixed point is an equivariant origin. Thus there is no hidden nonlinear workaround when the target is the underlying linear module.

**Tournament control.** For G<=S_n acting faithfully, a G-invariant tournament exists iff |G| is odd: an involution swaps some vertex pair and cannot preserve an orientation; for odd G no edge stabilizer can reverse its endpoints, so choose orientations independently on pair orbits. For odd G, choosing such a tournament gives an equivariant bijection between tournaments and **all** labelled graphs, respecting XOR/reversal. It is gauge-dependent, not canonical. For even G no equivariant set bijection to all labelled graphs exists. This differs from the type-bijection problem with the much smaller sign-even target.

## 5. Equal counts and even equal characters do not recover equivariance

For finite G-sets X,Y, an orbit-bijective equivariant map X->Y exists iff their orbit-indexed bipartite graph, with an edge G/H~G/K when H lies in a conjugate of K, has a perfect matching. A labelled equivariant bijection requires the stronger condition that corresponding stabilizers be conjugate. This is the inherited Hall mechanism.

Independent n=5 census: five tournament types with automorphism orders3,3,3,3,5 have only four eligible sign-even graph types. The maximum matching has size11 of12. For n<=4 a complete matching exists. The theorem excluding all n>=5 is inherited THM-4722, not proved merely by this finite check.

Even equality of full permutation characters is insufficient in general: take G=C2 x C2, its three order-two subgroups H1,H2,H3, and

    X = G/{1} disjoint union 2(G/G),
    Y = G/H1 disjoint union G/H2 disjoint union G/H3.

Each has six points; identity fixes six, every nonidentity element fixes two, and each has three orbits. Their complex permutation representations agree. But X has two G-fixed points and Y none, so no equivariant bijection exists. Full subgroup marks, not only element-fixed counts, determine finite G-sets.

## 6. Linear origins detect generators; affine quotients can hide them

**PROVED generator-saturation criterion.** Let G act linearly on an abelian group M, and let Q_S be the loop-deleted simple edge image on M/G of the translation edges x--x+s for s in S. In the undirected case replace S by S union(-S). Then

    Q_S=Q_T  iff  G(S union(-S)) minus{0} = G(T union(-T)) minus{0}.

Sufficiency follows by relabelling an edge. Necessity follows from the neighbors of the fixed origin [0], which are exactly the nonzero generator orbits. For F2 modules no minus signs are needed. This gives the common mechanism behind THM-4069: the generator-orbit spectrum is cycle length, and a tree basis realizes exactly lengths3 through diameter+1.

The fixed-origin hypothesis matters. For tournament types at n=3, the invariant kernel toggling one edge and the kernel toggling two edges are different generator orbits but induce the identical weighted quotient operator

    [[2,1],[3,0]]

in the order transitive,cyclic. Full reversal acts trivially on both types, and admissible sign-even dual graphs all have an even number of edges. The Fourier-visibility criterion explains the equality. Thus one cannot transplant the linear-origin recovery theorem to affine orbit spaces.

## 7. Independently audited new switching-class gap and explicit observable

The producer proposed the following corollary. Its proof is sound and independent exact checks pass.

Let A be the weighted single-arc reversal operator on **unlabelled tournament switching classes**, retaining each of the m=binomial(n,2) arc choices with multiplicity and retaining loops. Put L=mI-A. Its eigenvalues are

    2|E(X)|,

one per sign-even Eulerian graph type X. This is the affine-Fourier theorem with quotient by the cut space and dual equal to the Eulerian cycle space.

For every n>=4, gap(L)=8, with multiplicity one. A nonempty Eulerian graph has at least three edges. The sole three-edge possibility C3 is sign-odd. The sole four-edge possibility C4 is sign-even. Thus the minimal admissible nonzero support is uniquely C4 plus isolates. For n<=3 there is one source orbit and no positive gap.

For skew Seidel K, define

    F(T) = tr(K^4) - n(n-1)(2n-3).

This is invariant under relabelling K->PKP^-1 and switching K->DKD. Closed four-walks with repeated vertices contribute exactly

    2n(n-1)^2 - n(n-1) = n(n-1)(2n-3).

The remaining walks traverse simple four-cycles, eight oriented/rooted walks per cycle. Therefore F is eight times, up to a single reference phase, the signed C4 orbit Fourier sum. It is nonzero for n>=4 and L F=8F. This gives an explicit simple Fiedler observable.

Normalization: nonlazy **algebraic** gap is8/m (n=4:4/3, not an absolute-gap assertion); the half-lazy gap is4/m and its relaxation time m/4. This is not THM-4078's ordinary-Eulerian triangle-toggle operator or its gap2(n-2).

Exact quotient matrices and quartic values include:

- n=4: A=[[4,2],[6,0]], F=(8,-24), L spectrum(0,8).
- n=5: A=[[5,5],[3,7]], F=(40,-24), L spectrum(0,8).
- n=6: L spectrum(0,8,12,16,16,24); F=(120,24,-8,-40,-72,-72) in the recorded canonical order.

## 8. Reproduction and verification scope

Run from this directory:

```bash
python3 audit_even.py > normal.json
python3 -O audit_even.py > optimized.json
cmp normal.json optimized.json
python3 audit_switching.py > switching-normal.json
python3 -O audit_switching.py > switching-optimized.json
cmp switching-normal.json switching-optimized.json
```

`audit_even.py`: every labelled graph/tournament and vertex permutation n=0..5; no filters. **28,663** active exact gates: orbit counts, label counts, signed traces, graded edge census, transported Walsh invariance, every Fourier flip eigenfunction, exact quotient intertwiner, power-of-two determinant, full Hall matching and witness, XOR-nonclosure, full symmetric obstruction and odd cyclic positive controls.

`audit_switching.py`: every root-gauged labelled switching class and every Eulerian graph n=1..6; every vertex permutation. **19,530** active exact gates: type equinumerosity, signed transport, full Walsh basis, quotient eigenfunctions, simple gap8, quartic observable and exact global phase. No repository computation was imported.

Both scripts contain no assert nodes or floating literals. Normal and optimized outputs are byte-identical. The scripts write only output streams; they make no repository edits or network calls.

Source/outputs at this audit stage:

- audit_even.py SHA256 `97e1fa8837e8dd10e6dcd62397494543c78cd531716f93d143747ff97299e967`.
- normal.json = optimized.json SHA256 `cc37ef01079f1d31b3bc298921f689e629d94326c8c69afc798262110a9cbc7c`.
- audit_switching.py SHA256 `8acdfb9a9ce07db74d3edcccdf961acbaaaef0b787d6161a6be07cc993f4b043`.
- switching-normal.json = switching-optimized.json SHA256 `0bfb26315416c3f28a8b9047ed7192848031d0a19ae54e1247a8896c78b3e91b`.

Blockers: targeted repository message peek could not run because `.machine-id` is absent; no fabricated registration was made. This does not affect mathematical or computational evidence. Producer final report and root synthesis remain pending audit; this report is not a premature signoff on unseen files.

## 9. Additional independent abstract controls and producer proof audit

`audit_abstract.py` exhausts all30 subgroups of AGL(2,2), every invariant orbit-indicator difference kernel, and every one of the24 candidate point bijections for each group. It checks the signed invariant basis, full diagonal kernel algebra, the equivalence between a fixed origin and equivariant point-bijection/trivialization, and every generator subset pair for the linear-action saturation criterion. **1,306 exact gates pass**, normal/-O byte-identical. This independently checks the general theorem on affine actions not limited to tournaments.

- Script SHA256 `9373c7989867afdb9f19909311abff186dc21574582489aaaabb994e6c9c41f8`.
- Both output modes SHA256 `39e6f1fe8b2cc029d4e144cd24977b618393f6c464741936dd77cd0dde18a2ab`.

The producer's `EVEN_GRAPH_OPERATION_ENVELOPE.md` at SHA256 `b7c345c3bc4ed239fb1cef644832d17ec49d28c13fed10faf313d7bb44a071c1` received independent proof audit **PASS, no required theorem corrections**. Audited: Theorems A1--A3 including the pullback sign `(g chi)(c_g)`, origin independence, rational Fourier projection kernels, exact visibility iff; Theorem B including summed pushed kernel versus normalized fibre average; source-independent/translation-commuting iff; Corollary C and quartic observable; correct algebraic versus lazy normalization; inherited distinction between type maps, labelled maps and weighted operations.

The producer's verifier was independently rerun in normal and optimized modes. Both outputs are byte-identical to its supplied **81,177-gate** receipt. All15 source snapshot SHA256 values and every snapshot's exact bytes were checked against its pinned git blob at `4950589d056a64f8807774cf7da893dcafcd63a0`. The producer's script is not counted as an independent implementation; the three independent scripts above provide that evidence.

## 10. Further audited specialization: quotient by Eulerian reversals

This specialization was proposed by the producer after the preceding proof audit; its proof has also been independently checked.

Quotient tournaments by reversal along every Eulerian edge set Z. The resulting labelled torsor is canonically the outdegree-parity vectors p in F2^n satisfying sum(p)=m modulo2: an edge reversal toggles exactly its two endpoint coordinates, and the incidence map has kernel Z and image the even-weight subspace. Relabelling is ordinary coordinate permutation. The unlabelled states are therefore weights r with r=m modulo2.

Single-arc reversal is the exact weighted birth-death operator

    A[r,r-2]=binomial(r,2),
    A[r,r+2]=binomial(n-r,2),
    A[r,r]=r(n-r).

The dual objects are cuts K_(a,n-a), taking 0<=a<=floor(n/2). Every unequal-shore cut is sign-even. A balanced cut admits part exchange and is sign-even iff a is even, since that exchange reverses a^2 edges. Thus omit precisely the balanced odd-a cut. The full Laplacian spectrum is simple:

    {2a(n-a): 0<=a<=floor(n/2), not (2a=n and a odd)}.

For n>=3 its gap is2(n-1), simple; for n<=2 there is one state and no positive gap. The Krawtchouk vector

    K_a(r)=[z^a](1-z)^r(1+z)^(n-r)

is the corresponding eigenvector. The identity K_(n-a)(r)=(-1)^r K_a(r), restricted to the allowed parity r=m mod2, accounts exactly for the missing balanced character. The balanced odd-a vector vanishes **identically on the allowed r sector**; individual coordinates of other vectors can vanish, so an unqualified statement that only this vector 'vanishes' would be too loose.

`audit_scoreparity.py` independently checks every labelled tournament n=0..5, the actual outdegree-parity map, surjectivity/uniform cycle-space fibres, every native endpoint flip, and every Krawtchouk eigenvector entry/exact full rank for n=0..32. Normal/-O outputs byte-match.

- Script SHA256 `62dbe34801a9788f45475938be321b0a5e2655d3935d262b726eddd7d5f62264`.
- Both outputs SHA256 `f49bd5439fe3edc427b10b5c07d52748d845a508a46cac2a46f27deba3234154`.

This is an elementary coarsening consequence of the operation envelope, not a claim that Krawtchouk theory is new. It shows that the same source operation has three distinct observable gaps after three different quotients:4 on tournament types,8 on tournament switching types, and2(n-1) on score-parity types. These numbers cannot be compared without naming the quotient and retained measure.

## 11. Final audit disposition

The score-parity producer addendum `SCORE_PARITY_QUOTIENT.md` at SHA256 `94ce8274fbd5784a6af52e1ec23edf336bea05bf72b84ede69a1d9e2c88a00f8` passes independent proof audit with no required corrections. Its **116,047-gate** verifier was rerun normal/-O and both receipts byte-match the supplied output. The modular rank modulus1,000,000,007 was additionally verified prime by exact trial division. The independent `audit_scoreparity.py` contributes **16,406** active exact gates.

Root's `WALSH_REALIZABILITY.md` at SHA256 `8685f8db8b91de6ab5e913943ac9b987682c7a42ca229dbd792cd3dde025626f` passes independent audit: explicit n=4 matrix was reconstructed from the independent graph/permutation implementation using the root's bit1=up convention; its inverse was verified as an exact matrix identity; cone inequalities, converse diagonal, sign-gauge obstruction and positive-monomial orthant lemma all check. Root's normal/-O receipts match its supplied receipt. This signoff does not purport to audit a separate U2 positive-cone theorem; only the general comparison to a nonsimplicial cone is elementary here. The full receipt is `ROOT_WALSH_AUDIT.json`.

Total distinct gates in the **four independent implementations**:65,905. Producer rerun gates are additional reproducibility evidence, not added to that independent-implementation total. All four independent normal/-O pairs are byte-identical. All proof conclusions remain scoped as stated; finite audits are not substituted for all-order proofs.

Repository state inspected at close: branch `codex/session72-unifying-carriers...origin/main`, clean. This child made no shared-checkout edits or commits; publication remains the parent task's responsibility. No publication action was attempted or denied here.
