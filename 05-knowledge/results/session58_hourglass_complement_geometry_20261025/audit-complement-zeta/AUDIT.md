# Independent audit: complement charts and zeta bridges

Audit date: 2026-10-25. All work confined to this audit directory; no repository or remote changes. Original papers were consulted only as references, not copied into publication artifacts.

## Verdict

Mathematical claims PASS with the explicitly scoped inherited frame-counting theorem. Both corrected packets are publication-ready within their stated scopes. The absolute producer repaired all reported issues before final freeze: odd-prime field quantifier, full free-gauge argument, and explicit exception gates replacing assertions. Corrected normal and optimized replays are byte-identical; an injected false parity requirement is caught in both modes. The original assertion defect and its demonstration are retained only as correction lineage, not as a live blocker. The separately written independent verifier is also exception-based and passes both modes with injected failures caught.

The unrestricted rationality argument works for any finite base-field cardinality greater than one; the producer deliberately keeps the narrower odd-prime scope inherited from the frame theorem. Section 2 inherits this odd-prime hypothesis from the opening frame statement.

## Independent proof review: two charts

For a chosen total order, each unordered pair has exactly one coordinate, recording its backward direction. Reversing all arcs negates every bit, precisely graph complement. Relabelling vertices reverses the coordinate convention on each inversion, giving the XOR correction absent from ordinary graph relabelling. Thus simultaneous transport of graph and observer is lawful, while discarding the observer is not. One-edge graphs on 01 and 02 at order three encode respectively a transitive triple and a directed triangle. These are not isomorphic even modulo converse; orders at most two have only one tournament type. Reversing observer order changes every direction, so reversal classes of orders suffice precisely for the literal complement/converse-pair chart. The proposed minimality is limited to this decoder, not arbitrary encodings of isomorphism types.

For a marked Hamilton path, consecutive arcs are fixed forward and only nonconsecutive pairs are variable. The index map (i,j) -> {i,j-1} bijects these chords with pairs on n-1 positions. Complement negates chords only, giving K_P(T)=T^op XOR E(P). Reversing tournament and marked path together reflects the graph position set. These are different operations. At n=3 chart complement exchanges a transitive triple with a cycle, which converse cannot do. The supplied eight order-four rows independently replay to fibre multiplicities 1,1,1,5 and path-dependent K_P.

For permutation Burnside, an even cycle has an antipodal pair reversed by a suitable power, impossible for a fixed tournament. If every cycle is odd, the permutation has odd order, hence no unordered-pair stabilizer can reverse its two points. Each pair orbit then admits precisely two consistent orientations. This proves the stated fixed-count formula independently of enumeration. The class-level converse involution gives (t+s)/2; literal converse fixed points and self-converse isomorphism types are distinct notions.

Transporting all linear orders under an isomorphism permutes histogram summands; converse exchanges each complementary pair. This proves all-n invariance and mass n!. For n>=2 order reversal acts freely and pairs equal summands, proving even coefficients. Injectivity of the histogram is only FINITE-EXACT through n=5. The independent matrix implementation (rather than the producer's bitwise permutation rule) gives graph counts 1,2,4,11,34; graph-pair counts 1,1,2,6,18; tournament counts 1,1,2,4,12; self-converse counts 1,1,2,2,8; converse/profile counts 1,1,2,3,10. No all-n completeness follows.

## Independent proof review: frame chart and field quantifiers

Here p is an odd prime, U is a two-dimensional F_p-subspace, A=C+U, and a frame has p-element factors with unique sums covering A. C_p is inherited from the explicitly identified Session57 theorem; no larger-factor or nonprime-subfield extension is claimed.

For each frame, choose b0 in B and translate by -b0. Then B-b0 is contained in U and D+b0 is contained in A. If A is contained in k, this normalized representative lies in k even if the original frame is defined over an extension. Conversely an extension-field translation cannot identify two distinct base-field gauge classes: its parameter is a difference of base-field elements. Therefore classes neither appear nor coalesce after extension.

If translation t stabilizes the two factors without swapping, nonzero t lies in both difference sets, contradicting unique sums. If translation swaps them, the two factors are translates and have a common nonzero difference, again impossible. The parameter-to-unordered-frame map is therefore injective in both cases.

The difference map (b,d) -> d-b is injective: equality gives b+d'=b'+d. In the normalization 0 in B, B is contained in U and D is contained in A; therefore D-B is contained in A and has |A|=p² elements, so equals A. Since 2 is invertible, the entire forbidden translation locus is A/2, for every normalized shape. This proves the claimed identical puncture configuration, not merely equal cardinality. The independent p=3 exhaustion checks all 270 frames, 30 free gauge classes, and all normalized cases; a separate affine-offset embedding in F_3^3 checks each class has exactly 9 bad and 18 good translations.

Over a fixed base k=F_Q containing A, all p² punctures are k-rational and remain individually rational over every F_{Q^m}. Hence each affine-line chart has Q^m-p² points. These schemes are explicit chosen parameter spaces, not a canonical moduli functor. The other repair uses affine d-space over F_p with p² rational punctures and p^(dm)-p² extension counts; its base and dimension differ. The identical zeta denominator when Q=p^d does not identify these schemes or their bases.

## Independent proof review: zeta obstruction and finite etale construction

Assume a finite-type F_p-scheme has counts C(p^n-p²) for all n>=2 and put a=N1-C(p-p²). Its zeta is R(t) exp(a t), where R=((1-t)^(p²)/(1-pt))^C. Rationality forces exp(a t) to be rational. A rational function's logarithmic derivative has no nonzero constant polynomial part (at infinity it is O(1/t)). Thus a=0, forcing N1=C(p-p²)<0, contradiction. This is not merely the observation that the unspecified first count would be negative: rationality uniquely forces it. Reduction preserves all finite-field-valued points, and finite stratification into affine reduced pieces extends rationality to finite-type schemes, including nonseparated ones. The same logarithmic derivative argument rules out nonzero finite polynomial corrections to the logarithmic zeta.

For a finite set with involution, orbit lengths are one and two. Each one-cycle contributes (1-t)^(-1) and each two-cycle (1-t²)^(-1). Spec F_Q and Spec F_(Q²) realize these respectively, giving the asserted finite etale scheme. Arithmetic versus geometric Frobenius does not alter an involution, since it equals its inverse. The construction preserves exactly the involution-cycle information and does not recover tournament arcs or any other structure.

The deletion polynomial x^(p²)-x is squarefree and splits over F_(p²). Its rational roots over F_(p^n) are exactly F_(p^gcd(n,2)), so the deleted counts have the asserted parity. The n=3,p=3 hostile is 24 versus 18. In characteristic two the frame disjoint formula fails because 2t=0; this does not harm the finite-field puncture construction considered separately.

The formal absolute-zeta normalization is consistent: each factor contributes (Q-1)^(r-1) s^r/(s-1); the specified prefactor cancels this power. Holding r=p² fixed is essential. Negative N(1) excludes a literal finite-set cardinality, not every F1 framework. The integer-grid squaring example refutes only restriction of that particular lift: the outside section (-1,0) maps into the removed section (1,0).

## Sources and inheritance

Read AGENTS.md, ACTIVE-GUARDRAILS.md, RESEARCH-PROTOCOL.md, CORE-PAPERS.md, and the inherited Session57 prime_plane_frames.md. External primary checks on 2026-10-25:
- https://oeis.org/A000568 : definition, offset, small counts, and odd-part Burnside formula. Independent proof supplied above.
- https://stacks.math.columbia.edu/tag/04JI : Lemma 58.2.2 gives the Galois-set/etale-scheme equivalence; its finite-field involution instance is explicit above.
- https://www.mathnet.ru/php/archive.phtml?jrnid=mat&option_lang=eng&paperid=201&wshow=paper : primary record for Dwork's rationality theorem, Amer. J. Math. 82 (1960), 631–648. The direct short URL timed out; search retrieved the primary record. This audit uses the classical theorem, not a new verification of Dwork's proof.

## Reproduction and boundary

Run `python independent_audit.py`, `python -O independent_audit.py`, and `python negative_controls.py`. Normal and optimized independent outputs are identical. Negative controls perturb the order-five tournament count in a copied complement verifier and inject an explicit independent failure; both fail in both modes. The original assertion-based absolute failure is documented rather than hidden. Run `python corrected_negative_controls.py` for the repaired absolute verifier: the deliberately impossible odd nonfixed-count requirement is rejected under both modes. Replays are copies, not modifications of producer files. Hash manifests pin all reviewed sources and audit outputs. No source paper is placed in the artifact bundle, and no remote changes were made.
