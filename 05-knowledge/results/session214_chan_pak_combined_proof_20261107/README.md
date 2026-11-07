# Session214: forced shuffles, endpoint stability and two-apex surgery

**Status: PROVED structural results with independently checked finite certificates, in the scopes below. The unrestricted Chan–Pak Hamiltonian-path spectrum remains OPEN.** No literature-priority claim is made.

The target is to realize every odd Hamiltonian-path count at least 23. This session combines path merging, endpoint-defect compression, two-apex surgery and the earlier paired-source constructions. It proves three previously numerical extension statements, classifies a bounded-defect interface, and identifies exactly where two concrete construction strategies stop. It does not prove eventual completeness or extend the inherited unrestricted finite coverage range.

## 1. A constructive shuffle closes the general extension bound

[Proof and verifier](forced-shuffle/README.md); [independent audit](forced-shuffle-audit/FORCED_SHUFFLE_AUDIT.md).

Partition a tournament X into nonempty A,B with any arc A->B. Put h=H(X), p=H(X[A])H(X[B]), and let q count Hamilton words listing B before A. An explicit forced merge sends **every prescribed pair of A- and B-Hamilton words** to a nonordered Hamilton word, preserving both input subsequences. Restriction recovers the input pair, so the map is injective and

    h >= p+q,                 q <= (h-1)/2.

This is a positive path-gluing algorithm on the existing tournament, not an exact scalar recurrence for total H. Ordinary merging is inherited; forcing a nonordered image repairs the earlier collision problem.

The inherited insertion-slot count then gives, for a new vertex v->A, B->v,

    H(X+v) >= max(h+p,2h-p) >= (3h+1)/2.

If no A->B arc exists, H is unchanged. The same dichotomy holds for **any induced extension, with any number of added vertices**: preserve h or reach at least (3h+1)/2. The factor 3/2 is asymptotically optimal. These results settle the ordered-path lemma and the general growth/finiteness bounds in HYP-9469(1),(2),(6), not its unrelated one-centre density conjecture.

A positive deletion increment 2k forces the parent count h<=4k-1. Combining this with SCC factorization, Busch's strong-order bound and a complete core census proves the exact all-order parent sets:

- Increment 2: {1,3}, recovering the earlier all-order result.
- Increment 4: {1,5}.
- Increment 6: {3,5,9}.
- Increment 8: {1,5,9,11}.

The last three were previously finite observations. The finite cutoff is proved, rather than inferred from census stabilization.

## 2. Every fixed sole-exit endpoint defect has a finite dyadic description

[Proof, constructions and verifier](path-techniques/README.md).

At a positive endpoint count e_a, require the **output mark itself** to have exactly one outneighbor. For each fixed positive odd defect d=H-2e_a, an exact two-path-cover identity permits removal of certain singleton strong components. For strong outputs, a physical skeleton of order at most

    d+floor(log2(max(1,d-1)))+4

remains. Restoring the removed vertices gives finitely many affine dyadic rays. SCC factorization extends this finite description to arbitrary outputs. Thus the endpoint and total-H spectra at fixed d have only O_d(log X) values up to X, with an effective finite inverse.

At d=3, the total-path spectrum is exactly

    {5} union {6*2^q+3:q>=0}
        union {6*2^q-1:q>=1} union {10*2^q-1:q>=1}.

The packet supplies literal witnesses and physical constructions, including 23,47,95,... and 19,39,79,... . Novelty of those numerical families is not asserted; the contribution is the classification and finite-skeleton proof. Deleting a sole-exit vertex from a cycle tournament does **not** automatically leave a sole-exit marked endpoint. The result therefore does not cover arbitrary bounded-defect cycle-to-path conversions.

## 3. Moving the comb exit escapes one restriction, but leaves exact gaps

[Initial surgery proof](surgery/PROOF.md); [two-component addendum](surgery/ADDENDUM.md); [reproduction](surgery/README.md).

Let F consist of a transitive chain plus at most two apexes, one of which beats at most one chain vertex; the other cut and apex arc are arbitrary. This includes moving the binary comb's unique exit by reversing old arcs, beyond Session213's induced-retention hypothesis.

The exact F spectrum through 27,371 is all positive odd integers except

    7, 21, 7,477, 27,371.

The two new omissions are **all-order exclusions inside F**, not ambient forbidden counts. The proof combines hereditary reduction, Busch's strong-order bound, a finite exact census and a classification of the two-nontrivial-component case as B=>R_m, where B is one-apex and H(R_m)=2^(m-2)+1. Every other odd value in the stated interval has a finite-exact witness.

Allowing the restricted apex two chain exits produces explicit tournaments with H=7,477 on 11 vertices and H=27,371 on 13 vertices. These positive controls show exactly which class restriction fails. No completeness theorem for the two-exit class or minimum arbitrary edit distance is claimed.

Two closure consequences separate these gaps. Since 7,477 is prime, no finite dominance product of F members can realize it. In contrast, 27,371=101*271 is realized by such a product because both smaller factors occur in F. Also, swapping the apex names and taking tournament converse show that every two-apex presentation of either target must give **each** apex at least two incoming and two outgoing chain arcs. These are consequences of the proved exclusions and count-preserving symmetries, not a classification of all realizers.

## 4. The existing paired carriers have an exact synchronization language

[Proof and physical verification](paired-coverage/README.md).

For each of the four inherited carrier pairs at bases 677,759,1547 and 821, the common O/Z/J words preserving equal adjusted bases are exactly the words using only O and Z. The first J creates a signed separation that no later common suffix can repair. This is an all-length signed-cone theorem, not an extrapolation from short-word tests.

It rules out adding synchronized J macros to these particular pairs. Asynchronous words, source replacement, new pairs and other construction mechanisms remain outside its scope.

## Combined frontier and verification

The interfaces must stay distinct. The shuffle theorem controls induced growth but does not prescribe exact output H. Endpoint compression gives an exact finite description only with a sole-exit output mark. Two-apex surgery permits old-arc changes but has explicit class gaps. Common-word synchronization is a restriction on four marked carrier pairs, not on all numerical recursions.

Positive deliverables are the recoverable path-gluing algorithm, the classified physical dyadic rays, all finite witnesses in the stated surgery interval, and the two-exit realizers of its omissions. None supplies the missing unbounded all-odd construction. New cross-block operations, genuinely different marked states and noninduced rebuilding remain possible routes, not proved solutions.

Each packet contains its proof, reproducible code, deterministic receipts and source pins; independent reviews are in the adjacent audit directories. Normal and optimized Python checks agree. The all-order claims rely on the written injections, compression, invariant cones and proved finite cutoffs. Finite checks test those mechanisms and certify their explicitly bounded universes; they are not presented as evidence of eventual completeness.
