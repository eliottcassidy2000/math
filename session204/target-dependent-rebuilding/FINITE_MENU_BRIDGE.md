# A finite-menu bridge for target-dependent rebuilding

Status: RESEARCH CRITERION ONLY. No successful menu, closed guarded language, or all-target certificate is supplied. This note does not prove an eventual threshold. The underlying weighted-automaton constructions are elementary verification devices, not a claimed new theorem of tournament spectra.

## Physical representation

Let T(w) have apices0->1 and a transitive chain indexed by w in {0,1,2,3}*. A chain letter records which apices that vertex beats. The endpoint-forest construction specified in the appendix, truncated to at most one bulk-start component, gives fixed nonnegative integer28-by28 matrices A_a, initial column x, and nonnegative terminal row ell with

    H(w)=ell A_(last letter)...A_(first letter) x.

Later appended chain vertices cannot enter a bulk-start component, so the truncation is safe for this one-path readout. A rewrite p phi(w) q is an actual new tournament: it reorients the old apex incidences, can replace/delete chain vertices if phi maps letters to words or the empty word, and adds the stated end contexts. No preservation of the old induced tournament is claimed or needed.

## Guarded version

Choose a regular witness language L, finitely many regular guards G_i subset L, physical rewrites f_i(w)=p_i phi_i(w)q_i, and signs sigma_i in {-1,+1}. Require:

1. L has literal witnesses with H=23,25,...,45.
2. For each sign, the corresponding guards cover L.
3. f_i(G_i) is contained in L.
4. H(f_i(w))=2H(w)+sigma_i for every w in G_i.

These conditions would prove all odd counts at least23. Indeed induction on odd n>=47 uses the unique odd parent h=(n+1)/2 when n=1 mod4, and h=(n-1)/2 when n=3 mod4. It satisfies23<=h<n. A witness in L can use its covering guard, and image closure gives the next witness in L. The finite seeds cover the initial interval. No finite scan of n is a substitute for conditions2--4.

The language checks are finite automaton inclusion checks. Each count identity can be checked by a reachable-span calculation of dimension57 per guard-DFA state. Its state is

    z_w=(F_phi(w) F(p_i)x, F(w)x, 1),

with letter matrix diag(A_phi(a),A_a,1) and terminal row

    (ell F(q_i), -2ell, -sigma_i).

Here F(v) is the chronological matrix product for word v. For erasing letters, A_empty=I. Maintain exact rational spans for words reaching each DFA state. At most57 times the number of DFA states independent vectors can be inserted. Test the terminal row on the accepting-state spans. Thus every accepted word is checked, not only a finite sample.

## Selector-free version

There is a stronger sufficient test when no regular choice rule is known. Let L be every word with H(w)>=23. This is regular: saturate all nonnegative forest coordinates at23 after each letter. Saturation commutes with the nonnegative integer transfers, and the clipped H readout tells exactly whether H>=23. There are finitely many clipped states, however impractical their worst-case number.

For each sign sigma, choose a fixed finite menu of physical rewrites f_1,...,f_m. Put

    d_i(w)=H(f_i(w))-2H(w)-sigma.

Pointwise menu coverage is exactly

    product_i d_i(w)=0 for every w in L.

Each d_i is a rational weighted series of dimension57. Their pointwise product has an exact tensor-product realization of dimension at most57^m. A reachable-span zero-equivalence check, restricted to the eligibility DFA, therefore provides a finite, all-word certificate of the existential assertion: at least one physical rewrite succeeds on each input word. The successful index need not be supplied by a regular selector. This is only a finite verification criterion, and its naive dimension can be prohibitive.

Image closure here is automatic: a successful branch has H=2h+-1>=45, so its physical output lies back in L. Consequently the literal seeds plus successful product-zero certificates for both signs would prove threshold23.

Finite menus with additional conjunctive requirements can be handled by a sum of squares before taking the product. For example, a branch consisting of a count-preserving normalization g_i followed by f_i requires both H(g_i(w))=H(w) and H(f_i(g_i(w)))=2H(w)+sigma. Over integer evaluations, the sum of their squared defects is zero exactly when both equations hold. A product over branches again expresses existential success. This observation does not supply the missing normalizations.

## Current experimental boundary

- All6400 rewrites with arbitrary letter maps4->4 and at most one added letter at each end fail to produce a stateless binary compiler on tested seeds through chain length5, already at the second level.
- Full-word-dependent strategy choice improves matters; per-target existential paths improve it further. These quantifiers must not be interchanged.
- On parent fibers through chain length8, every odd h in23..511 has a representative admitting both children by the tested end-rebuild or free-interior-insertion menu. This is finite one-step readiness, not class closure.
- A freely chosen interior insertion is not one fixed morphism p phi(w)q. It must be encoded by an explicit guarded/transducer rule before the finite criteria above apply. No arbitrary position choice has silently been treated as a fixed finite menu.

No actual menu has passed an all-word check. The criterion is not the missing proof.

## Appendix: explicit 28-state convention

A state pattern is the unordered multiset of endpoint pairs of directed path components. An endpoint is port0, port1, or a bulk symbol*. A port may appear as the start or end of at most one component, with the singleton pair(i,i) permitted. A port internal to a path does not occur among the endpoints. Retain only patterns with at most one component starting at*. Empty patterns are included formally.

To introduce a chain vertex z of type a, its legal predecessor is either absent or the end of an existing component: every bulk end is legal, and port i is legal precisely when bit i of a is zero. Its legal successor is either absent or the start of an existing component at a port i whose bit is one. A bulk start is never a legal successor. Join the chosen components through z; choosing both endpoints of the same component is forbidden, as it would form a cycle. Absent predecessor/successor allows creating a singleton or extending only one end. Then replace endpoint label z by* and discard patterns with more than one bulk start. Sum multiplicities when distinct choices have the same resulting pattern.

This is a nonnegative integer transfer. Every resulting spanning path forest has a unique inverse after deleting the newest vertex and splitting its path there. A discarded bulk-start component can never be entered by a later newborn, so discarded patterns cannot contribute to the final one-path readout. Hence the transfers count H exactly at every word length.

The initial column has value1 on patterns {(0,1)} and {(0,0),(1,1)}, and0 elsewhere. The terminal row is1 on one-component patterns and0 elsewhere. The number N(w,b) of formal patterns on w active ports with at most b bulk starts obeys N(0,b)=b+1 and

    N(w,b)=3N(w-1,b)+N(w-1,b-1)+2(w-1)N(w-2,b),

with N(w,b)=0 for b<0. The cases are: the first port is internal, a singleton, a start ending in bulk, an end reached from bulk, or paired in either direction with another active endpoint. Thus N(2,1)=28. This specifies the matrices used by the research criterion without requiring neighboring research files. It does not instantiate a successful rewrite menu.
