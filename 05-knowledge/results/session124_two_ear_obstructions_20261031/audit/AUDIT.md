# Independent audit: two fixed ears, Session124

**Verdict:** the stated transfer, binary compiler obstruction (including the finite-phase corollary), ternary fixed-bank certificates, and common-spectrum theorem pass this audit within their explicit hypotheses. No claim about the full tournament Hamiltonian-path spectrum follows.

## Independent checks

`independent_audit.py` imports no producer code. It reconstructs the matrices as tensor sums of local append operators. Independently, it enumerates every allocation of spine vertices to the three increasing runs, and separately counts Hamiltonian paths by adjacency/subset DP.

* All 1,365 words of lengths0–5: exact36-coordinate transfer vector agrees with exhaustive partition enumeration; output agrees with adjacency DP. Both conserved prefix-unit slices equal1.
* 155 independent sole-outneighbor checks: decomposition `H=S-t+V`, three repeated-letter2 values, and the directed-triangle lower bound, for all31 seeds over{0,2} of length≤4. One-ear S,t,V are enumerated as two runs, not extracted from the36-state transfer.
* Both supplied rational certificates replay in all36 coordinates. Exact closure dimensions are `1,4,13,20,20` and `1,9,20,20`. Both conserved units and the output row lie in each final defect closure. All certificate suffixes are legal bank descendants; both coefficient sums are0.
* Observable closure dimensions are10,10,18 for{0,1},{0,2},{1,2}. Computing the invariant row closure is an exact future-output quotient construction, not a heuristic truncation.
* Four reconstructed generators are upper triangular in the claimed common basis, with identical diagonal multiplicities8,16,10,2. All16 length-two products have the squared diagonal. The explicit tensor filtration then proves the general-length characteristic polynomial.
* Exact annihilators and initial values verify both hostile controls: `H((01)^k)=(4^k+6k-1)/9` and `H(2330^k)=14*3^k-17`, k≥1. Polynomial Jordan terms occur; the2-mode can disappear on actual states.
* The producer's standalone Fraction verifier was separately executed normally and with `-O`; both report72 identities and4 rejected mutations.

`../receipts/independent_audit.normal.out` and `../receipts/independent_audit.optimized.out` are outputs of our own verifier, not copies of producer receipts. Their byte equality is checked. The certificate input hash is recorded there.

## Proof audit and quantifiers

### Transfer and corrected affine conditions

Every Hamiltonian path has a unique ear order and decomposition `P,a,Q,b,R`. Spine order forces each run to increase. The three choices for the next vertex are summed. Q/R first-vertex tests and P/Q last-vertex tests therefore give exactly the stated36 states and acceptance row.

For a ray `v(N)=N*u+w`, the required identities are `c(u)=1,c(w)=0,M_d*u=B*u,M_d*w=w+d*u`. The obsolete common-fixed-z identities do not imply these. On an actual prefix-unit slice the local operator is a reset with eigenvalues0,1; a2-eigenvector vanishes there. Hence a single full-state ray cannot synchronize two distinct one-letter updates. This does not alone exclude a quotient or phase-dependent representation.

The observable row closure W is stable under right multiplication by each allowed update, so its common kernel is update-invariant and invisible to every future output. Thus quotient identities certify scalar recursions. The finite seed searches reported by the producer remain finite negative searches, not universal nonexistence proofs.

### Binary fixed-length/finite-phase obstruction

A bounded-digit B=2^L scalar recursion on H>1 gives uniformly O(2^m) counts after m appended vertices. In a finite reachable bottom strongly connected phase component, a closed walk can visit every labeled edge. If its union alphabet varies both bits, its repetitions supply four distinct bounded-distance anchors: first Q, first R, last P, last Q. Repetition permits making anchors distinct even when a single letter supplies several properties. All old vertices can be assigned to P; all vertices between the early and late anchor windows have three independent choices. This gives3^(m-C) accepted paths with a constant C independent of repetitions. There is no random-word or typical-growth assumption.

If a bit is constant, at most two letters remain. The B distinct outgoing length-L words at every phase force all2^L words over the same binary edge. Hence the slow words needed by the case analysis remain available after every phase change. This cardinality step is what makes the finite-phase extension valid.

Edges containing3 offer a universal-sink word preserving H, impossible for H>1. Edge{0,1} either has an old y-outneighbor, which supplies the missing R anchor for3-growth, or has y a permanent universal sink and word1^L preserving H. Edge{0,2} is analogous until x beats only y. Then `H=S-t+V`; after at least one2, `t=0`, S stays constant, and V doubles. After the available prefix2^L0^L, S≥3. Two consecutive available2^L updates consequently require digit `-(B-1)S`, outside the allowed range, regardless of changing digit assignments. Thus no fixed digit assignment is needed in this final step.

Assumptions essential to the finite-phase proof: finite deterministic phase transitions; every reachable phase offers every digit; all physical words have the same fixed length; fixed ears and transitive spine. A phase encoding the entire unbounded count is not a finite phase system. Partial digit languages, changing ears, nontransitive additions, and unrestricted adaptive word choices are not excluded.

### Ternary restrictions

For radix3 and H>1, letter3 cannot implement an allowed digit, so the physical bank012 is forced even if digit assignments vary. Repeating1 gives at most O(2^m), whereas bounded ternary digits force at least `3^m(H0-1)+1`. Starting at H=1 is no escape for a complete bank: a positive digit first reaches H>1. The exact sum-defect certificate is stronger within this fixed bank and also excludes a constant nonzero defect via its zero coefficient sum.

For radix9 the seven slow length-two words can be excluded by repetition ONLY if each remains available. Therefore the forced nine-word bank and its certificate establish a **fixed physical bank** obstruction (allowing digit permutations), NOT a general changing-bank finite-phase obstruction. The explicit scope warning in the source is correct and indispensable.

### Spectrum and variable length

The common upper-triangular basis has letter-independent diagonal0,1,2,3, with multiplicities8,16,10,2. Products of L generators therefore have characteristic polynomial `z^8(z-1)^16(z-2^L)^10(z-3^L)^2`, whether diagonalizable or not. A repeatable fixed scalar action with B>1, H>1 and |d|≤B-1 has a strictly nonzero B^k coefficient; its annihilator consequently forces B=2^L or3^L. This claim assumes the fixed action, not merely a phase-dependent bounded choice of digits on successive repetitions.

For stateless fixed-radix banks, unique prime factorization then forces equal lengths. For finite phase cycles, the separate necessary law is `B^k in {2^L,3^L}` with k steps and total appended length L. Neither law is a global attainment theorem.

## Publication hygiene issue (not a theorem failure)

At audit time `two_ear/check.py` main still prints nullspaces for the old, incorrect fixed-z system before running its valid direct-DP checks. These outputs must not be described as corrected affine-accumulator evidence. Removing that legacy calculation or explicitly marking it obsolete would make the package less misleading. `actual_rays.py` and `observable.py` were inspected and do use the corrected common growth-eigenvector equations.

The earlier fixed-bank scope paragraph and later finite-phase strengthening are logically consistent but could be merged editorially to avoid readers mistaking the superseded escape-route list for the final scope.

No new excluded Hamiltonian-path count is proved. In particular, this audit supplies no evidence of a forbidden value beyond7 and21, and does not prove all other positive odd values are attained.


## Curated-copy hygiene resolution

The publication copy removes the obsolete fixed-z computation from `code/check.py`
and marks earlier open-length statements as superseded by the audited length3/4
results. The original research/audit history is preserved outside this package.
Run `../verify_all.py` to regenerate all included proof-checker and audit outputs;
no workspace paths or numerical LP dependencies are needed.
