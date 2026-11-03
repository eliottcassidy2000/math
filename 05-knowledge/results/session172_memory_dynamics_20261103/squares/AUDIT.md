# Independent audit of Session172 quadratic Hamilton memory

**PASS.** No mathematical corrections requested. Reviewed `NOTE.md` and `check.py`; independently replayed producer normal and Python -O, with byte-identical outputs equal to the declared certificate. Audited source hashes:

- NOTE.md: bef19f532dd0244163a315e5dc4f7d7df0941dd046ecabeef649c0de51c08257
- check.py: cf250701e222fd5984929797fc9a585dbd6d8f87c434e3bcc16232a91421cd12
- certificate.json: f45a5ec57555b3b8c0b4838e6268940654a9d9a95bfb99eb932d2f2ed6c3924e

## Proof review

1. **All-order response formula:** correct for n>=1 and the stated oriented-path convention. Endpoint insertion contributes twice the endpoint count. Interior insertion corresponds bijectively to ordered spanning two-path covers. Reversal and exchange give K_ij=K_ji, explaining the factor two on unordered pairs. Singleton paths cause no exception, including n=1.
2. **Recovery and edge influences:** singleton and pair contexts recover E,K exactly. A two-path cover never internally uses the joining edge, so toggling that edge leaves K_ij unchanged and changes H by precisely 2K_ij with the stated sign.
3. **Exact rank:** correct for every n>=3 over Q/R. Isolated i plus a spanning path with endpoints j,k yields exactly the two quadratic coordinates e_ij+e_ik, with zero E. There are n-2 edges, so a spanning two-path cover must use all edges and equal the two components. For n=3 the triangle wedge matrix has determinant two; for n>=4 the displayed three-wedge identity spans each coordinate. Spanning path graphs then yield endpoint vectors e_i+e_j, and subtracting quadratic records isolates the full endpoint span. The singleton/pair evaluation transform is invertible on this coefficient space. Minimum linear dimension is for the vector space formally spanned by all labelled graph states, not a nonlinear or arithmetic-restricted model.
4. **Equal-profile seed:** graph descriptions, degree sequences, full oriented component profile, endpoint vectors, and 0-versus-4 leaf response are correct. Slot polynomial expands correctly. Labels are necessary interface data; no square-sum realization has been proved.
5. **Arithmetic specialization:** the arriving-neighbor interval n+2<=s^2<=2n+1 is correct for new label n+1. Oriented counts 2,2,2,0,0,6,0,20 and joining coefficients 3,10 agree with replay. All supplied cut-cover witnesses were separately checked for ports, disjoint spanning vertex sets, and literal square sums.
6. **Failure of iterative closure:** both five-vertex sources have no cover with <=2 paths and hence zero one-step record. After two universal old-vertex connectors, the edgeless source has bipartition imbalance three. For the other graph, the three old components must remain intact because there are only two connectors; component orders, connector orders, and independent edge orientations give 3!*2!*4=48. This proves a genuine two-step separation.

## Independent computation

`independent.py` imports no producer code. It constructs coefficients by pairing full permutation path banks on complementary vertex subsets, then compares them with direct permutation Hamilton counts of extensions for all labelled graphs through order four (1,098 contexts). It also checks edge influences, computes the two profiles by permutation/cut enumeration divided by component-order factorials, checks the two-connector example directly among 7! orders, and verifies every published arithmetic joining certificate. Normal and -O receipts agree. Producer replay covers all 33,866 extension contexts through order five, exact rational ranks 1,2,6,10,15, and independent DFS/DP square counts.

## Scope

This is an exact one-new-vertex response theorem and a sharp all-graph linear-memory theorem. It is not a closed iterative state, a polynomial-time coefficient algorithm, a new all-n square-sum existence result, a physical arithmetic realization of the seeds, or a tournament-count spectrum advance. These limitations are stated correctly in the author's note. No forbidden repository path was read or scanned.
