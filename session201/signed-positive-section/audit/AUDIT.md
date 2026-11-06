# Independent audit: physical signed section

**PASS.** The frozen theorem and its scoped consequences are valid, relative to the inherited complete physical source grammar and canonical balanced inverse. No mathematical correction is required. Universal occupancy remains OPEN; in fact target 7 is excluded from this slice at every height.

Audited frozen PROOF.md SHA256: `782f12f77d57a0ad5bf395000ae9c2a5af38478b6f48cbd4fd968c0786cb7ed7`.
Audited verifier SHA256: `45824c80e200dd2e7aabce1b0915cb37c2b8336b32f64670e95ba324629184c5`.
The complete candidate manifest passed without modification.

## All-order proof review

1. The physical coordinates are used only for nontrivial sources. The singleton is treated as an allowed inverse precursor, not as a point in the nonnegative physical cone. The seed cannot be any of the positive-endpoint states where a terminal birth is inferred. This avoids a hidden seed exception.
2. The slice gives `(a,f,e,p)=(m−n,0,2m,m)` and m>0. Positive endpoint with E≠d forces terminal J. The predecessor's D,E,U are determined, but d is initially unknown.
3. When n=m, that predecessor is balanced and nontrivial; the inherited canonical graph inverse yields balanced v with D=U=m. For n<m, the Z-output inequality excludes Z at the next two predecessor stages. Only then is f=0 used to recover the previously forgotten d. No arbitrary cone coordinate is substituted for an actual graph predecessor.
4. The resulting balanced predecessor `(m,n,m,n)` is physical. Nonnegative E and odd n give n>0; the balanced graph inverse gives actual s with D=n, U=(m−n)/2. Its physical U≥D implies m≥3n. This proves both exclusions and disjointness.
5. Direct matrix multiplication verifies W=OJ^4, the displayed rank-two matrix, chronological JZ=OJ, and all adjoint generators. The exceptional JJ rebase is correctly restricted to its branch. W(singleton) is physical because its first letter O produces the edge seed.
6. The contractions are exact: U_q=6D_s+14U_s gives 14U_s<U_q, and the exceptional branch has U_q=6U_v. These are bounds conditional on an existing slice witness.
7. The scalar equality is valid in both directions: any source D count is 1 or a previous F count after removing terminal O letters; every F is achieved as a physical Z-output D. W maps all sources, including the singleton, into the slice. The inverse gives the reverse inclusion. This yields `{1} union F(S)`, not universal odd-target coverage.
8. The D=7 all-height reduction is sound. Strip terminal O letters; the remaining core has positive E. For Z, `2D_new−U_new≥3d≥3`. For J, `2D_new−U_new=D_old+(U_old−2E_old)≥D_old`, because `U_old−2E_old=a+f+2p≥0`. If D_old≥3 the bound follows. If D_old=1, the nontrivial physical cone and d≥1 force a=e=0,p=1, hence E_old=0,U_old≥2 and the same bound follows. Thus a D=7 core has 7≤U≤11. The complete finite bank contains no such core; this is an all-height proof after a proved bound, not an unbounded-search inference.

## Independent verification

- Full normal and optimized candidate runs passed and match each other and the frozen normal receipt byte-for-byte. At U≤20,000: 986,664 physical profiles; 11,394 slice profiles; 9,769 canonical and 1,625 exceptional; 1,636 distinct targets.
- `independent.py` uses positive-coordinate transitions, not candidate raw-state functions. At U≤1,000 it enumerates 11,721 source profiles and 191 slice profiles. Its entire slice agrees profile-for-profile with the candidate route. It tests actual precursor pairs and the D=7 core inventory.
- Matrix products are checked directly using separately entered integer matrices, rather than candidate state iteration.
- Five literal graph controls count paths by permutations and reserved covers by all subsets and permutations, rather than the candidate subset dynamic program. Singleton, edge, triangle, OZ, and ZJ precursors all give the predicted W profile, one root-start path, and one Hamilton cycle.
- Independent normal and optimized receipts also match byte-for-byte. Explicit exceptions enforce tests in both modes.

Artifacts are local audit outputs only. No candidate files or repository files were edited; no GitHub writes were performed. The parent owns integration and publication.
