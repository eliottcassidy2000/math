# Hostile referee report: distinguished-coloop round budget

Status: **ACCEPTED after corrections** (2026-07-28).

Audited candidate:

- `.scratch/ak_vandermonde_20260728/ROUND-RANK-THEOREM-CANDIDATE.md`
- canonical replay
  `04-computation/ak_distinguished_coloop_round_rank_audit.py`

The pin/component argument, coloop deletion identities, bicirculation/bridge
bound, two-round budget, and matroid-union topology theorem are valid over
`Q` under the now-explicit paid-row-multiset hypothesis.

Corrections incorporated during review:

1. independent-row two-forest sufficiency was separated from tied-label
   slope grammars;
2. the restricted-minor criterion and the exact
   `delta + gamma + kappa` topology/grammar/specialization decomposition were
   added;
3. the unrestricted generic-defect claim was supplied with the
   subset-minor argument and matroid-union rank formula;
4. paid strict/mode-three multiplicity and the failure of loose cost
   accounting were made explicit;
5. the complete slack-refined two-round window and the omitted `(g,u)=(5,3)`
   case were added;
6. rational-span, positivity, indexing, and firing/coloop quantifiers were
   made explicit.

Independent computation:

```bash
python3 .scratch/ak_rank_referee_20260728/independent_sympy_rank_defect_audit.py
python3 -O .scratch/ak_rank_referee_20260728/independent_sympy_rank_defect_audit.py
```

This uses SymPy ranks, rank-deletion coloop tests, and exhaustive enumeration
of every row subset for `delta`.  It agrees with the custom RREF on every
round and independently obtains `delta=sigma` throughout all four frozen
records.

Normal and optimized runs are byte-identical to
`independent_sympy_rank_defect_audit.out`.

```text
script_sha256: eb73fe86fc8fc3fca579c9f96db62d3f5dfbc17e7c1156c1199518faf2a6800c
output_sha256: eed5875acd62998e467052916012a24738c0af9a8197fb37b2fdd6941b2934e1
hash_basis: LF-normalized bytes
```

Administrative promotion conditions only:

- reserve a fresh theorem ID after fetching current `origin/main`
  (`THM-2841` was already occupied; `THM-2847` was free at fetch
  `f94d5aa794`, but must be rechecked);
- replace the reserved/proof-candidate status by the proved/finite-exact
  status;
- recompute hashes if the canonical replay is hardened or otherwise edited.
