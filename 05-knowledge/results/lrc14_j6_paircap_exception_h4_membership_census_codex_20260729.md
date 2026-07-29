# Pair-cap-exception H4 membership census

**Status: FINITE-EXACT SCOPED WORKLOAD CENSUS.**

On the `52` exact-pair-cap exceptions isolated by THM-2901, reconstruct
the literal carrier `C`, its mass `h`, and the attained allowed singleton
maximum `q_1`.  THM-2899 gives `q_1<3h/7`, so every hypothetical
five-cover contains at least two labels in

```text
H_4(C)={w allowed:c_C(w)>=(h-q_1)/4}.
```

The THM-735 discrepancy bound seals each core globally.  Exact membership,
not cutoff-universe cardinality, gives

```text
exception branches                         52
distinct bodies                            52
rank split                               51/1
actual H4 size range                    12..44
sum of actual H4 sizes                    1348
unordered actual H4 pair flags           18290
maximum analytic membership cutoff        2026.
```

This is the finite input universe for the `(k,s,ell)=(5,4,2)` heavy-link
child attack.  It does not resolve any pair flag, close any exception
branch, close a whole root, or prove LRC(14).

Ordinary and optimized outputs are byte-identical.  The semantic row
digest is

```text
f75c1440a58c8724cb6bcebd41343f81daa6d3bca302a6122ba3721373479cb6.
```

Artifacts:

```text
04-computation/lrc14_j6_paircap_exception_h4_membership_census_codex_20260729.py
SHA-256 63a80908a6380a877345f0cc4aba7a5e0ef2bb3d59b1b10d58367444ed406b75

05-knowledge/results/lrc14_j6_paircap_exception_h4_membership_census_codex_20260729.out
SHA-256 efe7e1138f4b4708ed47e526405b13b776d602c5e26b7319db06d7594036417e
```
