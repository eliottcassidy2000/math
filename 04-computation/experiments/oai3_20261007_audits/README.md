# Independent audits of mac-mini-2026-10-07-oaimath3 (2026-10-07)

Three adversarial auditors each wrote their own code. Their reports are kept verbatim; the corrections are applied in canon and recorded in MISTAKE-583.

| Report | Scope | Scripts |
|---|---|---|
| `A_coalescence.md` | THM-4581 and its corollaries (HYP-9213/9214/9217/9220, THM-4569 (7)) | `A/` (a1-a12; `a5_chain_sim.c` needs GMP: `cc -O2 a5_chain_sim.c -lgmp`) |
| `B_numbertheory.md` | THM-4566, THM-4567, and the THM-4555, THM-4512 and THM-4558 updates | `B/` (`exp2_audit.c` discriminant sieve to 2.1e11; SAT and PARI checks) |
| `C_ninefourths_tclock_zeroless.md` | THM-4568, THM-4569, THM-4580, HYP-9219 | `C/` (a1 = THM-4568, a2 = THM-4569, a3 = THM-4580) |

The reports cite their scripts by bare filename, as they ran in the session scratch directory; the same names are in the subfolders here. OEIS b-files, the #107 TeX and compiled binaries are not committed.
