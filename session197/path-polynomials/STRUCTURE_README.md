# Two exact structural bridges

These packets are separately frozen and independently audited. They do not assume that arbitrary positive algebraic states have graph realizations.

- `MINIMUM_COVER_POLYNOMIAL.md` proves p_nu(x)=x(x+1)...(x+nu-1)L_C(x)R_C(x) for every minimum root crossing-edge cover. The factors enumerate actual labelled prefix/suffix arrangements and residual bin profiles. This strengthens the scalar factorial-product identity and recovers the earlier matching lower bound.
- `ENDPOINT_RECIPROCITY.md` proves that the alternating sum over paths ending at a equals minus the number of Hamilton paths starting at a. A full squarefree matrix rank-one identity supplies the proof, with the singleton and diagonal-removal boundaries explicit.

Run `python3 audit_structure.py` to replay both verifiers under normal and optimized Python and compare frozen receipts byte-for-byte.

Independent reports:

- `session197-bridges/minimum-cover-audit/REPORT.md`, SHA256 `5267ec64e06b78cc5798219601ac46c23b2fe442ad7451d8ef3f7301ee74812d`.
- `session197-bridges/endpoint-reciprocity-audit/REPORT.md`, SHA256 `c1c4837b6f8fb41d713f148ba040af83f5d1283933f29fc684826db129a858bd`.

Neither packet claims the full tournament spectrum or a new universal unit-clock congruence. The physical primitive-cancellation construction and the binary deletion application are separate packets.
