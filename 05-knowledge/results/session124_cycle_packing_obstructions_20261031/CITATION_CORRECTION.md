# Correction lineage: OCF source numbering

The original cycle-realization proof and the inherited repository theorem
`01-canon/theorems/THM-002-ocf.md` identify the Irving–Omar Hamiltonian-path
formula as **Corollary 20**. The linked v1 actually labels it **Corollary 19**:

- [Irving–Omar v1 PDF, page 13](https://arxiv.org/pdf/2412.10572v1).
- [Irving–Omar v1 HTML, section 3.4](https://arxiv.org/html/2412.10572v1#S3.SS4).

The formula, its use of literal directed cycles modulo rotation, and its
attribution to Grinberg–Stanley are unchanged. This is a bibliographic repair,
not a retraction or alteration of the OCF theorem. No mathematical conclusion
or numerical receipt changes.

The original `cycle_realization/RESULTS.md` has SHA-256

```
f7ef4ba4b918f2dd6a3cb0c880403eca3975039a34ce4c464dff7994f07a5ac9
```

The independent audit reported the error without changing that source packet.
This curated publication copy replaces `Corollary20,` by `Corollary19,` in
`RESULTS.md` and records the original identities in `PROVENANCE.json`. Its
only other proof-file change is the illustrative command directory. The
inherited canon file is outside this package and has not been edited here.

Corrected source claim: the sum over permutations with nontrivial directed odd
cycles is equivalent to the weighted sum over disjoint literal odd-cycle
collections. Tournament path reversal identifies the complement-path count
with the original-path count.
