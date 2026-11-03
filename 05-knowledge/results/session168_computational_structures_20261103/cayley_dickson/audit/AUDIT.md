# Independent audit of S168 Cayley–Dickson transport

Audited `/tmp/s168_cd/NOTE.md` and `check.py`. No substantive mathematical defect found.

## Verified

- The recurrence, square identity, zero-translation defect, kernel equality with even bottom-half truth tables plus the top bit, and the all-level order proof are valid for the stated range n≥3. The associative n=2 exception is correctly separated.
- The degree-(n−1) homogeneous part is nonzero for every nonzero multiplier. The degree obstruction is correctly restricted to sign-target gates and unconditional index translations.
- Literal binary tournament encoding, ordinal profile multiplication, SCC-based induced cube, scheduled gate implementation, and fixed-labelled single-flip triangle obstruction are sound. The isomorphism-class clique barrier is necessary, not sufficient, as stated.
- The optional strong wrapper and relabelled seed 185 agree with the exact S167 ternary note at commit 596ae1ea47f640689304af59237e5f737c75c0b8 (read by targeted git show).

## Executed controls

`python3 /tmp/s168_cd/check.py` and `python3 -O /tmp/s168_cd/check.py` both passed and produced byte-identical results.json. They checked 1,364 basis products and 248 physical short-macro flips, with transport bits 7,12,21,38,71.

The supplemental `independent.py` in this audit directory passed normally and under -O, again with identical JSON. It:

- Enumerated groups from dense algebra-generated permutations without using sign() or full_permutation_group(): orders 4,8,128,4096 for n=1,2,3,4.
- Compared the entire actual zero-translation kernels at n=3,4 against all functions in the claimed kernel, not merely their dimensions.
- Recomputed ANF coefficients by direct subset summation and ran every full physical ANF macro for n=2,...,7: 43,176 macros and 539,424 literal single-arc updates.
- Classified all 264 single-arc reversals from the four two-block representatives via independently computed SCCs and complete six-vertex seed relabelling orbits, confirming the induced class cube.
- Checked the conversion from the six-vertex path partition coefficients to F(r)=r^6+12r^4+16r^2 at 16 integer inputs.

## Small editorial / attribution improvements

1. Section 6 calls the independent multiplication implementation “sparse-vector,” but multiply() uses dense integer tuples. Change this to “dense integer-vector.”
2. Add Ren and Zhao, *The twisted group algebra structure of the Cayley–Dickson algebra*, https://arxiv.org/pdf/2205.07986, Proposition 1, for the previously known twist recurrence. Their formula (3.3) is the same recurrence after substituting the commutation function; they explicitly credit Albuquerque and Majid for it. This keeps the explicit-coordinate presentation from appearing to claim a new general twist law. The sharp fixed-multiplier ANF-degree theorem was not found in this targeted literature check, but the check does not establish novelty.
3. Kirshtein, https://arxiv.org/pdf/1207.4230, uses the exact same doubling convention and Theorems 34–35 support the known one-sided group result. The printed Theorem 35 says n≥2, although its displayed elementary-abelian formula does not apply at the quaternion endpoint. NOTE correctly states n≥3 and independently handles n=2.

## Source hashes

- NOTE.md SHA256: 8d6fddfbcba42bd5bb9f45b347f6f59e1fa2169cb631aa2244627fd66252b421
- check.py SHA256: efb85d76d15875e64d12af1c4398c12648dc53c18b4a2acc043052756c1f7bc3
- results.json SHA256: cd8bb877fc901cdd0638e9cc09e7d95065a06674d95c0234cb88f76c9aa3b446
