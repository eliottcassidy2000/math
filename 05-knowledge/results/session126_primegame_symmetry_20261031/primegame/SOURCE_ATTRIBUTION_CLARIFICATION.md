# Exact source attribution: two published fraction lists

**The user's ordered list matches Conway's 1987 Theorem 1 exactly.** This does not mean it is the only published PRIMEGAME variant or the earliest historical version.

## Conway 1987: the rack audited here

Primary source: J. H. Conway, *FRACTRAN: A Simple Universal Programming Language for Arithmetic*, 1987, printed **p.4**, Theorem 1 (**PDF page1**). The labeled A–N fraction table is repeated on printed **p.9** (**PDF page6**).

- Numerators, in order: `(17,78,19,23,29,77,95,77,1,11,13,15,1,55)`.
- Denominators, in order: `(91,85,51,38,33,29,23,19,17,13,11,2,7,1)`.
- Hence the final three fractions are **15/2,1/7,55/1**.

The theorem is for starting integer 2. The proof is on printed pp.9–10. This page/table sequence was re-opened and directly checked against the primary PDF's extracted text on 2026-10-31.

[Primary chapter scan](https://gwern.net/doc/cs/computable/1987-conway.pdf) · [Publisher bibliographic record, DOI10.1007/978-1-4612-4808-8_2](https://link.springer.com/chapter/10.1007/978-1-4612-4808-8_2)

## A separately published variant

Jörg Endrullis, Clemens Grabmayer, and Dimitri Hendriks, *Complexity of Fractran and Productivity*, arXiv:0903.4366v3 (31 July 2009), printed/PDF **p.2**, gives the same first eleven fractions but final three **15/14,15/2,55/1**. The following page describes its prime-exponent output property.

[Primary 2009 paper](https://arxiv.org/pdf/0903.4366)

Avoid using “canonical” to silently select one of these two different lists. They are not equal one-step maps: at input 14, the 1987/user rack chooses 15/2 and gives 105, whereas the alternate rack chooses 15/14 and gives 15.

## Audit scope

Every script, divisor-loop proof, seed outcome, timing result, FUS comparison and CA implementation in this packet uses the **1987/user list**. Its all-time prime-exponent theorem is justified by the supplied direct divisor-loop proof for that exact list; it is not imported from the alternate program. No numerical result or arbitrary-seed classification is asserted to transfer unchanged between variants.
