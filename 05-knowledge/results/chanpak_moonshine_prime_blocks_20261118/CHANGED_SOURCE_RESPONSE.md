# Change the physical source, not merely the response coefficient

**PROVED source obstruction, exact general-response decoder, and replacement-source family.** None implies universal placement.

## An obstruction surviving every mask and every response A

For the original R profile (9,1,25,5), every retained-root output is

    H(W)=9+25W-5K(W)-(W mod2),

hence is always 3 or4 modulo5. This excludes the other three residue classes no matter how one changes A, the edit mask, word length, center knobs or decoder. Prepared descendants do not fix the problem: they are still words from the same source. Therefore response flexibility alone cannot make any recontact bank using this R source universally cover counts. This is stronger than a sparse-bank estimate and points to a specific missing physical resource.

## Actual replacement sources

The order-five outgoing masks

    L=(16,9,3,5,14), root1: (h,e,t,d)=(15,3,19,5);
    R=(24,17,3,6,12), root0: (h,e,t,d)=(15,3,21,5)

are independently checked by literal subset path counting. Their deletion count and weighted-start difference remain d=5,D=2, preserving every unit-response edit identity. But neither weighted start vanishes modulo5.

For the recursively fixed unit mask at r=4j+2, j>=1, put M=16^j, B=4M, u=(7M-2)/5. The addresses B+4BV and u+4BV have exact output difference

    C'_j-2V,
    C'_j=(247M+28-10j)/5.

Thus V_j=(247M+28-10j)/10 is a positive integer and gives a physical recontact at every j. Its common height is

    H_j=15+19B+(76B-10)V_j-5K(V_j).

It satisfies the all-order residue law

    H_j = 1-2j (mod5).

Therefore each of the five residues occurs infinitely often in this concrete recontact family. To prove the residue law, 16^j=26-10j (mod50) for j>=1 by induction; this gives V_j=2j (mod5). Also B=4 (mod5), so the displayed height reduces to 1+4V_j=1-2j. This is a physically proved change in local coverage, not an assumption that source profiles can be chosen freely.

## Stop privileging A=1: exact decoder for arbitrary response

For any positive odd A with A=1 (mod5), choose r=2 (mod4) large enough that s=(4B-A)/5 satisfies 0<s<2B. Take low words

    a=rev_(r+1)(x), b=rev_(r+1)(x+s), 0<=x<2B-s,

and optionally add the same bit at r+1 to each. With the replacement sources define

    C=19a-21b-5[K(a)-K(b)]-3(a mod2)+3(b mod2).

For their shared high address V the signed count difference is exactly C-2AV. Given a requested even gap g, accept precisely when

    V=(C-g)/(2A)

is a nonnegative integer. Then evaluate either actual output to test a prescribed absolute center. This is an exact physical simultaneous center/gap candidate test, not a theorem that a candidate exists for every center.

The identity follows from cancellation of high-high checksum terms; the remaining cross coefficient is 4B+5(rev(a)-rev(b))=A. Unit response removes divisibility restrictions on V. Nonunit responses reinstate those restrictions but change the absolute scale and can produce centers not in the unit-response bank. Thus they should be selected from the target, rather than either discarded or presumed universally helpful.

Replay: `python3 finite_edit_recontact/changed_source_replay.py`. It independently recounts the replacement source profiles, literal small returned graphs, the all-order family identities through j=10, and finite exact center banks for responses A=1,11,21,31,41 at r=6. These latter bank sizes and novelty counts are finite results only.
