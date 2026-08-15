---
id: THM-3394
title: "Twelve formerly missing Hadamard orders through 2000"
status: >
  PROVED + VERIFIED-EXACT + INERT RECONSTRUCTION + HOSTILE-AUDITED
  + INDEPENDENTLY AUDITED.
  Explicit sign matrices of orders
  668,716,892,1132,1244,1388,1436,1676,1772,1916,1948,1964
  are reconstructed from one checked-in 23,828-sign word and one frozen
  192-bit schedule.  A standard-library-only verifier checks every shape,
  sign, immutable matrix hash, and all 12,913,704 unordered pairs of
  distinct rows; every distance is exactly half the relevant order.
  This proves only these twelve finite existence statements and does not
  prove the Hadamard conjecture.
source: codex-2026-08-14 inert reconstruction of the supplied sign-puzzle bank
audit: literal fixed-template renderer, exact periodic/OA sidecars, normal/-O agreement, one-entry hostile, and an independent reconstruction/contract/Gram audit
depends_on: []
related:
  - THM-3393-hadamard-order-668-explicit-certificate
  - THM-2833-legendre-333-asymmetric-multiplier-exclusion
  - THM-3396-four-bit-pairwise-independent-fourier-cone
script: 04-computation/hadamard_twelve_order_bank_thm3394.py
data: 04-computation/hadamard_twelve_order_signword_thm3394.b85
output: 05-knowledge/results/hadamard_twelve_order_bank_thm3394.out
script_sha256: 7ae931b3cf268550287bd0621b9b85b8ea167126fadfb90d57b5106d0f82fb2d
data_file_sha256: 68f7ceebb67005bf1b968171f7e6897cc33bde68adbd63f14bd45edfeb7b3f06
output_sha256: d8efee90947015a7e6fc28a1685cc3d378357a85e1d4814953b32b17c5cd76a9
hash_basis: working-tree bytes (LF)
---

# THM-3394 -- twelve explicit Hadamard matrices

**PROVED + VERIFIED-EXACT + INERT RECONSTRUCTION + HOSTILE-AUDITED
+ INDEPENDENTLY AUDITED.**

For every

~~~text
n in {668,716,892,1132,1244,1388,1436,1676,1772,1916,1948,1964}
~~~

there is a matrix H_n in {+1,-1}^{n by n} satisfying

~~~text
H_n H_n^T = n I_n.                                           (1)
~~~

The proof is an explicit finite certificate bank.  It is deliberately
generator-independent at the final gate: after rendering each matrix, the
verifier checks all unordered pairs of distinct rows as bitsets.  If b_i is
the bitset of row i, with one bit for one of the two signs, then

~~~text
<H_i,H_j> = n - 2 popcount(b_i xor b_j).                     (2)
~~~

Every one of the n(n-1)/2 distances is n/2, and each diagonal norm is n.
Thus (2) proves (1) without floating point, a solver, or an appeal to the
construction templates.

## 1. Immutable input and exact schedule

The checked-in data file contains only base85 text.  Removing ASCII
whitespace gives 4,360 bytes with SHA-256

~~~text
7831fb93c91b82e70827345c2b72eee6702a30ffce3f5162b01107e5ce35fbc0.
~~~

Base85 decoding gives 3,488 zlib bytes with SHA-256

~~~text
1756297611d2bb403e9c4152ea91146428482983716aad165a3cc21396d5a61c.
~~~

Decompression gives exactly 23,828 ASCII signs, 12,188 plus and 11,640
minus, with SHA-256

~~~text
5b5fe8fa42f0d6a8b4e4c9926726d82a6aab8e1070c1ae4d1b430c1277e58db4.
                                                                    (3)
~~~

The renderer also embeds exactly these 192 schedule bits:

~~~text
111000101001100011100010110010001100001101111100
110001000110110011000100110111001000010101011000
011001011000000011100110100010001110011011101000
101001110111000011000111100111001100011110101100
~~~

Their SHA-256 is

~~~text
c18071311a78311819a20a1ba7823de6581bb907fd9b4b44a3a9b17860a370ff.
                                                                    (4)
~~~

They split into twelve 16-bit records.  The first three bits encode

~~~text
111 -> h, 110 -> g, 101 -> w, 100 -> v, 011 -> u,
~~~

and the last thirteen bits encode the unsigned primary length.  Decoding is
required to give

~~~text
(h,664), (h,712), (g,892), (g,1132), (g,1244), (v,1368),
(u,1408), (h,1672), (h,1768), (w,1904), (g,1948), (g,1964). (5)
~~~

The exact half-open consumption ledger is:

| # | code | primary slice | additional side slices | output | consumed |
|---:|:---:|:---|:---|---:|---:|
| 1 | h | [0,664) | none | 668 | 664 |
| 2 | h | [664,1376) | none | 716 | 1,376 |
| 3 | g | [1376,2268) | none | 892 | 2,268 |
| 4 | g | [2268,3400) | none | 1,132 | 3,400 |
| 5 | g | [3400,4644) | none | 1,244 | 4,644 |
| 6 | v | [4644,6012) | [6012,6892), [6892,7012), [7012,7132), [7132,7252), [7252,7372) | 1,388 | 7,372 |
| 7 | u | [7372,8780) | [8780,10460), [10460,11356), [11356,12252), [12252,13148), [13148,14044) | 1,436 | 14,044 |
| 8 | h | [14044,15716) | none | 1,676 | 15,716 |
| 9 | h | [15716,17484) | none | 1,772 | 17,484 |
| 10 | w | [17484,19388) | [19388,19724), [19724,19772), [19772,19820), [19820,19868), [19868,19916) | 1,916 | 19,916 |
| 11 | g | [19916,21864) | none | 1,948 | 21,864 |
| 12 | g | [21864,23828) | none | 1,964 | 23,828 |

The final cursor is required to equal both 23,828 and the decoded signword
length.  No sign is unused and no sign is read twice.

## 2. Construction contracts

The five codes are exact construction contracts, not claims that the
recorded low-order invariants alone characterize all valid inputs.

### The h contract

For an h record, split the primary word into four sign sequences of common
even period q.  Their row sums are

~~~text
(2,0,0,0),
~~~

and the sum of their four periodic autocorrelations is -4 at every nonzero
shift.  Equivalently, their negative supports form an SDS with parameters

~~~text
4-(q; q/2-1,q/2,q/2,q/2; q-2).                              (6)
~~~

The fixed renderer places the corresponding four-circulant core in the
same four-row/four-column bordered Goethals--Seidel template audited
independently for q=166 in THM-3393.  Its output order is 4(q+1).  The bank
uses q=166,178,418,442, yielding orders 668,716,1676,1772.

### The g contract

For a g record, split the primary word into four sequences of period q.
The sum of their periodic autocorrelations is zero at every nonzero shift.
The standard unbordered four-circulant Goethals--Seidel array therefore has
Gram matrix 4q I and output order 4q.  The exact row-sum sidecars are

| output | q | four row sums |
|---:|---:|:---|
| 892 | 223 | (11,11,11,23) |
| 1,132 | 283 | (19,19,19,7) |
| 1,244 | 311 | (21,19,21,1) |
| 1,948 | 487 | (17,1,17,37) |
| 1,964 | 491 | (29,27,15,13) |

In every row the sum of the four squared row sums is 4q, as required by
the zero-frequency identity.

### The w, v, and u contracts

These three composite records use their primary word, one unequal top-side
piece, four equal-length sign pieces, and a fixed literal block-rewrite
template:

| code | primary | top side | four equal sides | output increment |
|:---:|---:|---:|:---|---:|
| w | 1,904 | 336 | 4 by 48 | 12 |
| v | 1,368 | 880 | 4 by 120 | 20 |
| u | 1,408 | 1,680 | 4 by 896 | 28 |

Regarded as four binary factors across their equal side coordinates, each
four-row sidecar is balanced and pairwise orthogonal: it is an
OA(N,4,2,2).  In the triple order (012),(013),(023),(123), its exact
degree-three Walsh sums and degree-four Walsh sum are

| code | N | degree-three sums | degree-four sum |
|:---:|---:|:---|---:|
| w | 48 | (0,0,-8,-16) | -8 |
| v | 120 | (8,-8,0,0) | 24 |
| u | 896 | (-224,0,0,-224) | 448 |

Together with the zero first and second moments, these values are the full
four-bit Walsh transform and hence determine every sign-pattern
multiplicity by Fourier inversion; see THM-3396 for the general cone.
They are not asserted to be sufficient for the composite lift: the
unequal top-side piece and the complete fixed template remain load-bearing.
The proof of the three resulting Hadamard matrices is the full render
followed by the all-pairs certificate, not an inference from OA strength
two alone.

## 3. Exact matrix certificates

For each raw renderer output, the text hash is taken on n lines of n signs,
each terminated by LF.  Normalization multiplies rows and columns by signs
so that the first row and first column are all plus; its hash is computed in
the same text format.

| n | raw text SHA-256 | normalized text SHA-256 | pairs | distance |
|---:|:---|:---|---:|---:|
| 668 | bdeb5059d77e2703211082627b60441b8c888c928a55cc6f295e011941a387b0 | 73f1de1539849e1dc7e6085cc69c563fd2965c44970263e8203384bd1a46aa63 | 222,778 | 334 |
| 716 | 3adcb1bb2884467d9e34069a3b32950728adabcdb8b35a4503d20c3312664ee6 | 79ebae74dbadae11059aeabf77b61144ccf6fd9905c2bf8167ac98854651366c | 255,970 | 358 |
| 892 | e77fc79ab287f5f5ba5bbdc10191bdc7593839052fe1015c1fb6a2e974ab54de | 73d28c9e58d4c3bac1e2bc0cc0d6f0b504eb7719fa2370e9403692e4dfba7f6c | 397,386 | 446 |
| 1,132 | 7d1c1e892149e90330d58bb0cf9ef2c888078df1b35fb55f8724d580ebf7b743 | f0499c89d761d57d330471a1446a0f099723f663b5602285aa9dcc1070b0216b | 640,146 | 566 |
| 1,244 | 4cb747cf511eba1f203582b5121bdf6ab02671133e45579c1d023add8b2da143 | cffd1c77f09646b9b17e4c46ed9f6bf091bd56b2c2cc06c83c64054b44d74d33 | 773,146 | 622 |
| 1,388 | a6b92584eb803b87026709d64fe892dec8f7182a120e13de9edd3065cf05bf0b | 5ebac7a33165a39cb9e86e5b99eb42cc46e3c7006956b2c97885f8f6b4655db1 | 962,578 | 694 |
| 1,436 | e4d745a4d44f39a5671f9cd86f5c1d0aef93504dcfb2e253451cadc9e4086728 | 99984bf728bcf37a360ba0af57076a729653070b026271085f96dce35703a00e | 1,030,330 | 718 |
| 1,676 | 8e919c2bdb4d30c34817eb5650d2dd3d82d7c6504feccd96c5ca22a2191cdb99 | 787f9ec132cfa6900090e78bbebd8c1132af583793b4c65494ca52f7000a1e68 | 1,403,650 | 838 |
| 1,772 | 1852e951db69c44eb95b37ed741c3ff2e29691267eaf872d6a9da3a977236ba2 | 020f15819f78aab232da7ec9183caf172cc6ab1cdad191f4fabb9c63ec308af1 | 1,569,106 | 886 |
| 1,916 | be2073eeaa5399cfe104023829d2c6770b49dd2f07bf6347203f1cbd75577ae9 | 081197bb0c94e9e1aa628fa3170db50d7df45e2889494bf12ee5abcdfd4266da | 1,834,570 | 958 |
| 1,948 | fddc841ebf951f6e17e939551d058ea5df046251ea065b5f6e7ee2fd8d0f62ce | ee885b2ade31910d36fd7bc11321388c20aa64eaef081fecacf05f3b53d5d792 | 1,896,378 | 974 |
| 1,964 | 740b907cd442f1b7fd40dcc31f2b3aae9794842da6dc579f98dac1d0d9e1493d | 4a3e301084d016f6b11e32d425b9f8f54e79cd2b52ec8c6e64420635c13aaa53 | 1,927,666 | 982 |

The pair counts total 12,913,704.  The minimum and maximum Hamming distance
at every order are both the displayed distance, so this is an all-pairs
Gram verification, not sampling.  The twelve matrices contain 25,844,160
entries in total.

The normalized order-668 hash agrees exactly with the independently packed
matrix of THM-3393.  This is a second reconstruction route for that object.
The remaining normalized hashes freeze the other eleven representatives;
the theorem makes no uniqueness or equivalence-class claim.

## 4. Inertness, hostile control, and scope

The source puzzle ended in text resembling a shell decoder.  That text was
treated as untrusted data and was never executed.  The checked-in verifier
does not contain or invoke it.  In particular:

- the data file decodes only to the plus/minus word in (3);
- the Python program has no shell, subprocess, eval, dynamic import,
  network, or output-file-writing path;
- the small sed-like rewrite engine has fixed opcodes, parser, and templates
  embedded in the audited Python source; checked plus/minus data populate
  literal fields in the composite templates but cannot inject a delimiter,
  opcode, branch label, or Python text;
- the decoded alphabet, lengths, hashes, schedule, and cursor are checked
  before or during reconstruction; and
- the rewrite engine has an explicit finite step bound.

A one-entry flip in the order-668 matrix changes the first detected pair
distance from 334 to 333, and the verifier records that this hostile is
rejected.  Ordinary and optimized Python runs are byte-identical; the frozen
output records the exact-consumption, all-twelve Gram verdict, and finite
scope lines.

This theorem closes exactly the twelve displayed finite orders in the
repository's former through-2000 list.  It does **not** prove that a
Hadamard matrix exists at every positive multiple of four, does not
classify Hadamard matrices at any order, and does not establish a priority
claim.  In particular, the order-668 h construction bypasses but does not
solve the open length-333 Legendre-pair question of THM-2833.

Reproduce from the repository root with

~~~bash
python3 04-computation/hadamard_twelve_order_bank_thm3394.py
python3 -O 04-computation/hadamard_twelve_order_bank_thm3394.py
~~~
