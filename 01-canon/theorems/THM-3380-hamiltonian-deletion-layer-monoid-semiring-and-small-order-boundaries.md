---
id: THM-3380
title: "Hamiltonian deletion-layer monoid semiring and small-order boundaries"
status: >
  PROVED algebra + VERIFIED-EXACT witnesses + FINITE-EXACT complete
  small-order censuses.  The multiset of Hamiltonian counts in every vertex-
  deletion layer is a commutative multiplicative-monoid semiring compiler for
  ordered join; its augmentation is the shifted THM-3372 deletion transform.
  Equality of its fibres with (Gamma,Dham), observed through order seven,
  first fails at order eight in exactly three four-class fibres.  Its first
  skew-current collision is the structural self-converse palindrome ABBA
  versus BAAB.  The formal odd-cycle-count fugacity is determined by Dham
  through order eight, while at order nine it requires exactly the directed-
  triangle-factor coordinate b333; 17 of 91,780 Dham fibres split.  No
  all-order reconstruction, chronology, or completeness is asserted.
source: root/repository-archaeology-second-pass-2026-08-14
depends_on:
  - THM-002-ocf
  - THM-3324-tournament-deletion-response-gram-ordered-join-compiler
  - THM-3369-skew-deletion-response-and-ordered-join-orientation-current
  - THM-3372-multiaffine-deletion-transform-variance-and-skew-join-current
related:
  - THM-505-the-ocf-non-spectral-defect-H-equals-spectral-skeleton-plus-witt-defect
  - THM-506-permanental-companion
  - THM-3377-path-colour-deletion-compiler-and-skew-current
script: 04-computation/tournament_deletion_layer_fugacity_witness_thm3380.py
output: 05-knowledge/results/tournament_deletion_layer_fugacity_witness_thm3380.out
census_source: 04-computation/tournament_deletion_fugacity_n9_census_thm3380.c
census_output: 05-knowledge/results/tournament_deletion_fugacity_n9_census_thm3380.out
small_order_census_output: 05-knowledge/results/tournament_deletion_layer_census_n2_n8_thm3380.out
script_sha256: e3f318421c7bec2b20329099cdf5b0bda54a988fbd46778bfbcd153439a006c2
output_sha256: e7c5a0e069329c588053b92894e215965780639e1e3eb6f9f3e719d1ab615d25
census_source_sha256: 8ca2f4e6e98ec232796cf8813215b3e096d0f419227213c45fa281ae5e9f0644
census_output_sha256: 7a5a3605c6a0d8f47c9119564181f5ab603f9cf2eb0a13c244b8a13c78eda80f
small_order_census_output_sha256: 27de65f5a1f51b3d0405822829ed141da461c60018d48023dc007e7c0ebb24bb
semantic_sha256: f1e12cae71100526b0311e178306348dd4aeed3e356e595415d6cabe703087de
hash_basis: LF-normalized bytes
---

# THM-3380 -- deletion layers remember cards which moments forget

**PROVED algebra + VERIFIED-EXACT witnesses + FINITE-EXACT complete
small-order censuses.**

## 1. The Hamiltonian-card monoid semiring

For a tournament `T` on `n` vertices, write `H(empty)=1`.  Let
`N[N_(>0),times]` be the commutative monoid semiring with basis symbols `[h]`
and multiplication

```text
[a][b]=[ab].                                             (1)
```

Define the deletion-layer compiler

```text
M_T(t)=sum_(X subset V(T)) t^|X| [H(T-X)].              (2)
```

This is a histogram at every deletion size; it retains neither the labels of
the deleted vertices nor associations between different layers.  It is a
monoid semiring, not a group ring: Hamiltonian counts have no multiplicative
inverse in this coefficient system.

Every Hamiltonian path of an ordered join `X triangleright Y` consists of a
Hamiltonian path of `X` followed by one of `Y`.  Moreover deletion respects
the two factors.  Splitting a deleted set between them proves

```text
M_(X triangleright Y)(t)=M_X(t)M_Y(t).                  (3)
```

Let `epsilon([h])=h` be augmentation.  THM-3372's diagonal Hamiltonian
deletion transform `Dham_T(y)` satisfies

```text
epsilon(M_T(t))
 =sum_X H(T-X)t^|X|
 =Dham_T(1+t).                                          (4)
```

Thus `Dham` retains the sum of the card values in each layer, while `M`
retains their multiset.  Formula `(4)` is the exact controlled-forgetting map.

## 2. The first spectral/deletion separation is order eight

Use THM-3324's split spectral response `s=(P,zN)` and full marked-deletion
Gram `Gamma`.  Complete nauty representative censuses give

```text
order 7, 456 classes:
  M profiles=247,          (Gamma,Dham) profiles=247,
  the two equivalence relations are identical;

order 8, 6880 classes:
  M profiles=3328,         (Gamma,Dham) profiles=3325.  (5)
```

At order eight every `M` fibre still lies in a `(Gamma,Dham)` fibre, but
exactly three four-class `(Gamma,Dham)` fibres split into two `M` converse
pairs.  Hence the equality through order seven is a small-order coincidence,
first false at order eight.

The first split fibre has the following raw upper-triangle encodings.  The
identities of the witnesses are pinned independently by brute-force lexicographic
canonicalization; the generator row numbers are provenance only.

```text
row   raw bits
1049  1110100111011111011111111111
1133  1110100111011111111101111110
2620  1101000111111111011110111111
5234  1100011110100111111111101110.                      (6)
```

All four have

```text
P  =(1,8,56,224,704,1472,2176,1920,768),
zN =(0,8,56,352,1200,3200,5248,5248,2304),
Dham=(52,110,80,108,0,32,0,0,1).                       (7)
```

Stronger than Gram equality, their complete multisets of marked spectral
responses agree; the pinned serialization hashes are

```text
marked spectral deck  f707600d2804a502d7d91a2e4db6109f1aad26007f1679716194c451353b28c8,
Gamma                70c32ef0e0f4f5b53ae4f9b0dfa93a258cba051019fafd495877e52a59a22458.
                                                                    (8)
```

Nevertheless rows `1049/2620` have first-deletion Hamiltonian cards

```text
33,79,83,89,103,111,123,141,
```

whereas rows `1133/5234` have

```text
41,79,81,89,103,117,123,129.                            (9)
```

Both lists sum to `762`.  Their two-deletion histograms also differ while
both sum to `752`; all other `M` layers agree.  This localizes the missing
coordinate: an entire marked spectral deck need not determine the
Hamiltonian-card distribution, while augmentation `(4)` erases the
difference by construction.

## 3. A structural kernel for first skew currents

Let `A,B` be nonisomorphic self-converse strong tournaments and put

```text
T=A triangleright B triangleright B triangleright A,
U=B triangleright A triangleright A triangleright B.     (10)
```

Both words are self-converse: converse reverses the factor word and replaces
each factor by an isomorphic converse.  Unique ordered strong-component
factorization distinguishes the palindromes in `(10)`.  In contrast, `M`,
`Dham`, and THM-3324's `(s,Gamma)` interface are commutative under ordered
join and see only the common factor multiset `{A,A,B,B}`.  THM-3372's `xi`
and THM-3369's `Omega` are converse-odd, so both vanish on each self-converse
word.  Therefore even

```text
(M,Dham,s,Gamma,xi,Omega)                               (11)
```

cannot distinguish this family.

The least-order instance takes `A=K1`, `B=C3`.  Its two SCC-size words and
score sequences are

```text
(1,3,3,1), (0,2,2,2,5,5,5,7),
(3,1,1,3), (1,1,1,3,4,6,6,6).                          (12)
```

All `456` order-seven classes have distinct `(M,xi)` profiles.  Among all
`6880` order-eight classes, `(10)` with `K1,C3` is its unique collision.
Thus order eight is also the sharp first failure of this exteriorized
deletion-layer interface.

The mechanism suggests a nonabelian next response.  First exterior currents
behave like signed area and vanish on palindromes; a third-order iterated
marked-response tensor is a better hostile probe than another symmetric
moment.

## 4. Restore cycle count by a formal fugacity

Let `Omega(T)` now denote the odd-cycle conflict graph as in THM-002, not the
spectral current of §3.  For an independent cycle packing `S`, let `U(S)` be
its covered vertex union.  Define

```text
F_T(x,y)=sum_(S in Ind Omega(T)) x^|S| y^(n-|U(S)|).    (13)
```

Odd cycles cannot cross an ordered join, so packings split and

```text
F_(X triangleright Y)(x,y)=F_X(x,y)F_Y(x,y).            (14)
```

THM-3372 is the specialization

```text
Dham_T(y)=F_T(2,y).                                     (15)
```

For total covered size at most eight, the only partitions into odd parts at
least three are

```text
empty, 3, 5, 7, 3+3, 3+5.                              (16)
```

Consequently the covered size determines the number of packed cycles, and
`Dham` determines all of `F` for every tournament of order at most eight.

At order nine the sole ambiguity is full support:

```text
[y^0]F_T=a_9 x+b_333 x^3,
Dham_T(0)=2a_9+8b_333,                                 (17)
```

where `a_9` counts directed Hamiltonian cycles modulo rotation and `b_333`
counts unordered factors into three directed triangles.  All other covered
sizes remain governed by `(16)`.  Hence `(Dham,b_333)` determines `F`
exactly at order nine.

Two nonisomorphic order-nine witnesses have the common transform

```text
Dham=(394,656,664,304,262,0,50,0,0,1),                 (18)
```

but full-support fugacities

```text
161x+9x^3,                 157x+10x^3.                 (19)
```

Their complete polynomials agree away from `y^0`.  The trade of four
Hamiltonian 9-cycles for one triangle factor is invisible at `x=2`.  Anchored
Hamiltonian-cycle DP and direct odd-cycle-packing enumeration independently
give the two pairs in `(19)`.

The structured positive control `C3[C3]` has

```text
[y^0]F=207x+37x^3,       Dham(0)=710.                   (20)
```

Its `37` triangle factors are one internal factor plus the `36=(3!)^2`
transversal factors.

## 5. Exact census and information contract

The complete order-nine census has

```text
191536 tournament isomorphism classes,
 91780 distinct Dham profiles,
 91797 distinct (Dham,b_333) profiles,
    17 Dham fibres split by b_333.                       (21)
```

The typed contract is

```text
source:    all vertex-deletion Hamiltonian cards; odd-cycle packings
maps:      cards -> M -> augmentation Dham;
           packings -> F -> evaluation x=2
preserved: deletion-size histograms, ordered-join multiplication,
           covered size; cycle count only before x=2
destroyed: vertex labels, cross-layer card association, SCC order,
           palindrome order, cycle count at the first 9 versus 3+3+3 wall
sidecars:  Hamiltonian-card histograms for the spectral fibre;
           b_333 for order-nine fugacity
tests:     the four classes (6), ABBA/BAAB, and the pair (18)--(19). (22)
```

No all-order completeness, isomorphism reconstruction, chronology, owner,
Hamiltonian optimization, LRC, or dynamical conclusion is claimed.

## 6. Verification

The generator-free standard-library companion:

- independently brute-canonicalizes all frozen witnesses;
- recomputes `M,Dham,xi,s,Gamma,Omega` without importing another companion;
- checks the full marked-spectral-deck hostile `(6)--(9)`;
- proves the concrete `K1/C3` palindrome factorization and both product laws;
- directly enumerates every odd-cycle packing of the order-nine witnesses;
- independently computes `a_9` and `b_333`; and
- gives the positive `C3[C3]` control.

Its optional exact census mode consumes pinned default `gentourng` files for
orders two through eight (order one is trivial).  The two largest universe
counts and SHA-256 values are

```text
n=7: 456 rows,  164260b94960af0cc63faf3f178ceb95f4dd23bbca376ec23872c33b30d94261,
n=8: 6880 rows, fc96c6997724e54ccea3bd166f4117d9e27925d85f568e31b0623527e5139dad.
                                                                    (23)
```

The C census consumes the pinned `191536`-row default order-nine universe
with SHA-256

```text
4f7d6c43cfed87e1e5293dc751736efe2d7bc1554946cdc83f4026a575fbbbf8. (24)
```

It checks `(17)` class by class and produces `(21)`.  Builds with `-O3`,
`-O0`, and undefined-behaviour sanitization give the identical stored
transcript.  The inputs were generated by Debian nauty `2.7r3+ds-1` using
`nauty-gentourng 7`, `8`, and `9`; row indices are never called canonical.

Reproduce the generator-free audit with

```text
python3 04-computation/tournament_deletion_layer_fugacity_witness_thm3380.py
python3 -O 04-computation/tournament_deletion_layer_fugacity_witness_thm3380.py
```

Given the pinned representative files, reproduce the complete small-order
censuses with

```text
python3 04-computation/tournament_deletion_layer_fugacity_witness_thm3380.py \
  --workers 4 --census gentourng2.txt --census gentourng3.txt \
  --census gentourng4.txt --census gentourng5.txt --census gentourng6.txt \
  --census gentourng7.txt --census gentourng8.txt
gcc -std=c11 -Wall -Wextra -Werror -O3 \
  04-computation/tournament_deletion_fugacity_n9_census_thm3380.c -o /tmp/thm3380_n9
/tmp/thm3380_n9 gentourng9.txt
```

Ordinary/optimized Python and all three C builds byte-match their stored
outputs.  Artifact and semantic hashes are pinned in the frontmatter.

**QED.**
