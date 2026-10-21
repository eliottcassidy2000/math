# Order 35: a complete cyclic reduction and seven black ports per source

October 21, 2026. **PROVED**, with the explicitly **CITED** Feit–Thompson
solvability input and a **FINITE-EXACT** complete cyclic certificate bank:

\[
[z^{35}]R_6(z)=0,\qquad \deg R_6\le28. \tag{1}
\]

Together with this folder's [order-25 exclusion](README.md), this packet alone
leaves \(12,14,15,18,20,21,24,28\). Separate complete certificates at 12 or
14 may be combined later; they are not dependencies of this proof. The
inherited eventual recurrence order remains 18, and the uniform safe onset
improves to **41**: its numerator has degree at most \(28+12=40\).
No earliest-onset or full remaining-center census is claimed.

## 1. Source and imported group theorem

The exact source is a strong self-converse tournament with \(\Delta=6\).
The inherited [degree-35 parameter table](../deficit_six_centers_20261021b/COMBINED_CUTOFF.md)
and the [full-orbit run lemma in this folder](README.md) imply that an
order-35 candidate is vertex-transitive. This implication retains an actual
Hamilton cycle with exact type \((q,m)=(5,7)\) or \((7,5)\); it is not
regularity alone.

Every tournament automorphism group has odd order. We now additionally use
the following **CITED** published theorem, not an elementary lemma proved by
this experiment: every finite group of odd order is solvable. The primary
source is Feit–Thompson, *Solvability of groups of odd order*, Pacific Journal
of Mathematics 13(3), 1963, printed page 775,
[original article](https://msp.org/pjm/1963/13-3/pjm-v13-n3-p01-p.pdf).
The parent directly reread the theorem statement on October 21, 2026; its
proof is not re-audited here. The same imported input already appears in the
repository's [source audit](../bounded_black_defect_finiteness_20261020/SOURCE_AUDIT.md)
and in
[THM-4728 — pronormal Frattini lemma and the C15 boundary](../../../01-canon/theorems/THM-4728-hyp-9340-beyond-p-groups-defect-identity-pronormal-frattini-lemma-and-the-smallest-failure-c15-at-n-10.md).

Without that cited theorem, §2 is a proved statement for **solvable**
transitive automorphism groups. With it, it applies to every
vertex-transitive tournament of order 35.

## 2. The full automorphism group contains a regular cyclic subgroup

Let \(G\) be a faithful solvable transitive automorphism group on 35 vertices.
Choose a nontrivial minimal normal subgroup \(N\). It is elementary abelian
of order a power of a prime \(p\). Normality and transitivity make all its
vertex orbits have the same size, a power of \(p\) dividing 35. Faithfulness
excludes size one for every orbit. Consequently those orbits are \(q\)
blocks of size \(p\), where

\[
\{p,q\}=\{5,7\}. \tag{2}
\]

An abelian group transitive on a block acts regularly after factoring its
kernel: a point stabilizer fixes the whole block because it commutes with
the transitive action. Thus the image of \(N\) on each block is regular
\(C_p\). Choose an origin and coordinate on each block. Its translations
give a faithful linear embedding

\[
N\hookrightarrow\mathbb F_p^q,
\]

with each coordinate functional \(\ell_i\) nonzero. The embedding is
faithful because an element acting trivially on every block fixes every
vertex.

Put \(i\sim j\) when \(\ell_i,\ell_j\) are proportional. This is an
equivalence relation on the blocks. It is \(G\)-invariant: conjugation acts
linearly on \(N\), while changing a block coordinate rescales its nonzero
translation functional. The action of \(G\) on the \(q\) blocks is
transitive, so all equivalence classes have equal size. Since \(q\) is
prime, either all blocks are equivalent or all classes are singletons.

### Coupled case: all functionals are proportional

Faithfulness now forces \(\dim_{\mathbb F_p}N=1\), so \(N=C_p\) acts
regularly within each block. Choose a \(q\)-element \(g\in G\) whose action
on the blocks is a \(q\)-cycle. Such an element exists by choosing a Sylow
\(q\)-subgroup projecting to a Sylow subgroup of the transitive degree-\(q\)
block action; that image has order \(q\).

Neither 5 divides \(7-1\) nor 7 divides \(5-1\). Therefore the
\(q\)-power-order conjugation action of \(g\) on \(N=C_p\) is trivial:
it would have to have order dividing \(p-1\). The element \(g^q\) fixes
all blocks and centralizes \(N\). On each regular \(p\)-block, this
centralizer is \(C_p\), so

\[
g^q\in C_p^q.
\]

But \(g^q\) also has \(q\)-power order, hence \(g^q=1\).
The commuting group \(\langle N,g\rangle\cong C_p\times C_q\cong C_{35}\)
is regular: one factor moves freely inside blocks and the other cycles the
blocks. Here the regular subgroup lies in the chosen \(G\).

### Independent case: all pairs of functionals are independent

For every two distinct blocks the projection of \(N\) onto their two
translation coordinates has rank two. It can move any ordered pair of
vertices in those blocks to any other. Hence every arc between those blocks
has the same direction. The tournament is an actual lexicographic blow-up:
its quotient is a vertex-transitive tournament of order \(q\), and its
fibres are mutually isomorphic vertex-transitive tournaments of order \(p\).

Every transitive group of prime degree contains a regular prime cycle, by
the Sylow point-stabilizer index argument. Choose such cycles on the quotient
and on one fibre, and choose explicit isomorphisms from that fibre to all
others. These identify the **unchanged** source with a lexicographic product.
The quotient cycle and the simultaneous fibre cycle commute in its full
automorphism group, giving a regular \(C_q\times C_p=C_{35}\).

It is essential that this last subgroup is asserted in the **full source
automorphism group**, not necessarily in the initially chosen \(G\).
For example, the augmentation translation code
\(\{(x_0,\ldots,x_6)\in\mathbb F_5^7:\sum x_i=0\}\), semidirect by
cyclic block shift, is transitive on the lexicographic source with seven
five-vertex fibres but has no element of order 35. The shift-norm sum is zero
on this code, and its only constant vector is zero. Its source still has the
larger full-automorphism product action above. The checker retains this
precise boundary rather than presuming the smaller group has the needed
element. ∎

This proves that **every vertex-transitive tournament of order 35 is cyclic
Cayley**. More generally the argument works for two distinct odd primes when
neither divides the other minus one, with the same solvability input. It does
not claim the analogous result when that numerical guard fails.

## 3. Complete source universe and literal exclusion

Choose an origin for the regular cycle. Every candidate becomes

\[
x\to y\iff y-x\pmod {35}\in S,
\quad S\mathbin{\dot\cup}(-S)=(\mathbb Z/35)\setminus\{0\}. \tag{3}
\]

The seventeen inverse pairs \(\{d,35-d\}\), \(1\le d\le17\), are
independent binary choices. A mask bit zero selects \(d\), and one selects
\(35-d\). Thus **all \(2^{17}=131,072\) source masks** form a complete
necessary universe, with isomorphism duplication allowed.

For each mask, the bank supplies seven actual Hamilton words with distinct
complete ordered-adjacency keys. All are black by the positional-reversal
test. There are 917,504 word positions but only 707 distinct literal words.
This proves \(\Delta\ge7\) for every source; inversion is anti, so
\(\Delta\) is even and in fact at least eight. No deficit-six source exists.

The positive source certificates, not an unsuccessful search, establish
the exclusion. The separate deterministic producer
[produce35.cpp](produce35.cpp) found them in 14,094,183 DFS nodes total;
the largest frame used 262,205 nodes and nine completed words. These
observations do not claim an isomorphism census, full \(H\), full \(|A|\),
or a universal discovery complexity bound.

## 4. Exact affine key compilation accelerates the consumer

The 595-bit key of a fixed literal word \(W\) has the exact form

\[
K_W(m)=K_W(0)\mathbin\oplus\bigoplus_{i:m_i=1} C_{W,i}. \tag{4}
\]

Each ordered pair is controlled by exactly one inverse-pair mask bit. Thus
(4) is coefficientwise affine over \(\mathbb F_2\), not an integer
interpolation assumption. Let \(K_{W^{\rm rev}}\) be the key of the reverse
ordering. The word is blue exactly when

\[
K_W(m)\mathbin\oplus K_{W^{\rm rev}}(m)
=2^{595}-1. \tag{5}
\]

The consumer packs the key into the low 595 bits and the nonzero-blackness
syndrome from (5) into the next 595 bits. Every evaluation uses only XOR;
no arithmetic carries or field overlap are possible. Hamilton legality is
the conjunction of the word's 34 actual inverse-pair bit requirements. If
repeated requirements contradict, its legal source domain is empty.

All words are authenticated complete permutations before compilation. The
compiled columns are built internally from their actual cyclic differences;
the consumer accepts no external claim of affine coefficients or blackness.
For **every one of the 707 words**, a separate literal adjacency-matrix
implementation checks both complete keys and the Hamilton clauses at zero
and all seventeen single-bit masks. Since the key is proved affine, these
eighteen sources independently verify every coefficient. This is in addition
to, not in place of, the complete 131,072-frame legality, blackness, and
distinct-key replay. Extra whole-source checkpoints also use the literal
matrix consumer.

The full key equality is iff equality of source-automorphism word orbits,
as proved in the [order-25 note](README.md). It is not only a score invariant.
The compiled key preserves that full information. Path words remain
available for any later metagraph-port construction; this exclusion alone
does not infer the wiring between them.

## 5. Bank, API, checks and reproduction

[inputs35.json](inputs35.json) stores the 707 literal words and the complete
ordered stream of seven word indices per source mask. The latter uses two
bytes per index, then lossless compression; the whole file is about 54 KB.
[pack35.py](pack35.py) packages the independent producer output. A commitment
to the expanded 917,504 literal-word stream is retained. Verification
recomputes the commitment, but the proof rests on every source-bound path
and key test, not the hash.

[check35.py](check35.py) performs no discovery. The public source and
seven-word APIs require literal integer masks, complete immutable word
tuples, and exact vertex alphabets. Cache access is behind public literal
validation. Wrong-source, duplicate-key, malformed-word, float/bool, and
range hostiles reject. The bounded bank reader rejects missing/trailing
compressed data and out-of-range source-word references.

```text
python -B 05-knowledge/results/center_transitivity_20261021c/check35.py
python -B -O 05-knowledge/results/center_transitivity_20261021c/check35.py --verify
```

Normal and optimized replays pass **2,122,941 body checks**. Receipt LF SHA256:
`360ac9b31f78935cda6699c5dc08693f75d0a02079e777d1d31be720f5e47be1`.
The earlier order-25 receipt is unchanged.

Independent reviews before the finite replay: tiling_fold passed the group
reduction, including the centralizer and full-automorphism qualifications;
collatz_memory passed the affine-XOR compiler algebra and the basis-check
scope. The parent supplied and reread the primary solvability theorem.
Final independent full proof/API audit and optimized complete-bank replay:
tiling_fold, **PASS**, with the exact saved count and hash above. No correction
was requested. This packet is frozen.
