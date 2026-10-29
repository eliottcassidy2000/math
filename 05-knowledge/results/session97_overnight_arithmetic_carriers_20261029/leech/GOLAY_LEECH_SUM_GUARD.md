# From the 23-cycle Golay kernel to Leech: the indispensable sum guard

**Status:** classical construction, independently proved below; FINITE-EXACT
certificate for the explicit code and its complete minimum shell. Not an
external novelty claim. Date: 2026-10-29. Independent audit: PASS.

## Inheritance and scope

The input is the binary code constructed in
`../supersingular/check_cyclic_frobenius.py`: use the degree-12 polynomial

\[
g(X)=1+X^2+X^5+X^8+X^9+X^{10}+X^{11}+X^{12}
\]

(bitmask 7973). Embed the eleven shifts \(X^ig\), \(0\le i<11\), in
24 coordinates, with coordinate 23 zero, and adjoin the all-ones vector.
Our checker independently builds their span; no producer imports are used.
It checks dimension 12, self-duality, cyclic invariance on the first 23
coordinates, and weight enumerator

\[
1+759z^8+2576z^{12}+759z^{16}+z^{24}.
\]

This is the explicit extended binary Golay code \(C\). The genuine role of
23 here is a named order-23 coordinate permutation fixing the added coordinate.
Its membership in the Ogg-prime list was inherited from the separate arithmetic
lane; none of the proof below classifies that list or constructs the Monster.

Closest mechanism: the cyclic Frobenius kernel and Construction A in
`../supersingular/AFFINE_RECONSTRUCTION_AND_ISOSPECTRAL_LATTICES.md`.
Corrected near miss: Construction A of Golay is **not** Leech.
Hostile: the vector \(4e_1\) below. Niche bridge: exact parity gluing,
not equality of numerals or weight enumerators alone.

## Classical source and normalization

Daniel Allcock, *Orbits in the Leech Lattice*, Experimental Mathematics
14 (2005), 491–509, §2, gives exactly the parity/residue/sum construction
with inner product \(\frac18\sum x_i y_i\):
[author's primary research manuscript, p. 4](https://web.ma.utexas.edu/users/allcock/research/algor.pdf).
Its [errata](https://web.ma.utexas.edu/users/allcock/research/algor_errata.html)
should accompany literature use. We independently prove the required lattice
properties rather than rely on the orbit algorithm.
Robert Wilson's [notes, §5.3.1, equations (5.11)–(5.12)](https://webspace.maths.qmul.ac.uk/r.a.wilson/FSG/notes5.pdf)
corroborate the congruences and three minimum-vector families. These are
classical results; our contribution here is the repo-specific independently
checked bridge and its explicit information-loss guard.

## The construction and theorem

Let \(K\subset\mathbb Z^{24}\) consist of all \(x\) satisfying, for one
\(m\in\{0,1\}\),

1. \(x_i\equiv m\pmod2\) for every coordinate;
2. \(c(x)=\{i:x_i\equiv m+2\pmod4\}\in C\);
3. \(\sum_i x_i\equiv4m\pmod8\).

Set \(\Lambda=K/\sqrt8\), with the ordinary Euclidean inner product.
Then \(\Lambda\) is even unimodular, has minimum squared length 4, and has
196560 minimum vectors. It is the classical Leech lattice in the cited model.
The source's condition that *every* residue-class support be a codeword is
equivalent to condition 2: the only other nonempty class is its complement,
and the all-ones vector belongs to \(C\).

### Closure and index

For summands of parities \(m,n\), their new code bits are
\(c+d\) unless \(m=n=1\), when they are \(c+d+\mathbf1\).
Thus code closure and the all-ones word give addition closure. Sum guards add
correctly since \(4+4\equiv0\pmod8\). Negating an even vector keeps its code
bits; negating an odd vector complements them. So \(K\) is a subgroup.

Common coordinate parity costs index \(2^{23}\). The 24 residue bits then
cost index \(2^{24-12}=2^{12}\). To justify that every codeword and both
parities occur before imposing the sum guard, take \(x=m\mathbf1+2c\).
Its sum is divisible by 4, so exactly one binary choice remains modulo 8.
Adding \(4e_1\) preserves the first two guards and toggles that choice.
Hence the last index is exactly 2:

\[
[\mathbb Z^{24}:K]=2^{36}.
\]

Scaling by \(1/\sqrt8\) in dimension 24 changes covolume by
\(8^{-12}=2^{-36}\), so \(\Lambda\) has covolume 1.

### Evenness and integrality

If \(x=2z\), then \(z\bmod2\in C\). Double evenness gives
\(\sum z_i^2\equiv0\pmod4\), so \(\sum x_i^2\equiv0\pmod{16}\).
If \(x=\mathbf1+2z\), then \(z\bmod2\in C\) gives
\(\sum z_i^2\equiv0\pmod4\), while the sum guard gives
\(\sum z_i\equiv2\pmod4\). Therefore

\[
\sum x_i^2=24+4\sum(z_i^2+z_i)\equiv0\pmod{16}.
\]

Thus every squared norm in \(\Lambda\) is even. Polarization and closure
imply integral pairings. Integral covolume-one lattices are unimodular.

### Minimum and exhaustive shell

For even \(x\), a nonzero code support has at least eight coordinates of
absolute value at least 2, so \(\sum x_i^2\ge32\). With zero code support,
all coordinates are multiples of 4. A single \(\pm4\) violates the sum guard;
two are necessary at the minimum. For odd \(x\), all absolute values are at
least 1. If they are all 1, their negative-sign support is a codeword of weight
divisible by 4, so their sum is 0 modulo 8, violating the odd sum guard.
At least one absolute value is 3, giving \(\sum x_i^2\ge24+8=32\).

Equality allows exactly the following disjoint types:

| Integer-coordinate shape | Count | Exact sign condition |
|---|---:|---|
| \((\pm4)^2,0^{22}\) | \(4\binom{24}{2}=1104\) | all four sign choices |
| \((\pm2)^8,0^{16}\) | \(759\cdot2^7=97152\) | octad support; even number of negatives |
| \((\pm3),(\pm1)^{23}\) | \(24\cdot4096=98304\) | one vector per distinguished position and codeword |

For the last row, choose \(j\) and \(c\in C\). Set
\(x_i=1-2c_i\) except \(x_j=-3(1-2c_j)\).
The residue support is exactly \(c\), and the sum is
\(24-2\operatorname{wt}(c)-4(1-2c_j)\equiv4\pmod8\).
Conversely a minimum vector of this type uniquely recovers \(j,c\).
The counts add to **196560**.

## Why plain Construction A is different

\[
A(C)=\{z/\sqrt2:z\in\mathbb Z^{24},\ z\bmod2\in C\}
\]

is even unimodular, by self-duality and double evenness. A squared-norm-2
vector requires \(\sum z_i^2=4\): either \(\pm2e_i\), or four unit
coordinates. The latter would be a forbidden weight-four codeword. Thus
there are precisely 48 roots \(\pm\sqrt2e_i\), with root system \(A_1^{24}\):
the corresponding Niemeier lattice, **not Leech**.

At the Leech scaling, dropping only the sum guard admits \(4e_1/\sqrt8\),
a root of norm 2; it also admits \(\mathbf1/\sqrt8\), of norm 3.
That enlarged lattice has covolume 1/2 and is not the even unimodular Leech
lattice. The even-coordinate sector of \(K\) alone has index 2 in \(K\),
so retaining only it loses unimodularity. The odd glue restores volume while
the sum constraint removes the bad short vectors. This is a paired
constraint-and-glue mechanism, not a free deletion of roots.

## Carrier / loss ledger

| Source | Target / map | Preserved | Lost or added guard | Hostile/test |
|---|---|---|---|---|
| Named cyclic kernel | extended code span | coordinate action and exact code membership | one fixed extension coordinate | independent dimension and rotation checks |
| Golay code | plain Construction A | parity code, even self-duality | does not remove axis roots | 48 roots |
| integer residues | Leech lattice | named coordinate code and metric | retain sum modulo 8; include odd glue | \(4e_1\), all-ones rejected |
| Leech minimum shell | its cardinality | kissing number | loses incidence and isometries | no Monster identification follows |

This supports a reusable design principle: when lifting a finite combinatorial
carrier to a metric object, keep the extension/gluing cocycle as part of the
carrier. It does not prove a universal minimum-carrier theorem, identify
all symmetries of Leech, classify supersingular primes, or connect arbitrary
tournament exclusions \(\{7,21\}\) to Monster structure.

## Reproducibility

Run `python3 check_golay_leech_guard.py` and
`python3 -O check_golay_leech_guard.py`. Both return identical receipts with
**594001 active checks**. All 196560 minimum vectors are explicitly generated,
validated, and checked for duplicate tuples. The three-family exhaustiveness
is proved above, not inferred from a search box. The checker also validates
named order-23 action, roots of Construction A, guard counterexamples, and
deterministic closure/pairing samples. Universal closure/evenness/index are
proved above; the samples are supplemental, not substitutes for proofs.

Independent audit by the arithmetic-bridge auditor re-derived the proof and
counted the shell using compact integer-encoded vectors rather than our tuple
representation. Its `../audit/check_leech_independent.py` independently recovers
the code, all three shell counts and all guards; normal/-O receipts agree.
A separate residue count gives \( |K/8\mathbb Z^{24}|=2\cdot4096\cdot2^{23} \),
hence index \(8^{24}/2^{36}=2^{36}\), corroborating the index proof.
The cited Allcock errata affect later orbit tables, not the construction.
