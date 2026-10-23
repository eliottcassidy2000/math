# Fixed full Gram/theta collapse in the Paley order-eight conference switching class

Status: PROVED general algebra; FINITE-EXACT classification of the 128 specified labeled switches. Not a novelty claim for the inherited 189/333 split.

## Inheritance

- `01-canon/theorems/THM-472-tournament-hadamard-ceiling.md` already proves the skew-conference extremality criterion and records four order-eight isomorphism classes attaining it.
- `01-canon/theorems/THM-450-doubling-family-trichotomy.md`, section 6, and `01-canon/theorems/THM-454-blowup-product-and-lift-laws.md` already exhibit H=189 versus 333 under skew doubling.
- `01-canon/theorems/THM-505-the-ocf-non-spectral-defect-H-equals-spectral-skeleton-plus-witt-defect.md`, order-eight extension, already identifies the odd-cycle packing carriers used in section 4 below.
- The added conclusion here is that the **entire positive matrix Q itself**, not only its spectrum or determinant, and consequently its entire theta/Epstein lift, can remain identical across those different H-values. The three-vertex cycle/transitive switching example already shows theta-isometric collapse at H=3 versus 1; the conference example upgrades isometry to entrywise equality Q=8I.

## 1. General conference theorem

Let C be any skew tournament matrix of even order n with C²=−(n−1)I. Let d be any ±1 vector, D=diag(d), and C_d=DCD. Then C_d²=−(n−1)I, so

`Q_d=I−C_d²=nI`.

In particular every switch has exactly the same full finite heat trace, spectral zeta, lattice theta, Epstein zeta and determinant-normalized lattice shape:

`tr exp(−tQ_d)=n exp(−nt)`, `Z_Q(s)=n^(1−s)`,
`Θ_Q(t)=θ₃(0,exp(−πnt))^n`, `E_Q(s)=n^(−s/2) E_I(s)`.

The skew characteristic polynomial, resolvent numerator, and resolvent are respectively

`p(x)=(x²+n−1)^(n/2)`,
`q(x)=1ᵀ adj(xI−C_d)1=n x(x²+n−1)^(n/2−1)`,
`F(x)=q(x)/p(x)=n x/(x²+n−1)`.

Proof: `(xI−C_d)⁻¹=(xI+C_d)/(x²+n−1)` wherever defined, and `1ᵀC_d1=0` by skew-symmetry. Polynomial continuation supplies the numerator identity at singular x.

Since `A_d=(J−I+C_d)/2`, the matrix determinant lemma gives

`det(λI−A_d)=λ(λ−(n−2)/2)(λ²+λ+n/4)^(n/2−1)`.

Thus all these switches are adjacency-cospectral as well. At n=8 this is `λ(λ−3)(λ²+λ+2)^3`. In particular, neither Q nor the pair (p,F) recovers H even in this rigid family. If another result transports equality of (p,F) through a stated word algebra, that result may be applied; no unspecified word-algebra completeness is assumed here.

## 2. Exactly 2n universal-vertex switches for even n>2

Assume even n>2 throughout this section. Fix d_0=1 to remove the irrelevant sign d↦−d. For each vertex v and prescribed row sign ε∈{−1,+1}, there is exactly one switching making v a sink/source, namely the unique solution of

`d_v d_j C_vj=ε` for all j≠v, with d_0=1.

Distinct (v,ε) give distinct switches when n>2. Two sources or two sinks cannot coexist in any tournament. If one source and one sink coexisted, their rows would have inner product `−(n−2)`, contradicting CCᵀ=(n−1)I. Thus, for even n>2, the switching class contains exactly 2n members with a universal vertex. At n=2 there are only two switches, each with both a source and a sink, so the 2n count does not apply. Deleting it leaves a doubly regular tournament: writing the block with its uniform border, the conference-square equation gives B1=0 and B²=J−(n−1)I. This is structural, not a conclusion extrapolated from n=8.

At n=8, n−1=7, the exact classification below identifies these residual tournaments with Paley P7 or its reverse, so H=H(P7)=189 for all 16 universal members.

## 3. Complete specified 128-switch certificate

Use vertices 0,...,7, with C_0j=+1 for j>0, C_j0=−1. On vertices 1,...,7 put C_ij=+1 iff `(j−i) mod 7` belongs to `{1,2,4}`. Direct integer multiplication verifies C²=−7I. A mask m∈[0,127] specifies d_j=−1 iff bit j−1 is set, j=1,...,7; d_0=1.

The four isomorphism types are:

- Representative mask 0: 8 switches, scores (3,3,3,3,3,3,3,7), automorphism order 21, H=189.
- Representative mask 1: 56 switches, scores (2,2,2,4,4,4,4,6), automorphism order 3, H=333.
- Representative mask 3: 56 switches, scores (1,3,3,3,3,5,5,5), automorphism order 3, H=333.
- Representative mask 11: 8 switches, scores (0,4,4,4,4,4,4,4), automorphism order 21, H=189.

Reversal exchanges the first and fourth types, and the second and third types. Thus the 112 nonuniversal members are **two**, not one, isomorphism types. Complete mask lists are in the output certificate. The program verifies isomorphism by exact degree-compatible permutations, not a hash or spectral surrogate. H is independently computed by subset dynamic programming on Hamiltonian paths.

## 4. Where the extra 144 Hamiltonian paths reside

Write c_r for the number of directed simple r-cycles modulo cyclic rotation. Write c_33 for unordered disjoint pairs of directed triangles and c_35 for disjoint directed triangle–pentagon pairs. The tournament odd-cycle formula at n=8 specializes to

`H=1+2(c3+c5+c7)+4(c33+c35)`.

Direct cycle enumeration provides independent checks:

- Universal types: `(c3,c5,c7,c33,c35)=(14,42,24,7,0)`; H=189.
- Nonuniversal types: `(c3,c5,c7,c33,c35)=(14,42,48,19,12)`; H=333.

Therefore the 144 difference splits exactly into three equal contributions:

`ΔH=2·24 + 4·12 + 4·12 =48+48+48`.

They come from extra seven-cycles, extra disjoint triangle pairs, and extra triangle–pentagon pairs. Equality of the full adjacency characteristic polynomial does not imply equality of these simple-cycle/packing statistics; its signed cycle-cover combinations conceal this redistribution. No new odd-cycle formula is claimed: the experiment uses the inherited identity as a second computation path.

## 5. Fixed-gauge Fano codes, with a caution

The source masks are the binary [7,3,4] simplex code

`{0,23,46,57,75,92,101,114}`.

Its seven nonzero supports have size four, and their complements form a Fano plane. The complements (XOR 127) of the sink masks form **a different** simplex code

`{0,29,39,58,78,83,105,116}`.

Their complement supports form another labeled Fano plane. Closure under XOR, the weight enumerator `1+7z^4`, and exactly one line through every pair are verified directly. It would be incorrect in this fixed labeling to call the universal switches a code plus its complement: the source code and complemented sink code are different. This is a useful gauge-sensitive sidecar rather than a claim that the Fano/simplex connection itself is new.

## 6. Interpretation and reproduction

The positive lattice lift supplies a mathematically genuine theta/Mellin duality, but it forgets the orientation information that controls H. Here even supplying the full entrywise Q instead of its eigenvalues cannot repair that loss. A future bridge from theta zeros to H must include an orientation-sensitive sidecar; spectral flatness, all moments of Q, theta coefficients, and (p,F) are jointly insufficient in this example.

Run `python3 session35-zeta/CONFERENCE_COLLAPSE.py`. Uses Python standard library and SymPy. Output: `session35-zeta/CONFERENCE_COLLAPSE.out`. The complete universe is exactly the 128 switches of the explicitly given C, not all labeled order-eight tournaments. The classification and counts are FINITE-EXACT; the general conference algebra and the 2n universal-switch theorem for even n>2 are proved above.
