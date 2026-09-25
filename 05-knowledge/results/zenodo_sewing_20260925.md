# Exact symmetry and boundary transfer from the nine-loop symbol dataset

Status: **CITED / SOURCE-READ**, the deposited format and conventions; **PROVED**, the elementary representation transfer below; **FINITE-EXACT**, its bounded controls. This note does not recompute scattering amplitudes or prove Collatz convergence. Date: 2026-09-25.

## 1. What the linked source actually supplies

The primary source is [He, Jing, and Li, *The Symbols of Six-Gluon MHV Amplitudes through Nine Loops*, Zenodo 22800071](https://doi.org/10.5281/zenodo.22800071), deposited September 17, 2026. Its README, Wolfram loader, metadata tables, and manifest were read from the archive. No accompanying paper is named in its metadata, README, or loader; bounded title/author searches did not identify one. This is a limit of the source search, not a claim that no paper exists.

The release stores maximal-weight symbols of the BDS-like subtracted six-point MHV amplitude in planar N=4 SYM, through weight 18. Symbol-invisible constants are omitted. First-entry and last-entry recurrences are joined by a sewing matrix. Individual word coefficients can be recovered by sparse contractions, retaining the specified terminal coordinate. Basis blocks carry parity and D3 representations; their component factors are already included. These are the source's declared data and conventions, not an independent validation of integrability or the physics construction.

The useful methodological transfer is concrete: preserve both boundary bases, contract through the actual joining operator, and keep the terminal convention and symmetry factors. A compressed expression is exact only relative to those specified data. It does not itself imply that higher orders, all heights, or all trajectories have been covered.

## 2. Inheritance and the closest exact arithmetic mechanism

The nearest proved arithmetic mechanisms are the signed affine-word guards in [creation_fano_20260925.md, F3](creation_fano_20260925.md), the common-future and suffix certificates in [creative_decoder_20260925.md](creative_decoder_20260925.md), and the finite-state single-step carry machine in [creative_transducer_20260925.md](creative_transducer_20260925.md). They preserve the exact integer and legal parity condition rather than just the finite alphabet of branch labels.

The source does not provide a new proof of these inherited mechanisms. It suggests a useful organization: calculate the left and right boundary data independently, then join them using a small exact compatibility test. The parent [zenodo_bridge_20260925.md](zenodo_bridge_20260925.md) implements that proposal. Its guarded word matrix is

\[
 M_w=\begin{pmatrix}A_w&B_w\\0&D_w\end{pmatrix},\qquad
 A_w=3^{\#\text{odd symbols}},\quad D_w=2^{|w|}.
\]

For `w=uv`, a starting integer n, and a target r, its sewing residual is

\[
 K_{n,r}(u,v)
 =A_v(A_un+B_u)+D_u(B_v-D_vr).
\tag{1}
\]

The two boundary values are dyadic and triadic rationals:

\[
 \frac{A_un+B_u}{D_u},\qquad
 \frac{D_vr-B_v}{A_v}.
\tag{2}
\]

Their equality forces an integer, since `Z[1/2]∩Z[1/3]=Z`. Moreover the complete word's integral endpoint restores every intermediate parity guard. Thus this arithmetic seam preserves the needed condition, rather than just matching a dimension count.

**Independent audit of the parent statement.** At positive n,r and positive prefix/suffix lengths, (1) has exactly rank two as a rational matrix indexed by all words of those lengths, for either sign σ. Its row factors are `(A_un+B_u,D_u)` and column factors are `(A_v,B_v-D_vr)`. Comparing the all-even prefix with an odd-first prefix changes the first row coordinate by `2n+σ`, which is nonzero. The analogous two suffix columns have determinant `σ+2^(|v|+1)r`, also nonzero. Both factors therefore have rank two. Their integer registers can grow without bound; this is not a finite-state realization of the full parity-address map.

The parent's word-matrix peeling decoder supplies a second necessary distinction. It recognizes and reconstructs exactly which matrices were created by legal branch letters. That representation-membership certificate is not a certificate that every starting integer has a word reaching a terminal state. For the shortcut clock, the terminal bank must account for `1↔2`; literal arrival at 4 is not required from those two states.

Live board: boundary compatibility; legal-word reconstruction; symmetry blocks; parity sidecars; finite-order versus all-height coverage. Historical files containing “amplitude” or “hexagon” did not supply a proved scattering-symbol transfer in the targeted search. In particular, the speculative `amplituhedron_987_chem_s90br.out` is not a dependency. The least-used useful sidecar is the actual basis of the source's two-dimensional representation.

## 3. The source's D3 block has our exact three-color quotient

The source README gives the two-dimensional representation in a flip-even/odd basis, using the following column-action convention:

\[
 C=\begin{pmatrix}-1/2&-1/2\\3/2&-1/2\end{pmatrix},\qquad
 F=\begin{pmatrix}1&0\\0&-1\end{pmatrix}.
\tag{3}
\]

**Z1 (PROVED).** There is an explicit invariant integral lattice in this rational representation. With

\[
 B=\begin{pmatrix}1&-1\\1&1\end{pmatrix},
\quad
 \widehat C=B^{-1}CB=\begin{pmatrix}0&-1\\1&-1\end{pmatrix},
\quad
 \widehat F=B^{-1}FB=\begin{pmatrix}0&-1\\-1&0\end{pmatrix},
\tag{4}
\]

the lattice is `B Z²`. Both integral matrices preserve

\[
 N(a,b)=a^2-ab+b^2,
\qquad
 \widehat C^3=I,\quad\widehat F^2=I,
\quad\widehat F\widehat C\widehat F=\widehat C^{-1}.
\tag{5}
\]

The norm identity in the original coordinates is

\[
 3(a-b)^2+(a+b)^2=4N(a,b).
\]

All these statements follow by multiplication of the displayed two-by-two matrices. This is the integral hexagonal norm, with a specified lattice and basis; no arithmetic dynamics have yet been assigned to it.

Reducing this lattice representation modulo two gives

\[
 \overline C=\begin{pmatrix}0&1\\1&1\end{pmatrix},\qquad
 \overline F=\begin{pmatrix}0&1\\1&0\end{pmatrix}.
\tag{6}
\]

The first is exactly the matrix from [duck_zeckendorf_20260925.md](duck_zeckendorf_20260925.md): it rotates `(1,0)→(0,1)→(1,1)→(1,0)` and fixes zero. The second swaps the first two colors. These six distinct matrices are all of `GL₂(F₂)`, so reduction is faithful on D3. This is an actual map from the specified source representation to the recovered three-color action.

The basis step cannot be omitted: (3) has even denominators and cannot be reduced directly modulo two. Also, the construction defines a particular invariant lattice in the representation. It does not claim that arbitrary rational coefficients of the stored amplitude data belong to that lattice or have canonical residues modulo two.

## 4. Two exact losses: parity and arithmetic height

**Z2 (PROVED).** On an odd-parity E block, adjoining the source's parity generator `−I` gives twelve matrices. Reduction modulo two has kernel precisely `{I,−I}`. It therefore preserves the D3 color action but erases the parity factor. On an even-parity E block, parity already acts trivially before reduction.

Indeed the six D3 matrices already have distinct reductions, while g and −g have the same reduction. No D3 element equals −I, since its determinant-one subgroup has order three. Consequently the full group has twelve elements and exactly a two-element kernel.

The source distinguishes its order-three C from parity times C. On an odd-parity block, `−C` has order six because `(-C)^3=-I`; after reduction it has the same order-three color action as C. Keeping only the three colors cannot tell those source actions apart. Source parity, odd/even integer parity, the sign parameter in `3n±1`, and the neutral fourth color are distinct notions unless an additional map is supplied.

The norm is also lost: the vectors `(1,0)` and `(3,0)` have the same color but norms one and nine. The zero color contains nonzero vectors as well. Thus no height bound follows from the finite color action.

There is a further hostile lift worth retaining. The integer Fibonacci matrix

\[
 M=\begin{pmatrix}0&1\\1&1\end{pmatrix}
\]

has the same reduction as `widehat C`, but

\[
 M^3=\begin{pmatrix}1&2\\2&3\end{pmatrix}\ne I,
\quad \widehat C^3=I.
\]

An equality of color actions does not identify the integral Fibonacci growth operator with the finite-order source rotation. This is precisely where an integer or carry sidecar must remain.

## 5. A checkable sewing convention for the two-dimensional block

**Z3 (PROVED, with an explicit representation convention).** If both boundary vectors transform by the matrices in (3), a bilinear pairing `xᵀHy` is invariant under D3 exactly when

\[
 H=\lambda\begin{pmatrix}3&0\\0&1\end{pmatrix}.
\tag{7}
\]

Reflection invariance forces the off-diagonal entries to vanish. Rotation invariance then forces the first diagonal entry to be three times the second. In particular H=I is not invariant in the README's basis. If one boundary uses a dual rather than covariant basis, its transformation rule changes and (7) must be adjusted accordingly.

Likewise a parity-even scalar can pair covariant blocks of equal parity; opposite parities change the sign of the pairing. Their modulo-two color quotients cannot enforce that distinction. These are useful structural checks when transporting a sewing method. They do not instruct the reader to insert extra factors into the deposited matrices, whose README says those factors are already included.

The resulting transfer is deliberately exact in scope: a specified D3 representation reduces to the earlier three-color action; a guarded arithmetic seam uses two unbounded rational boundary registers; terminal and parity data remain explicit. Neither finite loop order nor finite representation dimension establishes all-height Collatz coverage. That is the remaining **OPEN** obligation in the parent construction.

## 6. Reproduction and audit universe

An independent proof audit of Z1--Z3 passed, including the invariant lattice,
all six/twelve group elements, the reduction kernel, and the pairing metric.

```text
python 04-computation/experiments/zenodo_sewing_20260925.py
python -O 04-computation/experiments/zenodo_sewing_20260925.py
```

The [script](../../04-computation/experiments/zenodo_sewing_20260925.py) verifies the rational change of basis, all six D3 elements, all twelve parity-times-D3 elements, the exact reduction kernel, the invariant metric, and 6,534 norm/intertwining cases on integer vectors in `[-16,16]²`. It also checks the six-to-three order collapse, the failed identity metric, distinct heights with the same color, and the incompatible integer Fibonacci lift. All checks remain active under optimization. The [output](zenodo_sewing_20260925.out) is deterministic.

The program does not execute the supplied Wolfram Language or parse the WXF tensors. Archive integrity, tensor metadata, and arithmetic sewing replays are handled in the parent lane. An unlocated accompanying paper and unrecomputed physical integrability conditions are not promoted to proved dependencies.
