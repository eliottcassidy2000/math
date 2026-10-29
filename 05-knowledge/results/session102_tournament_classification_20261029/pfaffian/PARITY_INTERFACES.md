# Rank-one Pfaffian cuts: parity, root covariance, and the second-port obstruction

**Status: PROVED algebraically; exact checker results in `results.json`.** This is a bounded interface result, not a classification of prime tournaments. No literature-priority claim is made.

## 1. Inheritance and scope

Pinned repository: `3b6e3ef82382da222fc59aa02914eea72d5a105b`.
Exact retrieved source blobs, not the older discovery checkout, are recorded in `sources/manifest.json`.

* Closest proved mechanism: **THM-3729**, `01-canon/theorems/THM-3729-rooted-pfaffian-response-and-sign-root-deletion-average.md`, identifies the sum over all odd principal rooted Pfaffian squares with `det(I+A) uᵀ(I+A)⁻¹u`. It already proves covariance `(A,u)→(DAD,Du)` and warns that resetting the root is invalid. Our top-degree linear rooted response is one term of its exact energy, not a replacement for that theorem.
* Prior algebraic source: **THM-174**, `01-canon/theorems/THM-174-pfaffian-determinant-identity.md`, uses the Pfaffian cofactor vector for odd tournament matrices. We independently prove all needed rank and sign facts below; we do not import its unrelated Hamiltonian-path consequences.
* Parity sidecar: **THM-1455**, `01-canon/theorems/THM-1455-skew-pfaffian-parity-and-a-mod-16-law.md`, recalls the classical vanishing of odd skew determinants and squared-Pfaffian principal minors.
* Canonical inherited hostile: THM-3729's five-vertex failure of unsigned square energy. Here the new operation-specific hostile is a four-vertex pair with indistinguishable separate scalar root energies and different joint attachments.
* Corrected near miss: an unrooted scalar is enough for an even/even rank-one cut's **closed Pfaffian**, but not for arbitrary future attachment. The missing coordinate is a mixed alternating bilinear response.
* Concept board: rank-one cut / parity; rooted cofactors / covariance; Schur transfer / mixed port coupling; scalar energy / information loss. All four are tested by the hostile and repaired by the multiport transfer.

The underlying Pfaffian block elimination is standard linear algebra, proved here for signs and self-containment. Retaining boundary Pfaffian signatures is also standard in matchgate theory: Jin-Yi Cai and Aaron Gorenstein, *Matchgates Revisited*, Theory of Computing 10(7) (2014), 167–197, https://theoryofcomputing.org/articles/v010a007/ . This supplies context, not a claim that our arbitrary tournament matrices are planar matchgates. No matchgate classification is imported.

## 2. Conventions and the parity law

Work over the real numbers (the polynomial identities also hold over any characteristic-zero field). Vertices are ordered, indices start at **zero**, `Pf([])=1`, and `Pf([[0,a],[-a,0]])=a`. Odd-dimensional Pfaffians are treated as zero only where explicitly convenient. For odd `m`, define

`p_i(A)=(-1)^i Pf(A with row/column i deleted)`, `0≤i<m`.

Expansion at a last root gives

`R(A;u):=Pf([[A,u],[-uᵀ,0]])=p(A)ᵀu`.                 (1)

With one-based indices, the same rooted convention is `(-1)^(i+1)`; using `(-1)^i` instead negates both cofactor vectors, so the odd/odd product below is unchanged, but a single rooted response changes sign.

Let `A` and `B` be skew matrices of orders `m,n`, with the A vertices first, and let

`S=[[A,uvᵀ],[-vuᵀ,B]]`.

When `m+n` is even,

* **m,n even:** `Pf(S)=Pf(A)Pf(B)`;
* **m,n odd:** `Pf(S)=R(A;u)R(B;v)`.

These identities do not require nonsingularity or tournament entries. For tournament cuts, require `u_i,v_j∈{−1,+1}` and internal off-diagonal entries `±1`; then the cross edges form a genuine signed rank-one cut. Neither the existence nor uniqueness of such cuts in arbitrary tournaments is asserted.

### Proof, including the cross sign

Let `α,β` be the internal two-forms and let `U=Σu_i e_i`, `V=Σv_j f_j`. The cut two-form is `U∧V`, with square zero. Hence

`exp(α+β+U∧V)=exp(α)exp(β)(1+U∧V)`.

A top coefficient uses zero or one cross pair. Parity forces zero for even/even and one for odd/odd. In the latter case choosing the pair `(i,m+j)` contributes deletion sign `(-1)^(i+m+j-1)=(-1)^(i+j)` because m is odd. Summing its internal Pfaffian products gives exactly `(p(A)ᵀu)(p(B)ᵀv)`. This is also a cancellation proof: terms with two or more cross pairs cancel because all 2×2 minors of `uvᵀ` vanish.

## 3. Gauge and reversal checks

For diagonal sign matrices `D_A,D_B`, let `D=diag(D_A,D_B)`. Covariantly transform

`A'=D_A A D_A`, `B'=D_B B D_B`, `u'=D_Au`, `v'=D_Bv`.

Then `S'=DSD`, `Pf(S')=det(D)Pf(S)`, and, in odd order,

`p(A')=det(D_A)D_A p(A)`, `R(A';u')=det(D_A)R(A;u)`.

Thus the parity formulas obey the same gauge. Root reset `u'→1` is **not** allowed without recomputing the response.

For even total order N, `Pf(-S)=(-1)^(N/2)Pf(S)`. To represent full reversal in the rank-one factors, use `A'=-A`, `B'=-B`, `u'=-u`, `v'=v`. For odd m,

`p(-A)=(-1)^((m-1)/2)p(A)`.

The extra minus in u gives precisely total factor `(-1)^(N/2)` in odd/odd. Merely replacing A and B by their negatives without reversing the cut is not full reversal. For odd rooted A, full reversal of the bordered matrix transforms **both** `A→-A` and `u→-u`, with factor `(-1)^((m+1)/2)`.

## 4. Mixed parity exposes the missing coordinate

Assume A has even order and is invertible, B odd. Give S a second external root `w=(x,y)` in A/B order. Then

`R(S;w)=Pf(A)[p(B)ᵀy+(p(B)ᵀv)uᵀA⁻¹x]`.             (2)

The term `uᵀA⁻¹x=-xᵀA⁻¹u` is alternating and generally nonzero. It is not determined by separate scalar root energies.

To prove (2), eliminate A in the bordered matrix. Its coupling columns are `[uvᵀ,x]`. The residual B block is unchanged because `uᵀA⁻¹u=0`, while its root column becomes `y+v(uᵀA⁻¹x)`. Equation (1) gives the result.

### Explicit smallest even-interior / two-port hostile

Let `A=[[0,1],[-1,0]]`, `u=(1,1)`, `B=[0]`, `v=y=1`. Use two alternatives `x₀=(1,1)`, `x₁=(1,-1)`. The parent S is identical in both cases; its final root attachment changes:

```
W0 = [ 0  1  1  1 ]    W1 = [ 0  1  1  1 ]
     [-1  0  1  1 ]         [-1  0  1 -1 ]
     [-1 -1  0  1 ]         [-1 -1  0  1 ]
     [-1 -1 -1  0 ]         [-1  1 -1  0 ]
```

All off-diagonal entries are ±1: both are honest tournaments. In both, `Pf(A)=1`, `det(I+A)=2`, and for either separate attachment vector z in `{u,x₀,x₁}`,

`zᵀ(I+A)⁻¹z=1`, `det(I+border(A,z))=4`, and the unnormalized THM-3729 odd square energy is 2.

Nevertheless `uᵀA⁻¹x₀=0`, `uᵀA⁻¹x₁=2`, so `Pf(W0)=1`, `Pf(W1)=3`, and determinants are 1 and 9. The first failed implication is “equal separate one-port scalar responses determine a joint two-port response.” This does **not** refute a grammar retaining the full vectors or the mixed bilinear term. Minimality is only in the stated even-interior/two-distinct-port model: a nonempty even interior needs two vertices, and two ports add two. No broader minimality is claimed.

## 5. Sufficient multiport skew transfer (exact repair)

For an even invertible interior A, k labelled ports, coupling matrix U, and skew port matrix E,

`Pf([[A,U],[-Uᵀ,E]])=Pf(A) Pf(E+UᵀA⁻¹U)`             (3)

whenever total order is even. Thus the pair `(Pf(A), UᵀA⁻¹U)`, retaining port labels and order, is sufficient for every external E on these same ports. This is a boundary-specific response: new coupling columns require recomputation, and it does not classify the interior tournament. If the specified port columns span an isotropic subspace for `A⁻¹` (equivalently `UᵀA⁻¹U=0`), the transfer reduces to the scalar prefactor. In particular, all columns proportional to one port vector satisfy this condition. For exactly two ports, a nonzero mixed entry changes `Pf(A)(E₁₂+K₁₂)`, so this zero-mixed-response condition is also necessary for the uncoupled scalar formula to hold for every external port edge. In (2), an attachment x proportional to u likewise kills the missing term; unrestricted new roots do not.

**Proof.** With `P=[[I,-A⁻¹U],[0,I]]`, direct multiplication yields `Pᵀ[[A,U],[-Uᵀ,E]]P=diag(A,E+UᵀA⁻¹U)`. Since `det(P)=1`, Pfaffian congruence proves (3).

Every even tournament skew matrix is invertible: its Pfaffian is a sum of `(2r−1)!!` terms, each ±1, hence is odd and nonzero. Consequently every odd tournament matrix has rank exactly n−1, since each one-vertex deletion is invertible and odd skew rank cannot be n.

For an **odd tournament interior** A, choose any pivot vertex a. Put `I={0,…,n−1}\{a}` in inherited order, `C=A[I,I]`, `b=A[I,a]`, `V=U[I,:]`, `z=U[a,:]`. Then C is invertible (including the empty n=1 case), and the complete exact transfer is

`γ=z+bᵀC⁻¹V`, `K=VᵀC⁻¹V`,

`Pf([[A,U],[-Uᵀ,E]])=(-1)^a Pf(C) Pf([[0,γ],[-γᵀ,E+K]])`. (4)

The pivot is moved to the end of the interior by `n−1−a` transpositions, whose parity is a since n is odd; eliminate C using (3). The result keeps one uneliminated odd mode **and** all mixed port responses. It is a sufficient finite multiport state for the specified attachment pattern, not a scalar one-port state. Formula (4) remains valid for a general odd skew A whenever the chosen C is invertible; arbitrary singular interiors need a different pivot rank or a full minor signature.

## 6. Reproduction and boundaries

Run `python3 check_pfaffian.py > results.json` from this directory, or run by its full path. The checker uses exact integers/Fractions, recursive Pfaffians, independent Gaussian determinant checks, explicit exceptions (not optimization-disabled assertions), exhaustive signed rank-one cuts with positive block sizes and even total ≤6, complete mixed rooted universes `(m,n)=(2,1),(2,3),(4,1)`, and seeded multiport checks for interior orders 1–7 and port counts 1–4 with even total. The cofactor sign and its direct gauge covariance, switching, reversal, every possible odd-interior pivot, and hostile are included. An additional 300 seeded weighted rank-one cases permit zero entries and singular blocks, checking the non-tournament polynomial scope. These finite checks support, but do not replace, the all-size proofs.

This says nothing about completeness of prime-tournament lists, Hamiltonian-path counts, fixed-root `H≥disc`, or LRC. The operation closure is exact only after retaining parity, order, gauges, specified port vectors, and mixed skew responses.

### Source-discovery limit

The three named theorem blobs were recovered exactly from the pinned commit. An optional broader `git grep` encountered unavailable lazy GitHub blobs (`Proxy CONNECT aborted`) and was canceled; that search is incomplete. No claim of novelty or exhaustive prior-art recovery is made.
