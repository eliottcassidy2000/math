# Gaussian annihilators and six-colour port obstructions

**Date:** 2026-10-14. **Recovery status:** reconstructed from retained research history after the earlier workspace disappeared. The accompanying self-contained script was rerun locally; no network access or repository modification was used for this recovery. These files have new hashes. They are not claimed to be byte-identical to the lost originals.

**Mathematical status:** PROVED below: the polynomial-density Gaussian-annihilator classification, its Keller consequence, the list-peeling criterion, and a fourteen-vertex necessary bound conditional on the repository's exact nine-apex census. FINITE-EXACT: thirty coefficient-kernel computations, six independent direct-Gaussian-moment matrix computations, four graph degeneracy/chromatic checks, forty-five residue checks, and the regular-pentagon calculation. The classification and combinatorial count were independently reviewed before recovery. The inherited geometric census has not been replayed here. No literature priority, new universal colouring bound, or resolution of planar JC is claimed.

## 1. Useful conclusions

1. A nonzero polynomial-weighted Gaussian functional cannot annihilate a planar Keller derivation. A direct identification of a Jacobian cokernel functional with a finite Gaussian/Hermite density is therefore impossible. A true cokernel functional, **if one exists**, must escape that finite-density model.
2. Four-chromaticity of an exterior graph is not enough for six-colour forcing at the current ports. EI17, EI19, EI21 and the Hochberg–O'Donnell fish are all three-degenerate. With at most two port neighbours per exterior vertex, each extends every six-colouring when the exterior has exactly its specified edges.
3. The more selective first test is a heterogeneous list core. A nonextendible context in the current `Q(zeta55)` geometry must contain such a core with at least fourteen exterior vertices. Any additional unit distances and seam edges must be retained before applying this filter.

## 2. Inheritance and scope

The audited source snapshot was repository commit `bd2d0d6660b3889179e226172a9572351bd122bc`.

- `01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md` concerns GMC(2)/NC2 and signed charge support. It does not imply planar JC. THM-2801 and MISTAKE-237 already repair that inference.
- `01-canon/theorems/THM-4620-keller-dual-field-module-kernel-ladder-and-de-rham-cokernel.md` and `THM-4626-surjective-dual-fields-are-locally-nilpotent-and-hyp-9255-is-the-planar-jacobian-conjecture.md` identify the planar obstruction with the relative de Rham cokernel. Surjectivity is an equivalent formulation, not a newly proved hypothesis for every Keller pair.
- `01-canon/theorems/THM-4614-three-colourable-contexts-cannot-change-the-six-colour-port-relation-in-q-zeta55.md` provides the doubling theorem and fixed-port nine-apex census.
- HYP-9247, that the entire `Q(zeta55)` unit-distance graph is three-colourable, remains open. The EI19/fish and EI21 notes from 2026-10-12 distinguish a proved EI19 field obstruction from still-open realization-curve questions.
- THM-4628's high-chromatic residue quotients do not furnish high-chromatic real-plane subgraphs. The direction of the map matters.

The new tests retain the actual functional, fibre periods, port labels, forbidden lists and geometric edges. Neither support alone nor an arbitrarily completed tournament replaces those objects.

## 3. Complete Gaussian-annihilator classification

Let `A=C[x,y]`, let `X,Y` be independent **real** `N(0,1)` variables, and extend expectation complex-linearly to `A`. Put

`Q=(x^2+y^2)/2`,

`{f,a}=f_x a_y-f_y a_x`,

`D_f={f,·}`, and `Lambda_h(a)=E[h(X,Y)a(X,Y)]`.

The pairing is complex-bilinear; no conjugate occurs in its definition. It is not the holomorphic moment functional for two independent circular complex Gaussians.

### Theorem

For nonconstant `f in A` and nonzero `h in A`,

`Lambda_h(D_f a)=0 for every a in A`

if and only if

`f in C[Q] and h in C[Q]`.

Equivalently, the space of polynomial densities annihilating `D_f A` is all of `A` when `f` is constant, `C[Q]` when `f` is nonconstant radial, and `{0}` when `f` is nonradial. Here radial means polynomial in the specified centered quadratic `Q`.

### Step 1: the Gaussian adjoint

The vector field `D_f=(-f_y,f_x)` has zero divergence. Gaussian integration by parts gives

`E[h D_f a]=E[a(h D_f Q-D_f h)]`.

Polynomial growth and Gaussian decay justify the absence of a boundary term. The polynomial Gaussian pairing is nondegenerate: if `E[a r]=0` for every polynomial `a`, take the polynomial obtained by conjugating the coefficients of `r`. On real points its product with `r` is `|r|^2`, so full support of the Gaussian forces `r=0`. Therefore annihilation is equivalent to

`D_f h=h D_f Q`.                                                    (1)

### Step 2: elementary algebraic exponential obstruction

The following more general algebraic lemma supplies the mechanism.

**Lemma.** In characteristic zero, for rational functions `f,V` and a nonzero rational function `h` in `k(x,y)`, the equation

`{f,h}=h{f,V}`

forces `{f,V}=0` and `{f,h}=0`.

**Proof.** Suppose `b={f,V}!=0`. The nonzero Jacobian makes `f,V` algebraically independent. Thus `L=k(x,y)` is a finite algebraic extension of `K=k(f,V)=k(f)(V)`. Set `delta=D_f/b`. Then `delta f=0`, `delta V=1`, and `delta h=h`. On `K`, the derivation is differentiation with respect to `V`.

Let the monic minimal polynomial of `h` over `K` be

`P(T)=T^r+sum_(j<r) a_j T^j`.

Because `h!=0`, its constant coefficient `a_0` is nonzero. Differentiate `P(h)=0` and subtract `rP(h)=0`. This leaves

`sum_(j<r) (delta a_j+(j-r)a_j)h^j=0`.

Minimality forces every coefficient to vanish, so `a_0'=r a_0`. Write the nonzero rational function `a_0=A(V)/B(V)`. Then

`A'B-AB'=rAB`.

The left side has degree at most `deg A+deg B-1`; the right side has degree `deg A+deg B`, since `r!=0` in characteristic zero. Contradiction. ∎

The equivalent norm argument is `R=N_(L/K)(h)`, with `R'/R=[L:K]`, impossible for a nonzero rational `R`. The minimal-polynomial proof above avoids invoking the norm identity and makes no assumptions about connected geometric fibres or special-fibre values of `h`.

Apply the lemma to (1) with `V=Q`. It gives `D_f Q=0` and `D_f h=0`.

### Step 3: the charge-zero kernel

Set `z=x+i y`, `w=x-i y`. Then

`Q=zw/2`, and `D_Q=i(z partial_z-w partial_w)`.

A monomial `z^a w^b` has eigenvalue `i(a-b)`, so

`ker D_Q=C[zw]=C[Q]`.

Since `D_f Q=-D_Q f`, the polynomial `f` is `p(Q)`. Nonconstancy gives `p'!=0`, and `D_f=p'(Q)D_Q`. Thus `D_f h=0` forces `D_Q h=0`, hence `h in C[Q]`. Conversely radial `f,h` make both terms of the adjoint residual zero. ∎

### Keller consequence and the precise missing information

If `f` has a polynomial mate `g` with `{f,g}=c!=0`, then `f` cannot be radial: every `p(Q)` has both first derivatives zero at the origin. Hence a Keller component admits no nonzero polynomial-density Gaussian annihilator.

This does **not** assert that a nonzero cokernel exists or that it vanishes. For `f=x`, the derivation is `partial_y`, which is surjective and has no nonzero algebraic-dual annihilator. The theorem concerns the intersection of the full annihilator with the subspace represented by polynomial Gaussian densities.

In the probabilists' Hermite basis, the pairing is diagonal with entries `i!j!`. Finite Hermite support is therefore equivalent to representation by a polynomial density. Consequently, **if** a nonzero algebraic-dual annihilator of a Keller derivation exists, it has infinitely many nonzero Hermite coordinates. No convergence, distribution class or integrability is implied by this formal statement.

The closest repository engine is the rational eigenvector obstruction in THM-4620(2), applied here after rational rescaling into `(f,Q)`. The algebraic exponential obstruction is classical Liouville/Rosenlicht material; see Rosenlicht, [Liouville's theorem on functions with elementary integrals](https://msp.org/pjm/1968/24-1/pjm-v24-n1-p16-s.pdf), Pacific J. Math. 24 (1968), 153–161. The self-contained argument above is what is proved here. A limited search not finding this exact Gaussian formulation does not establish novelty.

### Finite exact test bound

If `deg f=n` and `deg h<=m`, the adjoint residual has degree at most `n+m`. Thus testing `E[hD_f a]=0` on every monomial of degree at most `n+m` suffices for this fixed density cap, by nondegeneracy of the degree-bounded Hermite pairing.

The recovered script checks coefficient kernels for `d=0,2,4,6,8` and six polynomials:

- `x`, `x+y^2`, `x+(y+x^2)^2`, and `x+x^2y`: nullity zero;
- `Q` and `Q^2+Q`: nullity `floor(d/2)+1`.

For each of the same six examples at density degree six, an independent matrix of direct Gaussian moments agrees entry-by-entry with the integration-by-parts matrix and has the predicted rank. These are finite controls; the theorem is the algebraic proof.

## 4. An explicit period discarded by finite Gaussian densities

Take `f=x+x^2 y`. It is nonsingular, since `f_y=x^2` and `f_x=1+2xy` never vanish together. For `c!=0`, its fibre is `C*` with

`x=t`, `y=(c-t)/t^2`.

Here `D_f=-t^2 d/dt`. Define

`ell_c(a)=Res_(t=0)[-a(t,(c-t)/t^2)dt/t^2]`.

Then `ell_c(D_f a)=Res da=0`, but `ell_c(1)=0` and `ell_c(x)=-1`. This is a genuine nonzero residue/period functional despite the absence of any polynomial Gaussian-density annihilator. Its puncture is the missing information.

For `a,b>=0`, set `k=2b+1-a`. Expansion gives

`ell_c(x^a y^b)=(-1)^(k+1) binom(b,k)c^(a-b-1)` if `0<=k<=b`, and zero otherwise.

The formula is valid in all degrees; the script checks it and the annihilator recurrence for all forty-five monomials of total degree at most eight at `c=2`.

The derivation has rational slice `1/x`, but has no polynomial slice. Indeed the function field is `C(f)(x)` and `D_f=-x^2 partial_x`. Thus a hypothetical polynomial `g` with `D_f g=1` would be `g=1/x+r(f)`, with `r` rational. A finite pole of `r` at `c!=0` gives a pole on a component of `f=c`, where `1/x` is regular. A pole at zero gives a pole on `1+xy=0`, where `1/x` is also regular. Therefore `r` has no finite poles and is polynomial, contradicting polynomiality of `g-r(f)=1/x`.

This direct argument is necessary: positive generic first Betti number by itself would not rule out a polynomial mate for an arbitrary polynomial without resolving an open Jacobian question.

## 5. The six-colour heterogeneous list core

Let a finite graph `G` have fixed ports `B` with a proper six-colouring `tau`. Give exterior vertex `v` the list

`L_tau(v)={1,...,6} minus {tau(b): b in N(v) intersect B}`.

Repeatedly remove a vertex whose current exterior degree is strictly smaller than its list size. If all exterior vertices are removed, reverse the order and greedily colour, obtaining an extension of `tau`.

**Corollary.** A three-degenerate exterior in which every vertex has at most two port neighbours extends every proper six-colouring of the ports: all lists have at least four colours.

Independent exact checks give degeneracy three and chromatic number four for EI17 `(17,31)`, EI19 `(19,35)`, EI21 `(21,38)` and the fish `(23,42)`, where each pair denotes vertices and edges. The JSON stores explicit peel certificates and four-colourings. Exact three-colour backtracking is UNSAT, with respectively 94, 65, 307 and 116 recursive nodes under the recorded deterministic symmetry convention.

**Extra-edge caveat.** A unit-distance realization can have extra unit edges beyond an abstract edge list. The corollary applies only when the complete exterior unit-distance graph remains three-degenerate. Additional geometric edges, overlapping copies and seam edges must be included. No conclusion is asserted for every realization or every union of these graphs. A field-valued realization can refute HYP-9247 while still failing to provide six-colour forcing.

### Fourteen-vertex necessary bound at the current eight ports

Use precisely THM-4614's eight ports in `L=Q(zeta55)`. Its inherited exact census says every exterior vertex has at most two port neighbours; at most nine have two; those nine are mutually nonadjacent. The field unit graph is triangle-free, and two distinct planar points have at most two common unit neighbours.

If a six-colouring does not extend, take an inclusion-minimal non-list-colourable exterior induced graph `C`. Let `n=|C|` and let `s<=9` count the double-port vertices. Minimality gives exterior degrees at least four at those vertices and at least five elsewhere.

Because common-neighbour pairs have nonadjacent endpoints in a triangle-free graph,

`sum_v binom(d_v,2)<=2(binom(n,2)-e)`.

Equivalently,

`sum_v d_v(d_v+1)<=2n(n-1)`.

The degree bounds imply `30n-10s<=2n(n-1)`. With `s<=min(n,9)`, this excludes every `1<=n<=12`. At `n=13`, it requires `s>=8`. The marked vertices are independent, so their at-least-four neighbours belong to the `t=n-s<=5` unmarked vertices. Counting pairs of these neighbours gives

`6s<=sum_(marked v) binom(d_v,2)<=2 binom(t,2)<=20`,

contradicting `s>=8`. Hence `n>=14`. At fourteen, the two necessary inequalities allow only `s=6` or `s=7`; this is not an existence claim. The recovered script verifies that integer boundary separately.

The proof is conditional on the specified census. It is not a theorem about arbitrary ports or arbitrary triangle-free graphs.

## 6. Regular-pentagon Mycielski shortcut fails

Let `v_i` be a regular unit-edge pentagon with centre `O`. The common unit neighbours of `v_(i-1),v_(i+1)` are `v_i` and its reflected position. A faithful shadow must use the latter:

`u_i=O-phi^(-2)(v_i-O)`, where `phi=(1+sqrt(5))/2`.

These shadows have circumradius squared `(5-2sqrt(5))/5`, which is not one. A point equidistant from three noncollinear shadows is their centre `O`, which is not a common unit neighbour. Thus this particular regular-pentagon shadow construction has no common unit apex, even over the full real plane. This does not exclude flexible five-cycles or other Mycielski-like constructions.

## 7. Exact connection contract and next targets

**GMC to JC:** the map `h -> Lambda_h` takes polynomial Gaussian densities into the algebraic dual of `A`. Integration by parts preserves the precise annihilator equation. The theorem shows that this map contributes no nonzero Keller cokernel functional. Puncture periods and infinite Hermite tails are lost; the residue example exhibits the loss explicitly.

For `f=sum c_ij x^i y^j`, the full dual moments satisfy

`sum_(i,j) c_ij(i b-j a)m_(a+i-1,b+j-1)=0` for every `a,b>=0`,

with zero-coefficient negative-index terms omitted. A polynomial slice `g` also forces `ell(g^m)=0` for all `m>=0`. The next valid target is an explicitly justified tail or period condition on this recurrence. Finite truncations without a boundary rule can create spurious solutions.

**Tournament methods to HN:** a degeneracy order gives an acyclic orientation of actual unit edges. It is a partial orientation, not a tournament. It preserves a sufficient list-extension certificate while requiring the original graph, geometry and port labels. Arbitrarily orienting nonedges loses the target predicate. A proposed tournament invariant needs an exact operation on the boundary-colouring relation before it can help.

**First HN computational target:** compute all exact unit edges and heterogeneous list cores for the existing 570 port-colouring patterns. Empty cores prove extension cheaply. For a nonempty core, retain coordinates, port incidence and actual forbidden lists, then seek a certified eliminated pattern. Four-chromaticity or a nonempty core alone is intermediate.

The EI21/fish realization curves remain legitimate arithmetic targets. A successful field-valued point supplies a four-chromatic seed, but its complete unit graph and intended gluing must still pass the list-core test.

## 8. Recovered artifacts and reproducibility

The accompanying script is self-contained and uses no network or repository reads. The four adjacency lists were recovered from retained outputs of the source snapshot's `ei19_fish_placement_q_zeta55_20261012.py` and `ei21_placement_q_zeta55_20261012.py`. Counts, symmetry, triangle-freeness, colour witnesses and peeling certificates are checked afresh.

Run:

```sh
python3 liouville_planar_obstructions_20261014.py > result.json
python3 -O liouville_planar_obstructions_20261014.py > result-optimized.json
cmp result.json result-optimized.json
```

Recovery runs used SymPy 1.14.0; both outputs agree byte-for-byte. Explicit exceptions, rather than `assert`, enforce checks. Newly computed recovery hashes:

- Script SHA-256: `0725079d246e7ac5a7bb3d95710730af46bfd0f48b95df6014ec776d77d8e909`
- JSON SHA-256: `a00b023992fdbe4529d726f4e9b581f9bfb0b5a24608d7502c01fabd647446f9`

These replace no original artifacts until the parent reconciles them with the current repository. No public write, commit or push was made during recovery.
