# Session31: what the CM / honeycomb / tensor-category bridge preserves

Date: 2026-10-22. Scope: the user-defined family

\[
 A(m)=\frac{2}{\pi}m^{-1/m}\sin(\pi/m),\qquad
 C_{\rm MRB}=\sum_{r\ge1}(-1)^r(r^{1/r}-1),
\]

and the question whether the special geometry at `m=8` explains the actual
MRB defect `A(8)-C_MRB`, rather than only the nearby scale `A(8)`.

**Bottom line.** There is a precise scalar normalization through the modular
`S` matrix of `SU(2)_(m-2)`, and a precise Chebyshev relation to the honeycomb
loop parameter. Neither supplies an MRB character, a tensor functor, or a
formula for the defect. Several literal candidate identifications are
provably impossible. The remaining possibility is an *enriched analytic
construction* with a specified observable and summation/contour functional,
not an automatic consequence of Kazhdan–Lusztig equivalence.

Statuses used here: **PROVED** means the elementary arguments below;
**CITED** means the specified theorem in the primary literature;
**PHYSICS/CONDITIONAL** means the conformal interpretation, not a new rigorous
lattice scaling-limit theorem; **OPEN** means a missing construction. No
claim of novelty is made for the elementary observations.

## 1. Inheritance and the live concept board

Bounded startup read: repository `AGENTS.md`, startup packet, current frontier,
research protocol, guardrails, and the core-papers introduction. The startup
packet reported a clean cached checkout at `b44cbf65db89`. No Git mutation was
performed. The coordination status tool could not identify this worker because
`.machine-id` is absent; no registration was invented.

Closest existing positive mechanism:
`05-knowledge/results/factorization_vertex_supersymmetric_bridge_20261018.md`,
§4, supplies a *named change of variable* from the two Lee–Yang characters to
Beukers' Hauptmodul, then an identity for the Apéry generating function. That
is the appropriate standard for a real bridge: an entire function/series,
not an equality at one special number.

Hostile/corrected near miss: `MISTAKE-658` explicitly rejects using matching
arithmetic or modular vocabulary in place of a target-preserving map. The
same inherited factorization note has corrected overclaims and labels its
CM interpretation as analogy where no map has been proved.

Least-used relevant sidecar:
`05-knowledge/reference/CORE-PAPERS-VERTEX-ALGEBRA-CLASSICAL-LIMIT-2026-10-18.md`
records the distinction between a VOA, its associated variety/arc space, and
its character. Classical freeness or a Hilbert-series coincidence is not
an equivalence of representation categories. The available
`schur_index_vertex_algebra_mlde_20261018.{json,out}` records finite MLDE
checks, not an MRB construction. The corresponding `.py` at the parent-sent
path was absent, and these old computations were not rerun here.

Live board:
1. `A(m)`: algebraic trigonometric factor plus external `1/pi` normalization.
2. Honeycomb loop model: loop fugacity, edge fugacity, boundary sectors, phase.
3. Root-of-unity quantum category: fusion, trace, braiding, negligible objects.
4. Elliptic CM point: complex structure, nome convention, periods, missing
   Kähler modulus if interpreted as a torus CFT.
5. MRB functional: unbounded-index radicals and alternating accumulation.

## 2. Parameter dictionary: at least four different q's

For integer `m>=3`, put `theta=pi/m`.

| Symbol | Definition | Role |
|---|---|---|
| `q_WZW` | `exp(pi i/m)` | Quantum parameter for `SU(2)_(m-2)` |
| `q_loop` | `-exp(4 pi i/m)` | Convention with loop value `q_loop+q_loop^(-1)=n` |
| `n` | `-2 cos(4 pi/m)` | Weight of each closed loop in the O(n) model |
| `q_theta` | `exp(pi i tau)` | Theta-function elliptic nome |
| `Q` | `exp(2 pi i tau)=q_theta^2` | Usual VOA character/Dedekind-eta nome |

The alternative loop convention `n=-(q+q^(-1))` uses
`q=exp(4 pi i/m)` instead. Signs and inverse braidings must be stated, not
silently identified.

At `m=8`,

\[
 q_{\rm WZW}=e^{\pi i/8},\quad q_{\rm loop}=-i,\quad n=0,
 \quad \mu=2\cos(\pi/8)=\sqrt{2+\sqrt2}.
\]

At the separate CM point `tau=i sqrt(2)`,

\[
 q_\theta=e^{-\pi\sqrt2}=0.0117619805313891216891\ldots,
 \qquad Q=e^{-2\pi\sqrt2}=0.0001383441860207767254\ldots.
\]

Neither positive real nome is a root of unity. In fact `q_theta` is
transcendental by Gelfond–Schneider applied to the value
`(-1)^(i sqrt(2))` with logarithm `i pi`. Consequently one cannot identify it
with either algebraic root-of-unity quantum parameter. A generic quantum group
*can* be defined at this real parameter. Its usual finite-dimensional type-1
module category has infinitely many simple objects; it is not the finite
`SU(2)_6` modular category.

A claimed `q_theta^4` correction is a `Q^2` correction. More importantly, an
order of vanishing requires an analytic family near a cusp. From one number
`D=A(8)-C_MRB` at one fixed `q_theta`, writing `D=(D/q_theta^4)q_theta^4`
establishes no leading modular term.

## 3. The exact scalar bridge, and its loss ledger

**PROVED, using the standard SU(2) modular data.** Let `k=m-2`. With simple
labels `a,b=0,...,k`, the standard normalization is

\[
 S_{ab}=\sqrt{\frac2m}\sin\frac{(a+1)(b+1)\pi}{m},\qquad
 d_a=\frac{S_{a0}}{S_{00}},\qquad d_1=2\cos\frac\pi m.
\]

The root-of-unity quantum Weyl formula and `S00` normalization are documented
in [Coquereaux, §§2.1–2.2](https://arxiv.org/abs/1003.2589); the explicit
`SU(2)` total-dimension normalization also appears in
[Kitaev–Preskill, eq. (4)](https://arxiv.org/abs/hep-th/0510092).
Substitution, not a new representation theorem, gives

\[
 \boxed{A(m)=\frac{\sqrt{2m}}{\pi}\,m^{-1/m}S_{00}}.
\]

Equivalently, with `mu_m=2cos(pi/m)`,

\[
 A(m)=\frac{m^{-1/m}}\pi\sqrt{4-\mu_m^2},\qquad
 \boxed{n=-\mu_m^4+4\mu_m^2-2}.
\]

The second equality is `n=-2 T_4(mu_m/2)`.
For real `m>=4`, `cos(2pi/m)>=0`, so the Nienhuis candidate specializes to

\[
 x_{\rm N}(n)=\bigl(2+\sqrt{2-n}\bigr)^{-1/2}
            =\frac1{2\cos(\pi/m)}=\mu_m^{-1}.
\]

The algebraic identity is unconditional; interpreting `x_N` as the critical
point for every parameter is not proved by this algebra. At `m=8`, the exact
honeycomb SAW connective constant is **CITED**, namely Duminil-Copin–Smirnov,
[Theorem 1](https://arxiv.org/abs/1007.0575). Their key observable has marked
endpoints and winding phase, with spin `5/8` (Lemma 1). Their theorem is not
a construction of an MRB observable or a proof of the full SAW scaling limit.
For broader loop parameters, distinguish the exact scalar formula from
[Duminil-Copin–Glazman–Peled–Spinka's macroscopic-loop theorem for n in [1,2]](https://arxiv.org/abs/1707.09335).

Connection contract:
- Source: the explicitly indexed scalar family `m`.
- Target: `SU(2)_(m-2)` modular datum and the Nienhuis scalar parameter curve.
- Map: the displayed substitutions.
- Preserved: these scalar identities and the connective constant at `m=8`.
- Lost/not supplied: loop states, edge lengths, fusion products, marked
  endpoints, winding, the MRB index sequence, and its alternating functional.
- Required sidecar: a concrete observable and an intertwining identity for
  its whole generating function or contour integral.
- Cheapest hostile: compare the value of a closed loop, not just `mu`.

## 4. Exact no-go propositions

### 4.1 A(m) is not a fusion-category dimension

For every integer `m>=2`, `2m^(-1/m)sin(pi/m)` is a nonzero algebraic number.
Thus **A(m) is transcendental**, because a nonzero algebraic number divided
by transcendental `pi` cannot be algebraic.

In a pivotal fusion category, categorical dimensions give a ring character
of the finite Grothendieck ring. If `N_X` is the integer fusion matrix for
an object `X` and `d` is the nonzero vector of dimensions of simple objects,
then `N_X d=d(X)d` (up to the harmless transpose convention). Therefore
`d(X)` is an eigenvalue of an integer matrix, hence an **algebraic integer**.
Frobenius–Perron dimensions satisfy the same conclusion.

Consequently `A(m)` cannot literally be such a dimension. This does not
prohibit an analytic prefactor times a dimension, an evaluated character,
or a trace of an arbitrarily assigned complex scalar. The latter can encode
any number tautologically and provides no explanation.

### 4.2 Even pi*A(8) is not a fusion-category dimension

Put

\[
 \beta=\pi A(8)=2^{-3/8}\sqrt{2-\sqrt2}.
\]

Then `beta^4=3/sqrt(2)-2`, so

\[
 2\beta^8+8\beta^4-1=0.
\]

The reciprocal polynomial `t^8-8t^4-2` is Eisenstein at `2`.
Hence `beta` has degree `8`, with primitive irreducible polynomial
`2x^8+8x^4-1`; its monic minimal polynomial has a coefficient `-1/2`.
So **beta is not an algebraic integer**. Removing `pi` alone is insufficient
to produce a fusion-category dimension. This is an exact obstruction,
independently checked by symbolic minimal-polynomial computation.

### 4.3 The m=8 loop strand is not the SU(2)_6 fundamental object

In the complex linear Temperley–Lieb category `TL_delta`, the generating
self-dual strand `X` has categorical loop value `dim(X)=delta`. A complex
linear pivotal tensor functor preserves this value. At `n=delta=0`, it cannot
send `X` to the `SU(2)_6` fundamental object, whose dimension is
`sqrt(2+sqrt(2))`, not zero.

More generally, no such functor into a unitary fusion category equipped
with its canonical positive spherical pivotal structure can give `X` a
nonzero image: every nonzero object then has positive dimension. The
zero-image quotient is allowed, but destroys the polymer strand.

The failure is visible at two strands:

\[
 TL_2(0)=\mathbb C[e]/(e^2),\qquad e\ne0.
\]

This is a nonsemisimple algebra. In the usual diagrammatic negligible-morphism
quotient, `id_X` is negligible (`End(X)=C` and its trace is zero), so the
category generated by `X` collapses to the unit sector after semisimplification.
One cannot discard this loss and still claim to preserve the marked polymer
observable. Dilute models also retain vacancies; their full diagram algebra
is a dilute TL/Motzkin version, not just ordinary TL. See the explicit dilute
versus dense discussion in
[Jacobsen–Ribault–Saleur, §6, “Dilute O(n) model”](https://arxiv.org/abs/2208.14298).

A further simple warning: `SU(2)_6` has fusion graph `A7`. Its spectral radius
is indeed `mu`, but endpoint path counts start `1,1,2,3,6,10,...`; rooted
honeycomb SAWs start `1,3,6,12,...`. Equality of exponential growth rates is
not equality of the underlying counting problems.

### 4.4 Direct ordinary/finite-order twined character encoding fails

Let `g_r=r^(1/r)-1`. Consider a grading bounded below, with finite
multiplicity in each grade, and strip off the fixed initial character power.
No uniform bound on the dimensions across grades is needed. An ordinary
graded character then has nonnegative integral coefficients, whereas
`g_2=sqrt(2)-1` is not integral.

More strongly, let one grading-preserving automorphism act on the whole
graded space, with one fixed finite order `N` independent of the grade. Its
trace on every finite-dimensional grade lies in the one cyclotomic field
`Q(zeta_N)`. But for every prime `p`, `g_p` has degree exactly `p`: the
polynomial of `g_p+1=p^(1/p)` is Eisenstein `x^p-p`. Taking a prime exceeding
`[Q(zeta_N):Q]` is a contradiction. The same argument permits signs; for
supertraces it additionally requires that the twining automorphism preserve
parity, so its restrictions to the even and odd parts have order dividing
`N`. Thus the full raw MRB coefficient sequence cannot be the coefficients
of one such finite-order twined character.

The parity hypothesis is essential. On a single grade of dimension `1|1`,
the matrix `G(a)=[[a,1],[1-a^2,-a]]` satisfies `G(a)^2=I`, but its
supertrace is `2a`, which can be arbitrary. This matrix is not
parity-preserving. Ordinary traces are unaffected by this counterexample.

Scope: this rules out a *coefficientwise identification*, not an equality of
one final number with a special value, or an integral transform of characters.
An operator with freely prescribed infinite-order eigenvalues could reproduce
the coefficients by design; that would not yet explain them.

### 4.5 A finite-order nome change does not remove that obstruction

Replacing the character variable by `z^r`, shifting by a fixed power, or
interleaving finitely many characters still leaves all coefficients in one
fixed number field. The prime-degree argument persists. Any proposed
character encoding must therefore exhibit a genuinely nontrivial transform
or an unbounded family, not merely rename the nome or reindex finitely.

## 5. Which categorical theorems are actually available?

### 5.1 Generic Drinfeld/KZ equivalence

For a complex simple Lie algebra and `hbar` outside the nonzero rational
parameters, the Drinfeld category has finite-dimensional Lie-algebra modules,
ordinary tensor product, KZ associator, and KZ braiding. It is braided-tensor
equivalent to finite-dimensional type-1/admissible quantum-group modules at
`q=exp(pi i hbar)`. The functor preserves highest-weight labels and can be
compatible with the underlying vector-space functor.
[Neshveyev–Tuset, Theorem 3.3 and its following extension statement](https://arxiv.org/abs/0711.4302).
This is genuine transport of tensor structure and braid-group monodromy,
not transport of arbitrary scalar sequences.

In particular one can take `hbar=i sqrt(2)` to manufacture the real nome as
a generic quantum parameter. That is legal algebraically, but does not turn
the elliptic CM torus into the KZ local system, nor make this infinite
category the finite root-of-unity category.

### 5.2 Affine Kazhdan–Lusztig equivalence

For the present rank-one comparison, use `sl2`, affine level `k`, shifted
parameter `k+2`, and `q=exp(pi i/(k+2))`. The rank-one KL condition is
`k+2 not in Q_{>=0}`. This includes negative shifted levels and uses
appropriate KL/ordinary affine modules; its
quantum side uses Lusztig's **big** quantum group, not an unspecified small
quantum group or its semisimplification. See the parameter statement in
[Creutzig–Ridout–Rupert, §1.2](https://arxiv.org/abs/2112.13167), and the
`G(O)`-integrability formulation in
[Gaitsgory, §0.1.1](https://arxiv.org/abs/1810.09054).
The latter author's *extension* to a larger Iwahori-integrable setting is
explicitly conjectural in that paper; do not promote it with the original theorem.

### 5.3 Positive-integral WZW modular categories

The finite category `Rep L_k(sl2)` at positive integer `k` is the rational WZW
case. Its quantum model is the **semisimplified tilting** category at
`q=exp(pi i/(k+2))`, not the full category of root-of-unity modules.
For `k=6`, there are seven simple objects and affine central charge `9/4`.
This explains the formula in §3, but it is not the c=0 polymer theory.

A general VOA-to-modular-category theorem requires actual hypotheses:
simplicity, CFT-type/self-contragredience, complete reducibility, and
`C2`-cofiniteness in the standard setting. See
[Huang, Rigidity and modularity of vertex tensor categories](https://arxiv.org/abs/math/0502533).
“Modules of a vertex algebra form an MTC” without those restrictions is false.

### 5.4 The c=0 Virasoro obstruction and the correct logarithmic carrier

The simple Virasoro VOA `L(0,0)`, the `(2,3)` minimal model, is trivial.
Indeed `L_-2|0>` is singular in the universal c=0 vacuum module: the positive
Virasoro modes kill it, since `L_-1|0>=0` and `c=0`. The simple quotient
kills this conformal vector; the Virasoro VOA is generated by that vector,
so the quotient is just `C`. Its ordinary finite-dimensional module category
is `Vec`, not a nontrivial SAW category.

This does **not** say every VOA or CFT with c=0 is trivial. Retaining
nontrivial polymer fields requires a larger/nonsemisimple structure rather
than this simple Virasoro quotient. A rigorous available starting point is
`O_c^fin`, the finite-length modules of the universal Virasoro VOA with the
specified reducible-Verma simple factors, equivalently lower-bounded
`C1`-cofinite generalized modules. It has a braided tensor structure for all
`c`; the cited theorem alone does not assert finite semisimple modularity.
[Creutzig–Jiang–Orosz Hunziker–Ridout–Yang, Theorem 4.2.6](https://arxiv.org/abs/2002.03180).

For `m>=4`, the physical dilute branch has `g=2-4/m` and
`c=1-6(g-1)^2/g`, giving c=0 at m=8. The dense branch uses `g=4/m`,
giving c=-2 at the same m=8. At m=4 the branches meet at g=1; this labeling
is not asserted for m=3, where `2-4/m=2/3` and `4/m=4/3` exchange the
dense/dilute ranges. These branch formulas and their integrable-model
context are in [Zhou–Batchelor, §2.1](https://arxiv.org/abs/cond-mat/9611156).
They are not supplied by the bare complex number `q_loop`.

This distinction survives in explicit algebraic physics:
[Read–Saleur, §§3,5](https://arxiv.org/abs/hep-th/0701117) compare dense
c=-2 and dilute c=0 polymers at the same quantum-group parameter `q=i`,
with different vacuum/continuum data. Their dilute lattice Morita-equivalence
suggestion is explicitly a conjecture in §5, not a theorem. The familiar
triplet `W(2)` correspondence at q=i has c=-2; importing it as the c=0 SAW
category would change the problem. “Nonsemisimple modular category” also
uses a generalized notion, not the ordinary semisimple MTC hypothesis above.

## 6. A second loss: specializing the vacuum partition function at n=0

For a finite graph with no open strands, define the dilute loop partition sum
`Z_G(n,x)=sum_config n^(number of loops) x^(occupied edges)`.
At n=0, its only surviving configuration is the empty one, so `Z_G(0,x)=1`.
The derivative in n at zero retains single-loop data, while a two-point
sector with a marked open strand retains SAW data. Those are additional
observables, not recoverable from the bare constant vacuum value.

This is the concrete reason to demand endpoint/defect-insertion and trace
sidecars before trying to transport a statistical observable through a
category. It is also an elementary model of how semisimplification can erase
the information one hoped to explain.

## 7. What CM does and does not supply

The number `tau=i sqrt(2)` determines a CM elliptic curve and specialized
elliptic/theta quantities. It does not specify a quantum-group representation,
a module label, or an MRB summation functional.

There is a genuine but different CM-to-CFT theorem: for a torus sigma model,
rationality requires both complex structure `tau` and complexified Kähler
parameter `rho` in the same imaginary quadratic field.
[Gukov–Vafa, §5](https://arxiv.org/abs/hep-th/0203213).
A CM value of `tau` alone is insufficient. This was already correctly
flagged in the inherited repository's factorization bridge note.

Even after choosing `rho`, one must identify the relevant graded trace and
its evaluation. An equality involving `A(8)` does not identify the separate
constant `A(8)-C_MRB`. No source inspected here constructs that missing map.
This is a bounded research conclusion, not a theorem that no future map exists.

## 8. Concrete next tests, with stopping rules

1. **Normalized modular datum (closed scalar task).** The exact formula in
   §3 is the valid endpoint. Do not label it a functor or MRB evaluation.
   A proposal saying `A8` or `pi*A8` is a fusion dimension is refuted by §4.

2. **Whole-function MRB-to-character candidate (OPEN).** Specify a transform
   `T`, a VOA/module/inserted operator, and an equality
   `T[sum g_r z^r]=specified character or period` on a nonempty domain.
   The transform must not contain the unknown defect as a fitted constant.
   First test coefficient fields (§4.4), then singularities and monodromy,
   then independently held-out coefficients. Stop or repair on the first
   violated invariant. Agreement at a single CM point is insufficient.

3. **Logarithmic lattice candidate (OPEN, more faithful than SU(2)_6).**
   Work with dilute link patterns and their marked one-strand sector at
   loop value zero. Retain edge fugacity, boundary data, winding, and
   nonsemisimple extension/Jordan data. Any proposed morphism must intertwine
   the transfer operator and the chosen trace/marked observable. The
   two-strand nilpotent and the closed-loop test are cheap hostile controls.
   Merely transporting the abstract q=i quantum-group symmetry does not fix
   the dense/dilute branch.

4. **CM defect expansion (OPEN).** Before fitting powers of the nome,
   define an analytic `D(tau)` with `D(i sqrt2)=A8-C_MRB` from an independent
   construction. Determine its modular group, weight/multiplier, allowed
   poles, and coefficient field. Only then test a cusp expansion and CM
   evaluation. Without a family, a “q^4 scale” is just a numerical scale.

5. **Compare the inherited successful template.** The Lee–Yang/Apéry lane
   gives `t(tau)=(chi_0/chi_-1/5)^5` and an entire generating-function identity.
   An MRB bridge should reach that level of specification. Classical freeness,
   an MLDE fit, or matching a fusion spectral radius alone does not reach it.

## 9. Reproduction and limits

Run `python3 check_categorical_parameters.py > categorical_checks.json` in
this directory. It performs exact minimal-polynomial and Eisenstein checks,
checks the `A7` characteristic polynomial and a nonzero nilpotent in
`TL2(0)`, and checks the three scalar identities at 90-digit precision for
integer m=4,...,40. The general identities are proved above; the finite
numerical checks are controls, not their proof.

`A(8)=0.1878596613801383990492532045996126144757...`.
The MRB constant is not recomputed or certified here, and no arithmetic claim
about it or its defect is made. This note does not prove SAW conformal
invariance, a logarithmic categorical equivalence for dilute polymers, or an
MRB modular representation. It isolates exactly what remains to be built.
