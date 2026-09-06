# Audit of the actual cusp passport and its global consumer

**Status: ANALYTIC, PRIMARY-STATEMENT, AND EXACT SOURCE AUDIT PASS.**
September 6, 2026. This checks
[the cusp-passport note](planar_jc48_sep06_cusp_passport.md), the complete
producer, and normal/optimized/frozen replay. The audit supplied the
same-point re-access construction and the extension from five nodes to any
positive number of nodes; both are now explicit in the primary proof.
The final text and its exact scopes are accepted without further correction.

The conclusion excludes a **whole irreducible nonproperness support** with
normalization `A¹`, exactly one ordinary cusp, and `N>=1` ordinary nodes,
with no other singularities, for an actual polynomial Keller map over `C`.
It applies in arbitrary coordinate polynomial degree. It neither establishes
entry into this support class nor proves `JC(2)`.

## 1. The actual transfer, rather than an abstract conjugacy assumption

The primary local input was read in the downloaded author PDF of
[Shimada's lectures](https://www.math.sci.hiroshima-u.ac.jp/shimada/LectureNotes/LNZV.pdf),
Section 6, pp. 24–25, together with the braid presentation in Section 4.
The presentation uses actual fibre meridians. Substituting `p=3,q=2` in
`ell_(j+2)=m ell_j m^-1`, `ell_(j+3)=ell_j`, and `m=ell_1 ell_0`
gives the stated two-meridian braid relation. The local generators are not
being inferred solely from an abstract torus-knot group isomorphism.

Choose a common labelled fibre `Omega` over the local cusp complement.
The actual retained generic pages along a normal transversal give
`A subset Fix(sigma)`, with the globally constant cardinality `a` from the
smooth-support lemma. That cardinality counts actual affine points, not
all inertia-fixed labels. Let `g=sigma tau`; the braid relation gives
`g sigma g^-1=tau`.

The final convention in the primary note is internally consistent:
compose paths with the rightmost traversed first, let `c` access the normal
transversal, and replace it by `c'=c g^-1`. If `A_eta` is the actual endpoint
subset, then

```text
T_(c')^-1 A_eta = g T_c^-1 A_eta = gA,
T_((c')^-1 m c') = g sigma g^-1 = tau.
```

Thus the second subset is actual by **re-access of the same smooth point**.
There is no assertion that `g` is a preferred half-orbit along the cusp,
and no peripheral-centralizer identification is needed. The same operation
transports both the subset and its meridian, which is essential.

If a label lies in `A intersect gA`, it is fixed by `sigma` and `tau` and
hence by the entire local group. Its orbit is a singleton; uniqueness of
the finite normal extension gives a degree-one local sheet, isomorphic to
the smooth target neighborhood. Its generic cusp divisor is retained
because the label lies in the actual `A`. A deleted divisor in this copy
would have to be the whole irreducible cusp divisor; that is impossible.
The independently audited pure-boundary argument then prevents an isolated
deletion of the cusp point. This yields an actual affine point over the
cusp. Conversely, an actual cusp point is etale and supplies a singleton
local inverse sheet whose generic cusp page is retained. Therefore

```text
n_cusp = |A intersect gA|.
```

The purity hypothesis is load-bearing. A punctured degree-one local sheet
would otherwise defeat the reverse implication. The actual-subset
hypothesis is equally necessary: a deleted unramified sheet can be fixed by
both generators without being present in the affine source. Neither local
nor global transitivity is assumed in this argument.

## 2. Exact ledger and the all-degree sparse lemma

Write `delta=d-a`. The global source and target each have compactly
supported Euler characteristic one. Removing the unique cusp preimage and
the `2N` node preimages from the `A¹` normalization gives
`chi_c(S_smooth)=-2N`; the target complement has characteristic `N`.
At every ordinary node the already proved actual-subset formula gives
`d-2delta+omega_p` actual affine points, where
`omega_p>=max(0,d-2a)`. Hence

```text
1=dN-2aN+sum_p(d-2delta+omega_p)+n
 =n+sum_p omega_p.
```

No cusp specialization formula was used to obtain this ledger. In
particular `n` is zero or one. Choosing any one of the `N>=1` nodes gives
`d-2a<=omega_p<=1-n`. Consequently the actual outside set has size

```text
h=|Omega\(A union gA)|=d-2a+n<=1.
```

I independently derived the complete abstract classification. Put
`I=A intersect gA`, `U=A\gA`, `V=gA\A`, and `O=Omega\(A union gA)`.
The intersection `I` is fixed pointwise by both generators, and `gU=V`.

If `O` is empty, `sigma` is supported on `V` and `tau` on `U`. They
commute, and the braid relation forces them to be equal; disjoint supports
then make both identities. This gives `A=gA=Omega`, violating `a<d`.

If `O={o}`, `sigma` fixes `I union U` and `tau` fixes `I union V`.
For `x in U`, the image `tau x` lies in `U union {o}`. If it lay in `U`,
then `g x=sigma tau x=tau x` would lie in `U`, contradicting `gU=V`.
Thus `tau U subset {o}` and `|U|<=1`. This gives the same conclusion as
the primary note's equivalent braid-evaluation argument. Bijectivity
determines both permutations. The only possibilities with `a<d` are

| `d,a,n` | Complete local type |
|---|---|
| `2,1,1` | two identity generators; one common retained label |
| `3,1,0` | two adjacent transpositions on three labels |
| `4,2,1` | the three-label action plus one common fixed label |

All three types occur as abstract permutation passports. The four-sheet
control is therefore indispensable: dropping its fixed label by imposing
transitivity would produce a false strengthening of the local argument.
For `N>=2`, the degree-three passport already fails the ledger, because
each node would require an overlap of at least one. No finite enumeration
is needed for the all-degree classification.

## 3. Correctly typed classical finish and application

I independently opened the two primary PDFs supplied by root.
[Orevkov, Theorem 1.1](https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf),
p. 1, excludes generic mapping multiplicities two and three for a
polynomial map with nonzero constant Jacobian. The definition immediately
before the theorem is the generic number of preimages.
[Domrina's general-case theorem](https://www.mathnet.ru/php/getFT.phtml?jrnid=im&option_lang=eng&paperid=273&what=fullteng),
p. 1, excludes topological degree four and explicitly uses the same
generic-preimage definition. These are **CITED established theorems**;
this audit checked their precise statements, not their full proofs anew.

In characteristic zero this degree is exactly the finite function-field
degree `d` used by the passport. The result does not substitute coordinate
degree, curve degree, or intrinsic parametrization degree for it. Domrina's
general theorem has no additional one-dicritical assumption. The three
possibilities above are therefore excluded.

The zero-node case is outside the ledger-to-outside-set step because that
step chooses a node. Its separate no-self-intersection endpoint is already
available from the correctly versioned Jelonek theorem; the primary
statement here safely retains `N>=1`.

For `U=t^4+t`, `V=t^6+t^2+lambda t`, the separately audited classification
shows that the remaining three roots of
`128lambda^3-288lambda-283` have exactly one cusp and five nodes. This
consumer closes those three cases. Together with the existing nodal,
triple-point, and tacnode exclusions, every member of this **explicit**
one-parameter family is excluded as a sole support. The application does
not classify or exclude all intrinsic `(4,6)` curves.

## 4. Exact source and replay

The full source was read. Its permutation universe is all ordered pairs
in `S_d`, `2<=d<=5`, followed by all nonempty proper subsets of the first
generator's fixed set. Products act on the left, right factor first. The
declared filters are the braid relation, then `n<=1,h<=1`. There is no
unannounced quotient or transitivity filter. The structural labelled counts
are independently immediate: two identity choices, six ordered adjacent
transposition choices on three labels, and four choices of fixed fourth
label times six choices on the remaining three labels.

The second universe is only an integer-ledger relaxation:
`2<=d<=40`, `1<=a<d`, `1<=N<=8`, `n in {0,1}`. Its overlap feasibility
condition is exactly `N max(0,d-2a)<=1-n`. These rows need not have
nonnegative outside cardinality and are correctly kept separate from
permutation or actual geometric passports.

I ran both normal and optimized Python and compared both outputs to the
frozen file. All three are byte-identical: **16,815 always-active gates**,
**32 raw labelled sparse passports**, and **6,419 integer-ledger rows**.

```sh
python3 -B 04-computation/planar_jc48_sep06_cusp_passport.py
python3 -B -O 04-computation/planar_jc48_sep06_cusp_passport.py
```

```text
source SHA256 8a7d40c0f001a2063f976b949d49033087adc3ae1616072c1c644172918e3d6f
output SHA256 6c330f724ac39a6960e22d13bee1d54efcd3cd13ff891acb470de41eccd04d95
raw-passport SHA256 e79a983bbf01db369c0a9c8f6c49b3449c1940f38b1bddad54ed272b4b015d7a
```

The source certifies its finite universes. The actual-access transfer,
normal-extension retention, all-degree lemma, and classical consumer are
separate proof steps and have been checked as such. This audit is frozen;
root owns integration and any canonical promotion.
