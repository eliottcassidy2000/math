# Positive-divisor source flows: exact degrees and alternating/symmetric groups

2026-10-19. **PROVED; frozen after independent review.** This classifies the Galois groups
and individual coordinate degrees on the entire positive-divisor algebraic
locus of packet 128. It includes a proof that generic-target monodromy
survives the actual fixed nonzero scalar-time source substitution. The
full group classification explicitly uses Gareth Jones's theorem on
primitive groups containing a cycle, whose general proof depends on the
classification of finite simple groups. No historical novelty is claimed.

## 1. Statement and source parameter

Let K be any field of characteristic zero, let `m≥2`, and let `a≥1` divide
m. Put

    r=m/a,   k=2m+a,   n=m+a+1=a(r+1)+1,
    H=x^k t^m(1+x²t),   λ∈K*.

Let `(X,T)=exp(λ{−,H})(x,t)` be the source t-adic flow from packet 128.
Set

    u=x²t,   E=H=x^a u^m(1+u),   s=1/(x u^r).

Then `F=K(x,t)=K(E,s)` and

    u=Es^a−1,
    x=1/[s(Es^a−1)^r],
    t=s²(Es^a−1)^(2r+1).                            (1)

E and s are algebraically independent. Define the polynomial primitive
with zero constant term by

    P_E(v)=−Σ_(j=0)^r binom(r,j)(−1)^(r−j)
                    E^j v^(a(j+1)+1)/(a(j+1)+1).

Thus

    P'_E(v)=−v^a(Ev^a−1)^r,
    deg P_E=n,   leading coefficient ℓ=−E^r/n.       (2)

Packet 128 proves the actual source identity

    P_E(S)=P_E(s)+λ,   S=1/[X(X²T)^r],               (3)

with X and T recovered by (1), replacing s by S.

**Theorem.** Let `C=overline{K(E)}` and retain s as a variable. The
polynomial `P_E(V)−P_E(s)−λ` is irreducible of degree n over both F and
`C(s)`. Its splitting-field group over `C(s)` is

    G_geom = A_n  if a and r are both even,
             S_n  otherwise.                       (4)

If the geometric group is S_n, its group over F is also S_n. In the
remaining cases put

    δ=(−1)^(a/2)n ∈ K*.

The group over F is A_n exactly when δ is a square in K, and is S_n
otherwise. Moreover,

    F(X)=F(T)=F(X,T)=F(S),
    [F(X):F]=[F(T):F]=n.                            (5)

The individual coordinates are expressible by radicals over F exactly
in the quartic case `(m,a,k)=(2,1,5)`. In every other case in this positive
locus, both individual coordinates are not expressible by radicals.

The geometric convention in (4) treats E as part of the constant field.
The arithmetic statement is over the original source rational field F.
No numerical specialization of E or s is included.

## 2. The critical values are distinct

Choose `ρ∈C` with `ρ^a=1/E`. The critical points of P_E are zero, with
derivative multiplicity a, and the a points `ρζ`, `ζ^a=1`, each with
derivative multiplicity r. The corresponding critical values are

    0,   κρ^(a+1)ζ  (ζ^a=1),
    κ=(−1)^(r+1) a^r r! / ∏_(j=1)^(r+1)(aj+1) ≠0. (6)

Indeed substituting `ρζ` into the displayed primitive gives the factor
`ρ^(a+1)ζ` times

    −Σ_(j=0)^r binom(r,j)(−1)^(r−j)/(a(j+1)+1)=κ.

The last identity follows, for example, by integrating
`−v^a(v^a−1)^r` from zero to one and evaluating the elementary beta
integral; it is an identity of rational numbers and hence valid over K.
Consequently P_E has a+1 distinct critical points and a+1 distinct
critical values. Write their set as

    B={0} ∪ { b : b^a=κ^a E^(−a−1) }.               (7)

In particular there are no critical-value collisions in this family.

## 3. Indecomposability and primitive target monodromy

We give the elementary decomposition obstruction explicitly. Suppose
`P_E=A∘B_0`, with both polynomial degrees greater than one. If α is a
critical point of A, every root of `B_0−α` is critical for P_E and has
image `A(α)`. Critical injectivity forces all those roots to be one
point b. Thus

    B_0(v)−α=c(v−b)^d,   d=deg B_0.

Two different critical points α of A would give two different such
representations. Their derivatives force the unique root of `B'_0`
to be the same b, and evaluating there then forces the two α's to agree.
Thus A has just one critical point. In characteristic zero it is a
translated pure power, and the preceding expression for B_0 makes P_E
itself a translated pure power. That contradicts its at least two
distinct critical points. Therefore P_E is indecomposable over C.

Let Y be an independent target variable. The polynomial `P_E(V)−Y` is
irreducible over `C(Y)` and its group G acts transitively on its n roots.
Indecomposability implies that G is primitive. In field terms, a block
would yield a proper intermediate field of `C(v)/C(P_E(v))`; Lüroth's
theorem writes it as `C(B_0(v))`. Because infinity is the unique point
above infinity for a polynomial map, a generator can be chosen so that
B_0 and the resulting outer function A are polynomials. This would give
the decomposition just excluded.

For comparison with the literature, [Kreso–Tichy, *Diophantine equations
and the monodromy groups*, Lemma 3.4 and Proposition 3.9](https://arxiv.org/pdf/1601.07316)
prove the decomposition/primitivity equivalence and the stronger
2-transitivity conclusion for characteristic-zero polynomials with at
least two distinct critical points and all distinct critical values.
The argument above proves the indecomposability needed here directly.

## 4. Exact generic-target group

Since each finite critical value has just one critical point above it,
formal local lifting gives the following individual cycles in G:

- an `(a+1)`-cycle at the value zero, fixing `n−a−1=m` roots;
- an `(r+1)`-cycle at each of the other critical values, fixing all
  remaining roots.

For completeness, at a critical point c of ramification index e,
`P_E(c+w)−P_E(c)` starts with a nonzero multiple of `w^e`. Put the target
parameter equal to `η^e`. The e distinct nonzero leading coefficients
of w lift formally, and `η↦ζ_eη` cycles those e roots. All other roots
lift as power series in `η^e` and are fixed. Thus these are single
cycles, rather than products of several nontrivial cycles.

We now use exactly the following consequence of Jones's theorem: a
finite primitive permutation group containing a single cycle with at
least three fixed points contains the full alternating group. This is
[Gareth A. Jones, *Primitive permutation groups containing a cycle*,
Corollary 1.3, Bulletin of the Australian Mathematical Society 89 (2014),
159–165](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/BBF3021D8434AD2C0BB0D3D6CE637008/S000497271300049Xa.pdf/primitive-permutation-groups-containing-a-cycle.pdf).
No prime-length hypothesis is required. Its general proof depends on the
classification of finite simple groups; we do not replace it by an
unproved elementary extension of Jordan's prime-cycle theorem.

If `m≥3`, the `(a+1)`-cycle fixes exactly m points, so Jones's result gives
`A_n≤G`. If `m=2`, the only pairs `(a,r)` are `(1,2)` and `(2,1)`;
each has a simple critical point and hence a transposition. A primitive
group containing a transposition is S_n: the graph of its conjugate
transpositions has invariant connected components, so primitivity makes
it connected, and its edge transpositions generate S_n.

If either a or r is odd, one of the displayed cycles is odd. Thus
`G=S_n`. If both are even, (2) is minus a square polynomial; the
monic discriminant is a square over `C(Y)`, as computed in section 6.
Thus `G≤A_n`, and consequently `G=A_n`. This proves (4) for the generic
target equation `P_E(V)=Y` before any source substitution.

For the elementary subfamily a=1, Jones's theorem is unnecessary:
the cycles include a transposition and an `(n−1)`-cycle. The latter
makes a point stabilizer transitive on the other points, giving
2-transitivity, hence primitivity and S_n. More generally the usual
2-cycle/3-cycle arguments suffice when `min(a,r)≤2`.

## 5. Why the actual fixed scalar time preserves the group

Generic-target monodromy does not automatically survive replacement
of Y by a function of the source parameter. We prove the required
linear disjointness here.

Over `C(Y)`, let L be the splitting field of `P_E(V)−Y`, and let L_λ
be the splitting field of `P_E(U)+λ−Y`. Their finite branch loci are,
respectively, B and `B+λ`; infinity is the only other possible branch
point. These finite sets are disjoint. Indeed every difference of two
members of B has the form

    κρ^(a+1)(ζ_i−ζ_j),

where each ζ is either zero or an a-th root of unity. A nonzero such
difference cannot equal λ: taking a-th powers would make `E^(a+1)`
algebraic over K, whereas E is transcendental. A zero difference cannot
equal the stipulated nonzero λ either. Thus

    B ∩ (B+λ)=∅.                                    (8)

The intersection `I=L∩L_λ` is a finite Galois extension of `C(Y)`.
As a subextension of both splitting fields, it can be branched only
at infinity. A nontrivial finite cover of the projective line over an
algebraically closed characteristic-zero field cannot have just that
one possible branch point. If its degree is d and there are h points
above infinity, tame Riemann–Hurwitz gives

    2g−2=−2d+(d−h)=−d−h.

Since `g≥0` and `h≥1`, this forces `d=h=1`. See the
[Stacks Project, Riemann–Hurwitz, Section 53.12](https://stacks.math.columbia.edu/tag/0C1B)
for the characteristic-zero formula used here. Hence

    L∩L_λ=C(Y).

Finite Galois extensions with trivial intersection are linearly
disjoint. Embed the degree-n source field `C(s)` into L_λ by
`Y=P_E(s)+λ`. Then L and `C(s)` are also linearly disjoint, and

    Gal(L C(s)/C(s)) ≅ Gal(L/C(Y))=G.                (9)

The left side is the splitting-field group of the actual equation (3).
It acts transitively on n roots, so that equation remains irreducible
of degree n. Irreducibility over `C(s)` implies irreducibility over F.
This proves the fixed-time group and degree statements for every
`λ∈K*`, without specializing a generic-time Galois group by assumption.

The nonzero constant-time hypothesis is essential. At λ=0 the equation
has the factor `V−s` and the source branch is the identity. Nor does the
argument automatically apply to times allowed to depend on E, which
can cause the two finite branch sets to meet.

## 6. Arithmetic group and the exact square class

When `G_geom=S_n`, the arithmetic group over F contains this group and
is itself a subgroup of S_n, so it equals S_n.

Suppose now that a and r are both even. Then n is odd. For the generic
monic polynomial

    f_Y(V)=[P_E(V)−Y]/ℓ,
    h(V)=V^(a/2)(EV^a−1)^(r/2),

we have `f'_Y=(n/E^r)h²`. Therefore

    Disc(f_Y)= (−1)^(n(n−1)/2) (n/E^r)^n
                                      Res(f_Y,h)². (10)

Because r is even and n is odd, the square class of (10) in
`K(E,Y)*/K(E,Y)*²` is

    (−1)^((n−1)/2)n = (−1)^(a/2)n = δ.              (11)

Substituting `Y=P_E(s)+λ` preserves this expression with a nonzero
square factor: the substituted polynomial is separable. Thus its
arithmetic discriminant has square class δ in F. The geometric group
is A_n, so the arithmetic group is either A_n or S_n, and the
usual discriminant sign character distinguishes them. Finally a
constant δ is a square in `K(E,s)` exactly when it is a square in K,
since K is algebraically closed in its pure transcendental extension.
This proves the arithmetic assertion.

For example, `(m,a,k)=(4,2,10)` has `n=7`, geometric A7, and square
class `δ=−7`; over Q its arithmetic group is S7, whereas over
`Q(√−7)` it is A7. There are also genuine arithmetic alternating
examples over Q: `(m,a,k)=(16,8,40)` gives `n=25` and `δ=25`, so both
groups are A25. These show precisely why a universal full-symmetric
claim would be false.

## 7. Both individual coordinates generate the full extension

Let v be a chosen root of `P_E(V)−Y` in L, and define

    χ(v)=1/[v(Ev^a−1)^r],
    θ(v)=v²(Ev^a−1)^(2r+1).

Their degrees as rational maps of the v-line are

    deg χ=m+1<n,
    deg θ=2m+a+2=2n−a,
    n<2n−a<2n.                                     (12)

Neither degree is divisible by n. Consequently neither χ nor θ lies
in `C(P_E(v))=C(Y)`: membership would give a rational composition
`R∘P_E`, whose degree is n times the positive integer `deg R`.

After the source base change, χ(v) and θ(v) are the coordinate images
X and T of the selected algebraic branch. If either belonged to
`C(s)`, it would belong to

    L∩C(s)=C(Y),

contradicting (12). Hence both coordinates are nonrational even over
the extended source field `C(s)`.

The natural degree-n actions of A_n and S_n are primitive. Thus the
root field `C(s,v)` has no proper intermediate field over `C(s)`.
Each of its two nonrational coordinate elements therefore generates it.
Since the degree over `C(s)` is n and the degree over F is already n,
the tower-law bounds give (5) over F as well. This proves exact
individual degrees, rather than merely assigning a degree to a time
primitive or to an unselected algebraic equation.

## 8. Radical solvability and scope

The normal closure of each coordinate field is the splitting field in
(9), with arithmetic group A_n or S_n. If `n≥5`, either group is
nonsolvable, so neither coordinate is expressible by radicals over F.
The standard radical-solvability criterion is recalled, for example,
in [Sharifi, *Abstract Algebra*, Section 10.3](https://math.ucla.edu/~sharifi/notes/algebra-ch10.html).
The only case here with `n<5` is `(m,a,k)=(2,1,5)`, when the group is
S4 and the coordinates are expressible by radicals. This agrees with
the quartic case of packet 127; `(m,a,k)=(2,2,6)` gives its S5 quintic.

This packet concerns only the positive-divisor algebraic locus of
[packet 128](planar_two_exponent_flow_classification_20261019.md).
It does not classify Galois groups on the a=0 or negative-divisor loci,
or after numerical source specializations. It changes none of the
polynomiality obstructions and makes no claim about planar JC.

## Verification scope

The proof is universal. Its substantial external input is Jones's
primitive-cycle theorem; the source identity and birational coordinates
come from packet 128. The
[small exact checker](../../04-computation/experiments/planar_positive_divisor_galois_classification_20261019.py)
passed 31 checks on four specified witnesses `(a,r)=(1,2),(2,1),(2,2),(3,3)`
and the A25 integer square-class calculation. These test the primitive,
critical-value, coordinate, and discriminant identities. The
[normal output](planar_positive_divisor_galois_classification_20261019.json)
and [optimized output](planar_positive_divisor_galois_classification_20261019_optimized.json)
are byte-identical. Such checks do not prove the group classification or
replace the disjoint-branch field argument.

The [independent proof review](planar_positive_divisor_galois_classification_20261019_independent_review.md)
accepted the full mathematics without corrections. Its pinned draft
differs from this final version only in review status and verification
metadata, not in mathematical statements or proof.
