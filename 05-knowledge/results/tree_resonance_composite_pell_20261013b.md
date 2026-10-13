# Composite residue contexts, exact Pell genus, and a terminating golden control

2026-10-13. **PROVED, scoped and self-audited:** the prime-power and CRT
extension below (using the inherited CITED Weil/Hoffman bounds), the exact
genus formula, and universal entry for the explicitly defined golden-ring
control map. **FINITE-EXACT:** ring points, polynomial identities and the
complete golden box. **NUMERICAL:** floating-point spectra. No independent
peer audit or formalization. Six-colour forcing, planar Jacobian and
universal Collatz remain **OPEN**. The [Collatz companion](collatz_signed_resonance_pair_mass_20261013b.md)
contains the main computation and qualitative parity implication.

## 1. Inheritance, not another pass through closed constructions

[THM-4628, fixed-prime residue obstruction](../../01-canon/theorems/THM-4628-residue-colourings-of-the-unit-distance-graph-of-q-zeta55-are-kloosterman-forced-at-every-conjugation-fixed-prime.md)
is the nearest geometric mechanism. Its corrected fixed primes have
residue degree2 or10, not20, with quadratic residue half-field size
`q>=109`. Its scope explicitly left prime powers and composite ideals open.
The hostile is the reversed homomorphism implication: a graph mapping
onto a high-chromatic graph need not have that chromatic number. We use
the spectrum of the lifted graph itself, not that false shortcut.

For Jacobian, [THM-3765, normalized three-charge classification](../../01-canon/theorems/THM-3765-normalized-three-consecutive-charge-radial-keller-classification.md)
already excludes the entire smooth nonlinear polynomial ansatz containing
[THM-3757's Pell tower](../../01-canon/theorems/THM-3757-pell-chebyshev-three-charge-hyperelliptic-obstruction-tower.md).
We do not reopen it as a counterexample search. The useful refinement is
an exact formula for its genus and the differential obstruction it records.
The least-used sidecar is the mixed Chebyshev norm identity.

For square-triangular numbers, [THM-3335](../../01-canon/theorems/THM-3335-square-triangular-pell-markov-pythagorean-selector.md)
already supplies the orbit and its interfaces; [THM-3819](../../01-canon/theorems/THM-3819-chebyshev-pell-berggren-return-conductor-staircase.md)
already connects evaluation of this very Pell tower to integer states.
The [Pell/zeta note](pell_zeta_scale_20261011.md) already explains the common
quadratic unit in Apéry's limiting recurrence and why it supplies no new
irrationality theorem. Those stronger results are inherited, not re-proved.

## 2. Prime powers cannot evade the six-colour residue obstruction

Let `R` be the ring of integers of an unramified extension of `Q_p`, with
odd residue characteristic and residue field `F_q`. Let `S/R` be its
unramified quadratic extension, with involution bar. For r>=1 put

    G_(q,r)=Cay(S/p^r S, U_r),
    U_r={u: u*bar(u)=1 mod p^r}.

The graph has `q^(2r)` vertices and degree `d_r=(q+1)q^(r-1)`.
All norms and conductor depths refer to this specified local ring.

**R1 (PROVED).** Every nontrivial additive character of conductor `p^s`,
`1<=s<=r`, has eigenvalue bounded by

    |lambda| <= q^(r-s) 2 q^(s/2),
    |lambda|/d_r <= [2sqrt(q)/(q+1)] q^((1-s)/2).      (1)

At s1 this is the finite-field Kloosterman/Weil bound proved/referenced
in THM-4628. At s>=2 it follows from elementary local stationary phase,
as follows; it is not a claim that the finite-field formula can simply
be reused over a nonfield.

First reduce the character to its conductor: each norm-one point modulo
`p^s` has exactly `q^(r-s)` norm-one lifts, since the norm equation is
smooth. Its primitive sum is

    sum_(u in U_s) psi_s(Tr(a u)),                    (2)

where a is a unit of S and psi_s has exact conductor `p^s` on R. The trace
pairing is perfect because the extension is unramified. Choose an
anti-invariant unit delta with delta^2 a nonsquare unit in R. Around a
norm-one point u0, the local coordinate

    u(t)=u0 (1+delta*t)/(1-delta*t)

parametrizes its residue disk (`t in pR`). A critical residue of the phase
in(2) satisfies `Tr(a u0 delta)=0`, equivalently `a u0 in F_q`. Norm one
then says `(a u0)^2=N(a)` in F_q. There are zero or two such residues.
At either, the phase has second derivative `8 delta^2 a u0`, a unit;
there is a unique critical lift and an invertible local change of variable
that writes the phase as its critical value plus a unit times a square.
This follows by Hensel for the first derivative and the square root of
the remaining unit power series; only2 must be invertible.

For clarity, the resulting bound has a direct finite-sum interpretation.
At s=2m, summing in each fibre over the last m digits vanishes unless the
first m digits are critical. There are at most two such points and each
contributes modulus `q^m`. At s=2m+1, the last m digits impose criticality
to depthm and the next digit contributes a nondegenerate quadratic Gauss
sum of modulus sqrt(q), so each of at most two critical disks contributes
`q^m sqrt(q)`. The Gauss modulus follows by squaring its absolute value
and using additive-character orthogonality (2 invertible). This proves
`|(2)|<=2q^(s/2)`, hence(1), including p3. No archimedean approximation
or chromatic computation enters this proof.

Hoffman's bound now gives

    chi(G_(q,r)) >= 1+(q+1)/(2sqrt(q)).               (3)

Thus every power of a conjugation-fixed prime of `L=Q(zeta55)` still
requires at least7 residue colours. The local extensions there are
unramified: the fixed primes lie away from5 and11. Both degree2 and
degree10 primes are covered; q need not be a rational prime.

## 3. Composite ideals and the exact amount the spectral test retains

Let `I=product P_i^(r_i)` for a nonempty set of **distinct
conjugation-fixed primes** of L. Colour the points of the semilocal ring
`A={x in L: v_(P_i)(x)>=0 for all i}` by their residues moduloI.
The resulting finite unit graph is exactly

    G_I = tensor_product_i G_(q_i,r_i).              (4)

This requires a lifting assertion, not just CRT on vertices. Every local
norm-one residue u is `w/bar(w)` with w a unit: choose
`w=t+u*bar(t)` and choose the residue of t so w is nonzero. Such a choice
exists in the quadratic residue field; the same equation gives
`w=u*bar(w)` at every finite depth. CRT chooses one global integral w
with all these unit residues. Then `w/bar(w)` is an exact norm-one
element of L integral at every P_i, lifting all desired directions
simultaneously. Every residue vertex also lifts. Therefore a colouring
factoring through I is proper on A iff it is proper on(4).

This concerns the **semilocal domain**, not just the integral points
`O_L`, whose available exact integral unit directions can be smaller.
It also does not assert that the finite quotient graph embeds in the plane.

**R2 (PROVED).** For arbitrary nonempty regular graphs with degrees d_i
and normalized least eigenvalues `-eta_i`, their tensor product has
normalized least eigenvalue `-max eta_i`. Eigenvalues multiply; any
negative product has a negative factor of modulus at most max eta_i,
and every other normalized factor has modulus at most1. Equality is
obtained using that least factor and all other trivial eigenvalues.

Consequently the tensor product's Hoffman bound is exactly the minimum
of its factors' Hoffman bounds. Equations(1)--(4) give

    chi(G_I) >= 1+1/max_i[2sqrt(q_i)/(q_i+1)] >6.     (5)

Products and higher powers do not rescue a six-colour residue construction
in this class. They also do not automatically improve this spectral bound.
There is no conclusion here for ideals involving conjugation-swapped or
ramified primes, arbitrary additive periods, or geometric colourings that
retain information beyond the residue. The unrestricted plane problem is
not settled by an obstruction to this particular colouring rule.

**Controls.** Full FFT spectra were checked at `3^r (r<=4)`,
`7^r (r<=3)`, `11^r (r<=2)` in the nonsplit ring `u^2=-1`, and composite
moduli21,63,441,297. Their complete spectra agree with the CRT tensor
spectra. At the relevant first prime109, use `u^2=2`; the checker
constructs all11,990 norm-one points modulo109^2 exactly. At131 primitive
sample characters the even-depth stationary-phase formula agrees to
less than1.5e-13. The graph has141,158,161 vertices; no huge adjacency
matrix is needed. Its numerical Hoffman bound is6.599479, the same as
at109; the analytic Weil floor is6.268044, which already proves at least7.

## 4. The Pell tower's genus is exactly floor(n/2)

Retain THM-3757's definitions and put `t=2z-1`,

    a=1/[2(n+1)], b=1/(2n),
    P_n(t)=a C_(n+1)(t)-b C_n(t),
    S_n(t)=a U_n(t)-b U_(n-1)(t),
    c_n=P_n(-1)=(-1)^(n+1)(a+b).

The derivative identity `d_z P_n=U_n-U_(n-1)=G_n` gives
`chi_n=P_n-c_n`. The mixed Pell identity

    C_(n+1) C_n-(t^2-1) U_n U_(n-1)=t

and the two individual Pell identities give the exact cancellation

    P_n^2-(t^2-1) S_n^2 = a^2+b^2-2ab t.            (6)

Since `4F_n=(t^2-1)S_n^2`, the generic-fibre polynomial becomes

    Delta_n(z)=(Lambda+c_n)^2
                -2(Lambda+c_n)P_n(2z-1)
                +a^2+b^2-2ab(2z-1).                 (7)

For every n>=1 it has degree **exactly n+1** over `k(Lambda)`: the leading
coefficient of P_n is nonzero and `Lambda+c_n` cannot vanish there.
THM-3757 already proves generic squarefreeness. Therefore its smooth
projective genus is exactly `floor(n/2)`, improving that file's former
lower bound and finite check of equality through n8.

This gives genus1 at depths2,3 and genus9 at depths18,19. The earlier
no-rational-mate conclusion is unchanged: `dz/Y` is a nonzero holomorphic
time form in positive genus, so cannot be a rational derivative. Equation(7)
describes the closed family more accurately; it does not evade THM-3765.
The checker verifies the complete identities at every n1..20 independently
by polynomial expansion, including both requested genera.

The square-triangular connection is literal: evaluate the underlying norm
identity at t3, whose unit is `3+sqrt(8)`. The inherited Pell orbit obeys
`x_k^2-8y_k^2=1` and `T_((x_k-1)/2)=y_k^2`. Golden arithmetic uses the
different unit `phi^2`, of trace3 rather than6. Both are determinant-one
two-embedding recurrences, but their fields, integral lattices and growth
rates differ. A shared recurrence format does not identify their solutions.

## 5. A complete deterministic-tree control over the golden integers

The [earlier golden pair-chain experiment](../../04-computation/experiments/golden_20261007_readers/golden/golden_pair_chain.py)
introduced a distinct map on `Z[phi]`, with `phi^2=phi+1`:

    G(x)=x/2                         if x=0 mod2,
    G(x)=(phi*x+rho(phi*x))/2         otherwise,

where rho takes both coefficients modulo2 into `{0,1}`. This map is
**not asserted conjugate to Collatz**, nor is it the beta-map in
[THM-4528](../../01-canon/theorems/THM-4528-collatz-parity-in-base-phi-golden-beta-map-and-holonomy.md).

**G1 (PROVED).** Every algebraic integer in this map eventually reaches
one of the four fixed points `0, phi^2, phi^3, phi^4`.

Let `u=a+b phi`, `v=a+b phi'` be its two real embeddings. Every step obeys

    |u'| <= (phi/2)|u|+phi^2/2,
    |v'| <= max(|v|/2, (|phi'|/2)|v|+1/2).

The maximum in the second line is essential: the even branch has factor
1/2, larger than |phi'|/2. At x4 the odd-branch-only inequality is false.
The corresponding limiting radii are `phi^4<7` and `1/(2-|phi'|)<1`;
for the second, the positive excess above that radius contracts by at
most1/2 on either branch. More explicitly, for
`H=max(|u|/7,|v|)` and H>=1,

    H(Gx) <= kappa H(x),  kappa=(8phi+1)/14<1.       (8)

The second embedding's factor is `(|phi'|+1)/2=phi/2<kappa`.
So every orbit enters the invariant box `|u|<=7, |v|<=1` in a bounded
number of steps computable from its original algebraic integer. That box
has exactly15 lattice points: `|b|<=8/sqrt5<4`, `|a|<3`, and direct exact
inequalities select the points. Its complete transition graph has the
four stated fixed points and no other cycles. Every box point reaches a
fixed point in at most four further steps. The output lists every edge.

This is a genuine universal certificate mechanism in a nearby system:
an explicit norm bound followed by a complete finite transition table.
It explains the critical missing condition in attempts to import it into
Collatz. The classical odd branch expands real size by3/2, while the
golden control contracts at both embeddings. A representation of Collatz
would have to retain and pay the extra source-dependent carry rather than
silently replacing that expansion. The [ROOT-supported mean theorem](collatz_signed_resonance_pair_mass_20261013b.md)
remains an averaged statement and does not provide this pointwise norm.

### A family of certificate-bearing systems, and the linear-projection test

The control generalizes to a precise elementary design rule. On the
integer lattice of any number field K, take finitely many branches
`x -> (a_j x+b)/m`, with m>0 an integer, digits b in a fixed finite set,
and branch rules that keep the output integral. If every embedding sigma
and branch obeys `|sigma(a_j)|<m`, then each embedding has a strict affine
contraction bound. Choose radii strictly above all its limiting bounds.
The resulting Minkowski box is finite and eventually traps every source.
An explicit entry bound followed by its complete finite transition graph
is a certificate to the corresponding cycle. This proves eventual
periodicity for this whole family, not necessarily a unique cycle.

For Collatz the additional datum must be an authenticated projection
`pi` satisfying `pi(Gx)=T(pi(x))`, plus the source's actual lifted state.
Merely tagging an unrelated contracting orbit with the integer n does
not supply that equation.

There is a sharp cheap obstruction to the simplest proposed projection.
Suppose `pi(x)=Tr(c x)+d`, with c nonzero, and one fixed digit/branch agrees
with the Collatz **odd** rule on a set of integral states whose differences
span K over Q. Subtract the intertwining equations at two such states:

    Tr(c a_j h)/m = (3/2)Tr(c h)   for every h in K.

Nondegeneracy of the trace pairing forces `a_j=3m/2`. That branch therefore
expands every embedding by3/2, contradicting the all-contracting criterion.
This is a **PROVED scoped obstruction**, not a rejection of enriched
representations. It locates the necessary reroute: a nonlinear projection,
an invariant source subset with smaller affine span, infinitely many
authenticated branches, or an explicit accumulated carry whose decrease
is paid over several steps. The source subset and projection equation must
remain part of the proposed structure.

## 6. Genus1 to genus9: the information added by one more bit

The [marked fish cover](fish_marked_elliptic_cover_20261013.md) already
proves that its first three square-root coordinates define a degree-eight
genus-one cover, whereas the fourth gives a degree-sixteen genus-nine
cover. The latter double cover adds sixteen ramification points:

    2g_new-2=2(2g_old-2)+16.

Its character decomposition is seven genus-one pieces, one genus-two
piece and seven genus-zero pieces. Our finite branch-set calculation
reproduces this; it is not new evidence for an actual fish placement.

The relevant connection to a deterministic Collatz tree is the distinction
between an extra **label** and an extra **compatibility equation**. A
free binary prefix merely refines the source cylinder. A square-root
coordinate imposes a lift and creates ramification. Genus measures the
resulting differential space, not the number of bits or the Collatz
stopping time. Keeping all four roots is analogous to retaining `(k,e)`
and all source digits instead of retaining just the pair chain's level k.
There is presently no map identifying the fish's genus with a Collatz
invariant; none is used in a proof.

## Reusable maps and next obstructions

| Source -> target | Preserved | Lost unless retained | Decisive next test |
|---|---|---|---|
| Parity tree -> pair-state transfer | Absorption mass with integer edge weights | Original source address | A pointwise decreasing carry bound |
| Prime power -> residue field | Conductor1 character sums | Primitive higher conductors | R1 pays these by stationary phase |
| Composite ideal -> factors | All unit directions by simultaneous lifting | Geometric injectivity | Use R2 only for residue-factor colourings |
| Pell profiles -> scalar unit recurrence | Norm identity and growth | Generic time-form class | Equation(7) retains its exact genus |
| Golden integral map -> finite box | Both embeddings and all transitions | None for this specified map | Any proposed Collatz intertwiner must account for3/2 |
| Fish C4 -> elliptic C3 | Three square-root coordinates | Fourth root and its ramification | A point on C3 must lift to C4 with marked ports |

The next geometric route should retain within-residue information, or
examine conjugation-swapped/ramified factors with their actual norm image.
Repeating fixed-prime-power contexts is now excluded as a six-colour
residue construction. For Jacobian, changing the ansatz must change the
time-form class; simply taking another Pell depth is already ruled out.

Reproduce:

```sh
python3 04-computation/experiments/residue_pell_checks_20261013b.py --output /tmp/residue-pell.json
python3 -O 04-computation/experiments/residue_pell_checks_20261013b.py --output /tmp/residue-pell-O.json
```

[Saved output](residue_pell_checks_20261013b.json). Polynomial/norm/finite-box
claims are exact. Spectral comparisons use floating point and are labelled
NUMERICAL; the uniform colour lower bound rests on the proof above.
