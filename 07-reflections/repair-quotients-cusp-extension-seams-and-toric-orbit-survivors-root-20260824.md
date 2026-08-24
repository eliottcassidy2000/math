# Repair quotients, cusp extension seams, and toric orbit survivors

**Session synthesis, 2026-08-24.** The proved advances integrated here are
[THM-3989](../01-canon/theorems/THM-3989-cusp-log-laurent-conductor-and-nondividing-depth-reduction.md)'s
filtered depth-`2:3` supplement,
[THM-3990](../01-canon/theorems/THM-3990-componentwise-harmonic-obstruction-and-repair-quotient.md),
[THM-3991](../01-canon/theorems/THM-3991-periodic-unimodular-toric-cusp-factorial-euler-obstruction.md),
[THM-3992](../01-canon/theorems/THM-3992-reduced-two-three-cusp-jet-repair-and-first-node-residual.md),
[THM-3994](../01-canon/theorems/THM-3994-affine-p-double-collision-separates-two-addresses-from-one-a1-centre.md),
[THM-3995](../01-canon/theorems/THM-3995-scale-two-parity-hole-support-and-integer-variance-tariff.md),
and
[THM-3996](../01-canon/theorems/THM-3996-etale-node-address-balance-cycle-and-nonproperness-dichotomy.md).
They do not prove `JC(2)`, LRC(14), positive curvature on `S^2 times S^2`, or
a complex structure on `S^6`. The two global 2026 preprint claims remain under
their source audits.

## Inheritance and live board

The anchor was the first nondividing cusp depth `2:3`. Its closest proved
mechanism was THM-3989's Laurent conductor and Euclidean target-shear descent;
the hostile was `A=x^2+2u`, whose leading symbol is a square but whose next
coefficient does not lift. The corrected near miss was “surjective associated
graded symbols imply a filtered approximate root.” The least-used sidecars
were THM-3406/3412's filtration--observer failures.

The niche was the componentwise Poisson step in the incoming curvature source.
The source residual set is a connected two-torus, so its one average is the
complete final linear obstruction; a disconnected mixed-sign residual is the
canonical hostile. The wildcard was the incoming periodic toric cusp: its
Euler characteristic is a top-cell orbit count, with a nonunimodular
five-tetrahedron cube as the nearest failure.

The board stabilized at six concepts:

1. image of an allowed repair operator;
2. the quotient class which that image cannot change;
3. component averages and integral Smith torsion;
4. cusp normalization seams invisible in the depth symbol;
5. translation-orbit strata and compactly supported Euler characteristic;
6. LRC's nonlinear auxiliary-cover residue after the two-comb failure;
7. complete oriented node-address graphs versus one-line clutch samples.

## Typed comparison

| Lane | Source and allowed operation | Survivor / preserved predicate | Destroyed information | Sidecar and cheapest hostile |
|---|---|---|---|---|
| THM-3990 smooth | `w`, corrected by `Delta chi` | component averages; existence of a common strict sign | the zero-mean part on each component | component labels; averages `+1,-1` |
| THM-3990 integral graph | `Z^V`, corrected by `L Z^V` | free component sums plus Smith torsion | real rank forgets lattice index | Smith form; `C_3` class `(1,-1,0)` of order `3` |
| THM-3989 filtered cusp | depth-one symbol `h`, corrected by adding `P_0=k[p,y]` | `theta(h)=[s^3]h*[s]` in `k[s]/k[s^2,s^3]`; physical liftability | the associated graded forgets the constant seam and all higher seams | specialization `tau=-s^2`; `x^2+2u` |
| THM-3991 toric cusp | cells modulo the translation sublattice | top-simplex orbit count `d*n!` and vertex-orbit component count `d` | attaching maps, clutch orientation, period functions, analytic existence | `chi_c`; nonunimodular cube count `10` versus volume `12` |
| THM-4000 LRC method | danger combs covering the exact auxiliary set | fixed-radius containment certificate | a cover count has no physical time, arrival, or runner owner | exact set family; no one/two combs cover, `(8,9,10)` does |

The shared move is precise:

```text
derive the actual allowed operation
  -> quotient by its image/orbits
  -> test the surviving class in the consumer's category.
```

This is a research operation, not an object-level reduction. There is no map
from a graph Laplacian to the cusp Hamiltonian bracket, from a toric cell orbit
to an LRC safe arc, or from a component average to a Laurent coefficient. The
operators, coefficient categories, and consumers differ. What transfers is
the audit order and the demand to retain the correct sidecar.

## The new cusp coordinate

For every

```text
F=sum_i f_i(s) tau^i in B_2,
```

the exact color `p=0`, equivalently `tau=-s^2`, sends
`(x,u,p,y)` to `(-s^-1,-1,0,0)`. Therefore

```text
sum_i (-1)^i [s^(ell-2i)]f_i=0                    (ell>=1).
```

At depth one this computes the extension, not merely an obstruction:

```text
0 -> P_0 -> P_1 -> s k[s] -> 0,
theta(h)=[s^3]h*[s] in k[s]/k[s^2,s^3].
```

A two-term lift `(h,q)` exists exactly when `[s]q=[s^3]h`. Hence a depth-two
coefficient with leading `h^2` has a square approximate root modulo `P_0`
exactly when `h|a_-1` and this equality holds for `q=a_-1/(2h)`.

The first four normalized `2:3` rows still align over `k(s)` with a formal
base `R`: the first free term is `C=R^3+kappa R^-1+...`. The seam decides
whether that rational formal base is a polynomial element of the filtered
ring. In the hostile slice

```text
A=x^2+lambda*u+F(p,y),                     lambda!=0,
```

the first seam can be paid, but only if

```text
F_p(0,0)=lambda^2/12.
```

This condition is necessary, not sufficient. The next exact task is to keep
the integration scalar `K`, the positive coefficients `c_1,c_2,...`, the
scalar moment, and the entire seam tower until the forced `F` terminates,
becomes nonalgebraic, or hits a unit obstruction.

An intermediate calculation omitted the first positive contribution
`[s^0]c_2` at seam `ell=4` and incorrectly suggested `F_y(0,0)=0`. That claim
is **RETRACTED**. The free value `c_2(0)` absorbs the apparent condition; the
canonical supplement contains a regression gate for this exact failure.

## Extensions of the incoming problems

THM-3990 was independently replayed on all `1,099` labeled simple graphs with
at most five vertices, under both Laplacian signs, plus exact Smith controls
and positive rational-weighted examples. The replay exposed three wording
boundaries: both cokernel displays are isomorphisms rather than literal
equalities; `C_3` is the smallest *simple*-graph torsion hostile; and the
compact two-region proof must dispatch empty regions. A two-vertex double-edge
multigraph already has Smith form `(2,0)`.

For the curvature source, the proved consequence is narrowly favorable: its
displayed `Sigma` is connected, so nonzero integral really is the full
Poisson obstruction. The load-bearing symbolic notebook identities and the
global positivity argument remain separate audit obligations. Incoming exact
controls now verify the omitted `P1(L(h_a))=0` identity symbolically on both
generic parametrized components, relative to the reconstructed formulas, and
make the three nonzero terms of `V_bc^(2)` cancel at one generic exact point.
A typed background symmetry reduces global `V_bc=0` to the source's saved
`V_ac=0` block; because that antecedent has not been independently rebuilt,
the global identity remains conditional and the headline theorem remains
under audit.

The incoming reduced-cusp calculation THM-3992 makes the seam extension
strictly sharper. Once its bracket/moment argument has forced `h=gamma*s` and
centered the nodal line, put `eta=[s]q0`. The `A,ell=2` and `C,ell=1` seams
eliminate all unseen second-jet coefficients and give

```text
eta*(gamma*eta-a)=0.
```

Thus `eta=0` is not merely formally plausible: the exact `P_1/P_0` sequence
constructs a square-root lift of `A` modulo `k[p,y]`. The other branch is the
unique first-seam nonliftable one and pays `[s^2]b=2a`. A hostile `B_2` packet
survives the negative rows through `-2` and fails at `-1`, so this is a typed
branching theorem rather than a Keller construction. Further exact work now
separates both lanes. Every square lift in the first branch has constant nodal
remainder `a!=0`, so no simultaneous pure cube lift exists; the correct
replacement is the moving-node identity

```text
(C-L)^2=(A-K)(A+K/2)^2.
```

The second branch reaches its next unused seam exactly when
`a^3+2gamma=0`. Its scalar shadow `16I^3=27gamma^2` is weaker because it
forgets orientation. On that seam, `A2(0)=q2^2-32/(9a^9)` and `C3(0)` absorbs
the following seam; an arbitrary one-parameter exact jet survives through
`ell=6`. This is finite consistency only, so algebraization and positive rows
replace the obsolete row-`-1` test.

THM-3994 supplies the correct warning for both node addresses. Its two
double-resultant seams have the same total multiplicity but respectively
carry two reduced supports and one curvilinear length-two support with an
`A1` Rees chart. This is a preservation-law connection only: its basepoint
addresses are not identified with THM-3992's normalization clutches. What
transfers is the requirement to retain labelled incidence before taking a
class or Smith quotient.

THM-3996 supplies the missing ownership law without importing THM-3951's
completion-boundary forest. For an integral nodal target outside the
nonproperness locus, finite-etale degree makes every normalization owner
balanced, hence every complete address edge lies on a directed cycle. Thus
distinct companions in THM-3992 require an additional address or a nonproper
node; a complete two-address packet is a directed two-cycle. This implication
is one-way. Disconnected normalization, an incomplete census, deleted
properness, and ramified endpoint maps are exact hostile boundaries.

THM-3991 was independently checked through dimensions `1..6` with diagonal
index controls `d in {1,2,4}`. The proof must use compactly supported Euler characteristic on the
locally closed toric strata and then compactness to return to ordinary Euler
characteristic. In the exact one-nonzero-Euler-fibre grammar, sphere homology
forces `(n,d)=(1,2)` or `(2,1)`; irreducibility forces `(2,1)`. Higher-rank
work must therefore exhibit an explicit escape: another Euler fibre, a
verified finite quotient, a nonunimodular/singular correction, or a different
degeneration grammar. Merely changing period functions cannot change the
factorial count.

For LRC(14), THM-4000's exact three-versus-one obstruction says the current
unweighted auxiliary-cover operation has saturated. A lawful next probe must
change the operation—weighted cover, a coupled clock, or a carrier retaining
arrival—not reinterpret the cover number as a linear cokernel. LRC(14)
remains open. THM-3995 nevertheless gives the scale-two survivor the sharper
fixed trimmed carrier

```text
(2/21+1/(42t),8/63-1/(126t)) union
(55/63+1/(126t),19/21-1/(42t)),
```

support at most `4(t-1)/(63t)`, and an integer variance correction. This can
strengthen body-specific endpoint certificates. In THM-4000's notation it
gives the exact hybrid sufficient certificate

```text
disc_t(G)<meas(G)^2*(59t+4)/(4(t-1))
          +(s_t/t^2)*theta_t*(1-theta_t).
```

THM-4000 already closes its bounded `[1,20]` family; the hybrid closes no
arbitrary type and does not replace the missing cross-phase control.

## Reproduction anchors

- THM-3989 supplement: `13,007` gates; script/output/semantic hashes
  `ba52ebbb...d1b67`, `2a413b04...c312`, `2e57c7d6...a815`.
- THM-3990 independent audit: `33,517` gates; hashes
  `538fcf10...a9db`, `0e18bcf4...f5ea`, `2aa53145...5adc`.
- THM-3991 independent audit: `1,801` gates; hashes
  `2fcf3788...d385`, `21a65442...04fa`, `ea1b877b...4208`.
- THM-3992 combined cusp/seam audit: `119` exact pass rows; hashes
  `4c7fb202...cfa2b`, `1e44afac...dc5c8`.
- THM-3996 node-balance audit: `10` canonical pass rows plus an independent
  exhaustive binary-digraph hostile suite through four vertices.

The full paths and full hashes are pinned in the theorem frontmatter and the
results index. All hash-bearing companions pass under normal and optimized
Python after LF normalization where the host emits CRLF.
