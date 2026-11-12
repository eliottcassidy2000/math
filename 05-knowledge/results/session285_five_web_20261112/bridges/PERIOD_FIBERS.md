# Five local parameters, two to five global parameters

**Status: PROVED, with the hypotheses below.** This is an elementary global
consequence of the local stress normal form, not a claim of literature priority.
It is a bridge theorem, not a Chan–Pak result.

## 1. The exact global obstruction

Let Ω be a connected plane domain, and let h be a single-valued nonvanishing
holomorphic function. The net's principal directions are encoded by the phase
of h. Use the stress lane's convention

    p = σxx + σyy,     q = σxx − σyy − 2i σxy.

Require equilibrium without body forces, harmonic trace Δp=0, and alignment
q/h real. The local theorem in [stress](../stress/) gives

    R=q/h=a|z|²+bz+bar(b)bar(z)+c,
    p=d−2 Re ∫ h(z)(az+bar(b)) dz,                    (1)

where a,c,d are real and b=u+iv. The primitive starts at a fixed basepoint,
so d is the trace there.

**Global fiber theorem.** The real vector space of single-valued aligned
stress fields has dimension

    5 − rank_R P_h,                                  (2)

where the period map on three real coordinates is

    P_h(a,u,v)(γ) = Re ∮γ h(z)(az+u−iv) dz.           (3)

In particular its dimension is between two and five, regardless of the number
of holes. At most three independent real period constraints are needed.

**Proof.** Local forms glue with identical polynomial coefficients: a real
polynomial vanishing on an open overlap is zero. Thus a,b,c are global
constants. The real primitive defining p is path independent exactly when
(3) vanishes on every closed curve. Conversely this vanishing constructs p
and (1) gives equilibrium, alignment and harmonic trace. The local parameter
map is injective, so dimension is kernel dimension plus free c,d, proving
(2). The domain of P_h is R³ even when H_1(Ω) has infinite rank. Periods depend
only on homology because the holomorphic one-form is closed. □

**What the theorem does not say.** The set of all global displacement fields
can have further compatibility/period requirements. We classify stress fields
under harmonic-trace compatibility, not displacement single-valuedness or
traction boundary conditions. The five-dimensional space includes q=0, where
principal directions cease to be distinguished. A prescribed nondegenerate
net is realized on the open set where R≠0.

## 2. All four global dimensions occur even for rational isothermal coordinates

For distinct real poles j=0,…,r−1 choose complex c_j, and set

    f(z)=z−Σ_j c_j/(z−j),
    g=f'=1+Σ_j c_j/(z−j)²,       h=g².               (4)

Delete the poles and zeros of g from the plane. On the resulting connected
Ω, f is a global holomorphic local coordinate, and Re f, Im f give regular
orthogonal isothermal foliations locally. Global injectivity is not asserted.
Zeros of h add no period constraints. The only possible periods are residues
at the poles of (4). Put

    G_j=1+Σ_{k≠j} c_k/(j−k)²,
    D_j=−2Σ_{k≠j} c_k/(j−k)³,
    r_j=Res_j h=2c_jD_j,
    s_j=Res_j zh=j r_j+2c_jG_j.

Multiplication by the common nonzero scalar −2π shows that each period row
on (a,u,v) is

    [ Im(s_j), Im(r_j), −Re(r_j) ].                 (5)

These exact choices realize every rank:

| coefficient list c_j | real period matrix | rank | global stress dimension |
|---|---|---:|---:|
| [1] | [0,0,0] | 0 | 5 |
| [i] | [2,0,0] | 1 | 4 |
| [1,i] | [2,4,0]; [0,−4,0] | 2 | 3 |
| [1,i,1+i] | [5/2,9/2,−1/2]; [6,0,4]; [−9/2,−9/2,−7/2] | 3 | 2 |

The last determinant is nonzero, as the exact verifier checks. Thus neither
the number of holes nor the local five controls the global dimension alone:
the period *rank* is the sufficient statistic. To retain genuine distinguished
stress axes globally on a bounded region containing these generators, start
with R=1 and perturb within the kernel slightly; compactness keeps R positive.
This establishes the same dimensions for a local open family of nondegenerate
aligned stresses, after choosing boundary curves away from poles/zeros.

## 3. The full-state loop law

Choose the parameter vector θ=(a,u,v,c,d) at the basepoint. Analytic
continuation around γ fixes its first four entries and sends

    d ↦ d−2P_h(a,u,v)(γ).

Thus the monodromy matrix is M_γ=I+e_5 ℓ_γ, where the last two entries of
ℓ_γ are zero. Every nontrivial M_γ has rank-one deviation from identity,
and

    (M_γ−I)(M_η−I)=0,
    M_γ M_η=M_(γη)=I+e_5(ℓ_γ+ℓ_η).

All holonomy is therefore abelian (nilpotency class at most one), with a
square-zero deviation algebra in this representation,
even when the fundamental group of the punctured plane is nonabelian. The
stress-space dimension is the dimension of the common fixed subspace. This
is an exact finite-state loop law: it neither removes the topology nor
replaces it by merely counting holes. Distinct nonhomologous loops can be
invisible, and all commutators are invisible. To recover the underlying net
or fundamental group requires additional information.

## 4. The bridge to prior repository work

The mechanism is **local compatibility → finite carrier → global holonomy
annihilator**. It connects to, but does not identify, three different subjects:

- The [pentagon context note](../../pentagon_context_number_fields_20261011.md)
  keeps all five omitted colour pairs, rather than just five local list sizes.
  The joint diagonal condition is its missing global sidecar. Our period rows
  play that sidecar role, but no map between stress and colouring solutions is
  asserted.
- The ready-return pentagon keeps full continuation data around a loop; a scalar
  readout alone can forget what changed. Here a local primitive is valid but
  changes by a period on a loop. The common principle is exact loop transport,
  not the numeral five; see the session's holonomy lane.
- Session281's all-modulus but power-sparse language shows that satisfying every
  congruence need not place an absolute integer target. Our local stress law
  similarly requires a separate global gluing test. The distinction is exact;
  a theorem transporting one obstruction to the other is not supplied.

**Conditional design principle.** If a tournament construction admits a
finite-dimensional linear continuation carrier and lawful closed-word moves,
compute the subgroup of carrier defects first, then seek a physical family
inside its annihilator. This is a research prescription, not evidence that
all positive odd targets are covered. Physical binary choices and positivity
can make the resulting feasible set disconnected or empty even when the
linear annihilator has large dimension.

## 5. Reproduction

`python3 verify_periods.py` checks residue formulas independently against
symbolic Laurent residues, all ranks, determinant certificates, and the free
c,d modes. It does not need a numerical contour integrator or floating ranks.
