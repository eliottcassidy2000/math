---
id: THM-2355
title: "Component-deletion Gram reconstruction and twist-energy phase transport"
status: >
  PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT. A finite
  grouped complex current is reconstructed at the energy level by its
  singleton and complete pair-union energies; equivalently its Boolean
  energy has interaction degree at most two. For every cyclic twist of
  order p>=3, the first energy Fourier mode is the exact complex relative
  phase transport, so twisted pair energies on a spanning tree reconstruct
  all component phases up to one common phase. Untwisted tree energies do
  not. Target-character energy recovers the autocorrelation of the
  residue-grouped currents, not the currents themselves; a real
  full-support C_13 perfect-autocorrelation array is an exact hostile.
  If two arrays lie in phase cones whose widths sum to less than pi, their
  cross-correlation support is their full support difference. This gives a
  precise quadratic alternative to the THM-2303 phase tree and sharpens
  the THM-2344 inverse boundary, but no lawful canonical pair probe, cone
  confinement, LRC target landing, scalar-row exclusion, or LRC(14)
  closure is proved.
source: codex-2026-07-25-component-energy-phase-transport
depends_on: []
related:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2216-residual-capacity-hinge-gram-law
  - THM-2303-terminal-component-phase-current-and-defect-rank
  - THM-2344-correlation-inverse-rigidity-and-aligned-tooth-twist-hostile
script: 04-computation/component_deletion_gram_twist_energy_thm2355.py
output: 05-knowledge/results/component_deletion_gram_twist_energy_thm2355.out
script_sha256: 4a01e5da66f5927023b20b2953c37b20f562fc20b3f98724b85536eb5a0bc330
output_sha256: 1b8871c34b9f27a43bf5748ecf5bf360544b5eb5d5f4107e5260b42b56a345a0
hash_basis: working-tree bytes (LF)
---

# THM-2355 -- component-energy phase transport

**PROVED + VERIFIED-EXACT CANDIDATE UNDER INDEPENDENT AUDIT.**

THM-2303 identifies the current obstruction exactly: rooted energies retain
component sizes, while the desired Fourier coefficient is the complex sum
of the component currents.  Its phase-tree lemma repairs this by retaining
complex phase ratios.  There is a complementary repair which is often
easier to ask an exact finite packet to supply:

```text
retain energies after lawful component unions or cyclic component twists.
```

Polarization turns those quadratic observables back into the missing Gram
entries.  This is the complex-current analogue of THM-2216's cut-metric
polarization and the finite-query specialization of THM-2176's continuation
principle.

## 1. The complete pair-energy formula

Let `E` be a finite nonempty set, `n=|E|`, and let

```text
z_e in C,                         e in E.
```

For every subset `S` put

```text
Q(S)=|sum_(e in S) z_e|^2.                            (1)
```

For a nonempty `T subset E`, define its Boolean Mobius coefficient

```text
Delta_T Q
 =sum_(S subset T)(-1)^(|T|-|S|)Q(S).                (2)
```

Then exactly

```text
Delta_{e} Q=|z_e|^2,

Delta_{e,f} Q=2 Re(z_e conjugate(z_f)),

Delta_T Q=0                         when |T|>=3.      (3)
```

Indeed, expansion of (1) gives

```text
Q(S)
 =sum_(e in S)|z_e|^2
  +2 sum_({e,f} subset S)Re(z_e conjugate(z_f)),     (4)
```

and Boolean inversion gives (3).  Thus the energy table has interaction
degree at most two even though the component phases themselves are
continuous.

Summing (4) over the complete graph gives the exact reconstruction

```text
Q(E)
 =sum_({e,f} subset E) Q({e,f})
  -(n-2)sum_(e in E)Q({e}).                          (5)
```

Consequently all singleton energies and all pair-union energies determine
whether the grouped current vanishes.  They do not reconstruct a preferred
global phase, which is irrelevant to that verdict.

There is also a single-deletion identity.  Write

```text
Q_{-e}=Q(E minus {e}).
```

Then

```text
sum_e Q_{-e}
 =(n-2)Q(E)+sum_e Q({e}).                            (6)
```

For `n!=2`, (6) reconstructs `Q(E)`.  The exceptional value `n=2` is
sharp: the two configurations

```text
(z_1,z_2)=(1,1),                 (1,-1)              (7)
```

have identical singleton and deletion ledgers, while their grouped energies
are `4` and `0`.

Complete pair data in (5) cannot be replaced by untwisted pair energies on
a spanning tree.  On the path `1--2--3`, the currents

```text
(i,1,i),                         (i,1,-i)             (8)
```

have the same three singleton energies and the same two edge-union energies,
but total energies `5` and `1`.  The missing chord carries a real Gram
entry.  This is the quadratic version of the distinction in THM-2303
between a connected graph of phase transports and a graph of bare service
edges.

## 2. A cyclic twist recovers the complex phase transport

Let `p>=3`, let

```text
zeta=exp(2*pi*i/p),
```

and for two component currents define the real twist-energy response

```text
Q_(e,f)(t)=|z_e+zeta^t z_f|^2,       t in Z/pZ.      (9)
```

Use the normalized Fourier convention

```text
Qhat(k)=1/p sum_t Q_(e,f)(t)zeta^(-kt).              (10)
```

Expansion gives

```text
Q_(e,f)(t)
 =|z_e|^2+|z_f|^2
  +zeta^(-t)z_e conjugate(z_f)
  +zeta^t conjugate(z_e)z_f.                         (11)
```

Character orthogonality, and the fact that `1!=-1 mod p`, therefore give

```text
Qhat(0)=|z_e|^2+|z_f|^2,

Qhat(1)=conjugate(z_e)z_f,

Qhat(-1)=z_e conjugate(z_f),

Qhat(k)=0                         otherwise.         (12)
```

If both currents are nonzero, the normalized first mode

```text
Qhat(1)/(|z_e||z_f|)
```

is exactly THM-2303's phase transport from `e` to `f`.  Hence cyclic
pair-twist energy responses on the edges of any spanning tree reconstruct
all nonzero component currents up to one common element of `U(1)` and decide
whether their sum vanishes.  This is not an analogy: (12) literally
constructs every edge ratio used by the phase-tree lemma.

The condition `p>=3` is sharp.  At `p=2`, the `+1` and `-1` modes coincide,
so only the real cross term survives.  The pairs `(1,i)` and `(1,-i)` have
identical two-twist energy responses and opposite phase transports.

## 3. What a whole target-twist energy actually remembers

Let `G` be a finite abelian group and let

```text
Z:G->C
```

be the current already grouped by target residue.  For a character `chi`
put

```text
F(chi)=sum_(q in G)chi(q)Z(q),

E(chi)=|F(chi)|^2.                                   (13)
```

Define the periodic autocorrelation

```text
C(h)=sum_(q in G)Z(q+h)conjugate(Z(q)).              (14)
```

Then finite Fourier inversion gives

```text
E(chi)=sum_(h in G)chi(h)C(h),

C(h)=1/|G| sum_chi conjugate(chi(h))E(chi).          (15)
```

Thus the full twist-energy ledger reconstructs the residue
autocorrelation, not the individual grouped currents.  In particular

```text
E is constant
 iff C(h)=0 for every h!=0.                          (16)
```

Condition (16) is perfect periodic autocorrelation and does not force
singleton residue support.

There is an exact real full-support hostile already on `C_13`.  Put

```text
epsilon(1)=epsilon(-1)=-1,
epsilon(k)=1                     otherwise,

u(x)=1/13 sum_(k in C_13)epsilon(k)zeta^(kx).        (17)
```

Then

```text
u(0)=9/13,

u(x)=-4/13 cos(2*pi*x/13)          for x!=0.         (18)
```

Every entry is real and nonzero, because an odd-order root cannot have
zero real part.  Yet the Fourier transform of `u` is the sign vector
`epsilon`, so

```text
sum_x u(x+h)u(x)=delta_0(h),                         (19)
```

and all thirteen twist energies are equal.  Full target-residue support,
real-valuedness, and nonzero grouped currents therefore do not by themselves
force target variance.  This is the energy version of THM-2344's
perfect-autocorrelation boundary.

## 4. Acute phase cones make support honest

For arrays `U,V:G->C`, put

```text
R_(U,V)(h)
 =sum_q U(q+h)conjugate(V(q)).                       (20)
```

Assume that the arguments of the nonzero values of `U` lie in a circular
phase interval of width `alpha`, and those of `V` lie in one of width
`beta`, with

```text
alpha+beta<pi.                                      (21)
```

For a fixed `h`, every nonzero summand in (20) lies in one common open
half-plane: its phase belongs to the difference of the two phase intervals,
whose width is less than `pi`.  A nonempty sum of such vectors cannot
vanish.  Therefore

```text
support R_(U,V)
 =support(U)-support(V).                            (22)
```

This **cone-difference lemma** is sharp at the strict boundary.  On `C_3`,
take

```text
U(0)=V(0)=1,
U(1)=i,
V(1)=-i,
U(2)=V(2)=0.                                        (23)
```

The two cone widths are `pi/2`, their sum is `pi`, and

```text
0 belongs to support(U)-support(V),
R_(U,V)(0)=1+i*i=0.                                 (24)
```

For autocorrelation, one array in a phase cone of width less than `pi/2`
satisfies

```text
support C=support(Z)-support(Z).                    (25)
```

Hence a nonzero cone-confined residue current has constant twist energy if
and only if its support is a singleton.  This is the precise positive
replacement for the false rule "full residue support forces variance."

## 5. LRC connection and loss ledger

At a fixed frequency in THM-2303, take

```text
z_e=I_e(n).
```

Then (5) says that lawful singleton and complete pair-union Fourier energies
would decide the grouped coefficient without choosing phase branches.
Alternatively, (12) says that lawful thirteen-twist pair energies on a
spanning tree would manufacture the missing phase transports exactly.
Current canon supplies neither complete probe family; the theorem names the
minimal quadratic service to seek.

For THM-2334's target-character family, group the atomic current first by
its target residue `q`.  Equation (15) is then the exact energy-side
pushforward of the `169` twists.  It preserves residue-lag correlations and
forgets the individual components inside each residue.  The real hostile
(17)--(19) proves that even full grouped residue support is not a sufficient
replacement for coefficient-sensitive dispersion.

For THM-2344's endpoint arrays, (22) gives a sharp new conditional
separator.  If their nonzero values lie in phase cones with width sum less
than `pi`, a shifted-delta cross-correlation forces both endpoint supports
to be singletons.  THM-2323's acute-sector Galois argument suggests the
right geometry, but it acts on a fixed-colour physical cross-correlation;
no current theorem intertwines that cone with the grouped target-residue
arrays.  That root/target cone intertwiner, or the lawful pair-twist service
in (12), is the remaining sidecar.

The complete map is

```text
component currents
  -> singleton energies
     loses every relative phase;

component currents
  -> complete pair-union energies
     preserves the real Gram matrix and total vanishing;

component currents
  -> cyclic pair-twist energies
     preserves complex phase transport;

residue-grouped currents
  -> whole target-twist energies
     preserves only periodic autocorrelation;

acute phase cone
  -> support-exact autocorrelation
     rules out signed CAZAC cancellation.                        (26)
```

No canonical LRC row is proved to supply the required pair probes or cone.
No scalar row is excluded and LRC(14) remains open.

## 6. Exact companion

The dependency-free exact companion uses rational Gaussian currents and an
explicit `Q[zeta_p]` model.  It checks:

- `56` complete-pair polarization and deletion identities;
- `448` vanishing Boolean interactions of order at least three;
- the two-component deletion and untwisted-tree hostiles;
- `84` cyclic twist-DFT coefficients at `p=3,5,13`;
- `42` exact group-energy/autocorrelation inversions;
- every entry and lag of the real full-support `C_13` hostile;
- an acute-cone support-difference example; and
- the width-sum-`pi` cancellation boundary.

Reproduce with

```bash
python3 04-computation/component_deletion_gram_twist_energy_thm2355.py
python3 -O 04-computation/component_deletion_gram_twist_energy_thm2355.py
```

Both transcripts must match

```text
05-knowledge/results/component_deletion_gram_twist_energy_thm2355.out
```

byte-for-byte after LF normalization.  Independent audit is pending.
