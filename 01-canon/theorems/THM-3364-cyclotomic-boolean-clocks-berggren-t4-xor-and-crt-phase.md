---
id: THM-3364
title: "Cyclotomic Boolean clocks, Berggren T4 XOR, and CRT phase"
status: >
  PROVED analytic structure + VERIFIED-EXACT.  Ultimately periodic subsets of
  the positive integers modulo finite symmetric difference are exactly the
  clopen subsets of the profinite integers.  Their full finite cyclotomic
  Fourier transforms are convolution idempotents: they recover Boolean
  operations, cyclic phase, least eventual period, the complete periodic
  Dirichlet tail, and hence the harmonic finite part.  Harmonic density is
  only the zero mode.  The three Fibonacci/Berggren ancestry rays and binary
  Cassini clock factor by CRT into six phase-distinguished channel states.
  Separately, every Berggren parent with its three children is a labelled
  transitive T4 whose moving outer arc is a reset/XOR bit; its two depth-count
  divisibility clocks have equal harmonic residue but opposite quarter-phase.
  Finite startup, general subsets of N, actual triples, unrestricted ancestry,
  owners, global tournament edges, LRC phase, and physical current are not
  recovered.
source: codex-2026-08-14-cyclotomic-boolean-clocks
depends_on:
  - THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction
  - THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit
  - THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar
related:
  - THM-2438-poisson-newton-ternary-half-and-harmonic-divisor-incidence
  - THM-3345-berggren-gaussian-source-groupoid-prime-toggle-circuit-and-ancestry-cost
  - THM-3358-admissible-composite-parabolic-compiler-and-hensel-normal-offset-atlas
  - THM-3363-d14-complement-clock-small-lrc-terminal
script: 04-computation/cyclotomic_boolean_clocks_berggren_t4_xor_thm3364.py
output: 05-knowledge/results/cyclotomic_boolean_clocks_berggren_t4_xor_thm3364.out
script_sha256: 0bf40f7736ca65da2704f0f932b223a5946ca75275bc260b4aa8ee84718ad293
output_sha256: f4b06bd63509e23b390f1133c0bd2a1442c709b8f9cdbe4a04a5bdcff95cbecc
hash_basis: LF-normalized bytes
---

# THM-3364 -- cyclotomic Boolean clocks and the local Berggren `T4`

**PROVED analytic structure + VERIFIED-EXACT.**

## 1. Inheritance and the missing coordinate

[THM-3359](THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar.md)
proves that every modular support of an integer C-finite sequence is ultimately
periodic.  It assigns that support its density, harmonic pole, and deletion
scar, but explicitly shows that density loses tournament orientation and
Berggren ancestry.  The closest positive sidecars are
[THM-3339](THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md),
whose three ancestry rays and Cassini parity form a six-state clock, and
[THM-3357](THM-3357-berggren-three-branch-walsh-level-collapse-and-parent-circuit.md),
whose sibling comparison has one reset/flip bit.

The missing coordinate is not another scalar density.  It is the phase of
every finite cyclic Fourier mode.

| field | connection contract |
|---|---|
| source | an ultimately periodic support `H subset N_(>0)`, modulo finite symmetric difference |
| target | a finite-support convolution idempotent on `Q/Z` |
| map | Fourier transform of the corresponding clopen indicator on `Zhat` |
| preserved | eventual membership, Boolean operations, least period, density, cyclic phase, periodic Dirichlet tail and harmonic finite part |
| destroyed | finite startup unless retained separately; values, multiplicities, words, owners and geometry attached to indices |
| hostile | all residue cosets have the same Fourier magnitudes; arbitrary subsets need not extend to clopens; transpose-symmetric walk data still lose directed observers |

No literature-priority claim is made for this elementary finite Fourier
calculus.  The contribution here is the typed interface and its two exact
Berggren applications.

## 2. The cyclotomic Boolean idempotent

Let `H subset N_(>0)` be ultimately periodic.  Choose a period `p` and a
periodic extension

```text
f_H : Z/pZ -> {0,1}                                     (1)
```

that agrees with `1_H(n)` for all sufficiently large `n`.  Changing finitely
many terms of `H`, or replacing `p` by a multiple, does not change the clopen
indicator on the profinite integers `Zhat`.

For `alpha in Q/Z`, define

```text
Phi_H(alpha)=integral_(Zhat) f_H(x) conjugate(chi_alpha(x)) dx, (2)
```

with normalized Haar measure.  Equivalently,

```text
Phi_H(j/p)=(1/p) sum_(r=0)^(p-1) f_H(r) exp(-2 pi i j r/p), (3)
```

and `Phi_H` vanishes outside `(1/p)Z/Z`.  Thus `Phi_H` has finite cyclotomic
support.  Fourier inversion recovers the complete eventual indicator:

```text
f_H(r)=sum_(alpha in Q/Z) Phi_H(alpha) chi_alpha(r).     (4)
```

On finite-support functions on `Q/Z`, use convolution

```text
(Phi*Psi)(alpha)=sum_beta Phi(beta)Psi(alpha-beta).      (5)
```

Pointwise multiplication before Fourier transform becomes (5).  Therefore

```text
Phi_(H intersect K) = Phi_H*Phi_K,
Phi_(H xor K)       = Phi_H+Phi_K-2 Phi_H*Phi_K,
Phi_(H union K)     = Phi_H+Phi_K-Phi_H*Phi_K,
Phi_(complement H)  = delta_0-Phi_H,                    (6)
Phi_H*Phi_H         = Phi_H.                            (7)
```

Conversely, the inverse transform of any finite-support convolution
idempotent is pointwise idempotent and hence `{0,1}`-valued.  Equations
(2)--(7) therefore identify the Boolean algebra of clopens of `Zhat` with the
finite-support convolution idempotents on `Q/Z`.

The zero mode is precisely

```text
Phi_H(0)=delta(H),                                      (8)
```

the rational cycle density and harmonic residue in THM-3359.  It is only one
coordinate of the idempotent.

Let `ord(alpha)` be the additive order of `alpha in Q/Z`.  Fourier inversion
also gives the least eventual period:

```text
p_min=lcm{ord(alpha): Phi_H(alpha)!=0},                 (9)
```

with lcm of the empty set defined as one.  Indeed a translation fixes `f_H`
exactly when every supported character is trivial on that translation.

## 3. Dirichlet and harmonic profiles retain phase

Let `f_H(n)` now denote the periodic extension for all positive `n`, and put

```text
P_H(s)=sum_(n>=1) (1_H(n)-f_H(n))/n^s.                  (10)
```

The sum in (10) is finite.  For `Re(s)>1`, termwise Fourier inversion gives

```text
D_H(s)=sum_(n in H) n^(-s)
      =P_H(s)+sum_alpha Phi_H(alpha) Li_s(exp(2 pi i alpha)). (11)
```

The `alpha=0` term is `delta(H) zeta(s)`.  Every nonzero root-of-unity
polylogarithm in (11) is regular at `s=1`.  Consequently the harmonic partial
sums satisfy

```text
sum_(n<=x,n in H) 1/n = delta(H) log x+C_H+O(1/x),      (12)
```

where

```text
C_H=P_H(1)+delta(H) gamma
    +sum_(alpha!=0) Phi_H(alpha)[-log(1-exp(2 pi i alpha))]. (13)
```

Conjugate modes make (13) real.  Thus the zero mode decides logarithmic
divergence, while the nonzero modes and the finite transient decide the
renormalized harmonic constant.  This is the precise sense in which a subset
of the natural numbers is a selected subseries of the harmonic series: its
indicator, not merely the resulting scalar mass, is the object being
transported.

The restriction matters.  An arbitrary subset of `N` need not have an
eventual period or a canonical clopen extension to `Zhat`; its Fourier object
need not be finite.  THM-3359's square support (finite harmonic mass) and
`ceil(k log(k+1))` support (density zero but divergent harmonic mass) remain
decisive hostiles to extending (8)--(13) to all subsets.

## 4. The exact `2 x 3` Fibonacci/Berggren clock

THM-3339 classifies the three primitive Fibonacci-index ancestry rays as

```text
R_2=(BA)^*,       R_0=A(BA)^*,       R_1=C(BC)^*,       (14)
```

occupying indices `k=2,0,1 mod 3`, respectively.  These are three selected
ancestry languages, not one `C3` action on the Berggren tree.  Modulo finite
startup, their spectra satisfy, with `omega=exp(2 pi i/3)`,

```text
Phi_(R_a)(j/3)=omega^(-ja)/3.                           (15)
```

Cassini parity supplies

```text
P_b={k:k=b mod 2},       Phi_(P_b)(1/2)=(-1)^b/2.       (16)
```

For every `t mod 6`, put

```text
S_t=R_(t mod 3) intersect P_(t mod 2)={k:k=t mod 6}.    (17)
```

The coprime clocks in (17) give the exact convolution/CRT factorization

```text
Phi_(S_t)=Phi_(R_(t mod 3))*Phi_(P_(t mod 2)),          (18)
Phi_(S_t)(j/6)=zeta_6^(-jt)/6.                         (19)
```

THM-3339 obtains the fixed channel dictionary from consecutive oriented Farey
flanks (hence from reduced fractions), namely

```text
t:       0    1    2    3    4    5
order:  cab  cba  bca  bac  abc  acb.                  (20)
```

All six states have density `1/6`, the same Fourier magnitudes, and the same
THM-3359 scar coefficient `11/72`.  The phases in (19), together with the
fixed Fibonacci index carrier and dictionary (20), distinguish all six.
Discarding phase, the dictionary, or the carrier destroys that conclusion.
Nothing in (18) chooses a `K4` owner or extends these three languages to every
Berggren word.  Here six counts the `S3` order states of three channels, not
six tournament vertices; THM-3339's separate edge-product `T6` needs its own
weighted comparison data.

## 5. Every Berggren parent and its children form a local `T4`

Use THM-3357's positive Euclid parameter

```text
u=(m,n),       0<m<n,
(a,b,c)=(n^2-m^2,2mn,n^2+m^2),                         (21)
```

and its left, middle, and right children `Lu,Mu,Ru`.  Let

```text
N_L=||Lu||^2,       N_M=||Mu||^2,       N_R=||Ru||^2.  (22)
```

The parent hypotenuse is `c`, and direct subtraction gives

```text
N_L-c=4n(n-m)>0,
N_M-c=4n(n+m)>0,
N_R-c=4m(m+n)>0,                                      (23)
N_M-N_L=4b>0,
N_M-N_R=4a>0,
N_R-N_L=4(b-a).                                       (24)
```

The equality `a=b` would force `c^2=2a^2`, so it never occurs.  Orienting
toward larger hypotenuse makes the labelled quartet

```text
{parent,L,M,R}                                         (25)
```

a transitive tournament.  The parent is always lowest and `M` always highest;
only the labelled outer order moves:

```text
b>a:  parent < L < R < M,
b<a:  parent < R < L < M.                              (26)
```

This makes the user's local size-four tournament intuition exact.  It is
local, not a tournament on a whole level or on the entire ternary tree.

Let `sigma(u)=sign(b-a)` and encode the outer order by

```text
epsilon(u)=(1-sigma(u))/2 in {0,1}.                    (27)
```

THM-3357's sign automaton becomes the reset/XOR system

```text
epsilon(Lu)=0,
epsilon(Mu)=epsilon(u) xor 1,
epsilon(Ru)=1.                                         (28)
```

Thus the ternary tree is a gluing of labelled transitive `T4` chambers along
parent/child vertices, with one Boolean orientation sidecar.  If labels are
forgotten, the two tournaments in (26) are isomorphic.  If only the local
relation is retained, comparisons between unrelated chambers are missing.
If every pair is instead compared by hypotenuse, nonlocal ties remain: the
distinct primitive parameter nodes `(1,8)` and `(4,7)` both have hypotenuse
`65`.  Recording equality as a missing edge or as both directed arcs preserves
the tie but leaves the category of tournaments; breaking it creates a gauge.

## 6. Recurrence, mod-five supports, and opposite harmonic phase

Let `A_d` and `B_d` count depth-`d` parents with `epsilon=0` and `epsilon=1`.
The root `(1,2)` has `(A_0,B_0)=(1,0)`.  From (28),

```text
[A_(d+1)]   [1 2][A_d]
[B_(d+1)] = [2 1][B_d],                                (29)
```

so

```text
A_d=(3^d+(-1)^d)/2,
B_d=(3^d-(-1)^d)/2,                                    (30)
X_(d+2)=2X_(d+1)+3X_d.                                 (31)
```

Modulo five, (29) has the exact four-cycle

```text
(1,0),(1,2),(0,4),(3,4),(1,0).                         (32)
```

Shift depth to the positive index `q=d+1` and define

```text
H_A={q:5 divides A_(q-1)}={q:q=3 mod 4},
H_B={q:5 divides B_(q-1)}={q:q=1 mod 4}.                (33)
```

Both supports have harmonic residue `1/4` and THM-3359 scar coefficient
`7/32`.  Density and scar cannot distinguish the two labelled local-`T4`
types.  Their quarter-modes do:

```text
Phi_(H_A)(1/4)= i/4,
Phi_(H_B)(1/4)=-i/4.                                   (34)
```

Indeed their harmonic constants satisfy

```text
C_(H_B)-C_(H_A)=pi/4.                                  (35)
```

This follows directly from (13), or from the reflection identity
`psi(3/4)-psi(1/4)=pi` for the two harmonic progressions.

Their symmetric difference is the odd positive integers:

```text
H_A xor H_B={q:q=1 mod 2},                              (36)
```

with residue `1/2`, scar `3/8`, and sole nonzero nonconstant mode
`Phi(1/2)=-1/2`.  Equations (33)--(36) are depth-count divisibility clocks.
They do not select individual tree nodes or reconstruct a branch word.

## 7. Exact verification and stopping boundary

The companion constructs cyclotomic polynomials and exact finite DFTs without
floating point.  It exhausts every Boolean word of periods `1` through `10`
and checks inversion, (6)--(9), termwise Dirichlet decomposition, and minimal
period recovery.  It separately verifies all six CRT spectra (18)--(20), the
local `T4` formulas and automaton through a bounded injective Berggren tree,
the recurrence and mod-five clocks (29)--(36), and the hypotenuse-65 hostile.
Ordinary and optimized outputs match exactly; hashes and census are frozen in
the frontmatter and stored transcript.

The theorem proves a complete invariant of the **eventual periodic support**,
not of any object indexed by that support.  It cannot recover finite startup,
multiplicity, actual Pythagorean triples, unrestricted ancestry, owners, LRC
phase/current, JC flux, or global tournament arcs.  THM-3363's complement-clock
terminal is operation-compatible at the level of a finite residue word, but
its pointwise geometric carrier identity is an indispensable extra sidecar;
Fourier Boolean algebra alone does not prove an LRC cover or exclusion.
