---
id: THM-1226
title: THE GCD-PERIOD PROJECTIVE-CHARGE OBSTRUCTION — exact vertex loads, a factorial strict-high counterfamily, and the disconnected-spectrum rescue
status: PROVED (exact arithmetic obstruction and protected-needle embedding; conditional disconnected-G_gt transfer using THM-1221's complete projective branches)
source: codex-2026-07-19-S82
depends_on: [THM-1166, THM-1221]
related: [HYP-7678, HYP-7870]
script: 04-computation/lrc14_gcd_period_projective_charge_obstruction_referee_codex_S82.py
output: 05-knowledge/results/lrc14_gcd_period_projective_charge_obstruction_referee_codex_S82.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCGCDPeriodProjectiveCharge.lean
script_sha256: c61bc3536c4360737b7d796baa7e7f5c3695ee8cd80eb64a6461c09ca2501b46
output_sha256: 4b3d53a98eff6e924a56a8ef9104a317e751597ab82ae5f86e24f5dbf11df4ac
formalization_sha256: 41180f856d93579ee37053c44d8fb141c01cbe7ae7b329a0e879d873adcf5426
---

# THM-1226 — the gcd-period projective-charge obstruction

## Statement

For an edge `s_i=g_ij a_ij`, `s_j=g_ij b_ij`, with
`g_ij=gcd(s_i,s_j)` and `(a_ij,b_ij)=1`, put

```text
eta_ij=rho_ij(1-rho_ij),
kappa_ij=eta_ij ab/(a+b).
```

Then the periodic positioning error factors exactly as

```text
eta_ij/g=kappa_ij(1/s_i+1/s_j).                         (1)
```

For a tree `T`, define the projective load at vertex `i` by

```text
Lambda_i(T)=sum_(ij in T) kappa_ij.                      (2)
```

Writing `H=sum_i 1/s_i`, the total THM-1166 period error is therefore

```text
E_T=sum_(ij in T) eta_ij/g_ij
   =sum_i Lambda_i(T)/s_i.                               (3)
```

There is no absolute finite `C` such that every seven-speed packet admits a
tree satisfying both

```text
sum_(ij in T) rho_ij >=15/154,
E_T<=C H.                                                (4)
```

This remains false even when **every** spanning tree clears the first
inequality and the strict-high graph `G_gt={rho>1/63}` is complete.

The exact finite projective branches of THM-1221 do give a complementary
positive theorem.  If `G_gt` is disconnected, some THM-1221 floor tree obeys

```text
E_T<=C_disc H,        C_disc=448916/194775.              (5)
```

Consequently every interval `I` of length `L` satisfies the conditional
positioned-tree bound

```text
sum_(ij in T) mu(I intersect D_si intersect D_sj)
 >=(15/154)L-(448916/194775)H.                           (6)
```

Thus the connected `G_gt` branch is exactly the unbounded-projective-height
frontier of this gcd-period method.

## 1. Exact charge and oriented forms

Since `s_i=g a` and `s_j=g b`,

```text
1/s_i+1/s_j=(a+b)/(gab).
```

Multiplication by `eta ab/(a+b)` proves (1).  Summing it over a tree proves
(3).  Equivalently, after rooting `T` and charging every edge to its child,

```text
E_T=sum_(i != root)
 [eta_(i,p(i)) * (s_i/gcd(s_i,s_p(i)))]/s_i.             (7)
```

This gives a precise sufficient replacement for the false universal claim.
If a rho-heavy rooted tree has every tail height

```text
s_i/gcd(s_i,s_p(i))<=K,
```

then `rho<=1/14` gives `eta<=13/196`, and hence

```text
E_T<=(13K/196)H.                                        (8)
```

The relevant carrier is therefore not merely a low/high edge color.  It is a
reduced ratio channel together with its projective height and an oriented
vertex-load assignment.

## 2. The translated-factorial counterfamily

Let

```text
A=27720 N,       N>=1,
R={1,2,3,5,7,11,13},
S_A={A+r:r in R}.                                       (9)
```

Every difference of two elements of `R` divides `27720`, while the elements
of `R` are pairwise coprime.  If `r<s`, then

```text
gcd(A+r,A+s)
 =gcd(A+r,s-r)
 =gcd(r,s-r)
 =gcd(r,s)=1.                                           (10)
```

Thus all twenty-one pair periods are the whole circle.  THM-1166's folded
formula and defect bound give, for every edge,

```text
rho(A+r,A+s)
 =1/49+[F(r+s)-F(s-r)]/[196(A+r)(A+s)]
 >x_A:=1/49-1/(4A^2).                                  (11)
```

Here `A=0 mod 14`.  The exact twenty-one folded corrections, in lexicographic
pair order on `R`, are

```text
(20,16,8,0,-16,-24,32,16,0,-32,-20,
 24,0,-48,-16,0,-24,-8,0,0,16).                        (12)
```

Since

```text
1/49-5/308=9/2156,
```

the inequality `9A^2>539` implies `x_A>5/308`.  It already holds at
`A=27720`.  Hence `G_gt=K_7`, and every spanning tree satisfies

```text
sum_(ij in T)rho_ij>6*(5/308)=15/154.                   (13)
```

On the other hand every `g_ij=1`.  Because `rho<1/2`, the map
`x -> x(1-x)` is increasing on the relevant range, so

```text
E_T>6 x_A(1-x_A)>6(5/308)(303/308)=4545/47432.          (14)
```

But

```text
H_A=sum_(r in R)1/(A+r)<7/A.                            (15)
```

Consequently, uniformly over **all** spanning trees,

```text
E_T/H_A >(6A/7)x_A(1-x_A)
         >(4545/332024)A -> infinity.                   (16)
```

Given any proposed finite `C`, choosing
`A>332024 C/4545` contradicts the second inequality in (4).  More precisely,
the fixed correction bank (12) gives the uniform asymptotics

```text
rho_ij=1/49+O(A^-2),
E_T=288/2401+O(A^-2),
H_A=7/A-42/A^2+O(A^-3),
E_T/H_A=(288/16807)A+O(1).                              (17)
```

No choice of maximum tree, Pareto tree, root, or charge orientation repairs
this family: both reduced endpoint heights on every edge are order `A`.

## 3. A protected covered needle where the error loses the answer

The obstruction above is to the unsigned period-error proof, not to HYP-7678's
phase-aware F7 statement.  This can be seen inside the exact protected-needle
geometry.  Put

```text
q=A+1,
S_A=q+{0,1,2,4,6,10,12},
W={ (3q+1)/2+k : 0<=k<=5 },
m=(3q+11)/2,
I=[1/q-1/(14m), 1/q+1/(14m)].                           (18)
```

The number `q` is odd, so the six core speeds are integers, and
`|I|=1/(7m)`.  The core's uniform clearance beyond `1/14` is bounded below by

```text
(5q-77)/(14q)>0                                         (19)
```

for `q>=17`.  For every deleted speed `q+d`, `d<=12`, its distance from the
integer center throughout `I` is at most

```text
12/q+(q+12)/[7(3q+11)].                                 (20)
```

The clearance of `1/14` above (20) is exactly

```text
(q^2-517q-1848)/[14q(3q+11)].                           (21)
```

The numerator is `-288` at `q=520` and `236` at `q=521`, and increases
thereafter.  Our first value is `q=27721`.  Hence `W` is safe throughout `I`
while **all seven** deleted danger combs contain `I`.  Every one of the
twenty-one restricted pair intersections equals `I`, so every local tree has
weight `6|I|`.

Yet the separately summed period-error envelope remains order one and gives
no positive edge credit.  Indeed, when `g=1`,

```text
L rho-rho(1-rho)=rho(L-1+rho).
```

For `L<=139/154` and `rho<=1/14`, this is at most

```text
rho(139/154-13/14)=-2rho/77<0.                          (22)
```

Thus the gcd-period abstraction permits incompatible worst endpoint phases
on six tree edges simultaneously.  The actual wall events in (18) are instead
perfectly aligned.  This proves that the negative result does not refute F7;
it identifies information destroyed before F7 is asked.

## 4. The disconnected strict-spectrum rescue

THM-1221 classifies every disconnected `G_gt` packet into finite projective
branches.  Identity (1) turns those ratio banks directly into positioning
constants.

If the weak-high graph `G_ge={rho>=1/63}` is disconnected, normalization at
its singleton component places all vertices in

```text
V_<={1} union R_<,
```

where `R_<` is the fourteen oriented strict-low ratios.  Exact evaluation of
all pairs in this fifteen-vertex bank gives

```text
max kappa=85975/342804,                                 (23)
```

attained on reduced channel `20:33`.  Since a seven-vertex tree has degree at
most six, (3) gives

```text
E_T<=(85975/57134)H.                                    (24)
```

Now suppose `G_ge` is connected but `G_gt` is disconnected.  In the `1+6`
cut, all vertices lie in `{1}` plus the twenty-four oriented closed-low
ratios.  The complete twenty-five-vertex bank has

```text
max kappa=224458/584325                                 (25)
```

at the vertex pair `(5/9,9/5)`, whose mutual reduced channel is `25:81` and
whose pair mass is `97/4725`.  Therefore

```text
E_T<=(448916/194775)H.                                  (26)
```

THM-1221's `2+5` branch has only twelve normalized packets.  Their complete
pair bank has

```text
max kappa=43774/276507,
E_T<=(87548/92169)H.                                    (27)
```

The `3+4` branch is empty, and both constants in (24),(27) are smaller than
(26).  This proves (5).  Summing THM-1166's interval discrepancy lower bound
on the THM-1221 floor tree proves (6).

If the seven combs cover `I`, Hunter supplies the reverse local tree bound
`<=H/7`.  Thus the disconnected branch has the explicit harmonic crown

```text
H >= [ (15/154)/(448916/194775+1/7) ] L
  = (417375/10488302)L.                                 (28)
```

The connected strict-high family (9) shows exactly why (28) does not extend
through the last component branch using gcd periods alone.

## 5. Tournament and alternate-vertex audit

For (9), switching at `1/63` marks all twenty-one edges strict-high.  The
speed-order gauge is transitive: score histogram `(0,1,2,3,4,5,6)`, no
directed cycles, seven singleton SCCs, one Hamiltonian path, and zero low-edge
flips.  These fingerprints are constant while `E_T/H` diverges.  Thus the
runner tournament preserves global connectivity and the `15/154` Hunter
credit but destroys the projective height that controls localization.

Using gcd periods as vertices is worse on this family: every edge becomes the
same period-one vertex.  That quotient destroys the common wall alignment in
(18).  Ratio channels retain `a:b` and `ab/(a+b)` but still omit interval
phase.  The information-preserving local candidates are

```text
vertices = wall-crossing events or tooth addresses,
edge sidecar = signed endpoint cocycle
               delta_ij(I)=L rho_ij-mu(I intersect D_i intersect D_j),
hypergraph = simultaneous owner sets and the THM-1218 deletion circuit.      (29)
```

The challenged assumption is that a globally heavy runner tree can be
localized by summing absolute per-edge period oscillations.  Equation (18)
shows the errors are correlated signed events, not independent capacities.

## 6. Formalization and verification boundary

The dependency-free referee checks (1) in `2,691` integer pairs; freezes all
twenty-one corrections in (12); verifies four exact members of the infinite
family; proves the spectral, error, harmonic, and asymptotic-facing
inequalities; checks every core and deleted endpoint in (18); and recomputes
the strict, closed, and twelve-packet projective maxima in (23)--(27).

`LRCGCDPeriodProjectiveCharge.lean` kernel-checks the symmetric and oriented
charge identities, the seven-vertex load consumer, all displayed transfer
constants, the localized cleared-denominator consumer, the strict-high
quadratic threshold, and both protected-needle margin polynomials.  The
arbitrary-speed Haar formula, THM-1221 branch classification, and finite-bank
identification remain explicit external providers rather than being
overstated as formalized.

Frozen hashes:

```text
source         c61bc3536c4360737b7d796baa7e7f5c3695ee8cd80eb64a6461c09ca2501b46
output         4b3d53a98eff6e924a56a8ef9104a317e751597ab82ae5f86e24f5dbf11df4ac
formalization  41180f856d93579ee37053c44d8fb141c01cbe7ae7b329a0e879d873adcf5426
```
