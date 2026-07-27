---
id: THM-2460
title: "Idempotent semantic-word copy and word-block cosupport boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. A fixed
  THM-2305 pure/fork word already retained by
  a THM-2452 packet is a target-neutral Boolean coordinate and can be
  copied to the bare endpoint with zero additional loss. If only the
  unsplit return residual is retained, the three word strata enlarge
  the complete bank from 128 to 384 atoms, giving drift, eligible-
  energy, and coefficient-square floors D/147456, D/1916928, and
  D/3864526848. On a common C_13 root chart the future word is
  root-constant, so the THM-2457 cosupport graph is the disjoint union
  of three word blocks. A fixed block has edge floor M_sigma/16384;
  globally some block has edge floor M/49152. The copy preserves an
  already fixed owner, clock, and lawful source/word skew, but it is a
  redundant common filter rather than a second address-bearing word
  harmonic. A sharp rational split-base hostile separates selected
  drift from all same-word cosupport. The remaining sidecar is
  same-word root-image coupling or a word-conditioned exclusion of
  the THM-2456 uniform-offset locus. No scalar row is excluded.
source: codex-2026-07-26-semantic-word-copy
depends_on:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2349-first-depth-one-delayed-shallow-restart
  - THM-2365-lawful-target-coshift-and-h-drift-dichotomy
  - THM-2401-common-filter-endpoint-or-first-death-certificate
  - THM-2442-delayed-word-septimal-source-completion
  - THM-2452-indicator-idempotent-aggregate-endpoint-restoration
  - THM-2457-complete-atom-root-cosupport-graph-and-semantic-word-hostile
related:
  - THM-2445-twenty-four-cell-graft-owner-conditioning
  - THM-2448-right-endpoint-cospan-transition-atlas
  - THM-2456-two-root-replica-uniform-offset-boundary
script: 04-computation/lrc14_semantic_word_copy_thm2460.py
output: 05-knowledge/results/lrc14_semantic_word_copy_thm2460.out
script_sha256: 6d1bee5e09628407480fbf80930eb4299e602c58511df33c59aa4b96a16b6010
output_sha256: c73e99206b9faa631a2fa4e918f4136bf98572c2ecd5078c285b970340928270
hash_basis: working-tree bytes (LF)
---

# THM-2460 -- the future word is a free matched filter, not a new address

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2457 proves that the seven present local truth bits do not determine
the future THM-2305 word: one live strict complete atom contains both a
pure and a fork stratum. That hostile rules out inference. It does not
rule out retaining the actual word independently.

The actual word has a stronger property. At every positive
base-thirteen clock it is:

```text
Boolean,
target-neutral,
constant on each physical C_13 root fibre.                     (1)
```

Consequently, once a pure/fork word is already present in the packet,
one may copy that same word to the bare endpoint by idempotence without
another triangle-inequality loss. The operation preserves the semantic
label, but the copy is a redundant common filter. It must not be
reinterpreted as a second independent terminal-word harmonic.

## 1. The fixed-word idempotent slide

Fix a source owner `j`. Call the other two blocker labels `a,b`, and
put

```text
Sigma_j={{a},{b},{a,b}}.                              (2)
```

Use the intrinsic THM-2305 terminal strata

```text
Q_(j,sigma),                  sigma in Sigma_j.
```

They are pairwise disjoint and partition

```text
R_j=A_0 minus D_(c_j)                               (3)
```

up to null endpoints. Fix an integer clock

```text
k>=1,                 R=13^k,
```

and define the future-word filters

```text
W_sigma(x)=1_(Q_(j,sigma))(R x).                    (4)
```

Let `F_theta,E_theta,Delta_r` have the meaning of THM-2452 at one
lawful target/deep twist. Suppose a fixed semantic word `sigma_0` is
already retained:

```text
F_theta E_theta=F_theta,

F_theta W_(sigma_0)=F_theta.                        (5)
```

This is exactly the THM-2349/2365 word-retaining packet, and in the
source-completed branch it is the phasewise packet of THM-2442.

Let `{P_omega:omega in Omega}`, `|Omega|=128`, be the complete local
atom partition of THM-2452. Put

```text
H_omega(r,theta)
 =integral_T F_theta P_omega Delta_r,               (6)

H_((omega,sigma),(nu,tau))(r,theta)
 =integral_T
    F_theta P_omega W_sigma
    E_theta P_nu W_tau Delta_r.                     (7)
```

Boolean orthogonality gives the exact identity

```text
H_((omega,sigma),(nu,tau))
 =1_(omega=nu)1_(sigma=tau=sigma_0) H_omega.        (8)
```

Thus the already selected word may be placed on both endpoint filters
without changing the whole nonnegative table. There are still only
`128` effective atoms. The source owner `j`, clock `k`, and
pure-versus-fork label `sigma_0` are all retained.

The phrase **zero additional loss** is relative to the already
word-retaining table in (5). If one begins with THM-2452's bare-word
drift `D_0`, the existing THM-2452 BV invoice

```text
D_word>=mu(Q)^2 D_0/4                               (9)
```

still requires its stated large-clock threshold. Equation (8) imposes
no further factor beyond (9). Likewise, all coefficient-first
extension counts

```text
N_i=(1,16,8,4,2,1)                                 (10)
```

remain unchanged after a fixed word has been selected.

The proof is an indicator-boundary proof. Form both whole Boolean
products in (7), use idempotence, and only then apply a common
Poisson--Abel smoothing. Separately smoothed word factors do not square
to themselves at positive Abel radius.

## 2. The unsplit three-word invoice

There is a useful preselection form. Put

```text
W_R=sum_(sigma in Sigma_j)W_sigma
    =1_(R_j)(R x).                                  (11)
```

Suppose the present packet is already supported on the unsplit return
residual:

```text
F_theta W_R=F_theta.                                (12)
```

The complement of `R_j` is then a zero atom. The effective refined
bank is

```text
P_(omega,sigma)=P_omega W_sigma,

|Omega|*|Sigma_j|=128*3=384.                       (13)
```

Repeating the proof of (8) makes every off-diagonal local or word pair
vanish. If `T_R` is the return-supported table and

```text
D_R=||mathcal Q T_R||_2^2,
```

then

```text
T_R=sum_(omega,sigma)T_(omega,sigma).
```

The triangle inequality gives some semantic word atom with

```text
D(T_(omega,sigma))
 >=D_R/384^2
 =D_R/147456.                                      (14)
```

THM-2365's eligible-energy factor `1/13` and its `2016` eligible
colours give

```text
eligible energy
 >=D_R/1916928,                                    (15)

max eligible |B|^2
 >=D_R/3864526848.                                 (16)
```

For a coefficient-first THM-2452 extension of a first-failure cell,
the unsplit word bank has

```text
3N_i=(3,48,24,12,6,3)                              (17)
```

extensions, so one coefficient has modulus at least
`|B_alpha|/(3N_i)`.

Equations (14)--(17) apply only under (12). They do not split a wholly
unfiltered bare table, because the three semantic words partition
`R_j`, not the entire circle. The constants are sharp at the level of
linear conservation: `384` parallel equal projected components, or
`3N_i` equal complex coefficient components, attain the corresponding
triangle inequalities.

## 3. The word is a base coordinate on a root chart

Now add the common oriented physical chart required by THM-2457:

```text
phi(y,u)=(y+u)/13,              u in F_13.           (18)
```

For every `k>=1`,

```text
W_sigma(phi(y,u))
 =1_(Q_(j,sigma))(13^(k-1)y+13^(k-1)u)
 =1_(Q_(j,sigma))(13^(k-1)y).                      (19)
```

The last equality uses one-periodicity. Define

```text
w_sigma(y)=1_(Q_(j,sigma))(13^(k-1)y).              (20)
```

Thus the future word is constant in the root variable `u`. It is an
honest base partition coordinate, not another root mask.

Use THM-2457's disjoint semantic packets `A_0,F_0`, complete root
atoms `U_omega`, and root counts `a_omega,f_nu`. In this section they
are restricted to the return-supported packet:

```text
w_R A_0=A_0,                 w_R F_0=F_0,

w_R=sum_sigma w_sigma.                              (20a)
```

Without (20a), the complement of `R_j` is a fourth base block and the
three-word total-mass statements below do not apply. Refine the
return-supported packets by

```text
a_(omega,sigma)=w_sigma a_omega,

f_(nu,tau)=w_tau f_nu.                              (21)
```

The enhanced directed co-support matrix is

```text
M_((omega,sigma),(nu,tau))
 =integral a_(omega,sigma)(y)f_(nu,tau)(y)dy

 =1_(sigma=tau)
   integral w_sigma(y)a_omega(y)f_nu(y)dy.          (22)
```

Hence the graph is the disjoint union of three `128`-vertex word
blocks. Put

```text
M_sigma=sum_(omega,nu)
 integral w_sigma a_omega f_nu.                    (23)
```

If `M_sigma>0`, some directed edge in that same semantic word obeys

```text
M_((omega,sigma),(nu,sigma))
 >=M_sigma/128^2
 =M_sigma/16384.                                   (24)
```

Adjoining the two endpoints of this edge to any selected nonzero
matched atom in the same word gives THM-2457 service after a union of
at most three enhanced atoms, all carrying the same `sigma`.

If only the total

```text
M=sum_sigma M_sigma>0
```

is known, some block has `M_sigma>=M/3`, and consequently some
word-preserving edge has

```text
M_edge>=M/(3*128^2)=M/49152.                       (25)
```

The sharp THM-2457 energy bounds remain valid blockwise. For a
word-preserving service mass `M_sigma`,

```text
sum_(q!=0)|J_sigma(q)|^2>=M_sigma^2/342732,

max_(q!=0)|J_sigma(q)|>=M_sigma/2028.               (26)
```

The same word carries every nonzero root colour qualitatively.

## 4. The exact quantitative sidecar

Let the selected Boolean matched table be

```text
T_(kappa,sigma):G^3->[0,1],

bar(T_(kappa,sigma))
 =13^(-3)sum_z T_(kappa,sigma)(z).                 (26a)
```

As in THM-2457, orthogonal projection and `0<=T<=1` give

```text
D(T_(kappa,sigma))<=bar(T_(kappa,sigma)).
```

Choose a maximizing entry `z_*` with

```text
T_(kappa,sigma)(z_*)
 >=bar(T_(kappa,sigma)),
```

and let `S=S_(kappa,sigma,z_*)` be its underlying Boolean physical
support. Then

```text
mu(S)=T_(kappa,sigma)(z_*)
     >=bar(T_(kappa,sigma))
     >=D(T_(kappa,sigma)).                         (26b)
```

An arbitrary entry support need not satisfy (26b); the maximizing
selection is load-bearing. Because this table already contains
`W_sigma`, equation (19) automatically gives

```text
pi_13(S) subset {y:w_sigma(y)=1}.                   (27)
```

Let `I` be a set of complete atoms in the same word block and put

```text
Omega_(I,sigma)
 ={y:w_sigma(y)=1,
      sum_u A_I(y,u)>0,
      sum_u F_I(y,u)>0}.                            (28)

M_sigma(I)
 =integral_Y w_sigma(y)
    (sum_u A_I(y,u))(sum_u F_I(y,u))dy.
```

The word-preserving version of THM-2457's root-image hypothesis is

```text
pi_13(S) subset Omega_(I,sigma)                     (29)
```

up to null sets. Under (29), the same measure argument gives

```text
M_sigma(I)>=mu(pi_13(S))>=mu(S)
              >=D(T_(kappa,sigma)).                (30)
```

This is the exact sufficient sidecar which feeds selected drift into
same-word root service. The word part of (29) is free by (27); the
positive simultaneous semantic descendants and common orientation
are not.

## 5. Clock and target/source gauge typing

The operation preserves a clock; it does not identify clocks.

1. On a genuine THM-2305 prescribed-return arm,

   ```text
   k=k_j=lambda_j+1,
   ```

   and (8) preserves that exact first prescribed expiration.

2. On the all-row THM-2349/2365 route, `k` is a later row-dependent
   BV restart clock. Equations (8) and (19) preserve that delayed
   semantic handoff, but it must not be renamed the original
   prescribed expiration.

3. THM-2452's eventual delayed-word argument may increase `k` again.
   It preserves the chosen owner and word label at the resulting
   clock, not first-return minimality.

The target action is lawful because every transported word factor has
the form

```text
J_i(Rw_i x-R(s eta_i+t ell_i)/13).
```

Since `13|R`, the second term is an integer. Every `W_sigma` is
therefore target-neutral separately, not merely after summing the
three words.

The septimal source transform is different. Put

```text
epsilon=R mod 7=(-1)^k.
```

For a source phase `ell`, THM-2442 uses the lawful shifted word

```text
Q^(epsilon)_(sigma,ell)(y)
 =T_sigma(y)g(c_jy-epsilon ell/7),                  (31)

W^(epsilon)_(sigma,ell)(x)
 =Q^(epsilon)_(sigma,ell)(Rx).                     (32)
```

Apply (8) separately for every `ell`, copying the same
`W^(epsilon)_(sigma,ell)` to both endpoint filters, and only then take
the finite `C_7` transform. The equality is phasewise, so a previously
selected nonzero source colour is unchanged and there is no additional
partition loss.

Holding `W_sigma` fixed while translating the present source is not
lawful. It omits the `R ell` word phase and is exactly the shortcut
rejected by THM-2409 and repaired by THM-2442.

## 6. Relation-current typing and the fake second harmonic

The copied word in (7) is a boundary filter. Collapse

```text
W_sigma^2=W_sigma                                  (33)
```

before Fourier expansion. The original THM-2334 address is then still

```text
r=u+R beta+m e_c-v.                                (34)
```

The word harmonic `R beta` vanishes in the mod-`13` target quotient,
but remains active modulo seven. With a lawful source transform, its
septimal phase is retained exactly as in THM-2442.

If the two copies in (7) are expanded as independent physical word
factors, one introduces two harmonics `beta,gamma` and apparent
addresses containing `R(beta-gamma)`. Those allocations are not
canonical. They recombine only through the Fourier form of
idempotence,

```text
Q_hat(n)=sum_beta Q_hat(beta)Q_hat(n-beta).         (35)
```

For example, on `C_4` the singleton word `Q=1_{0}` has all four
normalized Fourier coefficients equal to `1/4`. Every output in
(35) has four nonzero summands `1/16`, although only their sum `1/4`
is intrinsic.

Accordingly, after the whole-table slide one expands the collapsed
original representation and uses THM-2365 to select a fresh exact
`X,m`. The semantic word, finite target/deep character, graft, and
lawful source colour remain. The deepest Fourier factor again makes
`m` a `91`-unit. No old exact `X`, old word harmonic `beta`, or
separate right-word address is preserved.

## 7. Two sharp stopping hostiles

### 7.1 Selected drift and semantic co-support can occupy different words

Let the rational base be the disjoint union of three intervals

```text
Y_0=[0,1/3), Y_1=[1/3,2/3), Y_2=[2/3,1),
```

carrying the three word labels. Put a fixed nonzero selected matched
table entirely on `Y_0`. On the common root chart, let the semantic
packet `A_0` occupy root `0` on `Y_0 union Y_1`, and let `F_0` occupy
root `1` only on `Y_1`.

The packets are rootwise disjoint. Their global simultaneous
same-parent mass is `1/3`, but

```text
M_(word 0)=0,

M_(word 1)=1/3.                                    (36)
```

Thus all selected table mass lies in word zero and all co-support lies
in word one. A THM-2457 edge exists globally, but adjoining it to the
selected table destroys the unique word. This rational common-chart
control proves that neither endpoint idempotence nor separate
co-support implies (29).

The control is abstract rather than a realized nine-comb scalar row.
It is sharp for the logical interface: any stronger theorem needs a
physical coupling hypothesis.

### 7.2 A word label does not remove the uniform-offset locus

Equation (19) makes the word constant on each root chart. Therefore
any of THM-2456's `8`-, `20`-, or `104`-chart uniform-offset controls
may be tensored with one fixed word label without changing a single
owner/replica density or convolution equation.

The actual canonical word can still forbid those mixtures by
restricting which base charts are physically realizable. That is a
new finite realizability theorem, not a consequence of Boolean
idempotence. The correct test is THM-2456's occupancy system refined
by the fixed data

```text
(owner j, clock k, word sigma).                     (37)
```

THM-2457's live strict-profile hostile supplies the complementary
boundary: before this refinement, one complete local atom genuinely
meets both a pure and a fork word on positive intervals. The extra
coordinate separates those intervals; the local bits alone do not.

## 8. Exact gain and remaining frontier

The theorem removes one possible ambiguity:

```text
an independently retained semantic THM-2305 word
can be carried through THM-2452 endpoint restoration
without another local-mask loss.                   (38)
```

It does not infer a word from the local atom. Nor does it prove:

- a common oriented physical root chart for the selected table;
- positive `A/F` descendant co-support in that same word block;
- the root-image inclusion (29);
- exclusion of the word-conditioned THM-2456 uniform-offset atlas;
- preservation of a prescribed first-expiration clock after a later
  BV restart;
- an old exact relation address.

The sharp next alternatives are therefore:

```text
prove same-word root-image coupling (29),

or

run the physical uniform-offset realizability atlas inside each
fixed (j,k,sigma) block and prove f=0.               (39)
```

The semantic word itself is no longer the missing Boolean coordinate.
The remaining obstruction is positive physical common-root service in
that word, together with the canonical relation-current interpretation.
No scalar row is excluded and LRC(14) remains open.

## 9. Exact companion

Run

```text
python 04-computation/lrc14_semantic_word_copy_thm2460.py
python -O 04-computation/lrc14_semantic_word_copy_thm2460.py
```

The dependency-free companion uses only integers and `Fraction`, with
explicit `require` checks active under optimized mode. It verifies:

- all fixed-word and unsplit-word diagonal censuses;
- the `384`, `147456`, `1916928`, and `3864526848`
  invoices and all adaptive counts;
- root constancy at six positive clocks;
- the three-block graph census and sharp `1/16384`, `1/49152`
  edge floors;
- target neutrality and septimal activity;
- the noncanonical four-allocation `C_4` word-harmonic control;
- the exact rational split-base hostile; and
- invariance of an averaged uniform-offset control under a constant
  word tensor.

Normal and optimized transcripts match the stored output
byte-for-byte after LF normalization. LF hashes are recorded in the
frontmatter.

An independent hostile audit rederived the fixed-word diagonal,
return-supported `384`-atom invoice, root-constant block
decomposition, clock and source-skew typing, relation-address
collapse, and both stopping hostiles. It caught and repaired the
missing maximizing-entry quantifier in Section 4 before promotion,
then independently replayed normal and optimized modes against the
stored transcript and reproduced both LF hashes.
