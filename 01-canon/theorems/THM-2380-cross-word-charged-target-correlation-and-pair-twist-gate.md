---
id: THM-2380
title: "Cross-word charged target correlation and pair-twist gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For two
  complex target currents on F_13^2, the complete Cayley difference
  energy of their cross-correlation is positive exactly when their
  target supports overlap away from zero. The cross-correlation
  inverts to the pointwise charged product d(q)conj(w(q));
  coordinate-line energies resolve the two target directions.
  Singleton energies and a full thirteen-phase pair-twist recover the
  product exactly, while one ordinary union and one nontrivial
  quadrature already recover its complex phase. On an affine target
  line q_1=b!=0, two t-slices of the ordinary-translation bank (for
  all s) provide that quadrature and recover the charged product
  without an artificial phase twist. The self-opposite line b=0 retains exactly the I+J
  ambiguity. The result is coefficient-theoretic unless the required
  endpoint-matched translations or pair twist are physically supplied. No
  scalar-row exclusion, ledger decrement, or LRC(14) consequence is
  asserted.
source: codex-2026-07-25-cross-word-charged-correlation
depends_on: []
related:
  - THM-2355-component-deletion-gram-and-twist-energy-phase-transport
  - THM-2356-finite-field-chirp-gram-tomography-and-bockstein-pairing
  - THM-2363-planar-graph-detector-dominates-word-support-energy
  - THM-2369-complete-line-target-dirichlet-and-balanced-observable-no-go
  - THM-2370-deletion-martingale-drift-conservation-and-sharp-clone-hostile
  - THM-2374-binary-allocation-complete-subcube-dirichlet-spectrum
  - THM-2375-gaussian-angular-complete-line-charge-tomography
script: 04-computation/lrc14_cross_word_charged_correlation_thm2380.py
output: 05-knowledge/results/lrc14_cross_word_charged_correlation_thm2380.out
script_sha256: 95ce85f9a4c9f71bb7ad82704ec375256c72d22129c0ebe1eaccd06b8bc89483
output_sha256: 74c4088a82d0adcb726246cfd70c358088bfcb0a91e069cd45ea1a8f6c8da532
hash_basis: working-tree bytes (LF)
---

# THM-2380 -- cross-word charged target correlation

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The target spectra of a deletion layer and a retained word need not meet:
THM-2370's deletion drift can live at a nonzero target which is dark for
the word. The missing invariant is not either spectrum separately, but
their charged pointwise product. A complete translated
cross-correlation recovers that product exactly:

```text
deletion current d(q)
  + retained-word current w(q)
  + endpoint-matched cross-correlation
  -> d(q)conj(w(q))
  -> exact nonzero-target overlap test.                         (1)
```

The algebraic reconstruction is unconditional. In the LRC application,
the endpoint-matched cross-word observable in the third line of (1)
remains an extra service.

## 1. Target Fourier conventions

Let

```text
G=F_13^2,                    N=|G|=169,

zeta=exp(2*pi*i/13),

chi_ell(q)=zeta^(ell_1 q_1+ell_2 q_2).                        (2)
```

Thus the two coordinate translation lines are dual to the displayed
coordinates `q_1,q_2`. An invariant formulation replaces `q_j!=0`
below by `q notin L_j^perp`.
For complex target currents `d,w:G->C`, write

```text
D(ell)=sum_(q in G) chi_ell(q)d(q),

W(ell)=sum_(q in G) chi_ell(q)w(q).                            (3)
```

Thus

```text
d(q)=1/N sum_(ell in G) conj(chi_ell(q))D(ell),                (4)
```

and similarly for `w`. Norms below use normalized counting measure:

```text
||D||_2^2=1/N sum_ell |D(ell)|^2
          =sum_q |d(q)|^2.                                    (5)
```

In an LRC use, (3) is imposed only after putting the deletion and word
currents in one common deepest gauge. Currents in different endpoint
gauges cannot be substituted into the theorem.

## 2. Cross-correlation is the charged product

Define the translated cross-word correlation

```text
R(delta)
 =1/N sum_(ell in G)D(ell+delta)conj(W(ell)).                  (6)
```

Expanding (3) and using character orthogonality gives

```text
R(delta)
 =sum_(q in G)chi_delta(q)d(q)conj(w(q)).                      (7)
```

Consequently,

```text
d(q)conj(w(q))
 =1/N sum_(delta in G)conj(chi_delta(q))R(delta).              (8)
```

There is no phase loss in (8): the entire complex same-target Gram
product is recovered, target by target.

## 3. Complete Cayley energy: an exact iff

Let the inner sum below range over all `e in G\{0}`. Parseval applied
to (7) yields

```text
E_Cayley(R)
 :=1/(N*338)
     sum_(delta in G) sum_(e!=0)
       |R(delta+e)-R(delta)|^2

 =sum_(q!=0)|d(q)|^2|w(q)|^2.                                 (9)
```

Here `338=2N`; hence the denominator is `2N^2`. It follows that

```text
E_Cayley(R)>0

iff

there exists q!=0 with d(q)!=0 and w(q)!=0.                   (10)
```

Thus complete Cayley energy is positive exactly when the deletion and
word target supports overlap at a nonzero target. In particular,
nonzero target energy in each current separately is not enough.

For a coordinate direction `j`, let `L_j<=G` be the corresponding
one-dimensional translation line. The directional version is

```text
1/(13N)
 sum_delta sum_(e in L_j\{0})
   |R(delta+e)-R(delta)|^2

 =2 sum_(q:q_j!=0)|d(q)|^2|w(q)|^2.                           (11)
```

Together with (9), the two coordinate-line energies separate pure
target overlap from the fork overlap where both target coordinates are
nonzero. Indeed, let `A,B,C` be the sums of
`|d(q)|^2|w(q)|^2` over, respectively,

```text
q_1!=0,q_2=0;       q_1=0,q_2!=0;       q_1 q_2!=0.
```

If `L_1,L_2` denote the two left sides of (11) and
`E=E_Cayley(R)`, then

```text
L_1=2(A+C),       L_2=2(B+C),       E=A+B+C,

C=L_1/2+L_2/2-E,  A=E-L_2/2,       B=E-L_1/2.     (12)
```

Equation (11) is a charged cross-word invariant; it is not the
uncharged one-current line energy of THM-2369.

## 4. Pair-twist reconstruction

Suppose one can form the endpoint-matched phase-twist intensities

```text
E_delta(t)
 =1/N sum_(ell in G)
    |D(ell+delta)+zeta^t W(ell)|^2,

                         delta in G, t in F_13.                (13)
```

Use the normalized `t`-transform

```text
Ehat_delta(k)
 =1/13 sum_(t in F_13)E_delta(t)zeta^(-kt).                   (14)
```

Expanding the square in (13) gives

```text
E_delta(t)
 =||D||_2^2+||W||_2^2
   +zeta^(-t)R(delta)+zeta^t conj(R(delta)),                   (15)

Ehat_delta(-1)=R(delta),

Ehat_delta(+1)=conj(R(delta)).                                (16)
```

Combining (8) and (13)--(16) gives the literal intensity inversion

```text
d(q)conj(w(q))

 =1/(13N^2)
   sum_(delta,t,ell)
     conj(chi_delta(q)) zeta^t
     |D(ell+delta)+zeta^t W(ell)|^2.                          (17)
```

At `N=169`, the denominator in the convention (2)--(5), (13)--(14) is

```text
13N^2=371293.                                                  (18)
```

The integer `371293` is not invariant under a change from normalized
to unnormalized transforms; it is asserted only with the displayed
normalization.

## 5. One quadrature suffices; one union does not

The full thirteen-phase bank in (13) is redundant once the two
singleton norms in (15) are known. Put

```text
U_delta(t)
 =E_delta(t)-||D||_2^2-||W||_2^2,

c=cos(2*pi/13),       s=sin(2*pi/13).                          (19)
```

Writing `R(delta)=x+iy`, (15) gives

```text
x=U_delta(0)/2,

y=(U_delta(1)-c U_delta(0))/(2s).                             (20)
```

Therefore one ordinary real union (`t=0`) together with one genuine
quadrature (`t=1`) recovers the full complex correlation.

One ordinary union alone is sharply insufficient. At a single
component, the two pairs

```text
(D,W)=(1,i),            (D,W)=(1,-i)                          (21)
```

have identical singleton energies and identical `t=0` union energy,
but their Gram phases are opposite. This is the same missing phase
that prevents squared Walsh or scalar angular spectra from identifying
a charged edge.

## 5a. Ordinary affine translations supply quadrature off the zero line

There is a sharper repair when the charged product is supported on a
known affine target line. Write target coordinates as

```text
q=(b,r),                    e_1=(1,0), e_2=(0,1),

z(b,r)=d(b,r)conj(w(b,r)).
```

Use only ordinary relative translations, with no scalar phase twist:

```text
U(s,t)
 =1/N sum_(ell in G)|D(ell+s e_2)+W(ell+t e_1)|^2
  -||D||_2^2-||W||_2^2.                              (21a)
```

The cross term and (7) give

```text
U(s,t)
 =R(s e_2-t e_1)+conj(R(s e_2-t e_1)).               (21b)
```

Taking the normalized two-variable transform yields the exact
symmetrized projector

```text
1/N sum_(s,t in F_13)
 U(s,t) zeta^(-sr+tb)

 =z(b,r)+conj(z(-b,-r)).                              (21c)
```

Thus real translated-union intensities lose only reflection across the
origin in target space. Suppose now that

```text
supp(w) subset {(b_0,r):r in F_13},                  (21d)
```

where `b_0!=0`. Since `13` is odd, the reflected line `-b_0` is
different from `b_0`, so (21c) recovers `z(b_0,r)` exactly.

In fact the complete `t`-bank is unnecessary. Put

```text
A_s=sum_r z(b_0,r)zeta^(sr).
```

The two ordinary relative shifts `t=0` and `t=b_0^(-1)` give

```text
U(s,0)=A_s+conj(A_s),

U(s,b_0^(-1))=zeta^(-1)A_s+zeta conj(A_s),           (21e)
```

and hence

```text
(zeta^(-1)-zeta)A_s
 =U(s,b_0^(-1))-zeta U(s,0).                         (21f)
```

The coefficient on the left is nonzero, and inverse Fourier transform
in `s` recovers every `z(b_0,r)`. Ordinary translation has supplied the
second quadrature because the nonzero affine charge rotates the common
line by a known phase.

For arbitrary complex coefficient currents, the hypothesis `b_0!=0`
is sharp. On the zero line, (21c) is the operator

```text
I+J,                 (Jz)(0,r)=conj(z(0,-r)).        (21g)
```

It has a nontrivial kernel. For example, take

```text
w(0,r_0)=w(0,-r_0)=1,

d(0,r_0)=d(0,-r_0)=i,                 r_0!=0.         (21h)
```

Then the charged energy is positive, but every `U(s,t)` in (21a)
vanishes. This is a coefficient-level hostile and deliberately violates
Hermitian symmetry. If both currents come from real tables, then

```text
d(-q)=conj(d(q)),              w(-q)=conj(w(q)),
```

so `z(-q)=conj(z(q))`, hence `Jz=z` and (21c) returns `2z` even on
`b=0`. In that Hermitian setting, support on one nonzero line `b_0`
also forces support on its reflected line `-b_0`; the single-line
hypothesis (21d) is therefore a complex coefficient grouping, not the
support of a full real table.

The promoted first-collision null-parent control lies on the same
`b=0` line but in the opposite eigenspace: on its twelve nonzero target
points one has `d=-w`, and the word is real-even, so (21c) returns
`-2|w|^2` there rather than zero. Its residual is therefore not
quadrature but canonical owner typing and restoration of the common
terminal endpoint.

## 6. Sharp disjoint-support hostile

Fix `q_0!=0` and take

```text
d(q)=1_(q=q_0),          w(q)=1_(q=0).                        (22)
```

This models a circulant retained word. To make both currents carry
separate nonzero targets, one may instead choose distinct
`q_0,q_1!=0` and put

```text
d(q)=1_(q=q_0),          w(q)=1_(q=q_1).                       (23)
```

In either example the deletion current has a live nonzero target and
the word current is nonzero, but

```text
d(q)conj(w(q))=0         for every q,

R(delta)=0               for every delta,

E_Cayley(R)=0.                                                (24)
```

Thus neither a deletion drift floor, nor a word-energy floor, nor the
existence of targets on both sides forces the overlap in (10).
The currents must be aligned at the same target.

For the `n=1` Boolean clone of THM-2370, the retained endpoint may be
circulant while the deletion layer drifts. Equation (22) is the finite
target model of that boundary: complete singleton and union data still
do not select a common charged target.

## 7. LRC typing and the remaining physical service

Equations (6)--(20) are lawful finite coefficient groupings once the
joint target coefficients are known. They do **not** by themselves
construct the physical observable (13).

At `t=0`, a same-endpoint sum can sometimes be interpreted as the
Fourier transform of a real union after disjointness or
inclusion--exclusion bookkeeping. For `t=1`, multiplication of one
current by `zeta` is not an indicator operation. Moreover, the current
word and the THM-2370 deletion layer can carry different right
endpoints. The cheapest missing sidecar is therefore:

```text
one common endpoint gauge
  + one lawful cross-word phase quadrature
  -> (20)
  -> R(delta)
  -> the exact overlap test (10).                              (25)
```

A coefficient-derived quadratic functional is enough for the algebra,
but calling it a physical LRC measurement requires an independent
realization. No such realization is asserted here.

Section 5a lowers that service on a known nonzero affine target line:
one common endpoint gauge and the two lawful `t`-slices in (21e), each
for all `s`, suffice. It does not construct that translation bank inside
the terminal owner word, and it does not repair the self-opposite `b=0`
line without an additional symmetry or reference.

## 8. Relation to nearby tomography

- THM-2355 reconstructs phases between labelled components from
  singleton and pair-twist energies. The present theorem groups two
  whole currents and identifies the exact targetwise product that such
  a twist must recover.
- THM-2356 and THM-2363 expose graph-channel energy but retain a
  singleton/phase-ratio boundary. Equation (8) names the missing
  cross-word ratio directly.
- THM-2369 proves that balanced one-current observables can miss target
  slope. The charged correlation (6) is deliberately not balanced:
  it retains the complex phase between two currents.
- THM-2374's squared Boolean Dirichlet spectrum and THM-2375's
  isotypic norms are Hermitian magnitude data. The pairs in (21) show
  why neither supplies (8).
- The independently audited finite-exact
  [`first-collision null-parent positive control`](../../07-reflections/lrc14-first-collision-null-parent-polarization-positive-control-codex-20260725.md)
  gives a complementary boundary. On THM-2377's reduced common-endpoint
  slice the parent is target-null, so `Ahat(q)=-What(q)` for every
  nonzero target and the real Gram product is already `-|What(q)|^2`;
  all twelve target and `144` deep/target colours survive. Thus
  disjoint support is not the immediate reduced-slice obstruction.
  Canonical owner typing and common-endpoint restoration to the
  terminal word remain open, so this is a positive control rather than
  a physical realization of (13).

No scalar row is excluded. The ledger remains `165`, and LRC(14)
remains open.

## 9. Exact companion and hostile audit

Run

```bash
python 04-computation/lrc14_cross_word_charged_correlation_thm2380.py
python -O 04-computation/lrc14_cross_word_charged_correlation_thm2380.py
```

and compare both transcripts byte-for-byte, after LF normalization,
with
`05-knowledge/results/lrc14_cross_word_charged_correlation_thm2380.out`.
The dependency-free exact companion checks:

- all `169` cross-correlation and inverse-target cells;
- `2,197` pair-twist orbit cells and six independent frequency-side
  Parseval controls;
- all `169` raw target projectors with denominator `371293`;
- the complete-Cayley denominator `57,122`, three direct translation
  controls, both coordinate-line energies, and the `(A,B,C)` recovery
  in (12);
- the exact quadratic-extension reconstruction (20);
- four raw frequency-side affine-union controls, all `169` symmetrized
  projectors in (21c), the two-shift determinant (21f), and both
  eigenspaces of the zero-line operator (21g);
- the ordinary-union phase ambiguity (21); and
- both disjoint-support hostiles (22)--(24), including the example in
  which both currents have nonzero-target energy.

Ordinary and optimized runs are byte-identical to the stored transcript;
the LF hashes are recorded in the frontmatter. An independent hostile
audit separately reconstructed every Fourier sign and normalization,
the factors `338`, `2`, and `371293`, the two-sample determinant, and
all three hostiles. It caught and repaired the need to fix the
coordinate pairing, use (9) together with the two line energies, and
include a hostile with both targets nonzero. QED.
