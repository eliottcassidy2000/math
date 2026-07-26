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
  quadrature already recover its complex phase. One real union alone
  is sharply insufficient. The result is coefficient-theoretic unless
  an endpoint-matched physical pair-twist service is supplied. No
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
script_sha256: da0e3b91be459fe65606bf5e728b25cc949741b90148590523b0f943b2986c23
output_sha256: 7a66fde3abbc9dfb7920185c7221a22a43d8ff77a9f47c6e5fad7c55ebc085f9
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
