---
source: codex-2026-07-25-zakharov-overlap-transfer-audit
status: >
  RESEARCH NOTE. The b-adic cylinder lemma below is proved, but its LRC
  application is quantitatively dominated by THM-2291. The itinerary and
  tooth-word paragraphs are stopping audits, not an LRC(14) reduction.
external:
  - https://arxiv.org/abs/2602.20143
---

# Zakharov overlap sees cylinders, not owner letters

Zakharov's sharp bipartite word-overlap theorem says that if
`A,B subset Omega^n` have no directed suffix--prefix overlap, then

```text
mu(A) mu(B)
 <= C_n
 :=(1/n)(n/(n+1))^(n+1).
```

The proof's first-hit decomposition and positive generating-function root
look close to the LRC expiration tree. The resemblance is real only after
the LRC sets are encoded by **whole dynamical cylinders**. An owner letter,
a blocker subset, or a binary danger itinerary does not retain enough phase
to make a suffix--prefix equality into a physical orbit collision.

Primary source: Dmitrii Zakharov,
[*An isoperimetric inequality for word overlap*,
arXiv:2602.20143v2](https://arxiv.org/abs/2602.20143).

## 1. The exact source-to-target map

Fix `b>=2`, put `Omega={0,...,b-1}`, and let

```text
T_b(x)=b x mod 1.
```

The word `w=(w_1,...,w_n)` names the half-open `b`-adic cylinder

```text
I_w=[m(w)/b^n,(m(w)+1)/b^n).
```

Let `E,R subset R/Z` be finite unions of intervals, considered modulo null
endpoints, and suppose their boundary sets have at most `J_E,J_R` points.
Define

```text
A_n={w:I_w subset E modulo null sets},
B_n={w:I_w subset R modulo null sets}.
```

Then

```text
mu(A_n)>=mu(E)-J_E/b^n,
mu(B_n)>=mu(R)-J_R/b^n.                    (1)
```

Indeed, every `b`-adic cell meeting both a set and its complement contains
a boundary point, so at most `J_E` source cells and `J_R` target cells are
lost.

Now suppose a suffix of `a in A_n` of length `j` equals a prefix of
`r in B_n` of length `j`. Concatenating the two words along that overlap
produces a word of length `2n-j`. Its cylinder lies, modulo endpoints, in

```text
E intersection T_b^(-(n-j))R.                       (2)
```

Thus the following is an exact contrapositive of Zakharov's theorem.

> **B-adic cylinder return lemma.** For every integer `n>=1`, if
>
> ```text
> max(0,mu(E)-J_E/b^n)
> *max(0,mu(R)-J_R/b^n)>C_n,                         (3)
> ```
>
> then there is an integer `k` with `0<=k<n` such that
>
> ```text
> E intersection T_b^(-k)R
> ```
>
> contains a positive-measure `b`-adic cylinder modulo endpoints.

The quantifiers and direction matter. The conclusion is one directed return
at one lag shorter than `n`. The hypothesis uses inner cylinders, not words
obtained by merely recording a factor of the state.

The connection ledger is therefore:

```text
source:
  two uniform full-shift word families;

map:
  word -> whole b-adic cylinder, suffix length j -> lag n-j;

preserved:
  direction, exact lag, and positive physical cylinder inclusion;

lost by an itinerary quotient:
  terminal phase, tail position inside the cylinder, marked root sheet,
  current blocker service, and sibling ancestry;

needed sidecar:
  whole-cylinder containment with a boundary count, or an endpoint-state
  automaton carrying the same information.
```

## 2. Why the repeated-first application is dominated

Apply the lemma with `b=13`,

```text
E=E_j
```

the exclusive-owner source, and

```text
R=R_j
```

the blocker-only target from
[THM-2291, repeated-owner BV mixing and delayed blocker handoff](../01-canon/theorems/THM-2291-repeated-owner-bv-mixing-and-delayed-blocker-handoff.md).
On every repeated-first row that theorem supplies

```text
mu(E)>=e_0:=5229541/593783190,
mu(R)>=eta_0:=2593/90090,
J_E,J_R<=2S.                                        (4)
```

The sufficient boundary-resolution condition

```text
13^n>=4S/e_0
      =(2375132760/5229541)S                         (5)
```

makes both inner-cylinder masses at least half their displayed floors.
Zakharov then gives a return whenever additionally

```text
C_n<e_0 eta_0/4.
```

The first such integer is `n=5805`; without the conservative half-mass loss,
the raw product `e_0 eta_0` first exceeds `C_n` at `n=1451`.

This does not improve the frontier. Condition (5) is exactly the BV scale in
THM-2291, while that theorem has no additional `n>=5805` tax. More strongly,
THM-2291 proves a positive return for **every** integer `k` after that scale
and gives the explicit mass floor

```text
13560199813/106987855174200.
```

The cylinder lemma gives only existence of one earlier-than-`n` return.
Retaining Zakharov's finer projection-drop generating function can reduce the
word-length tax because expiration expands the source image, but it cannot
remove the need to resolve the endpoint boundary bank. The existing BV
argument already uses that boundary information more efficiently.

## 3. The exact itinerary hostile

The sharp reason not to apply the theorem to binary owner itineraries is
already canonical. In
[THM-2269, marked expiration root spectrum and branch-state no-go](../01-canon/theorems/THM-2269-marked-expiration-root-spectrum-and-branch-state-no-go.md),
the exact branch

```text
x_0=79/338,
x_1=1/26,
x_t=1/2 for every t>=2
```

satisfies every displayed forward pointwise cover clause, while its
post-expiration normalized-core itinerary is

```text
00000...
```

That word has a suffix--prefix overlap of every length, yet the selected
owner supplies no future service on the branch. This does not contradict
Zakharov: the binary itinerary is only a factor of the root dynamics.
Equality of factor words has forgotten the terminal phase and sibling-root
gluing needed to produce (2).

The earlier
[THM-2261 expiration-image no-go](../01-canon/theorems/THM-2261-expiration-image-surjectivity-and-one-core-carrier-no-go.md)
already shows that raw guard-plus-owner data can have arbitrary expiration
image. THM-2269 strengthens the warning by retaining every forward cover
clause and still realizing the maximally self-overlapping zero word.

## 4. The two other meanings of “word”

The terminal blocker word in
[THM-2305](../01-canon/theorems/THM-2305-canonical-blocker-word-handoff-hypergraph.md)
is one of

```text
{a}, {b}, {a,b}.
```

It is an exact subset label, not a word under concatenation. The factor word
in
[THM-2337](../01-canon/theorems/THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar.md)
retains Fourier-factor indices, but its first jet is gauge-dependent after
the base/word split is forgotten. Neither object has Zakharov's
suffix--prefix operation without adding a new, proved composition law.

The chronological tooth word in
[THM-1253](../01-canon/theorems/THM-1253-irredundant-tooth-word-pays-every-handoff.md)
and
[THM-1256](../01-canon/theorems/THM-1256-binary-phase-word-landing-and-toothpick-run-break.md)
is a genuine ordered word, but a fixed LRC row supplies one positioned
interval chain rather than a positive-density family in `Omega^n`.
Its proof-bearing coordinates are tooth addresses, endpoint order, overlap
cells, and lcm weights. Projecting to owner letters can create many symbolic
prefix--suffix matches while destroying all four coordinates. Conversely,
THM-1256's `ABAB` exclusion is only a local metric restriction, not global
cross-bifix-freeness.

## 5. Stopping rule

Use Zakharov on LRC only after checking one of the following:

1. every word names a whole dynamical cylinder, with boundary loss charged;
2. the alphabet carries endpoint state, owner/service labels, and sibling
   ancestry so that a suffix--prefix match really concatenates orbit pieces.

Absent one of these sidecars, word overlap is telemetry. The all-zero hostile
shows that even maximal overlap need not create the missing owner collision.
With whole cylinders the transfer is rigorous, but THM-2291 already proves a
strictly stronger repeated-first return theorem. Consequently this connection
changes no frontier status and excludes no LRC(14) scalar row.
