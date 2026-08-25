---
id: THM-4064
title: "Rule 30 cyclotomic kernel character and C60 alias obstruction"
status: >
  PROVED universal cyclotomic fibre/Hasse-invariant theorem + PROVED
  characteristic-zero primitive-frequency theorem + FINITE-EXACT fixed
  Rule-30 left-front census + VERIFIED-EXACT + INDEPENDENTLY
  HOSTILE-AUDITED. Exponentiating the THM-4055 combined clock identifies
  mu_N with mu_60 x_(mu_g) mu_p. The C_60 projection has kernel mu_m, and
  the mu_m-invariants of K[P]/((P^N-1)^k) are exactly the pulled-back
  C_60 Hasse algebra. Every rational least-dyadic-period response has every
  odd Fourier frequency, so for p>4 it necessarily uses nontrivial kernel
  characters. In the THM-4047 fixed bank, labelled exact denominator data
  are lossless, while denominator, Stern depth/sign, and flank parity still
  alias ell_29(90) and ell_29(150). Moving-center behavior, temporal
  balance, every Rule-30 prize, and integral Smith transfer remain OPEN.
source: codex-frontier-synthesis-breakthrough-20260825 / Rule-30 observer lane
audit: >
  PASS after separating raw-bit and sign Fourier support, treating d=1 by
  an explicit None sentinel rather than a fake Stern depth, restricting the
  primitive-frequency lemma to characteristic zero, and weakening the
  consecutive/cyclotomic comparison to an integral Smith firewall. The
  primary recurrence/physical-front path performs 30,744,788 checks. The
  no-import audit performs 5,573,429 checks, conditional only on the seven
  frozen physical-query bits already audited by THM-4047. Both normal and
  optimized runs pass, with zero assert nodes and zero float literals; their
  independent factor and sign-sector semantic hashes agree.
depends_on:
  - THM-4010-confluent-consecutive-hasse-observer-kernel-index-and-smith-firewall
  - THM-4044-sixty-clock-hasse-alias-and-planar-jc-boundary-firewall
  - THM-4047-rule30-left-front-affine-monodromy-clock
  - THM-4055-sixty-dyadic-response-fibre-law
  - THM-4056-divisor-phase-compiler-and-duffin-schaeffer-lcm-clock
  - THM-4059-stern-brocot-depth-packet-character-and-divisor-star-convolution
related:
  - THM-4035-sixty-clock-separation-and-finite-kakeya-spine
  - THM-4048-rule30-periodicity-balance-and-model-firewalls
  - THM-4050-rule30-half-arc-marked-cylinder-and-radius-nine-hostile
script: 04-computation/rule30_cyclotomic_hasse_packet_thm4064.py
output: 05-knowledge/results/rule30_cyclotomic_hasse_packet_thm4064.out
independent_audit_script: 04-computation/rule30_cyclotomic_hasse_packet_thm4064_independent_audit.py
independent_audit_output: 05-knowledge/results/rule30_cyclotomic_hasse_packet_thm4064_independent_audit.out
script_sha256: 46f35b1081ab13d8a70bbcae35180eb7a00758fefaebfbee7b7cd83cb2f0d4a3
output_sha256: e753b1bfdf082d9f4871f116745e7f3bcd9aea4f097c5a523d51c746388c4dce
independent_audit_script_sha256: 11f559d7d41e3738b2942f0a5464a8eebfe5589f2f1d3e71fa211cbf9210ed03
independent_audit_output_sha256: 2baaeb90311e8b444979279e9e4006497b785850dac3c7627a33c074f3ab5e6d
factor_semantic_sha256: fc7679bc50f53005384a9f0c34f33ff651cf30b578a2828795495525adca427f
sign_sector_semantic_sha256: 41a87e724c878456bc94ae8deb5f65892562a1d689e280dc2ab53d77b87815d4
hash_basis: raw LF bytes
---

# THM-4064 -- the lost Rule-30 phase is a cyclotomic character

**PROVED in the universal scopes below; FINITE-EXACT only for the fixed
THM-4047 left-front bank.** THM-4055 found the missing cyclic torsor over a
sixty-phase address. Exponentiation makes that torsor the kernel of a map of
root groups. THM-4044 then turns the same kernel into an exact decomposition
of every confluent clock algebra. For rational dyadic responses, least period
forces primitive characters in the lost part.

## 1. Exponentiating the combined clock

Let

```text
p=2^e,       g=gcd(60,p),       m=p/g,
N=lcm(60,p)=60m.                                      (1)
```

Let `K` be a characteristic-zero field containing `mu_N`, and choose a
primitive `N`-th root `zeta_N`. The map

```text
C_N -> mu_N,          [t] |-> zeta_N^t                 (2)
```

is a group isomorphism. Moreover, if the representative `t/N` reduces to
`a/d`, with `(d,a)=(1,0)` at `t=0`, then

```text
order(zeta_N^t)=d=N/gcd(t,N).                          (3)
```

Thus THM-4056's labelled exact-denominator packet becomes literally the
point and order of a root of unity; no metaphorical transfer is needed.

Use the power maps from `mu_60` and `mu_p` to `mu_g`. Then

```text
Psi:mu_N -> mu_60 x_(mu_g) mu_p,
Psi(z)=(z^m,z^(N/p))                                  (4)
```

is an isomorphism, where

```text
mu_60 x_(mu_g) mu_p
 ={(x,y):x^(60/g)=y^(p/g)}.                           (5)
```

Indeed, both powers in `(5)` equal `z^(N/g)`. If both coordinates in `(4)`
are one, the order of `z` divides both `m=p/g` and `N/p=60/g`; these two
integers are coprime, so `z=1`. The fibre product has
`60p/g=N` points, proving surjectivity.

The first projection is

```text
mu_N -> mu_60,       z |-> z^m,       kernel=mu_m.    (6)
```

Every fixed `mu_60` fibre is therefore a `mu_m`-torsor. After choosing
primitive roots and a lift of the sixty-phase origin, this is exactly the
THM-4055 sidecar `q in C_m`; changing the lift translates `q`, so the torsor
still has no canonical origin.

## 2. The complete Hasse algebra splits by the same kernel

For `k>=1`, put

```text
A_(N,k)=K[P]/((P^N-1)^k).                             (7)
```

By THM-4044, this is the complete remainder algebra seen by depth-`k` Hasse
evaluation on `mu_N`. The group `mu_m` acts by `P |-> hP`. Since `m` divides
`N`, the ideal in `(7)` is stable. In the monomial basis of degree below
`Nk`, the weight of `P^j` is `j mod m`; hence

```text
A_(N,k)=direct_sum_(r in C_m) A_(N,k)^[r],
A_(N,k)^[r]=span_K{P^j:0<=j<Nk, j=r mod m}.           (8)
```

The invariant part is precisely the pulled-back sixty-clock algebra:

```text
K[Q]/((Q^60-1)^k)  ~=  A_(N,k)^(mu_m),
Q |-> P^m.                                           (9)
```

The map is well defined because `(Q^60-1)^k` maps to `(P^N-1)^k`. Its image
has basis `1,P^m,...,P^(m(60k-1))`, which is the complete invariant basis in
`(8)`; both sides have dimension `60k`. Thus retaining only observables
pulled back from `C_60` retains weight zero and destroys every nonzero
`mu_m`-isotypic component. The torsor loss of THM-4055 and the Hasse alias
are two representations of one kernel, not merely clocks of the same size.

## 3. Fourier interpolation identifies the destroyed characters

Let `w:C_p->K`, fix `zeta_p=zeta_N^(N/p)`, and define

```text
hat w(h)=sum_(a=0)^(p-1) w(a) zeta_p^(-ah).
```

The exact interpolant

```text
I_w(P)=1/p sum_(h=0)^(p-1) hat w(h) P^((N/p)h)        (10)
```

satisfies `I_w(zeta_N^t)=w(t mod p)`. Its frequency `h` has `mu_m` weight

```text
(N/p)h mod m.
```

For `e>=2`, `g=4` and `N/p=15`, so the weight is

```text
15h mod m.                                           (11)
```

There is a sharp characteristic-zero nonvanishing law. If `p>=2` and
`w:C_p->Q` has least period `p`, then

```text
hat w(h)!=0                  for every odd h.          (12)
```

To prove it, let `W(X)=sum_(a=0)^(p-1)w(a)X^a`. For odd `h`,
`zeta_p^(-h)` is primitive. If `hat w(h)=0`, its minimal polynomial

```text
Phi_p(X)=X^(p/2)+1                                   (13)
```

divides `W` over `Q`. Since `deg W<p`, comparison of the two halves gives
`w(a+p/2)=w(a)` for every `a`, contradicting least period. This also applies
to the sign response `s=1-2w` of a binary word, because `w |-> s` is
injective and preserves least period.

For `p>=8`, `m>=2`. Equations `(11)--(12)` show that every odd `mu_m` weight
occurs and is nonzero. In particular a rational least-period response with
`p>4` cannot live in the invariant algebra `(9)`.

Characteristic zero is load-bearing. In `Fbar_13`, take `zeta_4=5` and the
least-period-four word

```text
w=(5,1,0,0).
```

Then `zeta_4^(-1)=8` and

```text
hat w(1)=5+8=0 in F_13.                               (14)
```

Thus even characteristic prime to `p` does not suffice for `(12)`.

## 4. Exact census of the fixed Rule-30 bank

Now use only THM-4047's single-seed fixed left-front tails
`ell_r`, `0<=r<=100000`. For a tail of least period `p`, evaluate its response
on `C_N`, with `N` from `(1)`. Write

```text
r_60=t mod 60,
t/N=a/d in lowest terms,
D=D(a/d),       epsilon=(-1)^D,                       (15)
```

and let `lambda` be THM-4059's four-bit flank packet. At `d=1`, `D`,
`epsilon`, and `lambda` are the explicit sentinel `None`; no Stern depth is
assigned to `0/1`.

A tail **factors through** a displayed feature when its response is constant
on every feature fibre in `C_N`. The complete census is:

| retained feature | p=1 | p=2 | p=4 | p=8 | p=16 | p=32 | total |
|---|---:|---:|---:|---:|---:|---:|---:|
| `r_60` | 16 | 10 | 56 | 0 | 0 | 0 | 82 |
| `(r_60,d)` | 16 | 10 | 56 | 49 | 41 | 0 | 172 |
| `(r_60,epsilon)` | 16 | 10 | 56 | 0 | 0 | 0 | 82 |
| `(r_60,d,epsilon)` | 16 | 10 | 56 | 49 | 41 | 0 | 172 |
| `(r_60,d,D)` | 16 | 10 | 56 | 329 | 136 | 0 | 547 |
| `(r_60,d,lambda)` | 16 | 10 | 56 | 49 | 41 | 0 | 172 |
| `(r_60,d,a)` | 16 | 10 | 56 | 668 | 87118 | 12133 | 100001 |

Thus the unlabelled exact denominator sometimes recovers responses beyond
the bare clock, but depth sign and flank parity add nothing to it in this
bank. Full depth raises the count to `547`, while the labelled numerator is
lossless for all `100001` tails because `(d,a)` reconstructs the point of
`C_N`.

The kernel-character census below is for the **sign response** `1-2ell_r`,
not for the raw bit word. The injective bit/sign change leaves all factor
counts above unchanged. For periods exceeding four, the exact nonempty
`mu_m` weight sets and their multiplicities are:

| p | weight set(s) : multiplicity |
|---:|---|
| 8 | `{0,1}:625`; `{1}:43` |
| 16 | `{0,1,2,3}:84078`; `{0,1,3}:1360`; `{1,2,3}:1385`; `{1,3}:295` |
| 32 | `{0,1,2,3,4,5,6,7}:12015`; `{0,1,2,3,5,6,7}:61`; `{0,1,3,4,5,7}:4`; `{0,1,3,5,7}:3`; `{1,2,3,4,5,6,7}:40`; `{1,2,3,5,6,7}:7`; `{1,3,5,7}:3` |

Every set contains every odd weight, as `(12)` predicts. Exactly

```text
43+295+3=341                                         (16)
```

tails have only the odd weights. Equivalently, their sign word satisfies
`s(t+p/2)=-s(t)`, or their bit word has complementary halves.

## 5. A period-minimal alias with every coarse sidecar fixed

The canonical first useful hostile is the certified period-eight tail

```text
ell_29: start=90,
word=(1,1,0,1,0,0,1,0) indexed by absolute C_8 phase. (17)
```

Use `N=120`. The two times `90` and `150` have the same `C_60` phase `30`.
Their `C_120` residues and reduced fractions are

```text
90 -> 90 -> 3/4,             150 -> 30 -> 1/4.       (18)
```

They have packets

```text
(d,a,D,epsilon,lambda)
 =(4,3,3,-1,(0,1,1,1)),
 =(4,1,3,-1,(0,1,1,1)),                              (19)
```

but

```text
ell_29(90)=0,                 ell_29(150)=1.          (20)
```

Consequently `C_60` plus denominator, full depth, depth sign, or flank
parity still aliases the response. The labelled numerator `a`, equivalently
the point-torsor coordinate of `(6)`, is the first displayed separator.

The sign Fourier support of `(17)` is exactly

```text
{1,3,5,7} in C_8.                                    (21)
```

Under `(10)` these lift to exponents `{15,45,75,105}` in `C_120`; all have
the nontrivial `mu_2` weight. Thus `(18)--(21)` are the point-side and
character-side descriptions of the same lost bit. Period eight is minimal
for such a clock obstruction, since every dyadic period at most four divides
sixty.

## 6. Consecutive Hasse data do not inherit an integral Smith form

Exponentiation identifies finite cyclic point sets after choosing roots. It
does not identify the arithmetic evaluation orders used by THM-4010. Already
at three simple points,

```text
Q[X]/(X(X-1)(X-2)) ~= Q^3,
Q[P]/(P^3-1)       ~= Q x Q(zeta_3).                 (22)
```

These are not isomorphic as `Q`-algebras: the first has three split primitive
idempotent factors, while the second has one rational and one quadratic
field factor. Over a splitting field both become products of fields, but
that scalar extension still does not identify the integral lattices whose
Smith invariants THM-4010 computes.

Therefore no consecutive-node Smith exponent, factorial index, or modular
rank law transfers to the cyclotomic observer **by exponentiation alone**.
Any arithmetic comparison must specify a base order, a lattice, and an
integral comparison map. Equation `(22)` is a firewall, not a claim that
every cyclotomic arithmetic use requires a splitting field.

## 7. Verification and validity boundary

The primary script reconstructs the complete THM-4047 bank by the affine
recurrence with exact physical-front queries, checks all `100001` factor
flags, and computes the sign Fourier sectors by exact cyclotomic reduction.
The independent script does not import it: it reconstructs the bank by a
separate recurrence/label/subtractive-depth path, conditional on the seven
already-audited THM-4047 physical-query bits. It reproduces the two semantic
digests, hostile `(18)--(21)`, finite-character counterexample, fibre sizes,
and invariant dimensions. The universal statements `(2)--(14)` are proved
algebraically; they are not finite extrapolations.

The boundary is strict:

- the census concerns fixed left-front columns only; the Rule-30 center moves;
- the result proves no temporal balance, density, unbounded-growth statement,
  or other Rule-30 prize;
- denominator, Stern depth, and flank packets are observers, not dynamical
  causes;
- the Hasse decomposition gives no new Keller pair and no planar-JC result;
- `(12)` is a rational characteristic-zero theorem, with `(14)` its explicit
  finite-character failure;
- labelled numerator data are lossless here because they retain the full
  combined point, not because a coarser arithmetic statistic became complete;
- no integral Smith statement survives without the extra data named after
  `(22)`.

This theorem therefore closes the representation question left by THM-4055:
the missing sidecar is simultaneously a point of the `mu_m` torsor and a
nonzero `mu_m` character sector. It does not close any Rule-30 prize.
