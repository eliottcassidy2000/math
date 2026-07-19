---
id: THM-1283
title: A PROTRUDING CENTERED SURVIVOR EXPORTS A PROPER CARRIER-ENDPOINT-OWNER SEAM AND A STRICT RESIDUE/GCD TAIL TAX
status: PROVED (every protruding endpoint has a non-slowest fast owner; exact signed endpoint residual and outward tooth length; mirrored proper carrier-owner crossing; exterior seam disjoint from the full internal chronological invoice; survivor subtraction gives ell-eta>11/270; exact centered-error/residue and gcd integer cuts; center-to-endpoint winding congruence; terminal-word return-or-unique-owner corollary; the THM-1266 sharp five-rung local row is globally excluded; optimization-safe exact referee; sorry-free Lean arithmetic core). This strengthens the centered protrusion constraint and consumes THM-1274's endpoint branch, but does not prove six-comb noncoverage or LRC(14)
source: codex-2026-07-19-S82 terminal-endpoint continuation
depends_on: [THM-1198, THM-1237, THM-1250, THM-1253, THM-1264, THM-1267, THM-1274]
related: [THM-1199, THM-1252, THM-1266]
script: 04-computation/lrc14_terminal_endpoint_transfer_gcd_tax_thm1283.py
output: 05-knowledge/results/lrc14_terminal_endpoint_transfer_gcd_tax_thm1283.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCTerminalEndpointTransferGcdTax.lean
script_sha256: 357585477891fe1ed1e1fcaf970f7d73d94d860962e3c7961029d0f65e979f61
output_sha256: bfbc2c02aabefca2103acec72ea7bf134f26df63f1112249686c04c56fdeee57
formalization_sha256: 97adf0b28be6ef6abce4066ee47c7f048db851ce4ae93885e58e39b1b880de43
---

# THM-1283 — terminal endpoint transfer and gcd tax

## 1. The protruding endpoint has a fast owner

Use the setup and notation of THM-1267.  Thus

```text
G=[(14k+1)/(14c),(14k+13)/(14c)],
c<a=d_1<d_2<...<d_6,
```

the six strict fast danger combs cover `G`, and `S_a` is the complete closed
`a`-safe component through the centered `a`-spoke.  Its one-sided exterior
tail is nonempty.  Write

```text
sigma in {-1,+1}
```

for the protrusion side, `e_sigma` for the corresponding endpoint of `G`, and

```text
E_sigma=S_a minus G,
ell=|E_sigma|/|S_a|,
|S_a|=6/(7a).                                         (1)
```

THM-1267 also supplies the five-comb closed survivor `U` on `S_a` and the
six-bin density `f`, with

```text
U subset E_sigma,
integral_U f>11/360.                                  (2)
```

The endpoint `e_sigma` lies in the interior of `S_a`; hence `a` is not
strictly dangerous there.  Since the six open combs cover the closed gap
`G`, at least one owner

```text
x in {d_2,...,d_6}                                    (3)
```

has a strict danger tooth `T(x,n)` containing `e_sigma`.  This endpoint
owner exists before any chronological-word choice.  In THM-1274's terminal
branch it can be chosen to be the owner of the final selected tooth reached
by the protrusion-facing continuation.

## 2. The exact endpoint residual

Put

```text
e_+=(14k+13)/(14c),        e_-=(14k+1)/(14c).         (4)
```

Define the oriented integer residual `z_sigma` by

```text
z_+ =14cn-(14k+13)x,
z_- =(14k+1)x-14cn,                                  (5)
```

and define the positive endpoint numerator

```text
Q_sigma=c+z_sigma.                                    (6)
```

Strict containment of `e_sigma` in `T(x,n)` says exactly

```text
-c<z_sigma<c.                                        (7)
```

If `w_sigma` is the wall of `T(x,n)` lying outward from `G`, direct endpoint
subtraction gives

```text
q_sigma=|w_sigma-e_sigma|
       =Q_sigma/(14cx),
0<Q_sigma<2c,
0<q_sigma<1/(7x).                                    (8)
```

Let `g=gcd(c,x)`.  Equations (5)--(6) retain two pieces of arithmetic which
a bare endpoint incidence forgets:

```text
g divides Q_sigma,
Q_sigma == c+x (mod 14).                              (9)
```

Consequently

```text
Q_sigma>=g,
q_sigma>=g/[14cx]=1/[14 lcm(c,x)].                   (10)
```

There is a sharper residue version if desired.  Put

```text
h=gcd(g,14),       C=c/g,       X=x/g,
r_14=the least positive residue of C+X modulo 14/h.
```

Then `Q_sigma/g==C+X (mod 14/h)`, so

```text
Q_sigma>=g r_14.                                     (11)
```

Unlike an abstract lcm quantum, `Q_sigma` is attached to the actual endpoint,
side, and tooth address.

## 3. The adjacent carrier tooth makes a proper crossing seam

Immediately outside `G`, the carrier `c` becomes strictly dangerous.  On the
right the adjacent `c`-tooth is

```text
C_+=(e_+,e_++1/(7c)),
```

and on the left it is

```text
C_-=(e_--1/(7c),e_-).                                 (12)
```

Since `x>c`, equation (8) gives

```text
q_sigma<1/(7x)<1/(7c).                               (13)
```

The survivor `U` has positive mass and avoids the **closed** `x`-danger
comb.  Therefore the outward `x`-wall must occur strictly before the outer
endpoint of `S_a`; otherwise the entire exterior tail would lie in the one
closed `x`-tooth and `U` would be empty.  Thus the open interval

```text
N_sigma=(e_sigma,w_sigma)       in the outward order  (14)
```

lies in `E_sigma` and in both strict danger teeth `C_sigma` and `T(x,n)`.
It has the exact length `q_sigma`.

The two teeth cross properly.  In increasing real order their four walls are

```text
right protrusion:  L(T_x)<e_+=L(C_+)<R(T_x)<R(C_+),
left protrusion:   L(C_-)<L(T_x)<e_-=R(C_-)<R(T_x).   (15)
```

So reflection reverses the chronological order but preserves the outward
owner transfer `x -> c`.  At `w_sigma`, `x` is equality-safe while `c` is
still strictly dangerous.  This is a canonical neighboring-carrier wall,
not a sampled abstract pair edge.

Every raw handoff in THM-1253's full chronological seam invoice lies in
`int(G)`, whereas `N_sigma` lies strictly outside `G`.  Hence the endpoint
seam is disjoint from every already-paid internal handoff.  Its measure and
the internal seam invoice can be used simultaneously; no multiplicity or
Hunter credit has been counted twice.

There is a useful protected-needle corollary which does not require the
centered component.  Suppose a protected interval `I` contains `G` and

```text
|I|-|G|>=1/c.                                        (15a)
```

The two exterior margins of `G` in `I` sum to at least `1/c`, so one has
length at least `1/(2c)`.  At that endpoint choose any fast tooth which
covers it.  Its outward extension is less than

```text
1/(7x)<1/(7c)<1/(2c).
```

Hence its complete carrier-owner seam lies in `I`.  THM-1250 supplies a
located five-edge chronological spanning tree on the six fast owners inside
`G`; adjoining this exterior edge `{c,x}` gives a located six-edge spanning
tree on all seven owner vertices.  The exterior raw interval is disjoint
from the five internal tree seams.  Its scalar lcm term may already be
dominated by an existing protected-needle max-tree tariff, so this corollary
does not claim a second global scalar credit.  Its irreducible gain is that
the carrier vertex `c` is now attached at a named wall and address.

## 4. Subtracting the seam from the survivor tail

Normalize `S_a` affinely to `[0,1]`.  The endpoint seam (14) is the **inner**
prefix of the endpoint tail, of normalized length

```text
eta_sigma=q_sigma/|S_a|
         =(7a/6)q_sigma
         =a Q_sigma/(12cx)
         >=a g/(12cx).                               (16)
```

Because `U` avoids the closed `x` comb, it lies beyond `w_sigma`.  Therefore
it lies in the outer endpoint interval of normalized length

```text
y=ell-eta_sigma.                                      (17)
```

This observation removes the case split which lost the endpoint incidence
in THM-1267.  If `y<=1/6`, the density on that endpoint interval is exactly
`3/4`, so (2) gives

```text
11/360<integral_U f<=(3/4)y,
y>11/270.                                             (18)
```

If `y>1/6`, the same conclusion is automatic because `1/6>11/270`.
Consequently the uniform strict endpoint-tax law is

```text
ell-eta_sigma>11/270,
ell>11/270+a Q_sigma/(12cx)
   >=11/270+a gcd(c,x)/(12cx).                        (19)
```

This is genuinely stronger than `ell>11/270`: a fast tooth which covers the
carrier endpoint occupies a positive inner part of the protrusion where the
five-comb survivor cannot live.

The argument applies separately to every fast tooth containing the endpoint,
not only to the selected terminal occurrence.  Their outward segments are
nested at the common endpoint, and `U` lies beyond all of their walls.  Thus
one may replace `eta_sigma` in (19) by the maximum endpoint-owner value.  The
segments must not be summed: nested incidence supplies a maximum tax, not
several disjoint copies of measure.

## 5. Exact residue and gcd cuts

THM-1267's centered formula and nearest-integer bound are

```text
ell=1/2+7rho/6-a/(2c),             rho<=1/2.          (20)
```

Combining (19)--(20) yields

```text
270ax+45a Q_sigma<563cx.                             (21)
```

All entries are integers, so the exact phase-located cut is

```text
270ax+45a Q_sigma<=563cx-1.                          (22)
```

Using (10), every endpoint owner satisfies the coarser but phase-free gcd
cut

```text
270ax+45a gcd(c,x)<=563cx-1.                         (23)
```

Equivalently,

```text
(a/c)(1+gcd(c,x)/(6x))<563/270.                      (24)
```

Thus the old slowest-fast ratio ceiling is never attained even at its
former integer edge: the endpoint owner pays a located positive correction.
Equations (9), (11), and (22) are the stronger finite residue filter.

### 5.1 Retaining the centered error gives the sharp phase cut

The coarse constant `563` uses only `rho<=1/2`.  The actual centered phase
retains an integer which should not be discarded.  Put

```text
m=2k+1,
2cp-(c+a)m=sigma e,             e=2c rho in {1,...,c}. (24a)
```

Substituting `rho=e/(2c)` in (20) gives

```text
ell=(6c+7e-6a)/(12c).                              (24b)
```

Now (19) is equivalent to the strictly sharper exact phase cut

```text
270ax+45a Q_sigma<(248c+315e)x,                    (24c)
270ax+45a Q_sigma<=(248c+315e)x-1.                 (24d)
```

The earlier cut (22) is its corollary because `e<=c`.  Thus `e`, not merely
the sign `sigma`, is part of the terminal state.

There is also an exact center-to-endpoint congruence.  Equation (24a) gives

```text
(c+a)m == -sigma e                         (mod 2c). (24e)
```

The endpoint formulas (5)--(6) give

```text
Q_sigma-c+6x is divisible by 7,
x m == -sigma (Q_sigma-c+6x)/7             (mod 2c). (24f)
```

Eliminating `m` yields

```text
(c+a)(Q_sigma-c+6x) == 7ex                 (mod 14c), (24g)
gcd(c+a,14c) divides 7ex.                            (24h)
```

The congruence has a positive winding lift.  The center `n/x` of the endpoint
owner tooth remains on the outward side of the centered spoke `p/(c+a)`:
THM-1240's boundary clearance is greater than `5/(28c)`, whereas the endpoint
owner center is less than `1/(14x)<2/(28c)` from the endpoint.  Hence

```text
W=sigma[n(c+a)-px] is a positive integer,             (24i)

(c+a)(Q_sigma-c+6x)-7ex=14c W>0.                     (24j)
```

This is the exact neighboring-carrier/address descent payload.  Projecting
to the gcd cut (23) loses both `e` and the positive winding `W`.

### 5.2 The sharp local five-rung row fails the global endpoint condition

THM-1266's primitive positive control has

```text
(c,a,k,p)=(140,254,80,227),
sigma=+1,       e=126.
```

At the right endpoint, owner `x=256` with tooth address `n=148` has

```text
Q_+=172,
ell-eta_+=25/1536<11/270.                            (24k)
```

Numerically, the old coarse cut passes, but the exact phase cut fails:

```text
270ax+45aQ_+ =19,522,440,
(248c+315e)x-1=19,048,959.                           (24l)
```

Therefore the row which simultaneously realizes all five local toothpick
rungs and both aligned blocker cycles cannot be completed to a six-cover:
its endpoint-owned exterior segment leaves too little outer tail for the
compulsory five-comb survivor.  This is a global exclusion of the principal
local positive control, not a merely local return law.

## 6. The THM-1274 terminal word now has a typed child

In a deletion-minimal cover, at most one selected tooth can contain a fixed
endpoint of `G`.  At the right endpoint, if two did, the tooth with later
left endpoint would be contained in the earlier tooth after restriction to
`G`; at the left endpoint the tooth with earlier right endpoint would be
contained in the later tooth after restriction to `G`.  In either case one
would be deletable.  It follows that THM-1274's no-return continuation lands
on the unique terminal selected occurrence `T(x,n)`.

That occurrence has two exact outgoing payloads:

1. the exterior proper crossing `N_sigma`, ending at a wall where the old
   carrier `c` is strictly dangerous and `x` becomes safe; and
2. the residue/tax word `(sigma,n,Q_sigma,q_sigma)` in (5)--(10).

There is also a clean return/private-count split.  If another selected
`x`-tooth occurs, take the nearest such occurrence on the inward side of the
terminal tooth: before it at the right endpoint and after it at the left
endpoint.  In increasing chronological order, the intervening literal
consecutive subword is an endpoint-attached `x`-return with address jump
`R>=1`, so THM-1264 gives its exact holonomy

```text
sum omega=(1/7)sum 1/s_r-R/x>0.                      (25)
```

A closest repeated pair inside this subword is a six-owner return leaf of at
most six edges and carries the strict `6/5` ascent.  If no other selected
`x`-tooth exists, THM-1253's private count gives

```text
1>=ceil(x/(7c)),                 x<=7c.               (26)
```

Thus the terminal branch is no longer an untyped stop: it is either attached
to a literal return, or its endpoint owner is globally single-occurrence and
compact, and in both cases it exports the new exterior carrier seam and the
strict tax (19)--(23).  What is not proved is that the neighboring `x`-safe
component inherits a slowest-carrier six-fast cover; the owner set on that
component is mixed above and below `x`.

## 7. Exact mirrored guardrail

The THM-1248/1274 packet with

```text
(c,a,x)=(2,4,28)
```

gives an exact two-sided control.  On the left gap `k=0`, take the `28`-tooth
of address `1`; on the reflected right gap `k=1`, take address `27`.  In both
cases

```text
Q_sigma=gcd(2,28)=2,
q_sigma=1/392,
ell=1/12,
eta_sigma=1/84,
ell-eta_sigma=1/14>11/270.                           (27)
```

The carrier and endpoint-owner teeth cross in the two orders displayed in
(15).  This packet is globally lonely, so it is a guardrail rather than a
six-cover.  It shows that the gcd quantum can be sharp on both sides and that
the exterior seam is not itself a contradiction.

## 8. Tournament and alternate-vertex audit

Runner-speed order gives only the transitive edge `c<x`; it forgets the side,
tooth address, endpoint residual, and survivor exclusion.  Use instead the
two wall events as vertices:

```text
E_0: c changes safe -> danger at e_sigma,
E_1: x changes danger -> safe at w_sigma              (28)
```

in the outward gauge.  The pairwise observable is simultaneous strict
danger on the interval between the events.  Reflection is the switch, and
the outward order `E_0 -> E_1` is the tie Hamiltonian path.  The two-vertex
tournament has score histogram `(0,1)`, no directed cycles, two singleton
SCCs, and one Hamiltonian path.  Adding the outer survivor obligation gives
the forced event word

```text
carrier wall -> endpoint-owner wall -> survivor mass. (29)
```

The smallest faithful state is

```text
(sigma,e_sigma; x,n,Q_sigma; N_sigma;
 outer survivor interval and f-mass; terminal tooth occurrence).     (30)
```

It preserves proper crossing, exact address residue, invoice disjointness,
and survivor subtraction.  Projecting to runners destroys all four.  The
challenged assumption is that reaching the boundary of the chosen carrier
gap ends the chronological operation: the boundary itself begins a proper
adjacent-carrier seam.

## 9. Verification and formal boundary

The dependency-free exact referee enumerates both orientations of every
strict endpoint incidence in its finite bank.  It checks the signed residual,
the exact outward length, both congruences, the gcd/lcm quantum, proper
carrier-owner crossing, normalized tax identity, the endpoint quantile,
the integer cut, the private-count alternative, and the mirrored sharp
packet (27).  It contains no Python `assert` nodes; ordinary and optimized
outputs are byte-identical.

The sorry-free Lean module proves the endpoint-residual positivity and upper
bound, the outward-width comparisons, the endpoint-suffix quantile, the exact
normalization identity, the rational tax consumer, its gcd weakening, and
the integer strict-to-closed cut.  Selection of an endpoint owner from the
strict cover, identification of the adjacent teeth, containment of `U` in the
outer suffix, and extraction of the terminal word remain the named paper
topology providers.  There are no proof placeholders or `native_decide`
calls.

Frozen artifact hashes are

```text
source         357585477891fe1ed1e1fcaf970f7d73d94d860962e3c7961029d0f65e979f61
output         bfbc2c02aabefca2103acec72ea7bf134f26df63f1112249686c04c56fdeee57
formalization  97adf0b28be6ef6abce4066ee47c7f048db851ce4ae93885e58e39b1b880de43
```

THM-1283 supplies the phase-located endpoint tax and the first literal
neighboring-carrier seam requested after THM-1274.  It does not yet prove
that this child seam returns to a lower-rank centered spoke, nor does it add
the survivor mass and seam length as two copies of coverage excess.  The next
operation is to follow the `x`-safe side of `w_sigma` with the mixed owner set
while retaining `(sigma,Q_sigma)` and the disjoint internal seam invoice.
Global LRC(14) remains open.  ∎
