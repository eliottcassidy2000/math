# GMC(2): degree slopes beat atom separation, and the residue is a resonance band

**Session:** codex-2026-07-21-gmc2-degree-gap
**Theorem:** THM-2017
**Hypothesis:** HYP-8766

## Outcome first

The full two-dimensional nullcone conjecture is not proved here. What is proved
is a new infinite-dimensional slice: for

```text
P=Z^p a(s)+b(s)+Zbar^q c(s),       s=Z Zbar,
```

NC2 holds whenever the radial degree slopes are separated by more than one
primitive-return entropy unit. In the notation of THM-2017 this is

```text
|deg h-r deg b| >= r+1,
h=s^(pq/g)a^(q/g)c^(p/g),       r=(p+q)/g.
```

The open part of this entire three-weight family is compressed to the finite
band `-r<=deg h-r deg b<=r`. At both endpoints THM-2017 goes further: the
whole boundary layer converges to an explicit generalized hyper-Bessel
function, so only a discrete exceptional leading-coefficient locus remains.
For symmetric primitive charges and monomial `b,h`, even those exceptional
zeros are closed: the first `1/m` correction is a Bessel derivative, and ODE
uniqueness prevents the value and derivative from vanishing together.

## Why the first attractive proof was wrong

The scalar equation `E[P^m]=0` does not expose a coefficient polynomial in
independent atom variables. Distinct primitive-return atoms can land on the
same coefficient monomial and, more basically, different numerical summands
can cancel after the coefficients are specialized. The exact correction being
recorded concurrently as MISTAKE-211 is

```text
P=a Z^6+b Zbar^2+c Zbar^18,
E[P^4]=4*6!*a*b^3+4*18!*a^3*c,
```

which can vanish with `abc!=0`. Therefore "the least first-return atom must
vanish separately" is not a legal inference.

THM-2017 replaces atom separation by endpoint comparison in absolute value.
It controls the *sum of all other channels*, including channels proportional
to `m`, and proves that sum is `o(1)` relative to one nonzero factorial
endpoint. Cancellation is allowed but made too small to erase the endpoint.

## The creative coordinate: factorial slope, not charge support

The exact channel degree is affine:

```text
D(k)=d m+(e-rd)k.
```

This suggests a tournament whose vertices are channels, not monomials or
charges. Orient `k->l` when `D(k)>D(l)`, gauge by the sign of `e-rd`, and break
ties by increasing `k`. It is transitive, so the tournament itself is not the
theorem. Its value is diagnostic: it reveals the only two possible factorial
endpoints and isolates the tie/near-tie band.

Alternative vertex sets were challenged explicitly:

- monomials retain too much irrelevant phase data;
- charges forget radial factorial growth;
- primitive atoms invite the false separation inference;
- proof obligations are useful for workflow but not for asymptotics;
- radial-deficit channels preserve exactly the factorial slope and discard
  phases, so they are the right quotient for locating—not proving—dominance.

The quotient destroys coefficient phase and within-channel cancellation. The
analytic mixed-factorial bound supplies what the tournament cannot.

## Connection to incoming work

THM-2014 is the `p=q=1`, constant-endpoint model
`P=aZ+b(s)+cZbar`. It closes that full slice: for `deg b>=2` its uniform
all-channel estimate is the constant-endpoint specialization of the
`b`-dominant mechanism, while degrees zero and one are genuinely resonant and
need exact generating functions. This division is a useful hostile check:
for `P=Z+s+Zbar`, the `k=1` channel alone is `(m-1)m!`, so a naive
`E[P^m]~L(s^m)=m!` assertion is false. THM-2014 correctly excludes `d=1`
from that asymptotic and closes it differently.

HYP-8765's radial-channel tower and THM-2017 therefore meet at the right seam:
domination outside the band, exact/holonomic identities inside it.

## Challenging the apparent public-proof assumption

A current literature check found no public proof of GMC(2) or NC2. The July 20,
2026 preprint *Small Counterexamples to the Gaussian Moments Conjecture*
(arXiv:2607.18186) explicitly says dimension two remains unsettled and credits
a Sol instance for a counterexample in four variables, not for a proof in two.
The user's report may refer to a private/unpublished task result; until its
proof is available, the repository should not promote that report to a
theorem.

## The remaining attack

The degree-gap threshold is structurally sharp for ratio-to-endpoint
domination. At gap `Delta`, the first nontrivial channel has multinomial gain
`m^r` and factorial loss `m^Delta`; it is negligible only for `Delta>r`.
Inside the resonance band, trying harder to prove the same dominance statement
is the wrong task.

The proposed next move (HYP-8766) is to rescale the finite offset
`lambda=e-rd`. Fixed channels then form a generalized hypergeometric/Bessel
limit; proportional channels require a separate saddle. The target is not
"one channel wins" but "the resulting finite family of resonant entire
functions cannot be identically 1 unless the charged or neutral factor
vanishes." THM-2014 proves this program in its first nontrivial resonant model.

That is the honest frontier: full NC2 remains open, but an unbounded
three-weight region is closed and the obstruction is no longer diffuse—it is
a finite resonance band for each charge pair.
