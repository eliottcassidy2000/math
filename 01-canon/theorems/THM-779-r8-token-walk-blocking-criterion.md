---
id: THM-779
title: The r=8 token-walk criterion and redundancy cocycle at the prime-7 lens — exact event-word decidability; raw walls, active periods, and owner switches are noncompact
status: PROVED (blocking criterion, redundancy/event-word cocycle, and infinite refuters to raw-wall, active-period, owner-switch, and universal-extent targets) + VERIFIED (finite adversarial census, A8 state graph, and full-packet holonomy replay) + OPEN (normalized collision holonomy coupled to core incidence)
source: opus-2026-07-14-S302 criterion, corrected and strengthened by codex-2026-07-14-S10
depends_on:
  - THM-773   # the token k_a = -w_a^{-1} round(w_a x), the X^7 - X factorization criterion
  - THM-767   # zero-variance at the prime lens; chamber locking
  - THM-771   # the exact seven-owner defect frame
related: [HYP-2603, HYP-6835, HYP-6840, THM-777, THM-778, THM-783, THM-784, THM-786, THM-788, THM-794, MISTAKE-147, MISTAKE-148, MISTAKE-149]
verification:
  - 04-computation/lrc14_r8_token_walk_criterion_opus_S302.py
  - 05-knowledge/results/lrc14_r8_token_walk_criterion_opus_S302.out
  - 04-computation/lrc14_r8_raw_wall_refuter_codex_S10.py
  - 05-knowledge/results/lrc14_r8_raw_wall_refuter_codex_S10.out
  - 04-computation/lrc14_full_active_packet_holonomy_codex_S10.py
  - 05-knowledge/results/lrc14_full_active_packet_holonomy_codex_S10.out
---

# THM-779 — the r=8 token-walk blocking criterion

> **CORRECTION (S304/S10, MISTAKE-147/148):** the census constants (5, then
> 6 walls) are ratio-boxed artifacts; THM-784 gives unbounded runs. The
> criterion itself is unaffected. THM-786 proves the extent comparison only in
> its no-companion class; its advertised serving/sparse completion is withdrawn.
> THM-788 gives sound conditional inequalities through active f-periods, but
> THM-794/MISTAKE-149 makes active periods and genuine owner switches
> unbounded at `ceil(f/g)=2` and refutes the proposed universal extent.  The
> surviving carrier is collision return holonomy modulo diagonal sheet
> translation, with metric core incidence retained.

**Frame.** Lens c = 7, core P (|P| = 5), exceptions W = {w_1..w_8}, 7 ∤ w_a. By
THM-773's token algebra (refereed here exactly on 4,000 random rational points):
off its walls, owner a blocks exactly the sheet

`k_a(x) = -w_a^{-1} · round(w_a x)  (mod 7)`,

its walls sit at `x ∈ (1/2 + Z)/w_a` (mesh 1/w_a), crossing a wall steps the token
by `-w_a^{-1}`, and AT its own wall the owner blocks nothing (clearance exactly
1/14). The full-period structure: `x → x+1` carries every token by −1; the walk is
periodic in x with period 7.

## (1) The blocking criterion (PROVED; integer-decidable)

> The deck is fully blocked over an interval J (every sheet strictly blocked at
> every x ∈ J) **iff**
> 1. **(pieces)** on every wall-free piece of J, the eight tokens cover F_7
>    (eight values on seven letters — exactly one collision pair);
> 2. **(walls)** at every wall in J, the seven non-walling tokens are EXACTLY
>    F_7 — equivalently, the walling owner is a member of the collision pair
>    just before its wall;
> 3. **(no simultaneity)** no two owners wall at the same x in J (a double wall
>    leaves six tokens, which cannot cover F_7 — pierced instantly).

*Proof.* Zero-variance at the prime lens (THM-767(2)) gives each owner exactly
one blocked sheet off walls and none at them; coverage of Z_7 by the eight
tokens is (1); at a wall of a the other seven are constant and must cover all
seven sheets, which for seven values means exactly F_7 — and since the full
multiset on the adjacent pieces is F_7 plus one duplicate, removing a leaves F_7
iff a carried the duplicated value (2); (3) is counting. Conversely these
conditions clearly give coverage at every x. ∎

The whole check is an integer walk: sort the walls of the eight meshes, maintain
tokens, test surjectivity — no interval arithmetic. This replaces HYP-6840's
chamber/rainbow search (1,164 pieces of exact Fractions) with O(#walls) integer
steps, and reduces the r=8 escape-hatch question (Q1) to a symbolic-dynamics
question: **can the deterministic token walk remain inside the surjectivity
region SURJ ⊂ (F_7)^8 (density 28·7!/7^8 ≈ 2.45%) across many walls?**

## (2) The redundancy cocycle and exact event-word criterion

After a's wall the collision pair is (a, γ) with γ = the unique owner holding
`token_a − w_a^{-1}`; condition (2) at the NEXT wall (owner b) forces
**b ∈ {a, γ}** — the wall-owner schedule (fixed by the eight meshes) must agree
with the deterministic hop-target chain (fixed by the token algebra).  This is
an exact constraint on owner **switches**.  It does not bound repetitions of
the same owner.

At every covered generic state there is a unique duplicated sheet `z` and

```text
product_a(X-k_a)=(X^7-X)(X-z).                            (6)
```

Put `s_a=w_a^(-1)` and normalize `r_b=k_b-z`.  The offset multiset is

```text
{0,0,1,2,3,4,5,6}.                                      (7)
```

Owner `a` can cross a covered simple wall exactly when `r_a=0`; afterward

```text
r'_a=0,       r'_b=r_b+s_a       (b!=a).                 (8)
```

Indeed, at the wall the absent owner must have occupied one copy of the
duplicate root.  Its step moves it to `z-s_a`, duplicating the token already
there; renormalizing at that new root gives (8).

For an event word `a_1,...,a_m`, set

```text
S_j=sum_(ell<=j) s_(a_ell),
N_a(j)=#{ell<=j:a_ell=a}.
```

The word is supportable by a covered initial state if and only if every
occurrence `a_j=a` imposes the consistent requirement

```text
r_a^0=N_a(j-1)s_a-S_(j-1),                               (9)
```

and the fixed requirements use residue `0` at most twice and every nonzero
residue at most once.  If `h` owners remain unassigned, the number of supporting
states is `h!/product_q m_q!`, where `m_q` are the remaining multiplicities in
(7).  Formula (9) follows by iterating (8): just before event `j`, owner's
offset is `r_a^0+S_(j-1)-N_a(j-1)s_a`, and it must vanish.  Filling the unused
labels into the remaining multiset proves sufficiency and the count.

Between consecutive occurrences of one owner, the intervening inverse steps
must sum to zero.  In particular `ABA` is impossible.  This exact word test
combines directly with THM-778's centered Beatty schedule and replaces a raw
state search.

The rooted state space has

```text
binom(8,2) 6! = 20,160 = |A_8|                           (10)
```

states; there are `7*20,160=141,120` absolute covered configurations.  Ordering
the two root slots identifies the rooted set with an `A_8` torsor.  Unit anchor
moves are the seven-cycles `(1 2 ... 7)` and `(0 2 ... 7)`; their quotient is
the three-cycle `(0 1 2)`, so they generate `A_8`.  The exact transition graph
has `40,320` directed edges, indegree/outdegree two, one SCC, and `5,760`
monochromatic seven-cycles.  This also explains the old sector main term
`M_7(8)=20160/823543` in HYP-2603.

The strong connectivity is a guardrail: the normalized state automaton alone
cannot prove a universal exit.  Endpoint order and the metric core base are
indispensable.

## (3) The census (VERIFIED, adversarial)

- 120 random 8-tuples (w ≤ 500), quarter-period windows: **median maximal
  blocking run = 1 wall; 90th percentile 2; maximum 4.**
- Annealed to MAXIMIZE the run (400 steps): **K0 = 5 consecutive covered walls**
  at {8, 13, 19, 23, 92, 359, 372, 438} — the adversarial ceiling found.
- HYP-6840's one-moment survivor (P = {5,7,8,13,14}, W = {108,169,143,213,206,
  197,30,162}, x = 19/216) factors exactly as THM-773 predicted: tokens
  [None, 6, 5, 3, 1, 4, 2, 0] — absent owner 108 at its wall, the other seven a
  perfect rainbow — and its blocking run is EXACTLY 1 wall (piece-level failure
  immediately outside). The survivor was the algebra's minimal case, not the
  seed of a blocking family.

## (4) CORRECTION — every universal raw-wall K0 is false

The original version boxed the sentence “any core-safe component containing
more than `K0` walls cannot be fully blocked.” That did not follow from a finite
adversarial census, and its proposed universal form is **refuted** by
THM-784/MISTAKE-147. The fixed slow rainbow `{1,2,3,4,5,8,10}` blocks
`J=(5/16,7/20)` while `560N+1` inserts `21N` consecutive covered walls there.

The sound finite statement is only instance-relative: after the exact walk has
computed the maximal run `K(W,J)` for a specified tuple and interval, more than
that computed number forces a pierce for that instance. No uniform conclusion
can use the unnormalised number of walls.  THM-794 subsequently shows that
metric extent at the proposed scale and slow-owner switches are not uniformly
bounded either.  Direct intersection with a core-safe component must be
studied only after central packet return cycles are quotiented.

There is also a divisor-complete core-safe refuter, which puts the same failure
inside the arithmetic LRC14 residual rather than merely in an abstract deck.

Let

```text
P={1,2,11,12,13},       W_0={1,4,5,6,8,9,10},
I=[25/182,27/182].
```

At `x=1/7`, every core speed has clearance at least `1/7`; the
`13`-Lipschitz bound gives `I` inside the closed core-safe set.  Also
`I subset (1/8,3/20)`, and throughout that open chamber the seven `W_0`
tokens are the fixed permutation

```text
(0,5,4,1,6,3,2).                                        (11)
```

For any `m>=1`, put `A_m=182m+1` and `W=W_0 union {A_m}`.  At an `A_m` wall
the last owner is absent and the unchanged seven-owner stalk (11) still covers
all sheets; off its walls, that owner merely adds a duplicate.  The wall
indices inside `I` are exactly

```text
j=25m,...,27m-1,                                         (12)
```

because `(j+1/2)/A_m in I`; hence there are `2m` consecutive covered walls.
The thirteen-speed family

```text
7P union W
```

is distinct and contains a multiple of every modulus `2,...,14`.  Thus every
universal bound on raw covered-wall runs is false, even on the divisor-complete
LRC14 residual.  Taking the extra speed to be `1000` already gives eleven
covered walls in `I`.

## (5) What remains (corrected)

- **Visitor-cluster advance (PARTLY PROVED in THM-783/786):** the phi recurrence,
  period-sum law, single-visitor break, cluster balance, and a conditional
  metric-extent bound are valid.  They constrain genuine visitor clusters and
  owner switches.  They do **not** imply the sampled raw count `K0=6`: in the
  family above every complete fast-owner period is visitor-free, so the proved
  laws coexist with `2m` covered walls.  THM-794 realizes the opposite
  saturated mode: every period has all seven visitors, yet the same packet
  repeats arbitrarily often.
- **Normalized collision-holonomy exit (OPEN):** persistent-stalk contraction
  is only the first quotient.  THM-794 has `8H-8` genuine owner switches and
  `H-1` active periods at fixed ratio; one full packet returns the collision
  state by common sheet translation.  Divide out that central return action,
  then study SCC transitions using (9) and THM-778's centered mechanical
  schedule.
  THM-778 makes the schedule side exact: centered Beatty ranks merge the eight
  midpoint clocks, and their Euclidean parity cocycles recursively generate
  every owner/tie block.  Thus the lemma can be posed as non-synchronization of
  a centered mechanical schedule with this finite collision-hop transducer.
- Every specified finite component still falls to the direct piece check
  (criterion (1)).  THM-786's no-companion, transversal, and packet-polytope
  classes remain valid, but THM-794 refutes its universal extent conjecture.
  The dense residual is normalized packet persistence plus core incidence.
- **Simultaneous-wall stratum (PARTLY PROVED):** THM-778 shows that equal
  2-adic valuation gives a `gcd(u,v)` midpoint grid of double walls.  A
  core-safe interval of length at least `1/gcd(u,v)` pierces; for the five-core
  here, `gcd(u,v)>=21 max(P)/4` suffices.
- Every fixed family remains exactly decidable by the criterion and event-word
  test.  What is missing is a scale-uniform theorem on the diagonal-holonomy-
  quotiented collision coordinate and its metric embedding in the core-safe
  set.
- Lens generalization (7g | c strata, r = p+1 at other primes p ≤ 13 after
  rescaling δ) is mechanical but unwritten.
