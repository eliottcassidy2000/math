# Source-weighted inverse rays: bounded charge, sharp two-thirds debt

**Status: PROVED elementary, same-agent review and exact checks.** No independent
agent audit or universal Collatz conclusion is claimed. All cutoffs below are
positive integers. This is a result about the actual native charge, not a Haar
surrogate or a count of residue classes.

Inherited definitions: [paid native forcing and first-arrival frontier](../paid_root_cycle_budget_20261026/README.md),
[Eisenstein first-arrival energy](../root_eisenstein_first_arrival_20261018j.md),
and [moving-modulus renewal](../weighted_renewal_moving_modulus_20261026/README.md).
The direct-root rays existed in the earlier renewal obstruction. The new
contribution is their uniform cancellation and source-weighted refinement law.

## 1. Actual rays, with their integer clock

Let S be the positive odd integers other than 5 modulo8, and B(x)=4x+1.
For s in S, set a=2 when s=1 modulo8 and a=1 otherwise. Start at
(3s+1)/2^a and repeatedly undo B while the result is5 modulo8. The final
primitive value is P(s), and t=a+2j, where j is the number of undos, is its native
integer clock. In particular

    2^t (3P(s)+1) = 3(3s+1)+2^a.

Fix u in S and a in {1,2}. Write D=2^a+3, and choose the unique
tau in {a,a+2,a+4} with 2^tau(3u+1)=D modulo9. Put C=2^tau(3u+1).
Then precisely the a-preimages of u are

    s_k=(C 64^k-D)/9,  k=0,1,...;
    P(s_k)=u,  t_k=tau+6k.

Indeed t has parity a and the displayed congruence fixes its class modulo6.
The formula gives positive primitive sources of the required residue, and
(3s_k+1)/2^a=B^((t_k-a)/2)(u), so the stripping procedure is exact.
For a=1, C>=32; for a=2, C>=16. A contiguous index interval therefore also
represents an integer-clock interval, with no change of clock scale.

For first arrival at ROOT, omit the u=1,a=2,k=0 self-edge: source1 stops
at clock0. Its remaining ray is still an interval, so every bound below
survives. For a larger stopping frontier, already-absorbed sources must also
be removed according to their actual first-arrival status; this can make a
noninterval and is not silently included in the theorem.

## 2. The original weighted charge and the theorem

Let omega^2+omega+1=0 and Q(x+y omega)=x^2-xy+y^2. Set

    A_N(s) = sum_{j: B^j(s)<2N} omega^(B^j(s)) (omega-omega^2)
           + sum_{j: 2N<=B^j(s)<4N} omega^(B^j(s)+1).

These are precisely the inherited charges; each Q(A_N(s)) is0,1,or3.
The transported charges on a ray are (-1)^tau A_N(s_k), since its clock
increments are even. Their common sign does not change any assertion below.

For every u,a,N and contiguous index interval I (finite, a prefix, or a
tail), let z_k be these actual transported charges, omitting zero terms only
after evaluating them. Let T_k=sum_{j in I, j>=k} z_j, and put

    S_I = sum_{k in I} Q(z_k),
    D_I = sum_{k in I} max(0, Q(T_k)-Q(z_k)-Q(T_(k+1))).

Only finitely many charges are nonzero. Then

    Q(sum_{k in I} z_k) <= min(4,S_I),       D_I <= (2/3) S_I.       (R)

The amplitude constant4 is attained. The debt constant2/3 cannot be
uniformly decreased. These conclusions include every integer-clock interval
and every dyadic source bin intersected with such an interval.

## 3. Shell regrouping: why the bound is independent of depth on a ray

Reindex the original double sum by i=3k+j. With E=D-3=2^a,

    m_(i,k)=B^(i-3k)(s_k)=(C 4^i-E 4^(i-3k)-3)/9,
    0<=k<=floor(i/3).

For fixed i these m increase strictly with k. Their minimum is
((C-E)4^i-3)/9 and their maximum is less than (C4^i-3)/9.
The assumptions C>=16,E<=4 imply both

    max(shell i) < 2 min(shell i),
    2 max(shell i) < min(shell i+1).

For example the first difference is bounded below by
((C-2E)4^i+E4^(i mod3)-3)/9>0. The second follows from
((2C-4E)4^i+2E4^(i mod3)+3)/9>0.
Consequently the annulus [2N,4N) can interact with only one shell; the two
thresholds cannot both split shells. A threshold exactly at a lattice point
cannot occur here because m is odd and 2N,4N are even.

Modulo3, 64^k=1+9k modulo27 gives

    m_(i,k)=s_0+i+eta k modulo3,
    eta=-1 for a=1, eta=1 for a=2.

Write i=3ell+r. The complete shell sum is
omega^(s_0+r) sum_(k=0)^ell omega^(eta k). The three shells at
i=3ell,3ell+1,3ell+2 therefore cancel exactly. Truncating k at K simply
replaces ell by min(ell,K) in all three terms, preserving cancellation.
Starting at k=L replaces C by C64^L; the same proof applies.

This is an exact regrouping of the original source weights. It does not
reassign source labels or replace their clocks by shell number i.

## 4. Two-block normal form and the finite induction

Up to a common Eisenstein unit and the transported sign, the nonzero support
with possible zero padding has the form

    z_k = xi^k f (0<=k<m),   z_k=xi^k b (m<=k<m+n),
    xi=omega^eta.

Here r is0,1,or2, U_r is respectively0,1,1+omega, and the only types are

| threshold splitting the shell | f | b |
|---|---|---|
| 4N | omega | (omega-omega^2) U_r |
| 2N | omega-omega^(r+2) | omega |

Both thresholds in gaps are included by an empty or complete block. To see
this per source, the completed lower shells contribute U_r; the partial
shell contributes omega^r precisely for its initial k-segment. Use
(omega-omega^2)U_r+omega^(r+1)=omega. For a restricted source interval the
same word is clipped, and only its common phase and block lengths change.

Three consecutive terms of either constant-coefficient block sum to zero.
Their isolated positive debt is zero. Thus total charge depends only on
m,n modulo3, and

    D(m,n+3)=D(m,n),
    D(m+3,n)=D(m,n)+delta_(m mod3,n mod3).

The second increment repeats because the tail total is unchanged. This
reduces the all-length assertion to108 elementary cases: two eta, two
threshold types, three r, and nine residue pairs (m,n).
Direct Eisenstein arithmetic gives, in every case,

    Q(total) <= min(4, Q(f)m+Q(b)n),
    3D(m,n) <= 2(Q(f)m+Q(b)n),
    0<=delta<=2Q(f),                0<=m,n<3.

The possible (Q(f),delta) are (0,0),(1,0),(1,1),(1,2),(3,0).
Adding triples proves (R) for all lengths. `block_proof()` in
[the checker](check_root_rays.py) evaluates the entire finite certificate,
as well as separate recurrence controls. Its checks are exceptions and
remain active under `python -O`. The proof of the recurrence above supplies
the all-order step; bounded tests are not used as that step.

## 5. Why this is the actual dyadic separator debt

For k<l,

    v2((s_l-s_k)/2)=v2(C)-1+6k.

Thus modulo2^q in the true address (s-1)/2 there are singleton early
sources and at most one tail. Each nontrivial refinement splits off exactly
one source from that tail. Its energy increment is precisely the summand
defining D_I. Intermediate dyadic levels have zero increment. Intersecting
with a clock interval clips the same tail, rather than replacing it with
independent random samples.

The ray's minimal useful message can consequently retain
(u,a,N,index interval,integer clock offset,complex weighted sum). For an
operation that needs finer first-arrival membership, that membership must
also be retained. The message is sufficient for the declared interval/bin
grammar, not for arbitrary source deletions.

## 6. Sharp hostile family and the remaining cross-ray obligation

Take u=1,a=1, C=32, D=5. For h>=0 set

    N_h=(8*4^(9h+10)-8186)/9.

The unsigned charge word is exactly

    (omega,1,omega^2)^h, omega,1, -2-omega, 1+2omega.

To derive it, the active shell is i=9h+10 (r=1), and 4N crosses exactly
through k=3h+1: N_h=(m_(i,3h+1)+3)/4. The next term is beyond the threshold.
The two-block formula gives the word, including its two final norm3 terms.
One direct four-term calculation and the zero-sum triples yield

    S_h=3h+8,     D_h=2h+2,     Q(total)=4.

Hence D_h/S_h tends to2/3, and exceeds1/2 for h>=5. The first member of
this displayed family violating the half-budget has h=5 and
N=1153843746341072806340110295381334, with D=12,S=23.
This is not a counterexample to the global ROOT target; it refutes assigning
its half-debt budget independently to each ray. The uniformly bounded
amplitude coexists with linearly growing positive and negative refinement
increments, so amplitude alone does not control total positive variation.

Nor can four be applied after merging different rays. At u=1,N=1864022,
the signed sum of the two complete native a-rays is -1-4omega, of norm13.
The self-source has positive sign in either clock convention but must still
retain clock0 for first arrival. All larger multistep uses remain separate.

A constructive next coupling is exact: except for the minimal odd-clock
head when u=2 modulo3, the a=1 sources pair with a=2 sources by

    s_odd = 2 s_even + 1,    t_odd=t_even+1.

This follows from tau_odd=tau_even+1, or from that equality after shifting
the odd ray by one index when the two least positive clocks wrap modulo6.
Their signs are opposite. The task is to bound their **joint** dyadic
increments with their actual distinct source addresses, rather than summing
two positive-part invoices. Iterating to first-arrival ROOT or the627 frontier
also introduces different downstream phases and noninterval membership.
Uniform control of those cross-ray terms remains OPEN.

## 7. A paid pairing that removes the hostile family's growing debt

The second ray supplies a concrete repair of the sharp family, rather than
just a proposed cancellation. At the same N_h, its a=2 charge word is

    (-1-2omega,2+omega,-1+omega)^h,
    -1-2omega,2+omega,omega,omega^2.

It has source energy9h+8 and total -omega. The odd ray transports with
negative sign and has total -2omega. Including both gives total -3omega
of norm9, against original source energy12h+16.

Write the signed dyadic increment for splitting off source k as
Delta_k=Q(T_k)-Q(z_k)-Q(T_(k+1)). Direct evaluation of the displayed words
gives, from their small-source ends,

    odd ray:  (2,-3,-2)^h, 2,-3,-3,0;
    even ray: (-3)^(3h+2), -1,0.

For target1 the even ray splits at dyadic exponents3+6k and the odd ray
at4+6k, by the exact valuation identity in section5. They are consecutive
actual refinement stages. Pair these consecutive increments **before**
taking the positive part. Every paired increment is nonpositive. The
initial split between the two rays contributes4, so

    unpaired positive debt = 2h+6,
    paired positive debt = 4,       for every h>=0.

This is an exact all-h repair of this obstruction. It preserves both source
weights, both source addresses and all clocks. The source1 is already stopped
at ROOT with clock0; its positive charge is the same, and it is kept with that
tag rather than assigned a fictitious first-arrival loop. The dyadic pairing
is by the proved address exponents, not by overwriting that exceptional clock.
All other paired source clocks differ by1 as stated in section6.

The finite phase can be stored separately from two integer counters: total
charge depends on block lengths modulo3; source energy and positive variation
depend on the quotient counters. This is a small exact controller for the
declared ray grammar, not a finite-state representation of arbitrary Collatz
continuation. The remaining target is a uniform compatible pairing for other
cutoffs and targets, then for genuinely changing downstream controllers.

There is an explicit guard against overextending this repair. For target7,
take M_h=(22*4^(9h+11)-65530)/9. Its two signed words are

    (1+2omega,1-omega,-2-omega)^h, 1+2omega,-omega,-1,1+omega;
    (omega,omega^2,1)^h, omega,omega^2,2+omega,-1+omega.

The second ray's active upper shell is i=9h+11, cut through k=3h+1;
the first ray's lower threshold gives the first displayed word. The same
shell comparison used above verifies both for every h. Their paired
increments are (2,-9,-5)^h followed by2,-5,-4,0, with initial cross split6.
Consequently source energy is12h+14, endpoint energy13, but paired positive
debt is2h+8 and is unbounded. The pairing does not have a universal constant
invoice even for one fixed native target. Its *relative* cost tends to1/6,
so this is not a refutation of a suitably weighted budget.

Targets1 and7 have the same least clock choices (3,2), but different marked
Eisenstein phases. Their contrasting words show exactly why the terminal
source/phase cannot be discarded. The appropriate next target is a weighted,
phase-resolved paired bound, with endpoint correction and real integer
clock, rather than universal constant debt or a phase-blind controller.

**OPEN candidate, tested independently of the above proof:** for every
complete one-step inverse fibre and cutoff, let D_pair include the positive
initial odd/even-ray split, any unpaired minimal odd-clock head (when u=2
modulo3), and positive parts of the paired consecutive internal increments.
The candidate is

    D_pair <= S/2 + 3.

`check_paired_probe.py` tests217541 (u,N) pairs:98 targets, every N<=2048,
and adjacent cutoff events in37 shells for9 declared targets. It finds no
failure. The additive3 is attained at u=1,N=3, where D_pair=S=6. This is
FINITE-EXACT evidence only. It does not justify an all-parameter promotion,
nor summing a cost3 over uncontrolled terminal sets.

A concrete proof route is a coupled two-block automaton: retain target
phase modulo9, the relative positions of the two cutoffs, the block residues,
the exceptional head, and the exact integer clock offset. Price each allowed
transition by twice its positive paired increment minus its newly owned
source energy. To certify the candidate, one would need admissible transition
inequalities for the actual shared C and a potential bounding all path costs
by6. Arbitrarily combining two independent formal ray words changes the
source law and is not an admissible substitute.

## 8. Reproduction and exact scope

Run `python3 check_root_rays.py` from this directory (also with `-O`).
The receipt covers97 primitive targets through257, every N=1..2048, all
adjacent threshold events in24 shells, additional clipped-ray events through
shell240 and cutoffs of677 bits, and2080 dyadic-bin controls. The explicit
sharp family is separately checked through h=100. These are independent
literal B-chain versus shell-decoder calculations within one implementation
packet; no independent-agent audit is implied.
