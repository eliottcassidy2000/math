# Annular refinement, weighted controllers, and a logarithmic Gram target

2026-10-15. **PROVED:** the annular increment identity, deterministic fine-scale
envelope and conditional propagation theorem below. **FINITE-EXACT:** six
declared source windows, 32 fixed controllers, and two adaptive schedules.
**OPEN:** the uniform annular Gram premise, unconditional amplitude
nonvanishing, and universal ROOT entry. These are different targets.

[Program](../../04-computation/experiments/collatz_annular_overlap_20261015c.py),
[output](collatz_annular_overlap_20261015c.json),
[address/lattice companion](collatz_finite_address_lattice_20261015c.md).

## Inheritance and concept board

Anchor: weighted opposite-sign payment for a fixed integer source window.
Niche: finite witness detection in a lattice of observed edges.
Wildcard: the owner's 3–4–3, 4–6–4, and 4–8–4 address decompositions.

The closest proved mechanisms are the exact pair-fibre loss and full observed
graph closure in [the Gram/slack packet](collatz_gram_slack_20261015b.md),
and [THM-4658, two-octave-memory-energy-descent-critical-lines](../../01-canon/theorems/THM-4658-collatz-two-octave-memory-energy-descent-critical-lines.md).
The hostile is the native two-edge complementarity example: either edge alone
can pay nothing, while the pair pays the whole test row. The corrected near
miss is MISTAKE-618's omitted excess hypothesis, together with the false
inference from the sum of branch squared norms to contraction of their sum.
The least-used sidecar is the **literal annular source multiplicity**, rather
than a terminal occupancy count or an empirically estimated covariance.

Live board: integer clock; source multiplicities; joint cross blocks;
weighted pair cancellation; finite witness lattice; changing-generator metric.
Use the META-PATTERNS cards “Retain the product and shared labels before taking
scalar totals” and “Compose witness fibres before taking their scalar shadow.”

## 1. A three-label increment and its paid fine scales

Use the shortcut map T, ROOT-supported t, and w from THM-4658:

    T(n)=n/2 or (3n+1)/2,  t(n)=(-1)^tau(n) on ROOT, zero elsewhere,
    w(n)=3[n=2 mod3]-1,  H(N)=sum_(n<=N) w(n)t(n),
    Delta(N)=H(4N)-H(N),  D_X=sum_(X<=N<2X) Delta(N)^2.

No universal convergence is assumed. Since t(4m)=t(m) and w(4m)=w(m),

    Delta(N+1)-Delta(N) = sum_(j=1..3) w(4N+j)t(4N+j).       (A1)

Indeed, the fourth new source 4N+4 cancels the removed source N+1 exactly.
Three consecutive weights are a permutation of (-1,-1,2). Therefore
|Delta(N+1)-Delta(N)|<=4 at every integer N.

For any interval I=[a,a+L), let M_I=sum_I Delta and

    V(L)=(4/3)L(L^2-1).

The pairwise variance identity proves

    sum_I (Delta-mean_I Delta)^2 <= V(L).                     (A2)

Proof: |Delta(a+j)-Delta(a+i)|<=4(j-i), and
(1/L)sum_(i<j)(j-i)^2=L(L^2-1)/12. The bound is sharp for unrestricted
4-Lipschitz sequences; no claim of sharpness for Collatz labels is made.

For a partition P, define C_D(P)=sum_I M_I^2/L_I and
Z_D(P)=C_D(P)+sum_I V(L_I). Every interval refinement P' satisfies

    C_D(P) <= C_D(P') <= D_X <= Z_D(P') <= Z_D(P).             (A3)

For a split of lengths l,r, the mean difference has magnitude at most
2(l+r), the mean difference of the ramp 4n. Its ANOVA innovation is thus
at most V(l+r)-V(l)-V(r). Summing proves the upper monotonicity; the lower
monotonicity is ordinary ANOVA. This is the same mechanism as the sharper
mod-3 current envelope in [M1/M2](collatz_clock_mesh_20261015a.md), applied
to a different observable. It does not weaken or replace that result.

For X=2^k, b=2^floor(k/2), uniform bins of length b give

    0 <= D_X-C_D(P_b) <= (4/3)X(b^2-1) <= (4/3)X^2.          (A4)

Halving b reduces the leading fine allowance by four. Source interval
refinement and terminal-controller extension are distinct operations.
(A3) concerns actual means. An entrywise-absolute Gram relaxation need not
inherit its monotonicity under a bin split. Retain the old valid upper bound
and take the minimum with the new one; changing controllers on the same
rows separately has the proved nonincrease below.

## 2. The source row and its exact signed payment

For I=[a,a+b), count the N which use an original source n:

    c_I(n)=max(0,min(a+b-1,n-1)-max(a,ceil(n/4))+1),
    M_I=sum_n c_I(n)w(n)t(n).                                (A5)

Only a<n<=4(a+b-1) occurs. These are literal integer multiplicities. In
particular there is no dropped incoming prefix: it has already cancelled
in H(4N)-H(N). This does **not** justify deleting the prefix in the original
J-current rows.

Push all rows through a common authenticated native graph, retaining the
original source, signed step count, and terminal. If B_I(e) is the resulting
coefficient, with z_ROOT=1 and |z_e|<=1, put

    G_D=sum_I B_I B_I^T / b.

Then C_D(P_b)=z^T G_D z<=||G_D||_entrywise1. The ROOT coordinate is included.
All possible correlations between bins sharing a terminal are retained
before the absolute value is taken.

Under a common checked terminal extension A, G'_D=A^T G_D A. In each
ordered target-pair fibre let P and N be the positive and negative incoming
absolute mass, respectively. The inherited exact identity becomes

    ||G_D||_1-||G'_D||_1 = deleted_zero_mass + 2 sum min(P,N). (A6)

The zero term is allowed only with a witnessed non-ROOT cycle. Thus for a
fixed initial graph, a finite controller meets a budget B exactly when
its cumulative weighted payment is at least ||G_D(initial)||_1-B.
“Many collisions” or “many ROOT labels” cannot replace this invoice.

The companion proves finite attainment of the best value in any fixed
observed-edge union, even an infinite one. It supplies no uniform bound on
the size of the required finite witness. That distinction isolates the
missing arithmetic statement instead of hiding it in a limiting argument.

## 3. A sufficient logarithmic annular target

**Theorem A5 (PROVED implication; premise OPEN).** Suppose for every k>=32
there is a finite authenticated common controller on the rows (A5), with
X=2^k and b=2^floor(k/2), such that

    ||G_D,k||_1 <= (k^2/16) X^2.                              (AT)

Then E_(2^k)=sum_[2^k,2^(k+1)) H(n)^2 <= (k^2/16)2^(2k)
for every k>=32. In particular the logarithmic energy implication of
[G2](collatz_gram_slack_20261015b.md) applies, giving the same strictly
positive conditional residue interval as B3 of the
[amplitude boundary note](collatz_amplitude_boundary_20261013c.md).

**Proof.** From (A4), delta_k=D_X/X^2<=k^2/16+4/3. For k>=5,
sqrt(delta_k)<=k/4+1/2, since k/4+1/4>=4/3.
Write a_k=sqrt(E_X)/X. THM-4658 gives

    a_(k+2) <= (a_k+sqrt(delta_k))/2 + sqrt(3)*2^(-k/2).

For k>=6 the last term is <=1/4. If a_k<=k/4, the right side is
<=k/4+1/2=(k+2)/4. The two parity chains start from the inherited exact
census (with its (X,2X] endpoint convention converted to [X,2X)):

| k | E_(2^k) | E_(2^k)/2^(2k), displayed only |
|---|---:|---:|
| 32 | 6513622860706576771 | 0.353104203 |
| 33 | 7038100400636408653 | 0.095384047 |

Both are below k^2/16. These seeds come from the older independently
checked 2^34 census, not the new 2^37 run. Induct along both chains.
G2's tail comparison uses only sqrt(E_X)/X<=k/4 at k>=34; it is independent
of which Gram representation establishes that inequality. ∎

This is an alternative sufficient premise, not a proof that (AT) holds,
and not a claim that (AT) dominates the old J-row target. Different rows
have different source supports and observed-edge costs. The old constant
and incoming linear-log budgets remain available alternatives.

**Corollary A6 (PROVED implication; premise OPEN).** A directly linear-log
alternative is to require finite annular controllers with

    ||G_D,k||_1 <= (16k-482) X^2,  every k>=32.               (ALT)

This implies E_(2^k)<= (16k-480)2^(2k) for every k>=32. Thus for
k=34+j its annular target is 62+16j, and its resulting energy budget is
64+16j, exactly the premise of the incoming linear-log amplitude certificate
in [HYP-9283](../hypotheses/HYP-9283-collatz-energy-budgets-constant-fail-log-suffice.md).
The two preceding annular bounds at k=32,33 are part of (ALT), not silently
inferred from the measured energy seeds. We inherit the cited conditional
interval calculation; this session does not rerun its Arb certificate.

**Proof.** Put B_k=16k-480. By (A4), delta_k<=B_k-2+4/3<=B_k.
If a_k<=sqrt(B_k), then a_(k+2)<=sqrt(B_k)+sqrt(3)*2^(-k/2).
Squaring the latter adds at most

    2 sqrt(48k/2^k)+3/2^k < 3 < 32 = B_(k+2)-B_k

for k>=32: k/2^k decreases there and 48*32<2^32. The exact two seeds
above also satisfy E_(2^k)/2^(2k)<=B_k. Induct on the two parity chains. ∎

Neither numerical budget dominates the other over the entire stated range;
they cross. They give two explicit invoices for the same native
annular rows, without invoking a stochastic covariance model.

## 4. Exact census, including unsuccessful controllers

The six cases (k,depth) are (10,8),(10,16),(12,16),(14,16),(14,24),(16,32).
For each, observe the declared depth from **every** source 1..8X, stopping
at ROOT, and retain their edge union. Then extend every retained terminal
by q=1,2,4,8; also q=16 in the two still-failing cases (12,16),(14,16).
These are 32 fixed certificates, including the six initial ones.
Every extension uses the whole earlier graph and thus permits unequal-time
joins. The complete schedule is frozen before audit labels are calculated.

Selected scores below are ||G_D||_1/X^2, before the fine allowance:

| k, initial depth | initial | indicated extension | new score | budget k^2/16 | result |
|---|---:|---:|---:|---:|---|
| 10,8 | 23.289500 | q=8 | 3.234401 | 6.25 | passes |
| 10,16 | 3.234401 | q=0 | 3.234401 | 6.25 | passes |
| 12,16 | 22.841449 | q=16 | 9.096708 | 9 | fails |
| 14,16 | 46.766195 | q=16 | 5.035156 | 12.25 | passes |
| 14,24 | 13.512111 | q=1 | 11.278959 | 12.25 | passes |
| 16,32 | 15.274271 | q=0 | 15.274271 | 16 | passes |

The k=14,depth=24,q=1 success uses 320 new edges and leaves 2,244 of the
131,072 sources unresolved. It certifies an energy bound without finishing
every ROOT path. The nonzero scores are not convergence probabilities.

The cross block is useful. At k=14,depth=24, the common-terminal Gram for
the **sum of the bin means** H+Delta has bound 9.742687X^2. Combining the
separate H and Delta coarse bounds by Minkowski gives 22.752628X^2.
These numbers bound coarse sampled sums, not all of E_(4X); the integer
seams and within-bin variance still need payment. No full-window bound is
silently substituted for a coarse statistic.

### A weighted changing-terminal schedule

Keep the k=12,depth=16 graph and annular rows. In each round select the ten
non-ROOT columns of largest load sum_f |G(e,f)| (ties by terminal integer),
observe eight steps from each, retain the old graph, recompile, and update
the Gram. Stop at (AT), or at the declared cap of 30 rounds. The control
uses the same rule except that terminal integer order replaces load order.

| deterministic schedule | rounds | new distinct edges | final score / X^2 | unresolved sources |
|---|---:|---:|---:|---:|
| weighted load | 11 | 634 | 8.370896738 | 305 |
| terminal integer | 22 | 1181 | 8.832588121 | 288 |
| uniform q=16 | one batch | 1548 | 9.096708 | 298 |

The weighted schedule pays exactly 15,537,637,040 integer Gram units before
division by b=64. Its largest compiled native depth is 180. This is a
reproducible heuristic success, **not** an optimal scheduler or an all-scale
rate theorem. It uses neither ROOT audit labels nor an assumed white-noise
law. The control is not a claim of an equal-edge-budget comparison: both
adaptive rules run until their first success, and their costs are reported.

## 5. A hostile that isolates the genuinely missing ingredient

Let chi(n)=0,1,-1 on residues 0,1,2 modulo3. Set t(n)=chi(n), as a formal
bounded label sequence. It has t(1)=1 and t(2n)=-t(n) for **every** n, so
it satisfies the fourfold cancellation and every fine-scale rule above.
Nevertheless

    H(N)=-N-1_[N=2 mod3],  Delta(N)=-3N,
    D_X=(3/2) X(2X-1)(7X-1) ~ 21 X^3.                       (A7)

It fails the odd Collatz equation already at n=3. Thus source coverage by
halving identities and small local annular increments cannot prove (AT).
The required gain must exploit the **odd-step, source-weighted signed
joins**. This explains why the new universal refinement structure is useful
but insufficient, and why matched halving-only null models are a worthwhile
incoming control.

## 6. Incoming work, literature, and the next decisive questions

Incoming commit 93d1ee3885 extends the native census through n<2^37.
Its exact octave energies agree with both seeds above and report energies
17833025414500178295, 93862178945193578345, and 172802101944936995215
for k=34,35,36. We read its source and outputs; we did not repeat that large
run. Its statistical analysis is NUMERICAL, not a proof input here.
Three of five preregistered persistent-bias tests miss their bands; several
change sign. The old Gaussian surrogate also overpredicts the new window
energies. These observations weaken those particular predictions; they
neither disprove a future main term nor establish a uniform budget.

The incoming twisted split retains an exact positive halving memory after
weighting by the mod-3 character. Its averaged factor 2+sqrt(2) is a boundary
telescoping consequence of that memory, not independent evidence for a
stationary amplitude. Our hostile (A7) gives a deterministic reason to keep
the odd-step information. The four requested logic papers and the recovered
diagonal decompositions are integrated in the companion, with explicit maps
and limitations rather than new external proof assumptions.

The most useful remaining questions are now narrower:

1. Can one prove enough weighted opposite-sign overlap for the actual rows
   c_I(n)w(n), with both source clocks retained, to meet (AT)?
2. Can a block schedule retain the finite weighted improvement across k,
   with a proved invoice for the otherwise uncharged wait until a join?
3. Can a joint H/Delta bound use its observed negative cross contribution
   without assuming a fixed covariance or replacing its terminal cube?
4. Which extra odd-step constraints exclude (A7) quantitatively, rather
   than merely excluding that single example?

## Reproduction and audit

Run both scripts in normal and optimized Python; neither uses `assert` for
load-bearing checks. The annular program checks literal bin sums against an
independent descent computation **after** freezing every controller; checks
the multiplicity formula independently at k<=12; reconstructs transported
Grams from rows and pair fibres; bounds integer matrix products before using
int64; checks 3,280 short Lipschitz sequences and all their binary splits;
and verifies the two exact energy seeds with the correct endpoints.
The companion checks 750 toy graph/edge-order cases, address multiplicities,
Pell identities, the basis-change hostile, and (A7). Finite checks support
the stated finite universes; the universal statements have proofs above.

No claim of an independent peer audit is made in this session.

Reproduction commands, from the repository root:

```bash
python3 -B 04-computation/experiments/collatz_annular_overlap_20261015c.py
python3 -B -O 04-computation/experiments/collatz_annular_overlap_20261015c.py
python3 -B 04-computation/experiments/collatz_address_controls_20261015c.py
python3 -B -O 04-computation/experiments/collatz_address_controls_20261015c.py
```

Normal and optimized JSON outputs match byte for byte. SHA256:
annular `4d7eecb4a5c9b23b97cb32a0c72e089b70f604b81108e670e3aec4791fd4f113`;
address controls `2ec276d9536483ab41ef485b4514c41a71fe1151a6d45ebe35b0ef0be474c580`.
