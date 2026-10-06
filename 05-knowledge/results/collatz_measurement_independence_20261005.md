# Measurement independence: no certified bank certifies outside itself; the floor problem worked backwards; owed audits closed

2026-10-05, opus session `opus-2026-10-05-S7` (measurement-independence).
Owner's seed: "work what remains owed and new relevant tasks as they appear;
think about ways the required positive measurement bound could be proven or
computed, avoiding certificate recycling; reframe by assuming the positive
measurement bound with unique certificates and work backwards to fill the gap."

**Status.** PROVED (elementary): Proposition F (measurement independence),
Corollary F1 (bank-based kernel readouts), Lemma T4b sufficiency (constructive
Beatty words), the equivalence chain of section 3. FINITE-EXACT: the exhaustive
residue audit of the shadow theorem, the direct audit of the segment ledger,
cone classes to depth 16, the descent census to `2^28`, the first-descent
census to `2^24`, the mass split below `2^22`. VERIFIED numerical: the kernel
readouts of S7. CITED and typed: the four concurrent Codex notes of section 1.
**OPEN:** a certificate-free positive floor at any unresolved source; universal
Collatz. Not a canon promotion.

## 0. Answers

**How could the positive measurement bound be proven or computed?** Not from
measurements of the injection measure. Proposition F: every functional of
`lambda` computable from a certified bank, together with every proved global
identity, takes the same value whether an uncertified source is rooted or not,
because its forward orbit avoids the bank. Corollary F1 makes this concrete for
the Codex kernel readout `A_(m,d)=9H_(m,d+1)-8H_(m,d)`: computed from a bank
with uncertified mass `u`, the readout at any target outside the bank is at
most `-8u<0`, at every degree. S7 shows it on the bank of rooted leaves below
`2^18`: four targets just above the bank read negative at `d=0..40`; the
in-bank leaf 27 reads negative through `d=40` and equals its atom `8.87e-11`
at `d=80`. A readout certifies exactly the bank. The only non-recycling
content a floor can carry is a deadline, a bound on the odd stopping time
(T2 of the previous note, sharpened by Codex F1 to `W(n)<=2/((N+1)(N+2))`,
`N=L+K`); a proof of a deadline need not run the orbit, but every known
deadline is membership in a descent cylinder:
forward (Terras cylinders, the Codex controllers, their SF5 tree of
monotonically decreasing orbits) or backward (the shadows of the previous
note). Both are decidable families; neither reaches every integer.

**The reframing.** Assume (PMB): every source `n` has an independently
justified floor `epsilon(n)>0` on its weight, with the unique ROOT word as its
certificate. Working backwards gives the exact chain

    PMB  <=>  a computable deadline B(n) >= tau(n) for every n
         <=>  no positive integer lies in E_inf,

where `E_inf` is the set of 2-adic integers whose parity vector is rising at
every prefix, a compact set of Hausdorff dimension `h*=H(log_3 2)=0.9500`
(the Eggleston dimension of the information-dimension note). The deadline is
the depth at which `n` leaves the last rising cone that contains it. Any floor
function must therefore (i) have unbounded registers (the lookahead theorem),
(ii) respect the shadow order of Lemma B, and (iii) be logarithmic at best:
the records (`27` leaves at 96 T-steps, the maximum below `2^24` is 287)
force `B(n)>=c log n` with `c` at least the record constant. The gap is a
single counting statement: for `A` beyond about `log_2 X/(1-h*)`, no integer
below `X` lies in the depth-`A` rising cones. Section 3 measures how far the
exact census is from that statement: for `A<=23` the integers below `2^24` in
depth-`A` cones number exactly `N(A) 2^(24-A)` (ratio 1.0000 to the
rising-prefix count), and the thinning that must eventually reach zero is only
4% at `A=100`.

**Owed items closed.** The shadow theorem passes an exhaustive residue audit:
within every rising cone the chain exists with no exceptions. The converse,
that every smaller ancestor comes through a rising word, is false wordwise
(Codex integration audit, same day: 165 reaches 167 in 17 odd steps with
total valuation 27, slope below one, the carry producing the size increase)
and open setwise; the descent set contains the cone union, and my earlier
sentence claiming equality, as well as this session's first attempt to
"withdraw" the exception caveat, are both wrong and are corrected here. The
segment ledger passes a direct predecessor audit; the Beatty rule is proved in both directions by an explicit
word; the cone series reaches depth 16 (`0.46726`) against the census
`0.468669` at `2^28`; the basin minima are 29.7% of units and carry 6.2% of
the `W` head mass.

## 1. Inheritance, new tasks as they appeared, and the board

Four Codex notes landed during the previous session and were read as input:

- [source refinement floor](collatz_source_refinement_floor_20261005.md):
  SF1 `lambda(rho(n))^2 >= (160/441) W(n)^3` (sharp nonlinear transport through
  the leaf section); SF4 the evaluation-energy criterion
  `p_m >= 1/(c_0^(-1)+E_m)` with `E_m = sum (1/c_(i+1)-1/c_i)`, equivalent to
  atom positivity; SF5 a proper binary inverse tree of monotonically
  decreasing orbits with the source-only floor `W(n)>=w(3b,9b+1)`,
  `b=bitlength(n)`, recognized by a terminating descent. SF5 is a genuine
  non-recycling floor on an explicit infinite family; the family is defined by
  the very property (monotone descent to 5) that certifies it.
- [floor transport deadlines](collatz_floor_transport_deadlines_20261005.md):
  F1 `tau<=N` and `W(n)<=2/((N+1)(N+2))`, sharper than T2's `2/(L+2)`; F2 the
  superlevel set `{W>=epsilon}` is finite with largest element `S^B(1)`; F4
  conditional transport through a verified prefix; the hostile `S^j(3)` with
  `W=4/((j+2)(j+3)(j+4))->0` under inverse siblings (the sibling rays).
- [atom polynomial dual](collatz_atom_polynomial_dual_20261005.md) and
  [localized resolvent floor](collatz_localized_resolvent_floor_20261005.md):
  signed selectors whose exact readout increases to the atom, with the stated
  obligation "an independent positive bound for every source" and the
  admission that the numerical examples reuse completed words.
- [weight threshold receipts](collatz_weight_threshold_receipts_20261005.md):
  the finite superlevel compiler, and the Sun counterexample (THM-4026) as the
  hostile for witnesses escaping in ordinary size.

Concurrent commit `e0777a399` (Codex) integrated and audited the S6 note,
logging a mistakes entry (2026-10-05, refinement-shadow integration audit):
the wordwise rise/slope equivalence fails (165 -> 167), the T3 modulus is a
heuristic blanket comparison because aliases inside the residue cancel part
of the tail, T5 must count visits with multiplicity in general, T1 concerns
one-sided progressions with their lower endpoints, finite census proportions
are not densities, and no `tau = O(log n)` is proved. All accepted; the
statements of this note are made within those scopes.

Closest proved mechanism: Codex's criterion `(8)-(10)` and SF4. Canonical
hostile: the deleted geometric atom (every residue class positive, atom zero).
Corrected near miss: reading a bank-based readout as information about targets
outside the bank. Least-used sidecar: the forward orbit of the target as the
only object that can connect it to a bank. Board: **certified bank / forward
orbit avoidance / kernel readout / deadline / rising cone / Beatty depth /
basin minimum.**

## 2. Proposition F: measurement independence

Fix a finite set `C` of certified sources (orbits run to 1) and let `R(C)` be
all integers whose orbit meets `C`; these are exactly the sources certified by
`C` through transport, and `R(C)` is closed under `U`. Let `m` be a source
whose forward orbit avoids `C`, equivalently `m` is not in `R(C)`.

**Proposition F (PROVED).** Every quantity determined by the orbits of the
elements of `C`, and every proved identity among the weights (`sum lambda=1`,
the exact fibre rows, the counter formulas, the transport ratios), takes the
same value under the hypothesis "`m` is rooted" and under "`m` is unrooted".
Hence no lower bound on `W(m)` or on `lambda(m)` is derivable from `C` and
those identities.

*Proof.* The orbit of `m` is disjoint from `C` by assumption, so changing the
fate of `m` (rooted with its actual word, or unrooted with weight zero) changes
no orbit of an element of `C` and no counter of such an element. The identities
hold for every potential `G lambda'` of every nonnegative leaf measure
`lambda'` (Proposition C of the Pascal note, Codex P5), so they are satisfied
in both worlds. A derivation that used only these data would therefore prove
the same lower bound in a world where `W(m)=0`. QED.

**Corollary F1 (bank readouts; PROVED).** Let `p_j=lambda(6j+3)`, let `B` be
the set of bank indices, `u=1-sum_(j in B) p_j` the uncertified mass (exact,
since `sum p=1` is proved), and `h_m(j)=4t/(1+t)^2`, `t=2^(m-j)`. The best
bank-derived bounds are `H_(m,d+1) >= c := sum_B p_j h^(d+1)` and
`H_(m,d) <= b := sum_B p_j h^d + u sup_(j not in B) h_m(j)^d`. For a target
`m` outside the bank the supremum is attained at `j=m` and equals 1, so

    9c - 8b = sum_B p_j h^d (9h-8) - 8u <= -8u < 0,

because `9h-8<=0` at every bank index (`h<=8/9` off the target). For a target
inside the bank the supremum is `h_m(j_max+1)`, exponentially small, and the
readout converges to the atom as in Codex's (4). QED.

S7 (VERIFIED numerical): bank = rooted leaves below `2^18`, certified mass
`0.621951`, `u=0.378049`; targets at indices `43691, 43704, 43717, 43730` read
negative at `d in {0,5,10,20,40}`; the leaf 27 reads
`-6.2, -1.1e-4, -1.2e-6, -7.2e-11, 8.87e-11, 8.87e-11` at
`d=0,10,20,40,80,120`, the last two equal to `lambda(27)`; the threshold
`8(16/25)^d < lambda(27)` is `d>=57`.

Consequence. "Avoiding certificate recycling" has exactly one admissible
meaning: produce a deadline from the arithmetic of the source. A deadline is
a bound `tau(n)<=B(n)`; by F1 it gives the floor `2/((B+1)(B+2))`, and by F2
it is found by at most `B` forward steps. Every known deadline is a descent
cylinder: forward cylinders where the first `A` parity bits force a drop
(Terras; THM-4512's least-representative test; the controllers and SF5), or
backward cones where a smaller ancestor exists (the shadow theorem). Each
covers a decidable set of positive density, none covers every integer, and the
union over all depths of the forward cylinders is exactly the complement of
`E_inf`.

## 3. Working backwards from PMB with unique certificates

Let `E_inf` be the set of 2-adic integers all of whose parity prefixes are
rising (`3^(k_j) > 2^j` for every prefix length `j`, `k_j` ones). The rising
prefixes of length `A` number `N(A)`, computed by a two-dimensional lattice
count; `log_2 N(A)/A` is `0.7869, 0.8496, 0.8891` at `A=30,60,120` and tends
to `h*=0.9500` (the polynomial ballot factor is large). Each rising prefix is
one residue class of odd integers modulo `2^A`.

**Equivalences (PROVED).** (a) PMB gives, by F1, a computable `B(n)` with
`tau(n)<=B(n)`; conversely such a `B` gives the floor `2/((B+1)(B+2))`, which
is independently justified by the orbit replay of length `B`. (b) A computable
deadline exists for every `n` iff every positive odd integer leaves every
rising cone at a finite depth, i.e. `E_inf` contains no positive integer:
the integers in `E_inf` are exactly those whose orbit never drops below the
start, the divergent orbits and nontrivial cycles. (c) Unique certificates are
automatic: the ROOT word is the orbit.

**Where the gap sits (FINITE-EXACT profile).** The number of odd `n<X` with
first-descent time `sigma_T(n)>A` equals `N(A) X/2^A` whenever `A<=log_2 X`,
up to least representatives, because the event depends only on `n mod 2^A`.
The census below `2^24`:

| `A` | census `#{sigma_T>A}` | `N(A) 2^(24-A)` | ratio |
|---:|---:|---:|---:|
| 4 | 3,145,728 | 3,145,727.6 | 1.0000 |
| 12 | 925,696 | 925,695.9 | 1.0000 |
| 23 | 337,614 | 337,614.0 | 1.0000 |
| 30 | 199,643 | 199,551 | 1.0005 |
| 60 | 31,749 | 32,249 | 0.9845 |
| 100 | 3,845 | 4,004 | 0.9602 |

The maximum first-descent time below `2^24` is 287 T-steps. The statement to
be proved is that the census column reaches zero at a finite `A_0(X)` for
every `X`; the heuristic `N(A) X/2^A < 1` gives `A_0 ~ log_2 X/(1-h*)`, about
`20 log_2 X`, while the records grow like a smaller multiple. Between
`A=log_2 X` and `A_0` the integers in the cones are the least representatives
of their classes, and nothing but their archimedean size distinguishes them
from the 2-adic points of `E_inf` (THM-4027/4026 "archimedean alignment", in
the Sun thread, is the same shape). The previous note's T1 says that this
2-adic tower is, after the prefix price, the 3-adic residue tower of the
images; the 3-adic shadows of the same rising words (the descent set `D`) are
where the induction can be fed from below, and section 4 shows how much of the
mass they carry.

**What a floor function must look like.** Any `epsilon(n)` satisfying PMB
must (i) depend on unboundedly many digits of `n` (the finite-lookahead
theorem, in its specified-observer scope), (ii) satisfy
`epsilon(2*3^l-1)` consistent with `W(2*3^l-1)>=W(2^(l+1)-1)` for every `l`
(Lemma B), and (iii) be at least as large as the records along their sequence; no bound of
the form `tau = O(log n)` is proved anywhere in the thread, and none is
assumed here.
The Codex SF5 floor `w(3b,9b+1)` has all three properties on its tree; the
problem is its domain. No floor on all integers is proposed here.

## 4. Owed items

**Shadow theorem, exhaustive audit (FINITE-EXACT; scope: rising cones only).** For all 458 rising words
of length at most 7: every odd `m` below `4*2^(A+1)` follows `w` with exact
valuations iff `m = x_w mod 2^(A+1)` (2,281,464 checks), and every odd `n`
below `6*3^l` has an integral backward chain along `w` iff `n = x_w mod 3^l`
(2,368,206 checks), with the chain positive, returning to `m<n`. Inside a
rising cone there are no exceptions: integrality is exactly the class
condition and positivity is automatic since `2^a u-1>=1`. The class formula
`c_w 2^(-A) mod 3^l` agrees with `c_w (2^A-3^l)^(-1) mod 3^l`. This audit
does not address descents through non-rising words: for a word with
`2^A>3^l` the chain endpoint is below `n` exactly when `(2^A-3^l) n < c_w`,
a finite initial segment of the class, and such descents exist (165 -> 167).
Whether they ever occur outside every rising cone is open; the census of
section 4 counts all smaller ancestors, the cone series counts rising cones.

**Segment ledger, direct audit (FINITE-EXACT).** For the stopping and
single-rise rules with starts below `2^16`, the incoming sum at every odd
target below `4096`, computed by enumerating the recorded predecessors
directly (52,062 and 43,697 of them), equals the through mass plus the exit
mass, and the defect equals `nu(y)` minus the exit mass.

**Lemma T4b, both directions (PROVED).** For `l>=2` a primitive rising word of
length `l` exists iff `{l log_2 3} < log_2(3/2)`; depth 1 is the separate base
case `(1)` (its fractional part equals the threshold). In integer form the
test is `1 + ceil((l-1) log_2 3) <= floor(l log_2 3)`, which covers `l=1` as
well. Necessity was proved before. For
sufficiency take `a_1=1` and, for `j=1,...,l-1`,
`a_(l-j+1) = ceil(j log_2 3) - ceil((j-1) log_2 3)`, so that the suffix of
length `j` has sum exactly `ceil(j log_2 3)` and is not rising, while the
total `1+ceil((l-1) log_2 3) <= floor(l log_2 3)` makes the word rising. The
letters are 1 or 2, and `ceil(j log_2 3)` is the bit length of `3^j`. S3
builds the word at every admissible depth up to 60; the admissible depths
below 60 are `1,2,4,6,7,9,11,12,14,16,18,19,21,23,24,26,28,30,31,33,35,36,38,
40,42,43,45,47,48,50,52,53,55,57,59,60`, a Beatty set. S4 checks that the new
cone classes vanish exactly at the inadmissible depths through 16.

**Cone series and census (finite proportions, not proved densities).** New cone classes per depth through 16:
`1,1,0,1,0,2,8,0,28,0,124,602,0,2498,0,12319` (2,999,301 rising words
visited); the cone density series reaches `0.46726` at depth 16; the sieve
density of `D` below `2^28` is `0.468669` (units `0.703004`), constant to five
digits across the last four dyadic blocks. The remaining `0.0014` is carried by
deeper rising cones and, possibly, by non-rising descents outside every cone.

**Mass split below `2^22` (W head 2.850934).** Root 1.000000 (35.1%), leaves
0.656784 (23.0%), descent set 1.016542 (35.7%), basin minima 0.177608 (6.2%).
Basin minima are 415,230 of 1,398,101 units (29.7%) but carry only 6.2% of the
head mass: they are mass-poor because their atoms are sums over leaves above
them, all larger. This answers obligation (iii) of the previous note at head
level: the part of the weight that no induction from below can reach is small
in mass and large in count.

## 5. What remains

1. The counting statement of section 3 at the first unresolved scale: for a
   fixed `X` (say `2^24`), the census column is known to reach zero at
   `A=288`; a proof that it reaches zero for every `X` is the conjecture. The
   honest intermediate target is the growth of the record `A_max(X)` against
   `log_2 X/(1-h*)`; THM-4476/4499 (thin divergence, `o(X^(h*))`) are the
   current upper bounds on the census column at large `A`.
2. The exact value of `dens(D)`: the series over admissible Beatty depths
   converges; its terms at depths 14 and 16 (`2498/3^14`, `12319/3^16`)
   suggest a geometric tail, which would give a closed interval.
3. Codex SF4's energy `E_m` is the quantity a non-recycling argument would
   have to bound; Proposition F says it cannot be bounded from a bank. A
   source-only argument must bound `sum (1/c_(i+1)-1/c_i)` along the 3-adic
   tower of the images (T1), which is the deadline in disguise.

## 6. Reproduction and scope

[Script](../../04-computation/experiments/collatz_measurement_independence_20261005.py),
[output](collatz_measurement_independence_20261005.out),
[JSON](collatz_measurement_independence_20261005.json):

```text
python3 04-computation/experiments/collatz_measurement_independence_20261005.py --cone-depth 16 --sieve-bits 28 --head-bits 22 --json 05-knowledge/results/collatz_measurement_independence_20261005.json
python3 -O 04-computation/experiments/collatz_measurement_independence_20261005.py --cone-depth 10 --sieve-bits 20 --head-bits 18 --audit-len 5
```

4,658,591 explicit checks in 11 s (numba for the cone DFS, the sieves, the
counters and the first-descent census; exact integers and Fractions for the
audits). Hostiles: the inadmissible Beatty depths must give zero new classes;
the out-of-bank readouts must be negative at every tested degree; the census
and the prefix count must agree to four digits for `A<=23`. Limits: the
readouts are float64; the census is exact in its range; no finite statement is
read as an unbounded one, and nothing here certifies any new source.
