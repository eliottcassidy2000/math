# The first-untested-point law: a maturation retrospective

**kind-pasteur-2026-07-27-S135.** Provenance note, not truth source.
Companion card: META-PATTERNS "Interpolant death / mechanism
survival". Harvest: 24 break instances + 6 survivor programs,
gathered by a dedicated provenance sweep across MISTAKES.md, the
hypothesis indexes, and canon status lines.

## 1. The empirical law

Across every thread of this repo, conjectured closed forms fitted
to k data points and then extended have broken **at the very first
untested point** in the overwhelming majority of recorded cases:

```text
fit                              points  broke at   converted?
Leonardo strong floor            5       m=8        YES (Busch/THM-1370)
minH = m^2-5m+9 (its ancestor)   4       m=7        YES (same chain)
pure-blue interleave (HYP-4997)  6       n=9        YES (THM-2444)
its '+1 odd' repair              8       n=11       YES (THM-2454)
transitive-partner 2^{n-2}+1     3       n=7        partial (THM-646)
7*3^k gap tower                  2-3     63 at n=8  YES (THM-1370)
self-loops-only-on-mixed         2       n=6        YES (order parameter)
2 selfK = SC                     3       n=8        YES (THM-849)
kernel proportionality (C5)      2       n=7        YES (THM-646 iv)
A014574 immediate recurrence     ~19     term 21    YES (THM-2422)
blue self-loop 2^{n/2-2}         3       n=10       NO (honest open)
mod-9 bicycle fit                8(pts)  n=21       YES (mod-12, THM-2467)
diag-sum order-3 recurrence      3(!)    m=8        rule: need 2r+2 terms
diam(G_n)=n-2                    3-4     n=7        YES (FAS reading)
A215257 impersonation            5       k=6        YES (MISTAKE-062)
Walsh 2(n-2k)^k                  many*   n=8        YES (THM-848 EGF)
dim_nonspec = n-5                heavy*  n=9        YES (3-step chain)
Mersenne-controls-H              1       immediate  YES ({7,21} core)
c9/c7 ratio constancy            1       k=9 (50x!) YES (LEM-001)
Paley H = |Aut| 3^{(p-3)/2}      3       p=11       partial (THM-586)
tc(T_p) = 3^{(p-3)/2}            3       p=11       YES (THM-2454 note)
Carmichael tc                    1       p=19       YES (HYP-9028 kill)
LRC pairwise cap C(14-j,2)/91    3       j=5        YES (THM-576!)
Landau first bite                1       any        YES (double-zero form)
```

(*Walsh and dim_nonspec broke at the first point OUTSIDE an
accidental exactness boundary -- the subtler version of the same
event.) Notable outliers: the A014574 rule survived ~19 points
before term 21 (it was a *density* statement, not a closed form),
and the mod-9 bicycle fit survived 8 residue classes of data before
n=21 -- both still died on extension.

## 2. Why: selection, not bad luck

The fits that die are **minimal interpolants**: the shortest
formula consistent with the data. Minimal interpolants are exactly
the objects with a large family of equally short competitors that
agree on the fitted range and diverge immediately after it -- so
conditional on "we adopted the shortest fit", breaking at the first
untested point is the generic outcome, not misfortune. The repo has
independently converged on defensive rules that all say this:
MISTAKE-055 ("check Busch before refitting"), MISTAKE-198 ("an
order-r recurrence needs >= 2r+2 terms"), the S420 museum of
impersonations, THM-1370's "fifth-term-break motif" phrase.

## 3. The discriminator: what survives

Six programs made genuine out-of-sample predictions and survived,
several with pre-registered exact hits:

```text
survivor                        hits                mechanism lemma
alphabet law (THM-2453/2454)    3 pre-registered    center law + stack mult.
twin harmonics (THM-2447)       12 amplitudes       local factor lemma
D-graded gate tower (THM-1286+) 8 both-directional  composite-skip + cascade
leg/blue-parity law (THM-790)   n=8 full inventory  leg identity d_v+a_v=0
tower grammar (THM-2450)        8 census points     stack mult. + T5 control
LEM-035 survivor law            M=11 by hand        section-vector formula
```

The discriminator is sharp: **every survivor carried a proved
sub-lemma that generated the prediction; every casualty was a
pattern with no mechanism.** The two regimes are visibly different
BEFORE extension: a mechanism conjecture predicts structure
(inventories, class lists, amplitudes), a fit predicts only the
next number. THM-2454's n=13 prediction specified thirteen
(H, |Aut|, tc) triples; the Leonardo fit specified one integer.

Corollary for practice, stated as the card's trigger: when you
catch yourself excited by a formula that predicts only the next
number, either (a) find the sub-lemma that would make it predict
structure, or (b) spend the next hour trying to kill it at the
next point -- those are the only two moves with positive expected
value. The repo's own recent sessions ran this loop deliberately
three times in 48 hours (mod-9 bicycle, tc(T_p), Carmichael-tc:
all killed within the hour; alphabet law: mechanized, then three
straight hits).

## 4. Recurrences confirmed, welded, and honestly not found

- `91 = 7*13 = C(14,2)`: recorded four independent ways (THM-576's
  pair-avoidance caps, HYP-3092, HYP-2938's unital row, HYP-2445's
  PSL2(F13) bridge) -- but the `Z/91Z` stalk theorem family
  (THM-2430/2436 line) never states that its modulus equals the
  runner-pair count. One sentence would weld the analytic and
  arithmetic 91s; flagged for the next stalk-family session.
- `1001 = 7*11*13`: three threads; the single 1001-91 co-location
  is HYP-2938's `q=10: v=1001, h=91` row. Undeveloped.
- `{7,21}` vs `{12,24}`: death-star-S59w's trisection-vs-doubling
  reading, which honestly labels itself constructed-not-found. Its
  open question ("is the trisection axis a full ladder m*{1,3}?")
  remains open; my strong-floor work adds that 7 and 21 sit in
  Busch CREVASSES -- the {12,24} side has no crevasse reading yet.
- mod-12 constellation: the new bicycle law joins AP-12 protection
  (HYP-4356, Lean-formalized), the AP-completion killers = 0 mod 12
  (HYP-4152/4157), and the Erdos-Straus hard classes 1 mod 12.
  Different moduli mechanisms (lcm(4,6) recurrence periods vs
  1/12-time protection vs quadratic residues); no common cause
  claimed -- recorded to PREVENT a future numerological weld.
- The carry-family d=4 rung (A003269, parts {1,4}): the SEQUENCE is
  recorded five ways; the repo-native combinatorial OBJECT is not
  (OPEN-Q-098). My probe adds a falsifiable marker: if a d=4 world
  exists in a census, its palindromic counts should read
  pal14 = 1,2,1,2,1,3,2,4,3,6,4,8,5,11,... No current census
  matches; the natural candidate alphabet {1, S4} (the unique
  strong 4-tournament) fails purity because tc(S4) = 5 -- so a
  d=4 world, if it exists, lives under a DIFFERENT involution than
  sigma. That is a concrete hunting license, not a claim.

## 5. Stopping note

The meta-law itself was tested against its own standard: it is not
a fit (Section 2 gives the mechanism -- competitor divergence under
minimal-description selection), and it made one immediate
prediction (the mod-9 bicycle fit would die on extension) that was
confirmed the same session it was formulated. Counterindications
recorded on the card: laws with proved local mechanisms (THM-2447)
survive class-mean tests only AFTER detrending -- a true law can
LOOK broken if the carrier is left in (the dual failure mode).
