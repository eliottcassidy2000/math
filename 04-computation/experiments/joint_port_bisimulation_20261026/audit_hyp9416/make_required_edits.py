# Builds required_edits.json for the HYP-9416 audit.  Every "old" string is taken verbatim from the current file (as a
# full line containing a unique anchor) and is checked to occur exactly once; "new" is the audit's replacement.
import json, os
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..'))
HYP = '05-knowledge/hypotheses/HYP-9416-fresh-re-marking-two-port-quotient-joint-classes-are-minimal-and-line-responses-classify-strong-atoms.md'
NOTE = '05-knowledge/results/fresh_remarking_quotient_and_divisor_arrangements_20261026.md'
RIDX = '05-knowledge/results/INDEX.md'
HIDX = '05-knowledge/hypotheses/INDEX.md'
ALQ = '04-computation/experiments/cuboid_simplex_aliquot_20261026/aliquot_preimages.py'
AUD = '04-computation/experiments/joint_port_bisimulation_20261026/audit_hyp9416'
text = {f: open(os.path.join(ROOT, f), encoding='utf-8').read() for f in (HYP, NOTE, RIDX, HIDX, ALQ)}
edits = []

def line_with(f, anchor):
    lines = [l for l in text[f].split('\n') if anchor in l]
    assert len(lines) == 1, (f, anchor, len(lines))
    return lines[0]

def sub(f, anchor, old_part, new_part, why):
    """replace old_part inside the unique line containing anchor (old = that line, new = line with the replacement)"""
    line = line_with(f, anchor)
    assert line.count(old_part) == 1, (f, old_part)
    if text[f].count(old_part) == 1:          # the replaced part alone is unique in the file: use it as the old string
        edits.append({'file': f, 'old': old_part, 'new': new_part, 'why': why})
    else:
        edits.append({'file': f, 'old': line, 'new': line.replace(old_part, new_part), 'why': why})

def whole(f, anchor, new, why):
    edits.append({'file': f, 'old': line_with(f, anchor), 'new': new, 'why': why})

# ---------------- HYP-9416 ----------------
whole(HYP, 'title: "Fresh re-marking two-port quotient.',
      'title: "Fresh re-marking two-port quotient. Rooted two-port states [C; s, t] evolve by fresh (unnamed, multiset-counted) '
      're-marking along Hamilton paths from s to t followed by complement (V1), or along all Hamilton paths with ports reset (V2). '
      '(a) PROVED (audit reconstruction lemma): already the single sandwich e*h*e determines every strong joint two-port class for '
      'n >= 3, so the static observer of the two-sided gluings k*h*k\' with contexts of order <= 2 trivially admits no closed '
      'compression on strong states (scalar two-sided observers such as H(e*h*e) do compress). (b) One-sided gluings depend only on '
      'the one-sided classes (PROVED) and never separate the one-sided collision pairs; one V1 step, whose successor multiset is a '
      'joint-port statistic, separates all of them even from the bare class observation (FINITE-EXACT n <= 7); V2 separates none. '
      '(c) The port-forgotten line-response histogram of a strong atom determines the atom (FINITE-EXACT n <= 7; even the set of '
      'partner classes does). (d) If only the unrooted class must be determined at every step, the V1 quotient is a proper '
      'compression (FINITE-EXACT n = 4..7); the V2 quotient is the class partition (PROVED, trivial). This is an analogue, on '
      'two-port classes and with counted multisets, of the stochastic quotient mentioned in vitali-response C5, not that quotient."',
      'Statements (a) and part of (b), (d) are one-line theorems (audit); (b) is a joint-port effect of V1 by construction, not an '
      'upgrade of one-sided data; the identification with the C5 quotient is inexact (different carrier X_n vs two-port classes, '
      'counts vs distributions: the normalized-distribution quotient is coarser, 11031 vs 11033 strong cells for free at n = 7).')
whole(HYP, '  OPEN for all orders. FINITE-EXACT, exhaustively',
      '  (1) PROVED for every n >= 3 by the independent audit (reconstruction lemma: e*h*e alone determines a strong two-port state);',
      'Status must reflect that (1) is a theorem and which parts are trivial; the audit reproduced every finite number.')
whole(HYP, '  at n = 7). Validated against the published 289/253 hostile. Not yet independently audited.',
      '  the one-sidedness of glue1 and the V2/free statement are PROVED (trivial). (2) closure part, (3) and the V1 part of (4):\n'
      '  FINITE-EXACT, exhaustively over all tournaments, for n <= 7 (17824 rooted states and 14162 strong at n = 7), all numbers\n'
      '  independently reproduced; all-order versions OPEN. Validated against the published 289/253 hostile. Independently audited\n'
      f'  2026-10-26: {AUD}/.',
      'Second half of the status block (the two status lines are replaced by these four).')
sub(HYP, '- **Quotient.** For an observation `O`',
    'This is the bisimulation (operation-kernel) quotient of the vitali-response C5, with multisets in place of named guarded operations.',
    'It has the form of the operation-kernel quotient of vitali-response C5, but on a different carrier (two-port classes of all '
    'tournaments, not the presentations X_n), with V1 restricted to the current ports and with successor multisets counted rather '
    'than normalized: the normalized-distribution quotient is slightly coarser (V1 refining `free`: 11031 instead of 11033 strong '
    'cells at n = 7; refining `sand`: 5888 instead of 5889).',
    'Audit computed the distribution version (quotients.py, V1dist): it differs at n = 7, so the multiset quotient is not "the" stochastic quotient.')
whole(HYP, '1. **Joint classes are minimal (two-sided).**',
      '1. **Joint classes are minimal for the two-sided gluing observer: PROVED, trivially.** For every `n ≥ 3` the single two-sided '
      'gluing `e*h*e` determines a strong two-port state. In `e*h*e` the added entrance vertex `a` has out-degree 1 and the added exit '
      'vertex `b` has in-degree 1. Since `C` is strong, the only other vertex of out-degree 1 can be `s` (when `d⁺_C(s) = 1`), whose '
      'out-neighbour has out-degree at least 2, and the only other vertex of in-degree 1 can be `t`, whose in-neighbour has in-degree '
      'at least 2. As `a → s` and `t → b`, the vertices `a` and `b` are determined, and `s = N⁺(a)`, `t = N⁻(b)`, `C = e*h*e − {a, b}` '
      '(audit `sandwich_reconstruction.py`). Hence `glue2` separates any two distinct strong joint two-port classes in every order '
      '(through `n = 7` each of its nine contexts does so alone), and no closed compression refining it exists. This concerns this '
      'observer only: the two-sided scalar observer `sand = H(e*h*e)` does compress (5889 strong V1 cells of 14162 at `n = 7`).',
      'The former statement ("for every n >= 5", OPEN) is a one-line theorem valid from n = 3; the computational fact adds nothing beyond it.')
whole(HYP, '2. **Closure upgrades one-sided data.**', '2. **Fixed-port re-marking separates one-sided collisions.**',
      'V1 sees the joint ports by construction; nothing is "upgraded".')
whole(HYP, '   - `glue1` depends only on one-sided classes, so it separates none of these pairs.',
      '   - `glue1` depends only on one-sided classes (PROVED: an isomorphism of `(C, t)` extends by the identity to `h*k`, one of `(C, s)` to `k*h`), so it separates none of these pairs.',
      'Trivial lemma; verified on all states n <= 7 by the audit.')
whole(HYP, '   - The V1 quotient refining `glue1` separates all of them.',
      '   - The V1 quotient refining `glue1` separates all of them (FINITE-EXACT, `n ≤ 7`). This is not specific to `glue1`: the V1 '
      'quotient refining the bare class observation `free` separates every pair at its first refinement round, i.e. by the multiset of '
      'unrooted partner classes over the Hamilton paths from `s` to `t`, a joint-port statistic by construction. The number of `s → t` '
      'Hamilton paths alone (the size of that multiset) separates 1, 4, 16, 65 of the 1, 7, 20, 70 pairs at `n = 3, 5, 6, 7`. The V2 '
      'quotient separates none of them. One-sided observers that see only one side do not suffice: V1 refining only the free class of '
      '`h*e` leaves 12371 strong cells of 14162 at `n = 7`.',
      'Audit (quotients.py, collision_anatomy.py, extra_observers.py): the separation is a one-step joint statistic of V1, not a closure phenomenon.')
whole(HYP, '4. **Compression for the class observable.**',
      '4. **Compression for the class observable.** The V1 quotient refining `free` is a proper compression on strong states for '
      '`n ≥ 4` (5 strong cells of 12 at `n = 4`; FINITE-EXACT). The V2 quotient refining `free` is the class partition itself '
      '(PROVED, trivially: V2 successor multisets do not depend on the ports).',
      'Threshold is n >= 4 in the data; the V2 statement holds by construction in every order.')
sub(HYP, '| 3 | 8 | 2 | 4 (2) | 5 | 8 (2) |', '| 8 (2) | – | 2 |', '| 8 (2) | 2 | 2 |',
    'glue2 separates both strong states at n = 3 (indeed all 8 states).')
sub(HYP, '| 4 | 32 | 12 | 11 (5) | 14 | 30 (12) |', '| 30 (12) | – | 4 |', '| 30 (12) | 12 | 4 |',
    'glue2 separates all 12 strong states at n = 4 (indeed all 32 states).')
whole(HYP, '- at `n = 5, 6, 7` there are 7, 20, 70 pairs, in 2, 12, 39 fibres;',
      '- at `n = 3, 4, 5, 6, 7` there are 1, 0, 7, 20, 70 pairs, in 1, 0, 2, 12, 39 fibres (the `n = 3` pair is the cyclic triangle with `s → t` versus `t → s`);',
      'The n = 3 collision exists and is the smallest instance.')
whole(HYP, '- the singleton sandwich `H(e*h*e)` alone separates 6, 17, 64.',
      '- the singleton sandwich `H(e*h*e)` alone separates 6, 17, 64;\n'
      '- the number of `s → t` Hamilton paths alone separates 4, 16, 65 (together with `H(e*h*e)`: 7, 17, 68), and every pair splits at the first V1 refinement round of `free` (audit).',
      'Shows what actually separates the pairs.')
whole(HYP, '- there are 1, 1, 6, 35, 353 atoms, with 1, 5, 44, 547, 9985 distinct joint targets;',
      '- there are 1, 1, 6, 35, 353 atoms, with 1, 5, 44, 547, 9985 distinct joint targets (these counts, the blue counts and every coefficient of `F_A` and `B_A` are already in vitali-response `strong-atlas.json`; the audit matched all 396 sources coefficient by coefficient);',
      'Credit: the data are the vitali-response atlas; only the injectivity observation is new.')
whole(HYP, '- the histograms `F_A` are pairwise distinct in every order, also after forgetting the ports;',
      '- the histograms `F_A` are pairwise distinct in every order (vacuously at `n = 3, 4`), also after forgetting the ports; even the set of partner free classes, without multiplicities, separates all atoms through `n = 7`, while the multiset of partner score sequences gives 351 classes for the 353 atoms at `n = 7`;',
      'The full histogram is far from necessary for the classification (atoms.py).')
sub(HYP, '**Validation.** The rooted-lines R3 targets', 'and give `H(e*h*e) = 289` and `253`.',
    'and give `H(e*h*e) = 289` and `253`; they also have 4 and 2 Hamilton paths from entrance to exit, so their V1 successor multisets already differ in size.',
    'Same mechanism as item (2).')
sub(HYP, '- **vitali-response C5** (named guarded re-marking)',
    'That packet flags the fresh-re-marking quotient as a different, possibly smaller, question.',
    'That packet mentions a quotient under distributions of fresh uniform re-marking on its presentation carrier as a different, possibly smaller, question that must not stand in for the named protocol; the quotients here are a related construction on two-port classes (Session70 likewise records them as a different operation menu).',
    'Accurate description of what C5 flagged.')
whole(HYP, '  - For two-sided attachment observations there is no smaller quotient on strong states (1).',
      '  - For the observer `glue2` there is no smaller quotient on strong states (1), trivially, since `e*h*e` alone is a complete invariant of strong two-port states; the two-sided scalar observer `sand` compresses (5889 strong V1 cells of 14162 at `n = 7`).',
      'The former wording ("two-sided attachment observations") is contradicted by the V1/sand column of the table.')
whole(HYP, '  - Closure transfers joint-port information into one-sided observations (2).',
      '  - V1 (fixed-port) successor multisets are joint-port statistics, so one V1 step separates every one-sided collision pair (2); V2 adds nothing to one-sided data.',
      'Correct mechanism.')
whole(HYP, '  - only 74 partner classes are shared;',
      '  - only 74 partner classes are shared (already recorded in HYP-9414 (3), 74 between `T_000` and `T_100`); all 74 come from the phase-0 entrance family of the period-3 decoration;',
      'Credit HYP-9414 (3); audit refinement (control_audit.py).')
whole(HYP, "  - the period-3 decoration's three entrance phases carry equal path mass, 154305 each.",
      "  - the period-3 decoration's three entrance phases carry equal path mass, 154305 each. This is forced: A- = A+^op has the same ordered path-cover counts as A+ (fibre F4 of the tournament-composition packet), the T+ blocks carry the same masses, and the joint (entrance, exit) phase mass depends only on the phase difference (42930, 45036, 66339).",
      'The equal masses carry no information on the phase; credit the proved F4.')
sub(HYP, '  - This confirms HYP-9414 (1) at `m = 3`', 'This confirms HYP-9414 (1) at `m = 3`',
    'This re-checks HYP-9414 (1) and (3) at `m = 3`', 'Re-check, not new evidence.')
whole(HYP, '- All-order proofs of (1)–(3).', '- All-order proofs of (2) (closure part) and (3); (1) is proved (audit lemma).', '(1) is proved.')
whole(HYP, '- A structural reason why order-2 two-sided contexts suffice.',
      '- Answered: order-1 contexts suffice, by the degree reconstruction of `e*h*e`. Open: whether `A2(s→t)*h*A2(s→t)` (`k = k\'` the rooted arc entrance → exit), which is injective on all 8, 32, 176, 1408, 17824 states at `n = 3..7`, is a complete invariant of all two-port states in every order.',
      'Replaces an open question the audit answers with one it leaves open (single_context.out).')

# ---------------- results note ----------------
sub(NOTE, '2026-10-26, mac-mini-2026-10-25-pyth. Canon produced:', '(OPEN; FINITE-EXACT for n ≤ 7). Not yet independently audited.',
    f'(item (1) PROVED by the audit\'s reconstruction lemma; the rest FINITE-EXACT for n ≤ 7, all-order OPEN). Independently audited: `{AUD}/`.',
    'Status update.')
sub(NOTE, '**The open question.** vitali-response C5 answers', 'This section computes that quotient exhaustively for n ≤ 7.',
    'This section computes a related construction exhaustively for n ≤ 7: counted multisets (not distributions), on two-port classes of all tournaments rather than on the presentations X_n, with V1 restricted to the current ports. It is not literally the quotient C5 mentions; the normalized-distribution version is slightly coarser (11031 instead of 11033 strong cells for the class observable at n = 7).',
    'Inexact identification.')
whole(NOTE, '- **Two-sided attachments admit no compression on strong states.**',
      '- **The two-sided gluing observer admits no compression on strong states, trivially.** The nine two-sided gluings `k*h*k\'`, with contexts of order at most 2, separate all 2, 12, 92, 946 and 14162 strong joint classes at n = 3..7, and so does each of them alone. For `e*h*e` this holds in every order n ≥ 3: the added vertices are the out-degree-1 vertex `a` and the in-degree-1 vertex `b`; the only other candidates, `s` and `t`, satisfy `a → s` and `t → b` and are told apart by the degree of their neighbour, so `s`, `t` and `C` are read off (audit). A scalar two-sided observer such as `H(e*h*e)` does compress.',
      'Overgeneralized ("two-sided attachments") and trivial.')
whole(NOTE, '- **Closure turns one-sided observers into joint ones.**', '- **Fixed-port re-marking sees the joint ports directly.**', 'Correct mechanism.')
whole(NOTE, '  - A one-sided gluing `h*k` or `k*h` depends only on the exit-rooted or entrance-rooted class of h.',
      '  - A one-sided gluing `h*k` or `k*h` depends only on the exit-rooted or entrance-rooted class of h (an isomorphism of the rooted part extends by the identity). So no one-sided observation separates the 1, 0, 7, 20, 70 one-sided collision pairs at n = 3..7.',
      'Proof and the n = 3, 4 counts.')
whole(NOTE, '  - Closing under V1, i.e. repeated fresh re-marking and complement, separates every one of them.',
      '  - Closing under V1 separates every one of them, but no closure depth or attachment data is needed: a single V1 step from the bare class observation does it, because the V1 successor multiset is formed from the Hamilton paths from s to t and is a joint-port statistic by definition. The number of such paths alone separates 4, 16, 65 of the 7, 20, 70 pairs. Closing under V2 separates none.',
      'Audit: every pair splits at the first refinement round of V1/free (collision_anatomy.out).')
whole(NOTE, '  - This is a dynamic version of 289/253:',
      '  - So this is not a propagation phenomenon: the joint orbit is invisible to every one-sided observer and visible to any statistic built from s → t paths. The R3 pair itself has 4 and 2 Hamilton paths from entrance to exit.',
      'The former sentence misdescribes the mechanism.')
whole(NOTE, '  - So the "smaller stochastic quotient" exists for this observable, and only for observables that ignore attachment.',
      '  - So a smaller fresh quotient exists for this observable. It also exists for attachment observables that do not determine both one-sided classes: at n = 7, V1 refining `H(e*h*e)` has 5889 strong cells and V1 refining the free class of `h*e` alone has 12371; under V2, `glue1` stays at 14111.',
      '"only for observables that ignore attachment" is false (V1/sand in the lane\'s own table; extra_observers.out).')
whole(NOTE, '  - with fresh re-marking and any two-sided attachment observer, it is the joint two-port class, with no further compression;',
      '  - with fixed-port fresh re-marking (V1) and the observer `glue2` (indeed `e*h*e` alone), it is the joint two-port class, trivially, since that observer is a complete invariant of strong two-port states; weaker two-sided observers such as `H(e*h*e)` compress;',
      'False as stated ("any two-sided attachment observer").')
whole(NOTE, '  - with only one-sided attachment observers, closure forces the joint class anyway.',
      '  - with one-sided attachment observers that determine both one-sided classes (`glue1`, or just `h*e` and `e*h`), the V1 closure is again the joint class (FINITE-EXACT, n ≤ 7); with one side only, or under V2, it is not.',
      'False as a general statement (free class of h*e alone: 12371 of 14162; V2: 14111).')
whole(NOTE, '- **Beyond one run:** vitali-response C4',
      '- **Alignment memory:** vitali-response C4 shows that no constant number of extra bits repairs the paired-rooted observer for named guarded re-marking. Its family already uses a single global Hamilton path through k blocks (Session70), so it is not a statement about multi-run attachments only.',
      'Mischaracterizes the scope of C4 (also flagged by Session70).')
sub(NOTE, '**Classification.** In every order the histograms', 'So the line response determines the strong atom through order 7.',
    'So the line response determines the strong atom through order 7. Much less suffices: the set of partner free classes, without multiplicities, already separates all 1, 1, 6, 35, 353 atoms, and the multiset of partner score sequences gives 351 classes for the 353 atoms at order 7 (audit).',
    'Strength of the classification claim (atoms.out).')
whole(NOTE, '- Only 74 partner classes are shared.',
      '- Only 74 partner classes are shared, as already recorded in HYP-9414 (3); all 74 come from the phase-0 family below.',
      'Credit HYP-9414 (3).')
sub(NOTE, '**The Frobenius phase.** For the period-3 decoration', 'Injectivity makes the three families disjoint.',
    'Injectivity makes the three families disjoint. The equal masses are forced: A- = A+^op has the same ordered path-cover counts as A+, so every block-phase path count of T- equals that of T+ (fibre F4 of the tournament-composition packet); the T+ blocks also carry 154305 each, and the joint (entrance, exit) phase masses depend only on the phase difference (42930, 45036, 66339). The phase is invisible to path mass and visible only in the partner classes.',
    'Equal phase mass is a consequence of a proved result, not evidence about the phase.')
sub(NOTE, '**Why this pair is the right control.**', 'It therefore cannot predict a single partner. The datum must keep:',
    'It therefore cannot distinguish the two partner sets, which share only 74 classes. The data that suffice here (necessity is not proved) are:',
    'Overstatement ("cannot predict a single partner", "must").')
whole(NOTE, "vitali-response's HISTORY.md already threads the Vitali atom from March to October:",
      "vitali-response's HISTORY.md, with the October 11 Vitali recovery note it starts from, threads the Vitali atom from March to October:",
      'The quoted principle is in vitali_atom_recovery_20261011.md, not HISTORY.md.')
whole(NOTE, '- **THM-169 / THM-170:** the ±2 lambda switch',
      '- **THM-168 / THM-169 / THM-170:** the ±2 lambda switch; the labelled-lambda witness masks 82468 / 84707 with H = 141 vs 143 are in THM-168 and `collatz_effective_prefix_mass_20261005.md` §9 (HISTORY.md instead rechecks the twins 1354195 / 1351956 with H = 109 vs 111);',
      'Misattribution of the witness pair.')
whole(NOTE, '- **the recovery principle:**',
      '- **the recovery principle** (`vitali_atom_recovery_20261011.md`): "a target descends through an observation exactly when it is constant on each observation fibre";',
      'Source of the quote.')
sub(NOTE, '- **the endpoint, half-edge, joint-root, reflection, stabilizer and phase hostiles:**', '15/11, 11/13',
    '15/23 (vitali C2; the value 15/11 could not be traced in the cited packets), 11/13',
    'Untraceable value; the port-control witness in the cited packet is H = 15 vs 23.')
sub(NOTE, 'The one formal object behind all of these is the **operation-stable quotient**',
    'The one formal object behind all of these is', 'A common formal reading of all of these (HISTORY.md itself cautions that the historical uses of "Vitali atom" do not name one object) is',
    'Overclaim relative to the cited ledger.')
sub(NOTE, 'The new piece here is the **arrangement count dictionary**', 'The new piece here is the **arrangement count dictionary**',
    'The added piece here is the elementary **arrangement count dictionary** (classical region counts of generic affine and central arrangements, relabelled by divisors)',
    'Not a new result.')
whole(NOTE, '  - THM-1370: H ≠ 7, 21 at every order.',
      '  - THM-1370 (`THM-1370-h-spectrum-omits-7-21-all-n.md`; the number is also used by an unrelated Jacobian theorem): H ≠ 7, 21 at every order.',
      'Repository number collision; cite the file.')
whole(NOTE, '- **A second microcosm/macrocosm (new here).**',
      '- **A second microcosm/macrocosm (Goormaghtigh\'s 31 already appears in HYP-9362 and the Leech deep-hole note; new here is only its aliquot reading).**',
      'Repository prior art.')
sub(NOTE, '  - 31 = σ(2⁴) = σ(5²) is the only value', 'shared by two prime bases below 10¹² (p < 2000). This is Goormaghtigh\'s 31, so s⁻¹(31) contains 32 and 125.',
    'shared by two prime bases below 10¹² (complete: the audit checked every prime base p < 10⁶; the lane script only p < 2000). This is the classical Ratat–Goormaghtigh 31; Goormaghtigh\'s conjecture that 31 and 8191 are the only numbers that are repunits with at least three digits in two different bases is open, Bateman and Stemmler (with Horn) found 31 the only case below 1.275·10¹⁰ with both bases and the value prime, and Grantham (2024) found no new Goormaghtigh primes below 10⁷⁰⁰. So s⁻¹(31) contains 32 and 125.',
    'Unnecessary restriction and missing classical status.')
whole(NOTE, '- **The two sixes and abc (new remark).**',
      '- **The two sixes and abc (prior art: the equal-radical fibre 1+2=3, 1+8=9 with radical 6 is in THM-3833 (30), CORE-PAPERS-ABC-IUT (2) and abc_difference_families_20261011c; the Gersonides adjacencies {2, 4, 8} are in THM-4640).**',
      'Not a new remark in the repository.')
whole(NOTE, '  - 6 is the first n whose square lies in a later doubling row: 36 = 9·4.',
      '  - 6 is the first even n whose square lies in a later doubling row: 36 = 9·4. If odd n are allowed, 3 is first (9 is the leader of a doubling row after 3); aliquot124 chains §3 says exactly this.',
      'Wrong as worded; contradicts the general condition on the next line at n = 3 (a = 0).')
whole(NOTE, '  - The general condition is n = 2ᵃu with u^(k−1) > 2ᵃ.',
      '  - The general condition is n = 2ᵃu, u odd, with u^(k−1) > 2ᵃ (a ≥ 1 in the even case; Session69\'s map requires a ≥ 1).',
      'Parity qualifier.')
sub(NOTE, '  - All of this is PROVED in aliquot124 chains §2', 'aliquot124 chains §2', 'aliquot124 chains §§2–3',
    'The delay criterion is chains §3, eq. (7).')
whole(NOTE, '1. All-order proofs of HYP-9416 (1)–(3), and a structural reason why contexts of order 2 suffice.',
      '1. All-order proofs of HYP-9416 (2) and (3). Item (1) is proved by the audit\'s reconstruction lemma, which also explains why small contexts suffice: order 1 already does.',
      '(1) is proved.')

# ---------------- indexes ----------------
sub(RIDX, '[Fresh re-marking two-port quotients, strong-atom line responses, and divisor arrangements]',
    '**FINITE-EXACT (n <= 7); HYP-9416 OPEN; pre-audit at filing.** Companion to vitali-response C5 for fresh re-marking: two-sided order-2 gluings separate all strong joint two-port classes; closure turns one-sided observers into joint ones;',
    '**HYP-9416 (1) PROVED (audit lemma: e*h*e alone determines a strong two-port state); the rest FINITE-EXACT (n <= 7), all-order OPEN; independently audited (joint_port_bisimulation_20261026/audit_hyp9416).** Fresh-re-marking analogue (counted multisets on two-port classes) of the quotient vitali-response C5 mentions: two-sided gluings separate all strong joint two-port classes, trivially; one fixed-port re-marking step (a joint-port statistic) separates all one-sided collisions;',
    'Index must match the corrected statuses.')
sub(HIDX, '(2026-10-26; OPEN, FINITE-EXACT n <= 7, pre-audit at filing)',
    '(2026-10-26; OPEN, FINITE-EXACT n <= 7, pre-audit at filing) - fresh (unnamed) re-marking two-port quotient, companion to vitali-response C5: two-sided gluings of order <= 2 separate all strong joint classes (no closed compression); closure upgrades one-sided observers to joint;',
    f'(2026-10-26; (1) PROVED by the independent audit\'s reconstruction lemma, the rest FINITE-EXACT n <= 7 and OPEN in general; audited, {AUD}/) - fresh (unnamed) re-marking two-port quotient, an analogue on two-port classes of the quotient vitali-response C5 mentions: the sandwich e*h*e alone determines every strong joint class (trivially no closed compression); one fixed-port re-marking step, a joint-port statistic, separates all one-sided collisions;',
    'Index must match the corrected statuses.')
sub(HIDX, '- **THM-4780, HYP-9390, HYP-9391, HYP-9416 (2026-10-25, mac-mini-2026-10-25-pyth;',
    'HYP-9416 (2026-10-25, mac-mini-2026-10-25-pyth; PROVED / OPEN conditional on HYP-9276; independently audited, [audit]',
    'HYP-9416 (2026-10-25, mac-mini-2026-10-25-pyth; PROVED / OPEN conditional on HYP-9276 for THM-4780 and HYP-9390; HYP-9391 and HYP-9416 carry their own status below; THM-4780 and HYP-9390 independently audited, [audit]',
    'The combined header otherwise presents HYP-9416 as covered by the Pythagorean-hourglass audit.')

# ---------------- lane script ----------------
whole(ALQ, 'print(f"  equals {{p + 1 : (p-2, p) twin primes}}:',
      'print(f"  for v > 3 equals {{p + 1 : (p-2, p) twin primes}}: {[v for v in dh if v > 3] == [x for x in pred if x <= V]} (v = 3 also qualifies: 4 = 2^2 = 2*2)")',
      'The shipped aliquot_preimages.out prints a failing check (False) because v = 3 is a double hit via 4 = 2^2 = 2*2; regenerate aliquot_preimages.out with `aliquot_preimages.py 300` after this edit.')

# verify and write
for ed in edits:
    assert text[ed['file']].count(ed['old']) == 1, ed['old'][:80]
olds = [(ed['file'], ed['old']) for ed in edits]
assert len(set(olds)) == len(olds)
json.dump(edits, open(os.path.join(os.path.dirname(__file__), 'required_edits.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
# dry run: apply all edits in order to in-memory copies and check every old is still unique when reached
mem = dict(text)
for ed in edits:
    assert mem[ed['file']].count(ed['old']) == 1, ('sequential', ed['old'][:80])
    mem[ed['file']] = mem[ed['file']].replace(ed['old'], ed['new'])
print(f"{len(edits)} edits written; all old strings unique in the current files and applicable in sequence")
