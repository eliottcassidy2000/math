"""Re-target the audit edits to the current HEAD 1d167d914c (copies in latest_1d167/): drop the edits superseded by the
producer's rewrite of THM-4811 (C)/(E) and note sections 1-2, re-anchor three THM-4811 edits, update the reflection row,
validate (each old string unique, all edits compose) and write required_edits.json."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
old = json.load(open(os.path.join(HERE, 'required_edits_vs_eca666b44c.json')))
F11 = "01-canon/theorems/THM-4811-arc-strong-connectivity-is-the-exact-complement-distance-to-the-nonstrong-locus-black-refinement-and-the-dipole-law.md"
FR = "07-reflections/why-the-advances-work-ii-potentials-symmetry-budgets-positive-coordinates-macmini-20261027.md"
FN = "05-knowledge/results/dipole_calculus_black_distance_tc_gaps_collatz_numerator_20261027.md"
E = []
for e in old:
    if e['finding'] in ('m3', 'm5', 'm5-related'):
        continue                      # superseded: C2 proof text and member counts rewritten in 1d167d914c
    if e['finding'] == 'M4' and e['file'] == FN:
        continue                      # note table row fixed in 1d167d914c
    if e['file'] == F11 and e['finding'] in ('M2', 'm4', 'm6'):
        continue                      # re-anchored below
    E.append(e)
E.append({"file": F11,
  "old": "  matrix form K' = -K + 2A_P, klein's THM-4814 (D)",
  "new": "  matrix form K' = -K + 2A_P (also Session75 RESULT (6.1), on main and audited before this file, with the path spectrum\n  and the trace-norm tariff (6.2)-(6.3)), klein's THM-4814 (D)",
  "reason": "Session75 (d5f10e1cf2, 01:22) states K(tau_P T) = -K(T) + 2L as (6.1) before THM-4811's first commit (8ab6b3e358, 01:40).",
  "finding": "M2"})
E.append({"file": F11,
  "old": "Every black line at a non-strong class is such an edge (by (B) its other end has\n`lambda = 1`), and by (B)-(C) every non-pure-blue strong class of order at most 10 is joined to a `lambda = 1` class by\nblack moves (for every order, by Theorems A and B).",
  "new": "Every black line at a non-strong class is such an edge (its other end is strong by THM-4729 L3, so has\n`lambda = 1` by (B)), and every non-pure-blue strong class is joined by black moves to a `lambda = 1` class with a black\ndescent: at every order already by the inherited strong black `c3`-descent (black_descent_20261020c; the last strong\nclass of a black descent chain has `lambda = 1` and its last move is a black descent word), independently of Theorems A\nand B, and through order 10 also by (C1).",
  "reason": "(B) alone gives only lambda <= 1 at the other end; strongness of the other end is THM-4729 L3. The all-order reduction does not need the (pre-audit) classification: the inherited all-order black c3-descent suffices.",
  "finding": "m4"})
E.append({"file": F11,
  "old": "descent criterion of (B) verified on all legal words through order 8 (2,073,672 words at order 8, no failure);",
  "new": "descent criterion of (B) verified on all legal words of all strong classes through order 8 (2,073,672 words at order 8,\n  no failure; independent audit: every class, strong or not, through order 9, 2,143,304 words at order 8 and 270,447,284\n  at order 9);",
  "reason": "The producer's check ran on strong classes only (gentourng -c); 2,073,672 is the strong-class word count at order 8.",
  "finding": "m6"})
for e in E:
    if e['file'] == FR and e['finding'] == 'M4':
        e['new'] = "| same, black moves only (THM-4811 (C)) | `lambda` | 1 | a black word starting in a minimum out-cut and ending outside it | all but `C3`, `R5` and the family `E` of THM-4811 Theorem A (pre-audit); `E(M; Z)` is only its `t = 1` subfamily |"
        e['reason'] = "Stale row: it lists E(M; Z), which the current THM-4811 (Theorem A) and HYP-9421 (1) call an incomplete subfamily (E_2(S5) at order 11 is an exception outside it); in eca666b44c the same row stated the then-OPEN classification as fact."
    if e['file'] == FR and e['finding'] == 'm12':
        e['new'] = e['new'].replace("FINITE-EXACT only through order 9 (D2)", "FINITE-EXACT only through order 10 (D2)")
order = {"M1": 0, "M2": 1, "M3": 2, "M4": 3}
E.sort(key=lambda e: (order.get(e['finding'], 10 + int(e['finding'][1:])), e['file']))
D = os.path.join(HERE, 'latest_1d167')
texts, ok = {}, True
for e in E:
    texts.setdefault(e['file'], open(os.path.join(D, os.path.basename(e['file'])), encoding='utf-8').read())
    if texts[e['file']].count(e['old']) != 1:
        ok = False; print('BAD', e['finding'], e['old'][:80])
work = dict(texts)
for e in E:
    if work[e['file']].count(e['old']) != 1:
        ok = False; print('COMPOSE FAIL', e['finding'])
    work[e['file']] = work[e['file']].replace(e['old'], e['new'], 1)
print('valid against 1d167d914c:', ok, 'edits:', len(E))
json.dump(E, open(os.path.join(HERE, 'required_edits.json'), 'w'), indent=2, ensure_ascii=False)
