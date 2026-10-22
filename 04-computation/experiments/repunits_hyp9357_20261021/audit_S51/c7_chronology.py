"""S51 audit c7: chronology of the two pre-registrations against the data files (read-only: git log/reflog, file
times, and the public GitHub push events if reachable).  Local timestamps are self-reported by one machine; the
GitHub events give only server-side upper bounds for when a commit was on the remote.
Prints ALL CHECKS PASSED or the failures."""
import os, subprocess, datetime, json, time

FAIL = []
WT = r'C:\Users\Eliott\AppData\Local\Temp\math-wt-s51-repunits'
SP = r'C:\Users\Eliott\AppData\Local\Temp\claude\C--Users-Eliott-Documents-GitHub-ephrepos-math\e1fb236f-ec78-4211-aa52-c187130b32b0\scratchpad\s51'
EXP = os.path.join(WT, '04-computation', 'experiments', 'repunits_hyp9357_20261021')
TZ = datetime.timezone(datetime.timedelta(hours=-6))


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


def git(*a):
    return subprocess.run(['git', '-C', WT] + list(a), capture_output=True, text=True).stdout.strip()


def ts(s):
    return datetime.datetime.fromisoformat(s)


def ftime(p):
    st = os.stat(p)
    birth = getattr(st, 'st_birthtime', st.st_ctime)          # Windows: st_ctime is the creation time
    return datetime.datetime.fromtimestamp(birth, TZ), datetime.datetime.fromtimestamp(st.st_mtime, TZ)


print('(G) commits (author date = when the commit object was first made; committer date = after pull --rebase)')
C = {}
for h in ('a52d1d335c', '07d115c74e', '917fef912a', 'eaf70543a9', 'be7c76115f', '07c46ff196', 'ce5e055a0a', 'b1c598b6fa', '4712463a7b'):
    out = git('show', '-s', '--format=%aI|%cI|%s', h)
    if not out:
        print(f'     {h}: not found'); continue
    a, c, s = out.split('|', 2)
    C[h] = (ts(a), ts(c))
    print(f'     {h}: author {a}  committer {c}  {s[:90]}')
check(C['07d115c74e'][0] == C['917fef912a'][0] and C['eaf70543a9'][0] == C['be7c76115f'][0],
      'the rebased pre-registration commits keep the author dates of the originals (07d115c74e -> 917fef912a, eaf70543a9 -> be7c76115f)')
print('(F) data files (birth = creation time; NTFS tunnelling preserves it across git rebase rewrites)')
F = {}
for name in ('rep1.py', 'census_H_n9.pkl', 'inv.py', 'inv_1100.pkl', 'inv_2100.log', 'inv_2100.pkl', 'regress_train.py', 'evaluate_R.py',
             'evaluate_R7.py', 'N_seq_2100.pkl', 'decomp2.py', 'decomp2_1100.pkl', 'autcheck.py', 'aut_1100.pkl'):
    p = os.path.join(SP, name)
    if os.path.exists(p):
        F[name] = ftime(p)
        print(f'     scratch {name:18s} birth {F[name][0]:%H:%M:%S.%f}  mtime {F[name][1]:%H:%M:%S.%f}')
for name in ('repunit_rarity.py', 'repunit_rarity.out', 'two_arc_formulas.out', 'PREREGISTRATION_out_of_sample.md'):
    F[name] = ftime(os.path.join(EXP, name))
    print(f'     committed {name:32s} birth {F[name][0]:%H:%M:%S.%f}  mtime {F[name][1]:%H:%M:%S.%f}')
pre1 = C['07d115c74e'][0]; pre2 = C['eaf70543a9'][0]
print(f'     R1-R7 commit object created {pre1:%H:%M:%S}; inventory-2100 log created {F["inv_2100.log"][0]:%H:%M:%S} (run start); '
      f'pickle created {F["inv_2100.pkl"][0]:%H:%M:%S.%f}')
check(F['inv_2100.log'][0] < pre1 < F['inv_2100.pkl'][0], 'the 2100 inventory was launched before the R1-R7 commit and its result file appeared after it')
gap = (F['inv_2100.pkl'][0] - pre1).total_seconds()
print(f'     margin: the result file appeared {gap:.1f} s after the R1-R7 commit object (and {(F["inv_2100.pkl"][0] - C["917fef912a"][1]).total_seconds():.1f} s after its rebase)')
check(F['evaluate_R.py'][0] > F['inv_2100.pkl'][0] and F['N_seq_2100.pkl'][0] > F['inv_2100.pkl'][0], 'the R1-R7 evaluation scripts and N(n,h) tables were created after the inventory finished')
check(F['repunit_rarity.py'][0] > pre2 and F['repunit_rarity.out'][0] > pre2, 'repunit_rarity.py (first test-range decomposition) was created after the G1-G4 commit object')
check(F['decomp2_1100.pkl'][0] < pre2 and F['aut_1100.pkl'][0] < pre2, 'the training-range decomposition files (cap 1100) predate G1-G4, as stated')
txt = open(os.path.join(SP, 'rep1.py'), encoding='utf-8').read()
leak = '0.7 * top' in txt and 'for n in (8, 9)' in txt and "'3 | h'" in txt
check(leak and F['rep1.py'][1] < pre1, 'caveat: rep1.py (mtime before R1-R7) correlated the n = 9 class counts N(9,h) = S_9(h) (h > 661) over h <= 0.7 max H = 2349 with 3 | h, 5 | h, primality: the order-9 slice of the test range had been examined')
print('(S) server-side: public GitHub push events (upper bounds only)')
try:
    ev = []
    for pg in (1, 2, 3):
        r = subprocess.run(['gh', 'api', f'repos/eliottcassidy2000/math/events?per_page=100&page={pg}'], capture_output=True, text=True, timeout=60)
        ev += json.loads(r.stdout)
    pushes = [(e['created_at'], e['payload'].get('before', '')[:10], e['payload'].get('head', '')[:10]) for e in ev if e['type'] == 'PushEvent']
    rel = [p for p in pushes if '917fef912a' in p[1:] or 'be7c76115f' in p[1:] or 'ce5e055a0a' in p[1:]]
    for p in sorted(rel):
        print(f'     push at {p[0]} UTC: before {p[1]} head {p[2]}')
    first = min((p[0] for p in pushes if p[1] == '917fef912a'), default=None)
    if first:
        t = datetime.datetime.fromisoformat(first.replace('Z', '+00:00')).astimezone(TZ)
        print(f'     => 917fef912a was on the remote no later than {t:%H:%M:%S} local, i.e. {(t - F["inv_2100.pkl"][0]).total_seconds():.0f} s AFTER the 2100 inventory existed:'
              ' the order prereg < data is supported only by local timestamps')
except Exception as e:
    print('     (GitHub events not reachable:', e, ')')
print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
