"""make_tables.py -- markdown tables for q7_report.md from the run/validation JSON files (no transcription by hand).

usage: python3 make_tables.py D RUNDIR [k1 k2 ...]     (default: every k with a validation_k<k>.json)
Prints a per-k table (leaves C / DP, resolving leaves, orbits, CPU, ns/leaf, chunks, validation status)
and a per-(k, a) table.
"""
import sys, os, json, glob, re

def main():
    D = int(sys.argv[1]); rundir = sys.argv[2]
    ks = [int(x) for x in sys.argv[3:]] or sorted(int(re.search(r'_k(\d+)\.json', p).group(1)) for p in glob.glob(os.path.join(rundir, 'validation_k*.json')))
    # wall time per k from the driver log: last 'start' to 'COMPLETE' (after a restart only the last session)
    import datetime
    wall = {}; st = {}
    lp = os.path.join(rundir, 'log.txt')
    if os.path.exists(lp):
        for line in open(lp):
            m = re.match(r'(\S+ \S+) (start|COMPLETE) D=%d K=(\d+)' % D, line)
            if m:
                t = datetime.datetime.strptime(m.group(1), '%Y-%m-%d %H:%M:%S'); kk = int(m.group(3))
                if m.group(2) == 'start': st[kk] = t
                elif kk in st: wall[kk] = (t - st[kk]).total_seconds()
    print('| k | leaves (C search) | leaves (DP) | equal | resolving leaves | orbits | CPU s | wall s | ns/leaf | validation |')
    print('|---|---:|---:|:-:|---:|---:|---:|---:|---:|:-:|')
    rows = []
    for k in ks:
        v = json.load(open(os.path.join(rundir, 'validation_k%d.json' % k)))
        s = json.load(open(os.path.join(rundir, 'summary_k%d.json' % k)))
        c, dp = v['leaves_total_C'], v['leaves_total_DP']
        ns = 1e9 * v['cpu_sec'] / c if c else 0.0
        print('| %d | %s | %s | %s | %d | %s | %.1f | %s | %s | %s |' % (k, format(c, ','), format(dp, ','), 'yes' if c == dp else 'NO',
              v['resolving_leaves'], v.get('orbits', 0), v['cpu_sec'], ('%.0f' % wall[k]) if k in wall else '-', ('%.2f' % ns) if v['cpu_sec'] >= 1.0 else '-', 'ok' if v['ok'] else 'FAILED %s' % v['problems']))
        rows.append((k, v, s))
    print()
    print('| k | a | b | delta | reps | chunks | leaves (C) | leaves (DP) | resolving | CPU s |')
    print('|---|---|---|---|---:|---:|---:|---:|---:|---:|')
    for k, v, s in rows:
        for a in sorted(v['per_a'], key=int):
            p = v['per_a'][a]; ai = int(a); b = k - ai
            if 'nrep' in p:
                print('| %d | %d | %d | %d | %s | %d | %s | %s | %d | %.1f |' % (k, ai, b, ai - b, format(p['nrep'], ','), p['chunks'],
                      format(p['leaves_C'], ','), format(p['leaves_DP'], ','), p['found'], p['sec']))
            else:
                print('| %d | %d | %d | %d | - | - | 0 (weight lemma: W_min=%d > %d) | 0 | 0 | 0 |' % (k, ai, b, ai - b, p['wmin'], p['bound']))

if __name__ == '__main__':
    main()
