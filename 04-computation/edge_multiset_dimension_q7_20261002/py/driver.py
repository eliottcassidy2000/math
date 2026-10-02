"""driver.py -- resumable, checkpointed parallel driver for src/edimsearch (normal form B search).

usage: python3 driver.py D K [--a a1,a2,...] [--workers W] [--target SEC] [--engine E] [--hot H] [--run DIR]

For every a = ceil(K/2)..K (or the given list) the representatives [0, nrep) of a-subsets of Q_{D-1}
are split into chunks; each chunk runs `nice -n 10 edimsearch D K a repfile first last ...`.
A completed chunk is appended (with fsync) to RUN/chunks.jsonl together with its SUMMARY fields and the
RESOLVING lines it printed; on restart, only the uncovered parts of [0, nrep) are scheduled again.
(K, a) pairs without a representative file are allowed only if the weight lemma excludes them:
an a-subset A with |beta_c(A)| >= a - 2b > 0 for all c has at most b minority entries per coordinate,
so after translating the majority point to 0 its weight sum is <= (D-1) b, while a distinct vertices of
Q_{D-1} have weight sum >= W_min(a) (lightest vertices first); W_min(a) > (D-1) b excludes the pair.
At the end RUN/summary_k<K>.json is written after checking that the chunks of each a tile [0, nrep).
"""
import sys, os, json, time, subprocess, threading, argparse, math
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
EXE = os.path.join(ROOT, 'src', 'edimsearch')

def repfile(n, a): return os.path.join(ROOT, 'data', 'q%d' % n, 'reps_n%d_a%d.bin' % (n, a))

def wmin(n, a):
    """minimum weight sum of a distinct vertices of Q_n"""
    s, r = 0, a
    for w in range(n + 1):
        t = min(r, comb(n, w)); s += t * w; r -= t
        if r == 0: return s
    return None

def lemma_excludes(D, K, a):
    n, b = D - 1, K - a
    t = a - 2 * b
    return t > 0 and wmin(n, a) > n * b

def parse_summary(line):
    f = {}
    for tok in line.split()[1:]:
        k, v = tok.split('=', 1)
        try: f[k] = int(v)
        except ValueError:
            try: f[k] = float(v)
            except ValueError: f[k] = v
    return f

class Run:
    def __init__(self, args):
        self.a = args; self.D, self.K = args.D, args.K
        self.dir = args.run or os.path.join(ROOT, 'runs', 'd%d_k%d' % (self.D, self.K))
        os.makedirs(self.dir, exist_ok=True)
        self.chunkfile = os.path.join(self.dir, 'chunks.jsonl')
        self.logfile = os.path.join(self.dir, 'log.txt')
        self.lock = threading.Lock()
        self.done = []
        if os.path.exists(self.chunkfile):
            for line in open(self.chunkfile):
                line = line.strip()
                if line:
                    try: rec = json.loads(line)
                    except json.JSONDecodeError: continue   # torn last line after a crash: chunk is redone
                    if rec.get('D') == self.D and rec.get('K') == self.K: self.done.append(rec)
        self.t0 = time.time()

    def log(self, msg):
        line = '%s %s' % (time.strftime('%Y-%m-%d %H:%M:%S'), msg)
        with self.lock:
            with open(self.logfile, 'a') as f: f.write(line + '\n')
        print(line, flush=True)

    def record(self, rec):
        with self.lock:
            rec = dict(rec, D=self.D, K=self.K)
            with open(self.chunkfile, 'a') as f:
                f.write(json.dumps(rec, sort_keys=True) + '\n'); f.flush(); os.fsync(f.fileno())
            self.done.append(rec)

    def covered(self, a):
        return sorted((r['first'], r['last']) for r in self.done if r.get('a') == a and 'first' in r)

def uncovered(nrep, cov):
    gaps, pos = [], 0
    for f, l in cov:
        if f > pos: gaps.append((pos, f))
        pos = max(pos, l)
    if pos < nrep: gaps.append((pos, nrep))
    return gaps

def run_chunk(run, a, first, last):
    n = run.D - 1
    cmd = ['nice', '-n', '10', EXE, str(run.D), str(run.K), str(a), repfile(n, a), str(first), str(last),
           '-e', str(run.a.engine), '-H', str(run.a.hot)]
    t = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    out = p.stdout.splitlines()
    summ = [l for l in out if l.startswith('SUMMARY')]
    res = [l for l in out if l.startswith('RESOLVING')]
    if p.returncode != 0 or len(summ) != 1 or any(l.startswith('FATAL') for l in out):
        raise RuntimeError('chunk failed a=%d [%d,%d): rc=%d\n%s\n%s' % (a, first, last, p.returncode, p.stdout[-2000:], p.stderr[-2000:]))
    f = parse_summary(summ[0])
    if f['reps'] != last - first or f['found'] != len(res) or f['first'] != first or f['last'] != last:
        raise RuntimeError('inconsistent chunk summary: %s' % summ[0])
    rec = dict(a=a, first=first, last=last, wall=round(time.time() - t, 3), resolving=res,
               **{k: f[k] for k in ('reps', 'feasible_reps', 'leaves', 'nodes', 'fullchecks', 'found', 'scan', 'sec', 'engine', 'hotcap')})
    run.record(rec)
    return rec

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('D', type=int); ap.add_argument('K', type=int)
    ap.add_argument('--a', default=''); ap.add_argument('--workers', type=int, default=2)
    ap.add_argument('--target', type=float, default=120.0, help='target seconds per chunk')
    ap.add_argument('--engine', type=int, default=0); ap.add_argument('--hot', type=int, default=1024)
    ap.add_argument('--run', default='')
    args = ap.parse_args()
    run = Run(args); D, K = args.D, args.K; n = D - 1
    lockf = os.path.join(run.dir, 'driver_k%d.pid' % K)
    if os.path.exists(lockf):
        try:
            pid = int(open(lockf).read().strip()); os.kill(pid, 0)
            print('another driver (pid %d) is running on %s' % (pid, run.dir)); sys.exit(2)
        except (ProcessLookupError, ValueError): pass
    open(lockf, 'w').write(str(os.getpid()))
    avals = [int(x) for x in args.a.split(',')] if args.a else list(range((K + 1) // 2, K + 1))
    plan = []
    for a in avals:
        if not os.path.exists(repfile(n, a)):
            if lemma_excludes(D, K, a):
                if not any(r.get('a') == a and r.get('excluded') for r in run.done):
                    run.record(dict(a=a, excluded='weight lemma: W_min(%d)=%d > %d*%d' % (a, wmin(n, a), n, K - a), leaves=0, found=0))
                run.log('k=%d a=%d excluded by the weight lemma (W_min=%d > %d)' % (K, a, wmin(n, a), n * (K - a)))
                continue
            raise SystemExit('missing representative file for a=%d and the weight lemma does not exclude it' % a)
        nrep = os.path.getsize(repfile(n, a)) // 8
        plan.append((a, nrep))
    run.log('start D=%d K=%d plan=%s workers=%d target=%.0fs engine=%d hot=%d' % (D, K, plan, args.workers, args.target, args.engine, args.hot))

    # task generator with adaptive chunk sizes (per-a cost estimate from completed chunks)
    state = {a: dict(nrep=nrep, gaps=uncovered(nrep, run.covered(a))) for a, nrep in plan}
    glock = threading.Lock()
    def cost_per_rep(a):
        recs = [r for r in run.done if r.get('a') == a and 'first' in r and r['last'] > r['first']]
        if not recs: return None
        return max(1e-7, sum(r['sec'] for r in recs) / sum(r['last'] - r['first'] for r in recs))
    def next_task():
        with glock:
            for a, nrep in plan:
                g = state[a]['gaps']
                if not g: continue
                f, l = g[0]
                cpr = cost_per_rep(a)
                size = 1 if cpr is None else max(1, int(args.target / cpr))
                size = min(size, l - f)
                if f + size >= l: g.pop(0)
                else: g[0] = (f + size, l)
                return (a, f, f + size)
            return None
    errors = []
    def worker():
        while not errors:
            t = next_task()
            if t is None: return
            try:
                rec = run_chunk(run, *t)
                rem = sum(l - f for a2, _ in plan for f, l in state[a2]['gaps'])
                run.log('done a=%d [%d,%d) leaves=%d found=%d sec=%.1f ns/leaf=%.2f  (reps left in queue: %d)' % (
                    rec['a'], rec['first'], rec['last'], rec['leaves'], rec['found'], rec['sec'],
                    1e9 * rec['sec'] / max(1, rec['leaves']), rem))
                for l in rec['resolving']: run.log('FOUND ' + l)
            except Exception as ex:
                errors.append(str(ex)); run.log('ERROR ' + str(ex))
    th = [threading.Thread(target=worker) for _ in range(args.workers)]
    for x in th: x.start()
    for x in th: x.join()
    if errors: raise SystemExit('errors: %s' % errors)

    # final consistency check and summary
    summ = dict(D=D, K=K, per_a={}, excluded={})
    for r in run.done:
        if r.get('excluded'): summ['excluded'][str(r['a'])] = r['excluded']
    tot = dict(leaves=0, nodes=0, fullchecks=0, found=0, sec=0.0, reps=0)
    resolving = []
    for a, nrep in plan:
        cov = run.covered(a)
        pos = 0
        for f, l in cov:
            if f != pos: raise SystemExit('coverage error a=%d at %d (chunk starts %d)' % (a, pos, f))
            pos = l
        if pos != nrep: raise SystemExit('coverage error a=%d: ends at %d of %d' % (a, pos, nrep))
        recs = [r for r in run.done if r.get('a') == a and 'first' in r]
        s = {k: sum(r[k] for r in recs) for k in ('leaves', 'nodes', 'fullchecks', 'found', 'reps', 'feasible_reps')}
        s['sec'] = round(sum(r['sec'] for r in recs), 2); s['nrep'] = nrep; s['chunks'] = len(recs)
        summ['per_a'][str(a)] = s
        for k in ('leaves', 'nodes', 'fullchecks', 'found', 'reps'): tot[k] += s[k]
        tot['sec'] += s['sec']
        for r in recs: resolving += r['resolving']
    tot['sec'] = round(tot['sec'], 2); summ['total'] = tot; summ['resolving'] = resolving
    summ['ns_per_leaf'] = round(1e9 * tot['sec'] / max(1, tot['leaves']), 3)
    with open(os.path.join(run.dir, 'summary_k%d.json' % K), 'w') as f: json.dump(summ, f, indent=1, sort_keys=True)
    run.log('COMPLETE D=%d K=%d leaves=%d found=%d cpu=%.1fs ns/leaf=%.2f' % (D, K, tot['leaves'], tot['found'], tot['sec'], summ['ns_per_leaf']))
    os.remove(lockf)

if __name__ == '__main__':
    main()
