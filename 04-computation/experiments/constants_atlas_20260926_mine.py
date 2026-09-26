#!/usr/bin/env python3
"""constants_atlas_20260926_mine.py -- mine the repository for numerical constants that recur
across DIFFERENT research threads (session collatz-exponent-atlas-20260926, opus).

Universe: every .md file under 01-canon/theorems, 05-knowledge/results, 05-knowledge/hypotheses,
07-reflections and 06-writeups. Navigation/atlas files and agent inboxes are excluded: they
aggregate every thread and would give every constant every label.
Thread label: keyword vote over the file text; a file carries its top two labels only (>= 3 hits).
Constants kept: decimals with >= 4 significant digits (normalised to 4 s.f.), fractions a/b with
b >= 9 or a >= 3 (normalised, and also mapped to their 4 s.f. decimal so 19/12 and 1.583 meet),
and a curated list of integers of interest. Ubiquitous small numbers are dropped by these filters.
Ranking: number of distinct threads in which the constant occurs in >= 2 files each.
This is a lexical census, not a theorem: a shared number is a lead, never a connection.
Usage: python3 constants_atlas_20260926_mine.py <repo root> [min_threads]
"""
import os, re, sys, math, json
from fractions import Fraction
from collections import defaultdict

ROOT = sys.argv[1] if len(sys.argv) > 1 else '.'
MIN_T = int(sys.argv[2]) if len(sys.argv) > 2 else 3
DIRS = ['01-canon/theorems', '05-knowledge/results', '05-knowledge/hypotheses', '07-reflections', '06-writeups']

THREADS = {
    'collatz': [r'\bcollatz\b', r'\bsyracuse\b', r'3n\+1', r'3x\+1', r'3n-1'],
    'lrc': [r'lonely runner', r'\bLRC\b', r'loneliness'],
    'jc': [r'jacobian conjecture', r'\bJC\(2\)', r'dixmier', r'danielewski'],
    'nc2': [r'\bNC2\b', r'\bGMC\b', r'gaussian moment'],
    'tournament': [r'\btournament', r'\bpaley\b', r'kirkman'],
    'hadamard': [r'\bhadamard\b', r'crouzeix'],
    'amm': [r'\bAMM\b', r'monthly problem', r'12592', r'fair coin'],
    'pell': [r'\bpell\b', r'pythagorean', r'berggren'],
    'kakeya': [r'\bkakeya\b'],
    'erdos': [r'erd[oő]s'],
    'additive': [r'additive bas', r'2-4-6-8', r'binomial additive'],
    'abc': [r'\babc\b', r'\bIUT\b', r'mochizuki'],
    'bernoulli': [r'bernoulli'],
    'transcendence': [r'transcend', r'e-function', r'siegel'],
    'mahler': [r'\bmahler\b', r'z-number'],
    'graph': [r'kuratowski', r'square-sum', r'graceful'],
    'zeta': [r'riemann', r'\bzeta\b'],
    'sequences': [r'\bOEIS\b', r'integer sequence', r'\bsumset', r'sidon'],
}
THREAD_RE = {k: [re.compile(p, re.I) for p in v] for k, v in THREADS.items()}

DEC_RE = re.compile(r'(?<![\w./^])(\d+\.\d{4,})(?![\w/])')
FRAC_RE = re.compile(r'(?<![\w.^])(\d{1,4})/(\d{1,4})(?![\w.])')
INT_INTEREST = {41, 183, 189, 65537, 12592, 2187, 2048, 6561, 507, 1729, 496, 8128, 2310, 30030, 5040, 1296, 4096, 65536, 257, 641, 6700417,
                139, 6592, 2457, 1093, 3511, 691, 24, 168, 240, 336, 1024}
INT_RE = re.compile(r'(?<![\w.^])(\d{2,8})(?![\w.])')
TRIVIAL_DEC = {'0.5', '0.25', '0.75', '0.125', '1', '2', '3', '4', '10', '100', '1000', '0.1', '0.01', '0.001', '0.2', '0.3', '0.4', '0.6', '0.7', '0.8', '0.9', '1.5', '2.5', '0.05'}


def sig4(x):
    if x == 0:
        return '0'
    e = math.floor(math.log10(abs(x)))
    return '%.4g' % (round(x, 3 - e))


def labels_for(text):
    votes = {}
    for k, pats in THREAD_RE.items():
        n = sum(len(p.findall(text)) for p in pats)
        if n >= 3:
            votes[k] = n
    top = sorted(votes.items(), key=lambda kv: -kv[1])[:2]
    return [k for k, v in top]


def main():
    const_files = defaultdict(lambda: defaultdict(set))
    const_kind = {}
    nfiles = 0
    for d in DIRS:
        base = os.path.join(ROOT, d)
        if not os.path.isdir(base):
            continue
        for fn in sorted(os.listdir(base)):
            if not fn.endswith('.md'):
                continue
            p = os.path.join(base, fn)
            try:
                text = open(p, encoding='utf-8', errors='ignore').read()
            except Exception:
                continue
            nfiles += 1
            labs = labels_for(text)
            if not labs:
                continue
            rel = d + '/' + fn
            seen = set()
            for m in DEC_RE.finditer(text):
                try:
                    x = float(m.group(1))
                except ValueError:
                    continue
                if x > 1e7 or x == 0:
                    continue
                key = 'd:' + sig4(x)
                if key[2:] in TRIVIAL_DEC:
                    continue
                seen.add(key)
                const_kind[key] = 'decimal'
            for m in FRAC_RE.finditer(text):
                a, b = int(m.group(1)), int(m.group(2))
                if b == 0 or a == 0 or a >= b * 20:
                    continue
                if b < 9 and a < 3:
                    continue
                fr = Fraction(a, b)
                if fr.denominator == 1 or (fr.denominator < 9 and fr.numerator < 3):
                    continue
                key = 'f:%d/%d' % (fr.numerator, fr.denominator)
                seen.add(key)
                const_kind[key] = 'fraction'
                key2 = 'd:' + sig4(float(fr))
                if key2[2:] not in TRIVIAL_DEC:
                    seen.add(key2)
                    const_kind.setdefault(key2, 'decimal(from fraction)')
            for m in INT_RE.finditer(text):
                n = int(m.group(1))
                if n in INT_INTEREST:
                    key = 'i:%d' % n
                    seen.add(key)
                    const_kind[key] = 'integer'
            for key in seen:
                for lab in labs:
                    const_files[key][lab].add(rel)
    rows = []
    for key, th in const_files.items():
        threads = [t for t, fs in th.items() if len(fs) >= 2]
        files = set().union(*th.values())
        rows.append((len(threads), len(files), key, th))
    rows.sort(key=lambda r: (-r[0], -r[1], r[2]))
    print("files scanned: %d (labelled files only are counted)" % nfiles)
    print("constants occurring in >= 2 files of each of >= %d distinct threads; one example file per thread" % MIN_T)
    for nt, nf, key, th in rows:
        if nt < MIN_T:
            break
        ex = '; '.join('%s(%d):%s' % (t, len(fs), sorted(fs)[0].split('/')[-1][:48]) for t, fs in sorted(th.items(), key=lambda kv: -len(kv[1])) if len(fs) >= 2)
        print("%-14s %-22s threads=%2d files=%4d | %s" % (key, const_kind.get(key, ''), nt, nf, ex))
    out = {key: {t: sorted(fs) for t, fs in th.items()} for nt, nf, key, th in rows if nt >= MIN_T}
    with open(os.path.join(ROOT, '05-knowledge/results/constants_atlas_20260926_mine.json'), 'w') as f:
        json.dump(out, f, indent=1)


if __name__ == '__main__':
    main()
