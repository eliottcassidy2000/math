"""Independent complete top-twenty audit using the inherited native referee.

No mathematical producer is imported or executed. The native engine rebuilds
the atlas, every full-profile domain, every pair table and every Kruskal tree.
"""
from pathlib import Path
from hashlib import sha256
from itertools import combinations
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

sys.stdout.reconfigure(encoding='utf-8', newline='\n')
HERE = Path(__file__).resolve().parent
STEM = 'continuing10_20260908_lrc_upper_twenty_audit'
GATES = 0

def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise ArithmeticError(why)

def canonical(x):
    return json.dumps(x, sort_keys=True, separators=(',', ':')).encode()

def main():
    ap = argparse.ArgumentParser()
    filed = HERE.name == '04-computation'
    ap.add_argument('--root', type=Path, default=HERE.parent if filed else Path('C:/w/s0905'))
    ap.add_argument('--producer', type=Path, default=HERE.parent / '05-knowledge/results' if filed else Path('C:/w/continuing10_20260908_upper'))
    ap.add_argument('--work-dir', type=Path, default=Path(tempfile.gettempdir()) / STEM)
    args = ap.parse_args()
    profile = args.root / '04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json'
    raw = profile.read_bytes()
    need(sha256(raw).hexdigest() == '935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f', 'complete inherited profile pin')
    P = json.loads(raw)['levels']
    states = sorted({row[0] for row in P['6']['profiles']})
    need(states == P['6']['gcds'] and len(states) == 42, 'complete sheet alphabet')
    old_raw = (args.root / '05-knowledge/results/continuing8_20260906_lrc_minimum_tree_certificate.json').read_bytes()
    need(sha256(old_raw).hexdigest() == '580a7c930103aab3bea867ad463a90b0e0208323a90ee95a685ff811a761582d', 'old minimum-tree certificate pin')
    old = json.loads(old_raw)['new_scales']
    need(len(old) == 7646 and max(old) == 11995 and old.count(7200) == 1, 'old exact array cardinality and membership')
    need(sha256(canonical(old)).hexdigest() == '8ffc6d14b3883cf7e02c3ab02ddca5339d909a8051411def9096dee83b0aaed7', 'old semantic array pin')
    baseline = [t for t in old if t != 7200]
    need(sha256(canonical(baseline)).hexdigest() == 'f8c42793c4a5081d40a3937c9a79d9cd307e6c08c8b9ea3e94bac5388e47d16c', 'audited single-clock subtraction')
    upper = list(reversed(baseline[-20:]))
    need(len(upper) == 20 and upper[0] == 11995 and upper[-1] == 11940, 'entire twenty highest remaining clocks')
    cert_raw = (args.producer / 'continuing10_20260908_lrc_upper_twenty_certificate.json').read_bytes()
    certificate_sha = sha256(cert_raw).hexdigest()
    need(certificate_sha == '8a4ea5f6ec5669d0b687e3b52519ad8df9a43b0e9d64fa09a1565cfd3cec29ab', 'frozen complete producer-certificate pin')
    J = json.loads(cert_raw)
    records = J['clocks']
    need([r['t'] for r in records] == upper and J['removed_scales'] == upper, 'all and only declared clocks')
    need(J['profile_sha256'] == sha256(raw).hexdigest(), 'declared full-profile dependency')
    need(J['baseline_semantic_sha256'] == sha256(canonical(old)).hexdigest(), 'declared old baseline dependency')
    need(J['baseline_after_7200_semantic_sha256'] == sha256(canonical(baseline)).hexdigest(), 'declared after7200 baseline dependency')
    new = [t for t in baseline if t not in set(upper)]
    need(J['new_scales'] == new and len(new) == 7625 and max(new) == 11935, 'literal retained array')
    need(J['new_scales_sha256'] == sha256(canonical(new)).hexdigest() == '51c00db52c18732ec6128e4d9e1ed197813ab107bfd3dec7dc7ace9c5c3c67c7', 'retained semantic array digest')
    domains = sorted({tuple(d for d in states if t % d == 0) for t in upper})
    bank_pins = {tuple(d['domain']): d['words_sha256'] for d in J['domains']}
    need(len(domains) == 16 and set(domains) == set(bank_pins), 'complete divisor-domain bank')
    index = {d: i for i, d in enumerate(domains)}
    profiles = [(int(k), c, w) for k, row in P.items() for c, w in row['profiles']]
    lines = [str(len(profiles))]
    lines.extend(' '.join(map(str, [k, c] + w)) for k, c, w in profiles)
    lines.append(str(len(domains)))
    lines.extend(' '.join(map(str, [len(d)] + list(d))) for d in domains)
    lines.append(str(len(records)))
    for r in records:
        D = tuple(d for d in states if r['t'] % d == 0)
        need(list(D) == r['domain'], 'native clock divisibility determines domain')
        table = {tuple(ab): v[0] for ab, v in r['weights']}
        need(set(table) == set(combinations(D, 2)) | {(d, d) for d in D}, 'every unordered pair, including repeated sheets')
        parents = list(range(7))
        def find(i):
            while parents[i] != i:
                i = parents[i]
            return i
        total = 0
        need(len(r['owner_edges']) == 6, 'six owner-tree edges')
        for i, j, weight in r['owner_edges']:
            need(0 <= i < j < 7 and find(i) != find(j), 'owner edges are a positional tree')
            parents[find(i)] = find(j)
            need(weight == table[tuple(sorted((r['owner'][i], r['owner'][j])))], 'owner edge uses its margin-pair table')
            total += weight
        need(total == r['owner_tree'], 'owner total')
        lines.append(' '.join(map(str, [r['t'], index[D], r['word_count'], r['minimum_margin'], r['owner_E'], r['owner_tree'], r['event_evaluations'], r['compatible_atlas_edges']] + r['owner'])))
        lines.append(str(len(r['survivors'])))
        lines.extend(' '.join(map(str, word + [excess, credit])) for word, excess, credit in r['survivors'])
        lines.append(str(len(r['weights'])))
        lines.extend(' '.join(map(str, ab + [len(v)] + v)) for ab, v in r['weights'])
    for tree in [[(i, i + 1) for i in range(6)], [(0, i) for i in range(1, 7)]]:
        for mask in range(128):
            edges = sum(bool(mask & (1 << i)) and bool(mask & (1 << j)) for i, j in tree)
            need(edges <= max(mask.bit_count() - 1, 0), 'actual forest pointwise multiplicity bound')
    need(3 > 3 - 1, 'triangle sum fails the forest predicate')
    need(sum(0 if j == 1 else 10 for j in range(1, 7)) > 0, 'invented maximum star does not lower-bound zero actual path')
    work = args.work_dir / ('optimized' if sys.flags.optimize else 'normal')
    work.mkdir(parents=True, exist_ok=True)
    native_input = work / 'input.txt'
    native_input.write_bytes(('\n'.join(lines) + '\n').encode())
    compiler = shutil.which('g++') or 'C:/Users/Eliott/scoop/apps/gcc/current/bin/g++.exe'
    exe = work / ('referee.exe' if os.name == 'nt' else 'referee')
    options = ['-std=c++17', '-O3', '-DNDEBUG'] if sys.flags.optimize else ['-std=c++17', '-O2']
    subprocess.run([compiler, *options, str(HERE / (STEM + '.cpp')), '-o', str(exe)], check=True, capture_output=True)
    env = os.environ.copy()
    env['PATH'] = str(Path(compiler).resolve().parent) + os.pathsep + env.get('PATH', '')
    with (work / 'progress.txt').open('wb') as err:
        proc = subprocess.run([str(exe), str(native_input), str(work)], stdout=subprocess.PIPE, stderr=err, env=env)
    if proc.returncode:
        raise ArithmeticError((work / 'progress.txt').read_text()[-4000:])
    native = proc.stdout.replace(b'\r\n', b'\n')
    for i, D in enumerate(domains):
        need(sha256((work / ('words_' + str(i) + '.json')).read_bytes()).hexdigest() == bank_pins[D], 'complete unpruned word-list digest')
    print('INDEPENDENT_ENGINE raw spatial cells; signed ceil-floor phase sweep; unpruned126-profile census; Kruskal')
    print(native.decode().strip())
    print('RETAINED_SCALES', len(new), 'MAXIMUM', max(new))
    print('RETAINED_SEMANTIC_SHA256', sha256(canonical(new)).hexdigest())
    print('PRODUCER_CERTIFICATE_SHA256', certificate_sha)
    print('PASS', GATES, 'always-active Python input/completeness gates; output LF')

if __name__ == '__main__':
    main()
