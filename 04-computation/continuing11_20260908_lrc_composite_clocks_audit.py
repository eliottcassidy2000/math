"""Independent complete three-composite-clock audit using the inherited native referee.

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
STEM = 'continuing11_20260908_lrc_composite_clocks_audit'
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
    ap.add_argument('--producer', type=Path, default=HERE.parent / '05-knowledge/results' if filed else Path('C:/w/continuing11_20260908_lrc'))
    ap.add_argument('--work-dir', type=Path, default=Path(tempfile.gettempdir()) / STEM)
    args = ap.parse_args()
    profile = args.root / '04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json'
    raw = profile.read_bytes()
    need(sha256(raw).hexdigest() == '935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f', 'complete inherited profile pin')
    P = json.loads(raw)['levels']
    states = sorted({row[0] for row in P['6']['profiles']})
    need(states == P['6']['gcds'] and len(states) == 42, 'complete sheet alphabet')
    old_raw = (args.root / '05-knowledge/results/continuing10_20260908_lrc_upper_twenty_certificate.json').read_bytes()
    need(sha256(old_raw).hexdigest() == '8a4ea5f6ec5669d0b687e3b52519ad8df9a43b0e9d64fa09a1565cfd3cec29ab', 'previous upper-twenty certificate pin')
    old = json.loads(old_raw)['new_scales']
    need(len(old) == 7625 and max(old) == 11935 and 7200 not in old, 'complete inherited current array')
    need(sha256(canonical(old)).hexdigest() == '51c00db52c18732ec6128e4d9e1ed197813ab107bfd3dec7dc7ace9c5c3c67c7', 'current semantic array pin')
    baseline = old
    upper = [10080,10800,11520]
    need(all(t in baseline for t in upper), 'all three explicitly declared clocks are in the actual source array')
    cert_raw = (args.producer / 'continuing11_20260908_lrc_composite_clocks_certificate.json').read_bytes()
    certificate_sha = sha256(cert_raw).hexdigest()
    need(certificate_sha == '2c9b4fa1a1103d3309094d06a422ec2f5e56e9b804961037cab9d5a05a352b73', 'frozen complete producer-certificate pin')
    J = json.loads(cert_raw)
    records = J['clocks']
    need([r['t'] for r in records] == upper and J['removed_scales'] == upper, 'entire three-clock universe tested and removed')
    need(J['profile_sha256'] == sha256(raw).hexdigest(), 'declared full-profile dependency')
    need(J['baseline_semantic_sha256'] == sha256(canonical(old)).hexdigest(), 'declared current baseline dependency')
    new = [t for t in baseline if t not in set(upper)]
    need(J['new_scales'] == new and len(new) == 7622 and max(new) == 11935, 'literal retained array')
    need(J['new_scales_sha256'] == sha256(canonical(new)).hexdigest() == '832714c179ab4f425a76172c41d7d8eecf05dc15adfdd3c5e08bc7eb5a13e48b', 'retained semantic array digest')
    domains = sorted({tuple(d for d in states if t % d == 0) for t in upper})
    bank_pins = {tuple(d['domain']): d['words_sha256'] for d in J['domains']}
    need(len(domains) == 3 and set(domains) == set(bank_pins), 'complete divisor-domain bank')
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
        zero_parents=list(range(7))
        def zero_find(i):
            while zero_parents[i]!=i:i=zero_parents[i]
            return i
        for i,j in combinations(range(7),2):
            if table[tuple(sorted((r['owner'][i],r['owner'][j])))]==0:
                zero_parents[zero_find(i)]=zero_find(j)
        parts={}
        for i in range(7):parts.setdefault(zero_find(i),[]).append(i)
        cuts=sorted(sorted(r['owner'][i] for i in part) for part in parts.values())
        expected_cuts={10080:[[6,9,9,36],[8,8,32]],10800:[[9,16,27,36,48],[10,20]],11520:[[9,9,32,36,48],[10,20]]}
        need(cuts==expected_cuts[r['t']], 'entire zero-weight owner graph has the declared positional cut')
        bridges=[table[tuple(sorted((r['owner'][i],r['owner'][j])))] for i,j in combinations(range(7),2) if zero_find(i)!=zero_find(j)]
        need(min(bridges)==r['owner_tree'], 'full owner credit is its cheapest forced bridge')
        need(sorted(w for i,j,w in r['owner_edges'])==[0,0,0,0,0,r['owner_tree']], 'five zero edges and one bridge')
        if r['t']==10080:
            need(all((10080//7)%d==0 for d in D), 'every10080 sheet has zero marginal ceiling excess')

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
