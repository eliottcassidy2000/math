#!/usr/bin/env python3
"""Independent complete native referee: raw walls, unpruned words, Kruskal.

No mathematical producer is imported or executed. The companion C++ engine
reconstructs every domain and all pair costs; the JSON certificate is only an
expected answer. Python validates inherited pins and serializes the native input.
"""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

sys.stdout.reconfigure(encoding='utf-8', newline='\n')
HERE = Path(__file__).resolve().parent
STEM = 'continuing8_20260906_lrc_minimum_tree_audit'
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
    ap.add_argument('--producer', type=Path, default=HERE.parent / '05-knowledge/results' if filed else Path('C:/w/continuing8_20260906_root'))
    ap.add_argument('--work-dir', type=Path, default=Path(tempfile.gettempdir()) / STEM)
    ap.add_argument('--scout', action='store_true', help='historical development comparison only')
    args = ap.parse_args()
    # Pointwise graph controls: different actual trees and every danger subset.
    trees = [[(i, i + 1) for i in range(6)], [(0, i) for i in range(1, 7)],
             [(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)]]
    for tree in trees:
        for mask in range(128):
            edges = sum(bool(mask & (1 << i)) and bool(mask & (1 << j)) for i, j in tree)
            need(edges <= max(mask.bit_count() - 1, 0), 'induced actual forest multiplicity control')
    need(3 > 3 - 1, 'triangle overlap cannot be added as a forest')
    # Invented nonactual edges may have positive weights while the actual path
    # has weight zero. Their possible maximum is not a lower bound for it.
    zero_path = set(trees[0])
    invented_star_weight = sum(0 if (0, j) in zero_path else 10 for j in range(1, 7))
    need(invented_star_weight == 50 and sum(0 for _ in zero_path) == 0,
         'maximum possible tree is not actual-tree credit')
    profile = args.root / '04-computation/overnight12_20260906_lrc_decoder_descent_inherited_profiles.json'
    baseline = args.root / '05-knowledge/results/continuing7_20260906_lrc_residue_domains_certificate.json'
    raw = profile.read_bytes()
    need(sha256(raw).hexdigest() == '935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f', 'raw complete-profile pin')
    P = json.loads(raw)['levels']
    states = sorted({row[0] for row in P['6']['profiles']})
    need(states == P['6']['gcds'] and len(states) == 42, 'complete sheet alphabet')
    old = json.loads(baseline.read_bytes())['new_scales']
    need(sha256(canonical(old)).hexdigest() == 'c67f5e98eff3fe406a4416057c6063095290330a50f039e731ced0fc2ca4657a', 'complete baseline semantic pin')
    upper = sorted((t for t in old if t >= 12000), reverse=True)
    need(len(old) == 8195 and len(upper) == 549 and 7200 in old, 'complete inherited target selection')
    if args.scout:
        records = [json.loads(s) for s in (args.producer / 'range_12000.jsonl').read_text().splitlines()]
        hostile = json.loads((args.producer / 'forest_7200.json').read_bytes())
        worst = hostile['worst']
        records.append(dict(t=7200, domain=hostile['domain'], word_count=hostile['layers'][-1],
                            minimum_margin=worst[0], owner=worst[1], owner_E=worst[2], owner_tree=worst[3],
                            survivors=hostile['survivors'], weights=hostile['weights']))
        bank_pins = {}
        for f in args.producer.glob('domain_*.json'):
            J = json.loads(f.read_bytes())
            bank_pins[tuple(J['domain'])] = sha256(canonical(J['words'])).hexdigest()
        certificate_sha = 'SCOUT_ONLY'
    else:
        cert_path = args.producer / 'continuing8_20260906_lrc_minimum_tree_certificate.json'
        cert_raw = cert_path.read_bytes()
        certificate_sha = sha256(cert_raw).hexdigest()
        need(certificate_sha == '580a7c930103aab3bea867ad463a90b0e0208323a90ee95a685ff811a761582d', 'frozen complete producer-certificate pin')
        J = json.loads(cert_raw)
        records = J['clocks']
        bank_pins = {tuple(d['domain']): d['words_sha256'] for d in J['domains']}
        need(J['profile_sha256'] == sha256(raw).hexdigest(), 'declared profile dependency')
        need(J['baseline_semantic_sha256'] == sha256(canonical(old)).hexdigest(), 'declared baseline dependency')
        new = [t for t in old if t not in set(upper)]
        need(J['removed_scales'] == upper and J['new_scales'] == new, 'literal complete deletion and retained array')
        need(J['new_scales_sha256'] == sha256(canonical(new)).hexdigest(), 'retained array semantic digest')
    need([r['t'] for r in records] == upper + [7200], 'exact clock order, no missing selected clock')
    domains = sorted({tuple(d for d in states if t % d == 0) for t in upper + [7200]})
    if args.scout:
        bank_pins = {d: bank_pins.get(d) for d in domains}
    need(set(domains) == set(bank_pins), 'every domain and no extraneous domain')
    index = {d: i for i, d in enumerate(domains)}
    lines = []
    profiles = [(int(k), c, w) for k, row in P.items() for c, w in row['profiles']]
    lines.append(str(len(profiles)))
    lines.extend(' '.join(map(str, [k, c] + w)) for k, c, w in profiles)
    lines.append(str(len(domains)))
    lines.extend(' '.join(map(str, [len(d)] + list(d))) for d in domains)
    lines.append(str(len(records)))
    for r in records:
        D = tuple(d for d in states if r['t'] % d == 0)
        need(list(D) == r['domain'], 'native divisibility determines full domain')
        if not args.scout:
            table = {tuple(ab): v[0] for ab, v in r['weights']}
            parents = list(range(7))
            def find(i):
                while parents[i] != i:
                    i = parents[i]
                return i
            total = 0
            need(len(r['owner_edges']) == 6, 'six positional owner-tree edges')
            for i, j, weight in r['owner_edges']:
                need(0 <= i < j < 7 and find(i) != find(j), 'owner edges form an actual positional tree')
                parents[find(i)] = find(j)
                need(weight == table[tuple(sorted((r['owner'][i], r['owner'][j])))], 'owner edge references reconstructed pair table')
                total += weight
            need(total == r['owner_tree'], 'owner edge total')
        lines.append(' '.join(map(str, [r['t'], index[D], r['word_count'], r['minimum_margin'], r['owner_E'], r['owner_tree'], r.get('event_evaluations', -1), r.get('compatible_atlas_edges', -1)] + r['owner'])))
        lines.append(str(len(r['survivors'])))
        lines.extend(' '.join(map(str, word + [excess, credit])) for word, excess, credit in r['survivors'])
        lines.append(str(len(r['weights'])))
        lines.extend(' '.join(map(str, ab + [len(v)] + v)) for ab, v in r['weights'])
    work = args.work_dir / ('optimized' if sys.flags.optimize else 'normal')
    work.mkdir(parents=True, exist_ok=True)
    native_input = work / 'input.txt'
    native_input.write_bytes(('\n'.join(lines) + '\n').encode())
    cpp = HERE / (STEM + '.cpp')
    compiler = shutil.which('g++') or 'C:/Users/Eliott/scoop/apps/gcc/current/bin/g++.exe'
    exe = work / ('referee.exe' if os.name == 'nt' else 'referee')
    options = ['-std=c++17', '-O3', '-DNDEBUG'] if sys.flags.optimize else ['-std=c++17', '-O2']
    subprocess.run([compiler, *options, str(cpp), '-o', str(exe)], check=True, capture_output=True)
    env = os.environ.copy()
    env['PATH'] = str(Path(compiler).resolve().parent) + os.pathsep + env.get('PATH', '')
    with (work / 'progress.txt').open('wb') as err:
        proc = subprocess.run([str(exe), str(native_input), str(work)], stdout=subprocess.PIPE, stderr=err, env=env)
    if proc.returncode:
        raise ArithmeticError((work / 'progress.txt').read_text()[-4000:])
    native = proc.stdout.replace(b'\r\n', b'\n')
    for i, D in enumerate(domains):
        words_raw = (work / ('words_' + str(i) + '.json')).read_bytes()
        if bank_pins[D] is not None:
            need(sha256(words_raw).hexdigest() == bank_pins[D], 'complete unpruned ordered word-list digest')
    new = [t for t in old if t not in set(upper)]
    print('INDEPENDENT_ENGINE raw spatial cells; signed ceil-floor phase sweep; unpruned126-profile census; Kruskal')
    print(native.decode().strip())
    print('RETAINED_SCALES', len(new), 'MAXIMUM', max(new))
    print('RETAINED_SEMANTIC_SHA256', sha256(canonical(new)).hexdigest())
    print('PRODUCER_CERTIFICATE_SHA256', certificate_sha)
    print('PASS', GATES, 'always-active Python input/completeness gates; output LF')


if __name__ == '__main__':
    main()
