"""Verify sealed receipts; --replay reruns Python, --full also repeats C++ n<=8."""
from pathlib import Path
import argparse, hashlib, json, shutil, subprocess, sys, tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]

def need(value, message):
    if not value:
        raise RuntimeError(message)

def digest(path):
    # Git text files may be checked out with CRLF outside this packet.
    return hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()

def read(name):
    return json.loads((HERE / name).read_text(encoding='utf-8-sig'))

def run_json(command):
    result = subprocess.run(command, cwd=ROOT, check=True, capture_output=True)
    return json.loads(result.stdout.decode('utf-8-sig'))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay', action='store_true')
    parser.add_argument('--full', action='store_true')
    args = parser.parse_args()
    expected = {}
    for line in (HERE / 'SHA256SUMS').read_text().splitlines():
        h, name = line.split('  ', 1)
        need(name not in expected and Path(name).name == name, 'manifest path')
        expected[name] = h
    actual = {p.name for p in HERE.iterdir() if p.is_file() and p.name != 'SHA256SUMS'}
    need(actual == set(expected), 'manifest coverage')
    for name, h in expected.items():
        need(digest(HERE / name) == h, 'packet hash: ' + name)
    sources = read('SOURCES.json')
    for item in sources['executable_inputs']:
        need(digest(ROOT / item['path']) == item['sha256_lf'], 'input drift: ' + item['path'])
    for stem in ['audit', 'literal', 'phase']:
        need(read(stem + '-normal.json') == read(stem + '-optimized.json'), stem + ' normal/-O')
    census = {r['n']: r for r in read('census.json')}
    for n, r in census.items():
        need(r['presentations'] == 1 << ((n-1)*(n-2)//2), 'presentation universe')
        need(all(r[k] == r['strong_classes'] for k in ['joint_signatures', 'free_signatures', 'joint_support_signatures', 'free_support_signatures']), 'finite atom classifier')
    for r in read('audit-normal.json')['atom_pushforward']:
        need(r['joint_signatures'] == census[r['n']]['joint_signatures'] and r['free_signatures'] == census[r['n']]['free_signatures'], 'independent atom overlap')
    for r in read('literal-normal.json')['literal_reset_quotients']:
        if r['n'] < 3:
            continue
        for key, field in [('free/free_reset', 'free_pair_cells'), ('free/colour/free_reset', 'free_pair_colour_cells'), ('joint/free_reset', 'rooted_pair_cells'), ('joint/colour/free_reset', 'rooted_pair_colour_cells')]:
            need(r['protocols'][key]['cells'] == census[r['n']][field], 'independent pair-cell overlap')
    need(census[8]['classes'] == 6880 and census[8]['strong_classes'] == 6008, 'order-eight census')
    need(read('phase-normal.json')['unnamed_target_law_total_variation'] == '17071/17145', 'hostile law')
    print('Seals, executable inputs, normal/-O receipts and overlapping implementations agree.', flush=True)
    if args.replay or args.full:
        for stem in ['audit', 'literal', 'phase']:
            for optimized in [False, True]:
                command = [sys.executable, '-B'] + (['-O'] if optimized else []) + [str(HERE / (stem + '.py'))]
                need(run_json(command) == read(stem + '-normal.json'), 'replay: ' + stem)
                print(stem + (' optimized' if optimized else ' normal') + ': reproduced', flush=True)
    if args.full:
        compiler = shutil.which('g++')
        need(compiler, 'g++ is required for --full')
        with tempfile.TemporaryDirectory(prefix='unnamed-carrier-') as directory:
            exe = str(Path(directory) / ('census.exe' if sys.platform == 'win32' else 'census'))
            subprocess.run([compiler, '-std=c++17', '-O3', str(HERE / 'census.cpp'), '-o', exe], check=True)
            need(run_json([exe, '8']) == read('census.json'), 'C++ order-eight replay')
        print('C++ exhaustive order-eight census: reproduced', flush=True)

if __name__ == '__main__':
    main()
