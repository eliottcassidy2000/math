"""Corrected legacy census plus independent, streamed literal certificate audit.

Standard library only. Build products and stderr captures live in a temporary
directory. No producer Python modules or tournament engines are imported.
Fast stored-receipt comparison: python check.py
Full replay: python check.py --census [--output normal.json] [--orders 8 9 10]
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time


HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
SOURCE = HERE / 'blue_line_black_detours_corrected.cpp'
SOURCE_SHA256 = '4499981f6b76688d9477e40bb1aed3312a4bd226495cbab729d84d31212ff232'


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def arc(rows, u, v):
    return (rows[u] >> v) & 1


def edges(word):
    return frozenset(tuple(sorted((u, v))) for u, v in zip(word, word[1:]))


def blue(rows, word):
    n = len(rows)
    reversal = [0] * n
    for i, v in enumerate(word):
        reversal[v] = word[n - 1 - i]
    return all(arc(rows, u, v) == arc(rows, reversal[v], reversal[u])
               for u in range(n) for v in range(u + 1, n))


def move(rows, word):
    # Pair-by-pair definition, independent of producer's row complement/restore.
    keep = edges(word)
    n = len(rows)
    result = [0] * n
    for u in range(n):
        for v in range(u + 1, n):
            forward = arc(rows, u, v) ^ ((u, v) not in keep)
            a, b = (u, v) if forward else (v, u)
            result[a] |= 1 << b
    return tuple(result)


def validate_source(rows):
    n = len(rows)
    require(all(0 <= row < (1 << n) for row in rows), 'row bounds')
    require(all(not arc(rows, v, v) for v in range(n)), 'source loop')
    require(all(arc(rows, u, v) + arc(rows, v, u) == 1
                for u in range(n) for v in range(u + 1, n)), 'not a tournament')
    require(all(arc(rows, v, v + 1) for v in range(n - 1)), 'base not directed')
    require(arc(rows, n - 1, 0), 'base not closed')
    require(rows[0].bit_count() == n - 5 and rows[-1].bit_count() == 4,
            'wrong endpoint layer')
    require(blue(rows, tuple(range(n))), 'base not blue')


def validate_certificate(n, rows, words):
    validate_source(rows)
    current = rows
    parity = set()
    for word in words:
        require(sorted(word) == list(range(n)), 'word is not a permutation')
        require(all(arc(current, u, v) for u, v in zip(word, word[1:])),
                'word not directed at actual source')
        require(not blue(current, word), 'word blue at actual source')
        parity.symmetric_difference_update(edges(word))
        current = move(current, word)
    base = tuple(range(n))
    require(parity == set(edges(base)), 'path-edge XOR mismatch')
    require(current == move(rows, base), 'literal target mismatch')


def independent_universe(n):
    # Freely choose the three incoming ports of 0 in 2..n-2; reflection fixes
    # the other endpoint. Count nonpath pair-orbits entirely in the middle.
    free_pairs = {(u, v) for u in range(1, n - 1)
                  for v in range(u + 2, n - 1)}
    remaining = set(free_pairs)
    orbit_count = 0
    while remaining:
        u, v = min(remaining)
        orbit = {(u, v), (n - 1 - v, n - 1 - u)}
        require(orbit <= free_pairs, 'middle reflection not closed')
        remaining.difference_update(orbit)
        orbit_count += 1
    require(orbit_count == (n - 3) ** 2 // 4, 'orbit formula discrepancy')
    ports = math.comb(n - 3, 3)
    return {'port_sets': ports, 'free_middle_pair_orbits': orbit_count,
            'sources': ports * (1 << orbit_count)}


def hostile_controls():
    # An actual pure-blue target, so an unfiltered all-order promotion is false.
    source = (6, 4, 8, 19, 103, 79, 15)
    validate_source(source)
    target = move(source, tuple(range(7)))
    require(target == (122, 124, 121, 112, 32, 64, 16), 'n7 target')
    paths = []
    def visit(word, unused):
        if not unused:
            paths.append(tuple(word))
            return
        for v in sorted(unused):
            if not word or arc(target, word[-1], v):
                visit(word + [v], unused - {v})
    visit([], set(range(7)))
    require(len(paths) == 9 and all(blue(target, q) for q in paths), 'n7 hostile')
    # Check the independent replay itself rejects a blue base word as a repair.
    try:
        validate_certificate(7, source, (tuple(range(7)),) * 3)
    except RuntimeError as error:
        require(str(error) == 'word blue at actual source', 'wrong hostile rejection')
    else:
        raise RuntimeError('replay accepted a blue word')
    return {'n7_target_hamiltonian_paths': len(paths), 'all_blue': True,
            'blue_word_rejected': True}


def run_order(executable, n, temporary):
    universe = independent_universe(n)
    seen = set()
    digest = hashlib.sha256()
    lines = []
    command = [str(executable), str(n), '4']
    stderr_path = Path(temporary) / ('stderr-%d.txt' % n)
    began = time.monotonic()
    with stderr_path.open('wb') as errors:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=errors)
        try:
            for raw in process.stdout:
                if not raw.startswith(b'CERT '):
                    lines.append(raw.decode().rstrip())
                    continue
                raw = raw.replace(b'\r\n', b'\n')
                digest.update(raw)
                fields = list(map(int, raw.split()[1:]))
                require(len(fields) == 1 + 4 * n and fields[0] == n, 'certificate format')
                rows = tuple(fields[1:1 + n])
                words = tuple(tuple(fields[1 + j * n:1 + (j + 1) * n])
                              for j in range(1, 4))
                validate_certificate(n, rows, words)
                packed = sum(row << (n * u) for u, row in enumerate(rows))
                require(packed not in seen, 'duplicate marked source')
                seen.add(packed)
                if len(seen) % 16384 == 0:
                    print('n=%d: %d/%d literal certificates checked' %
                          (n, len(seen), universe['sources']), flush=True)
        except BaseException:
            process.kill()
            process.wait()
            raise
        require(process.wait() == 0, 'producer exit failure')
    require(len(seen) == universe['sources'], 'incomplete source universe')
    require(any(line.startswith('# decoder controls:') for line in lines), 'missing decoder controls')
    expected = 'k=4: %d / 0 / %d / 0 = 0 + 0' % (len(seen), len(seen))
    require(any(expected in line for line in lines), 'producer census mismatch')
    require(not any(line.startswith('NO ') for line in lines), 'unresolved source')
    result = {'n': n, 'universe': universe, 'unique_literal_certificates': len(seen),
              'certificate_stream_lf_sha256': digest.hexdigest(),
              'producer_arguments': [str(n), '4'], 'producer_stdout_without_certificates': lines,
              'producer_stderr': stderr_path.read_text().splitlines(),
              'elapsed_seconds': round(time.monotonic() - began, 3), 'status': 'PASS'}
    print('n=%d PASS: %d unique sources, all three moves black and literal target exact' %
          (n, len(seen)), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--orders', nargs='+', type=int, default=[8, 9, 10])
    parser.add_argument('--output', type=Path)
    parser.add_argument('--census', action='store_true', help='regenerate the full bounded census')
    parser.add_argument('--compiler', default='g++', help='compiler executable name on PATH or absolute path')
    parser.add_argument('--decoder-only', action='store_true',
                        help='exercise the minimal disconnected edge plus triangle at n=5')
    args = parser.parse_args()
    if not args.census and not args.decoder_only:
        from compare_receipts import paired_summary
        payload = json.dumps(paired_summary(HERE), indent=2) + '\n'
        if args.output:
            args.output.write_bytes(payload.encode())
        else:
            print(payload, end='')
        return
    require(all(n in (8, 9, 10) for n in args.orders), 'orders restricted to audited boundary')
    require(hashlib.sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256, 'corrected source hash')
    compiler = shutil.which(args.compiler)
    require(compiler is not None, 'compiler unavailable: ' + args.compiler)
    version = subprocess.check_output([compiler, '--version'], text=True).splitlines()[0]
    with tempfile.TemporaryDirectory(prefix='math-finite-head-') as temporary:
        executable = Path(temporary) / 'finite-head.exe'
        build = [compiler, '-O2', '-std=c++17', '-Wall', '-Wextra', str(SOURCE), '-o', str(executable)]
        built = subprocess.run(build, capture_output=True, text=True)
        require(built.returncode == 0, 'compile failed: ' + built.stderr)
        result = {'status': 'PASS', 'compiler': version, 'build_options': build[1:5],
                  'corrected_source_lf_sha256': SOURCE_SHA256,
                  'compiler_stderr': built.stderr.splitlines()}
        if args.decoder_only:
            tested = subprocess.run([str(executable), '5', '4', '0', '1'],
                                    capture_output=True, text=True)
            require(tested.returncode == 0, 'minimal decoder control failed')
            require('# decoder controls:' in tested.stdout, 'minimal control missing')
            result.update({'producer_arguments': ['5', '4', '0', '1'],
                           'disconnected_edges': [[0, 1], [2, 3], [3, 4], [2, 4]],
                           'producer_stdout': tested.stdout.splitlines(),
                           'producer_stderr': tested.stderr.splitlines()})
        else:
            result.update({'independent_hostiles': hostile_controls(),
                           'orders': [run_order(executable, n, temporary) for n in args.orders]})
    payload = json.dumps(result, indent=2) + '\n'
    if args.output:
        args.output.write_bytes(payload.encode())
        print('receipt:', args.output, flush=True)
    else:
        print(payload, end='')


if __name__ == '__main__':
    main()
