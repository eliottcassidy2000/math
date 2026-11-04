"""Independent integer-DP checks of the actual compiled parity functions."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import check

HERE = Path(__file__).resolve().parent


def core(family):
    a = [0]*11
    for u in range(11):
        for v in range(u+1, 11):
            d = (v-u) % 11
            edge = family == 0 or (d <= 5 if family == 1 else d in (1, 3, 4, 5, 9))
            if edge:
                a[u] |= 1 << v
            else:
                a[v] |= 1 << u
    return a


def birth(a, mask):
    n = len(a)
    out = a[:] + [0]
    for u in range(n):
        if mask >> u & 1:
            out[u] |= 1 << n
        else:
            out[n] |= 1 << u
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    files = [HERE/n for n in ('audit.py', 'audit.cpp', 'twelve.cpp', 'check.py')]
    pins = {p.name: sha256(p.read_bytes()).hexdigest() for p in files}
    controls = [check.tournament(n, code) for n in (2, 4)
                for code in range(1 << (n*(n-1)//2))]
    controls += [check.tournament(6, 408), [4, 17, 42, 131, 45, 75, 31, 119]]
    controls += [birth(core(family), mask) for family in range(3)
                 for mask in (0, 1, 2, 3, 31, 341, 682, 1023, 2047)]
    hostile8 = [4, 17, 42, 131, 45, 75, 31, 119]
    controls.append(check.inflate([2, 4, 1], check.inflate([2, 4, 1], hostile8)))
    compiler = shutil.which('g++')
    check.need(compiler is not None, 'existing g++')
    with tempfile.TemporaryDirectory(prefix='twelve-independent-') as temp:
        temp = Path(temp)
        inp, exe = temp/'input.txt', temp/'audit.exe'
        inp.write_text(str(len(controls))+'\n'+'\n'.join(
            str(len(a))+' '+' '.join(map(str, a)) for a in controls)+'\n', encoding='ascii')
        subprocess.run([compiler, '-O2', '-std=c++17', str(HERE/'audit.cpp'), '-o', str(exe)],
                       check=True, capture_output=True)
        actual = json.loads(subprocess.run([str(exe), str(inp)], check=True,
                                          capture_output=True, text=True).stdout)
    check.need(len(actual) == len(controls), 'output count')
    columns, true_cover, true_cubic = 0, 0, 0
    witnesses = []
    for a, observed in zip(controls, actual):
        n = len(a)
        s, e, c = check.sidecar(a)
        ep = sum((x & 1) << u for u, x in enumerate(e[-1]))
        sp = sum((x & 1) << u for u, x in enumerate(s[-1]))
        cover = all((c[u][v]-c[v][u]) % 4 == 0 for u in range(n) for v in range(n))
        # The endpoint-sum formula (3), not the truncated formula used in C++.
        cubic = [sum(check.cubic(a, u, v, (s, e)) << u for u in range(n) if u != v)
                 for v in range(n)]
        check.need(observed == [ep, sp, int(cover)]+cubic, 'independent integer parity comparison')
        columns += n
        true_cover += cover
        true_cubic += not any(cubic)
        if n == 12 and not any(cubic):
            witnesses.append(dict(rows=a, H=sum(e[-1]), S=s[-1], E=e[-1], cubic_columns=cubic))
    quartics = []
    for code in (1023, 1021):
        a = check.tournament(5, code)
        coefficient = 0
        for mask in range(16):
            x = (mask & 1) | (mask & 2)
            y = (((mask >> 2) & 1) << 1) | (((mask >> 3) & 1) << 4)
            coefficient += (-1)**(4-mask.bit_count())*check.h(check.extend(a, x, y))
        quartics.append(coefficient)
    check.need(quartics == [4, 2], 'general pairs retain a quartic difference')
    check.need(pins == {p.name: sha256(p.read_bytes()).hexdigest() for p in files}, 'source changed')
    result = dict(status='FINITE-EXACT independent implementation controls',
                  controls=len(controls), cubic_columns_checked=columns,
                  true_cover_controls=true_cover, zero_cubic_controls=true_cubic,
                  twelve_cubic_blind_witnesses=witnesses,
                  general_pair_quartic_hostile=dict(n=5, codes=[1023, 1021],
                      monomial='a_0 a_1 c_1 c_4', coefficients=quartics, actual_queries=32),
                  universe='All labelled orders2,4; order6 cover-symmetric and order8 endpoint-balanced hostiles; '
                           '27 declared order12 extensions (9 named masks for each of 3 cores); '
                           'one endpoint-parity-balanced order12 cyclic tower over the order8 hostile',
                  source_sha256=pins)
    text = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(text, encoding='utf8', newline='\n')
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
