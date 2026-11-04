"""T5 integer skeleton checks and exact identification of the six-vertex survivor."""
import argparse
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path
import check

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    files = [HERE/'structure.py', HERE/'check.py']
    pins = {p.name: sha256(p.read_bytes()).hexdigest() for p in files}
    controls = [check.tournament(n, code) for n in range(1, 6)
                for code in range(1 << (n*(n-1)//2))]
    hostile = [4, 17, 42, 131, 45, 75, 31, 119]
    controls += [check.tournament(6, 408), hostile, check.inflate([2, 4, 1], hostile)]
    entries = 0
    for a in controls:
        n, full = len(a), (1 << len(a))-1
        e, s = check.ends(a), check.ends(check.opposite(a))
        target = check.inflate([2, 4, 1], a)
        te, ts = check.ends(target)[-1], check.ends(check.opposite(target))[-1]
        for u in range(n):
            e2 = sum(e[m][u]*sum(e[full ^ m]) for m in range(1, full))
            b2 = sum(s[m][u]*sum(e[full ^ m]) for m in range(1, full))
            check.need(te[u] == e[-1][u]+e2 and ts[u] == s[-1][u]+b2, 'T5 integer skeleton')
            check.need((te[u]-s[-1][u]) % 2 == 0 and (ts[u]-e[-1][u]) % 2 == 0,
                       'T5 parity swap')
            entries += 1
        check.need(te[-2:] == ts[-2:] == [sum(e[-1])]*2, 'T5 exterior endpoints')
    standard = [14, 12, 56, 48, 35, 3]  # cyclic quotient, each vertex replaced by 0->1
    survivor = check.tournament(6, 22693)
    isomorphism = next((p for p in permutations(range(6)) if all(
        ((standard[u] >> v) & 1) == ((survivor[p[u]] >> p[v]) & 1)
        for u in range(6) for v in range(6))), None)
    check.need(isomorphism is not None, 'six-vertex survivor identification')
    check.need(pins == {p.name: sha256(p.read_bytes()).hexdigest() for p in files}, 'source changed')
    result = dict(status='FINITE-EXACT controls for proved T5', cores=len(controls),
                  internal_endpoint_pairs=entries, six_vertex_standard=standard,
                  survivor_code=22693, standard_to_survivor=isomorphism,
                  universe='All1099 labelled cores orders1..5, plus code408 at6, '
                           'the parity-balanced order8 hostile and its first cyclic extension at10',
                  source_sha256=pins)
    text = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(text, encoding='utf8', newline='\n')
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
