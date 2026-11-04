"""Complete inherited even-order <=8 C-skew sieve, without endpoint filtering."""
import argparse
import gzip
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path
import check

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BANK_SHA = 'b4e559ad6cb31980852034f9c7f61e366a7cb479fee4d9a641ad3356745200d3'


def literal_cover(a):
    n = len(a)
    c = [[0]*n for _ in a]
    for word in permutations(range(n)):
        for cut in range(1, n):
            if all(a[word[i]] >> word[i+1] & 1
                   for i in range(n-1) if i != cut-1):
                c[word[cut-1]][word[cut]] += 1
    return c


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    files = [HERE/'scout.py', HERE/'check.py']
    pins = {p.name: sha256(p.read_bytes()).hexdigest() for p in files}
    raw = (ROOT/'05-knowledge/results/hn_metagraph_information_20261026/enumeration/representatives9.json.gz').read_bytes()
    check.need(sha256(raw).hexdigest() == BANK_SHA, 'bank pin')
    bank = json.loads(gzip.decompress(raw))['members']
    rows = []
    for n, expected in ((2, 1), (4, 4), (6, 56), (8, 6880)):
        check.need(len(bank[str(n)]) == expected, 'bank cardinality')
        survivors = []
        for item in bank[str(n)]:
            a = check.tournament(n, item['code'])
            s, e, c = check.sidecar(a)
            if all((c[u][v]-c[v][u]) % 4 == 0 for u in range(n) for v in range(n)):
                check.need(n <= 6, 'unexpected survivor needs literal independent audit')
                check.need(c == literal_cover(a), 'literal ordered cover')
                survivors.append(dict(code=item['code'], rows=a, H=sum(e[-1]), S=s[-1], E=e[-1]))
        rows.append(dict(n=n, inherited_classes=expected, symmetric_mod4=len(survivors), survivors=survivors))
    check.need(pins == {p.name: sha256(p.read_bytes()).hexdigest() for p in files}, 'source changed')
    result = dict(status='FINITE-EXACT complete inherited even-order banks through 8',
                  completeness='Inherited authenticated representative bank, not re-proved here',
                  filters='None: no endpoint-parity, strongness or symmetry filter',
                  orders=rows, bank_sha256=BANK_SHA, source_sha256=pins)
    text = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(text, encoding='utf8', newline='\n')
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
