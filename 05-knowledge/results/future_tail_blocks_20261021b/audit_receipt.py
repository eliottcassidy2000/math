"""Independent literal-pair envelope consumer; not a substitute for source replay."""
from itertools import combinations
from pathlib import Path
import json


def q(a, b):
    return a*a - a*b + b*b


def f(v):
    return q(v[0]+v[2], v[1]+v[3]) - q(v[0], v[1])


def k(a, b):
    return f(tuple(x+y for x, y in zip(a, b))) - f(a) - f(b)


def main():
    data = json.loads(Path(__file__).with_name('receipt.json').read_text())
    for case in data['cases']:
        ports = [tuple(row[1:5]) for row in case['arrived_ports']]
        pending = [tuple(row[1:5]) for row in case['pending']]
        current = sum(map(f, ports))
        local = sum(max(0, f(b)+max([0]+[abs(k(a, b)) for a in ports]))
                    for b in pending)
        pairs = sum(abs(k(a, b)) for a, b in combinations(pending, 2))
        actual = [current, local, pairs, current+local+pairs]
        if actual != case['envelope']:
            raise ValueError(('literal-pair envelope disagreement', case['N'], actual))
        print(json.dumps({'N': case['N'], 'D': case['D'], 'envelope': actual,
                          'literal_pending_pairs': len(pending)*(len(pending)-1)//2}))


if __name__ == '__main__':
    main()
