"""Complete fixed-blue-path probes for singleton-padded strong centers."""
from collections import Counter
from itertools import permutations
import json
from pathlib import Path
from moves import blue, blue_tilings, induced


def strong(a):
    n = len(a)
    full = (1 << n) - 1
    for v in range(n):
        seen = 1 << v
        while True:
            new = seen
            for u in range(n):
                if seen >> u & 1:
                    new |= a[u]
            if new == seen:
                break
            seen = new
        if seen != full:
            return False
    return True


def paths(a, first=None):
    n = len(a)
    full = (1 << n) - 1
    def visit(p, mask):
        if mask == full:
            yield p
            return
        todo = a[p[-1]] & ~mask
        while todo:
            bit = todo & -todo
            todo ^= bit
            v = bit.bit_length() - 1
            yield from visit(p + (v,), mask | bit)
    for u in (range(n) if first is None else [first]):
        yield from visit((u,), 1 << u)


def run():
    report = []
    for n in range(3,10):
        count = Counter()
        hostile = []
        for a in blue_tilings(n):
            count['fixed_path_blue'] += 1
            if not strong(a):
                continue
            count['strong'] += 1
            regular = len({r.bit_count() for r in a}) == 1
            if regular:
                count['regular_strong'] += 1
            if any(not blue(a,p) for p in paths(a,0)):
                count['black_same_start'] += 1
                continue
            start_paths = list(paths(a,0))
            count['all_same_start_blue'] += 1
            witness = next((p for p in paths(a) if not blue(a,p)),None)
            if witness is None:
                count['pure_blue'] += 1
                continue
            count['mixed_hostile'] += 1
            if regular:
                raise RuntimeError('regular mixed fixed-start hostile')
            middle = induced(a,tuple(range(1,n-1)))
            cyclic_blowup = (a[0] == sum(1 << j for j in range(1,n-1))
                            and a[-1] == 1
                            and all(a[j] >> (n-1) & 1 for j in range(1,n-1)))
            reverse_transitive_path = all(
                bool(a[u] >> v & 1) == (v == u+1 or v < u-1)
                for u in range(n) for v in range(n) if u != v)
            if cyclic_blowup:
                if any(not blue(middle,p) for p in paths(middle)):
                    raise RuntimeError('cyclic hostile middle is mixed')
                family = 'cyclic_blowup_pure_blue_middle'
            elif reverse_transitive_path:
                family = 'reverse_transitive_path'
            else:
                raise RuntimeError('new endpoint-hostile family')
            count[family] += 1
            hostile.append(dict(a=a, scores=[r.bit_count() for r in a],
                                starting_paths=start_paths, black_path=witness,
                                family=family))
        record = dict(n=n,counts=dict(count),hostiles=hostile)
        report.append(record)
        print(n,dict(count),flush=True)
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__ == '__main__':
    run()
