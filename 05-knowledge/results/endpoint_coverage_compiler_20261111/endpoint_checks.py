"""Value content, difference content, and exact one-vertex response support."""
from collections import Counter
from itertools import combinations
from functools import reduce
from math import gcd
from pathlib import Path
import json
import inflation_bank as ib

HERE = Path(__file__).resolve().parent


def main():
    bank, levels, witnesses = [[]], [], []
    literal = 0
    for n in range(1,8):
        codes = {ib.cc.canonical(ib.cc.child(r,m)) for r in bank for m in range(1 << (n-1))}
        bank = [ib.cc.decode(n,c) for c in sorted(codes)]
        tally = Counter()
        for code, rows in zip(sorted(codes), bank):
            h,u,w,_,_,ends = ib.tc.path_tables(rows)
            a = [u[i]-ends[i] for i in range(n)]
            d = reduce(gcd, a + [w[i][j] for i,j in combinations(range(n),2)], 0)
            values = []
            for mask in range(1 << n):
                selected = [v for v in range(n) if mask >> v & 1]
                value = h + sum(a[v] for v in selected) - sum(w[i][j] for i,j in combinations(selected,2))
                values.append(value)
                if n <= 6:
                    ib.tc.require(ib.cc.hp(ib.cc.child(rows,mask)) == value, 'literal response polynomial')
                    literal += 1
            ib.tc.require(d == reduce(gcd, (v-h for v in values),0), 'difference content equality')
            c = reduce(gcd, values,0)
            ib.tc.require(c == gcd(h,d), 'value content equals gcd of anchor and difference content')
            tally[(c,d)] += 1
            if c == 1 and d > 2:
                witnesses.append(dict(n=n,code=code,rows=rows,H=h,C1=c,D1=d,values=sorted(set(values))))
        levels.append(dict(n=n,ordinary=len(bank),content_pairs=[[c,d,count] for (c,d),count in sorted(tally.items())]))
    report = dict(levels=levels,literal_child_controls=literal,gcd_one_with_nonprimitive_difference=witnesses)
    (HERE / 'endpoint_checks.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'levels':levels,'literal_controls':literal,'first_hostiles':witnesses[:2]},indent=2))


if __name__ == '__main__':
    main()
