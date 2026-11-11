"""Recount every new first point and both uses of all inherited module seeds."""
import json
from pathlib import Path
from construct import construct, ROOT
from verify import H, require

HERE=Path(__file__).resolve().parent
r=json.loads((HERE/'verification.json').read_text())
bank=ROOT/'05-knowledge/results/kps_s21_two_adic_lines_20261110/line_seeds.txt'
modules=[int(l.split()[0]) for l in bank.read_text().splitlines()
         if not l.startswith('#') and l.split()[2]=='module' and int(l.split()[0])<=4097]
cs=sorted({x['c'] for x in r['new_seeds']}|{1,3,5,7,11,13,19,27,1039}|set(modules))
checks=[]
for c in cs:
    x=construct(c);require(H(x['rows'])==x['expected_H'],f'constructor c={c}')
    checks.append([c,x['exponent'],x['order'],x['expected_H']])
for c in modules:
    first=construct(c);x=construct(c,first['start_exponent']+1)
    require(H(x['rows'])==x['expected_H'],f'module inflation c={c}')
    checks.append([c,x['exponent'],x['order'],x['expected_H']])
for c,s in [(3,1),(5,2),(-1,None),(2,None),(4099,None)]:
    try:construct(c,s)
    except ValueError:pass
    else:raise RuntimeError(f'bad request accepted {c},{s}')
(HERE/'constructor_checks.json').write_text(json.dumps(dict(checked_points=checks,count=len(checks),
    rejected_requests=[[3,1],[5,2],[-1,None],[2,None],[4099,None]]),indent=2)+'\n')
print(f'PASS {len(checks)} literal constructor recounts and five rejected requests')
