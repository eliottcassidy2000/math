from pathlib import Path
from collections import Counter
from math import comb,gcd
import json,re,ast,csv
P=Path(__file__).resolve().parent
checks=0
def ck(v):
 global checks
 if not v: raise RuntimeError('check failed')
 checks+=1
inv=json.loads((P/'inventory_375.json').read_text())['orders']
expected={17:{6:2},31:{6:2,7:2},89:{7:14,8:7,9:2}}
for h,want in expected.items():
 got={int(n):sum(row[1]==h for row in rows) for n,rows in inv.items()}
 ck({n:c for n,c in got.items() if c}==want)
cs={int(n):sum(row[1]==89 and row[2] for row in rows) for n,rows in inv.items()}
ck({n:c for n,c in cs.items() if c}=={8:1})
a=ast.literal_eval(re.search(r'P_d\(1\).*?: (\[.*\])',(P/'strong_ears_625.out').read_text())[1])
ck(len(a)==125 and a[(527-375)//2]==184)
rows=list(csv.DictReader((P/'n8_classes.csv').open()))
ck(sum(int(r['H'])==527 and int(r['strong']) for r in rows)==22)
h7=ast.literal_eval(re.search(r'n = 7: N\(n,h\) = (\{[^\n]*\})',(P/'h_multiplicities.out').read_text())[1])
ck(max(h7)==189)
for n in range(11,180):
 pair=sum(8*comb(n-m+2,2) if n>=m else 0 for m in (12,13))
 ck(pair==8*(n-11)**2)
 for m in range(8,13):ck(max(0,n-m+1)==n-m+1)
 for m in range(9,13):
  W=22*8+162*m
  direct=22*max(0,n-7)+162*max(0,n-m+1)+pair
  ck(direct==8*n*n+8*n+1152-W)
  ck((8*(n+9)**2+8*(n+9)+1152-W-direct)%72==0)
for n in range(8,180):
 direct=sum(c*max(0,n-m+1) for m,c in expected[89].items())
 ck(direct==23*n-149)
 ck((direct+(n%2==0))%2==0)
ck(23*7-149==12 and sum(c*max(0,7-m+1) for m,c in expected[89].items())==14)
for d in range(1,73):
 works=all((16*d*n+8*d*(d+1))%72==0 for n in range(11,83))
 ck(works==(d%9==0))
 ck(((23*d)%72==0)==(d%72==0))
ck(gcd(23,72)==1)
print(json.dumps({'status':'PASS','checks':checks,'strong_counts':expected,'S527':184,'W527':'unrecovered','N527':'8*n*n+8*n+1152-W527; n>=11','N89':'23*n-149; n>=8','periods_mod72':{'527':9,'89':72}},indent=2))
