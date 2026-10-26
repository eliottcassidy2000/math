import sys as _output_sys
_output_sys.stdout.reconfigure(newline="\n")
from pathlib import Path
import itertools, json, subprocess, sys
from math import gcd

p=Path(__file__).resolve().parent
def need(c,m):
    if not c: raise ValueError(m)

rows=[]
for m in range(2,16):
    total=good=0
    for a,b,c,d in itertools.product(range(m),repeat=4):
        x=(a*d-b*c)%m
        if gcd(x,m)==1:
            total+=1
            good+=x in {1,m-1}
    phi=sum(gcd(a,m)==1 for a in range(m))
    need(good*phi==total*len({1,m-1}),'determinant proportion')
    rows.append([m,total,good])
need(next(r for r in rows if r[0]==7)==[7,2016,672],'seventh root proportion')
for m,q,qi in ((7,((2,7),(7,25)),((25,-7),(-7,2))),
               (5,((2,5),(5,13)),((13,-5),(-5,2)))):
    need([[sum(q[i][k]*qi[k][j] for k in range(2)) for j in range(2)] for i in range(2)]==[[1,0],[0,1]],'displayed Laurent inverse')
    need(q[0][0]%m==2 and q[1][1]%m==pow(2,-1,m),'compensating character')
replay=subprocess.check_output([sys.executable,'-B',str(p/'check.py')])
need(json.loads(replay)==json.loads((p/'normal.json').read_text()),'producer replay')
need((p/'normal.json').read_bytes()==(p/'optimized.json').read_bytes(),'normal optimized bytes')
(p/'parent_audit.json').write_text(json.dumps({'status':'PASS','independent_rank2_determinant_counts':rows,'producer_replay_matches':True,'normal_optimized_bytes_equal':True},indent=2)+'\n', newline='\n')
print('PASS: F1 independent determinant census, exact inverse identities and producer replay')
