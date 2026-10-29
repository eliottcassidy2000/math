"""A second computational path: prime-power Smith pivots, not inverses.

Reads the exact saved matrices and checks invariant valuations independently
of the producer's determinant/exponent/rank reconstruction. Same agent.
"""
from pathlib import Path
import hashlib, json, subprocess, sys

ROOT=Path(__file__).resolve().parent
checks=0
def need(ok,msg):
    global checks
    checks+=1
    if not ok: raise RuntimeError(msg)

def valuations(a,p,precision=3):
    modulus=p**precision
    b=[[x%modulus for x in row] for row in a]; n=len(b); out=[]
    def val(x):
        if not x:return precision
        k=0
        while x%p==0:x//=p;k+=1
        return k
    for k in range(n):
        z,i,j=min((val(b[i][j]),i,j) for i in range(k,n) for j in range(k,n))
        if z==precision:return out+[precision]*(n-k)
        b[i],b[k]=b[k],b[i]
        for row in b: row[j],row[k]=row[k],row[j]
        unit=b[k][k]//p**z; inv=pow(unit,-1,p**(precision-z))
        for i in range(k+1,n):
            need(b[i][k]%p**z==0,'pivot divisibility')
            t=(b[i][k]//p**z)*inv
            b[i]=[(x-t*y)%modulus for x,y in zip(b[i],b[k])]
        for j in range(k+1,n):
            need(b[k][j]%p**z==0,'row pivot divisibility')
            t=(b[k][j]//p**z)*inv
            for i in range(n):b[i][j]=(b[i][j]-t*b[i][k])%modulus
        out.append(z)
    return out

def run():
    normal=subprocess.check_output([sys.executable,'-B',str(ROOT/'check.py')])
    optimized=subprocess.check_output([sys.executable,'-O','-B',str(ROOT/'check.py')])
    need(normal==optimized,'normal and optimized byte equality')
    data=json.loads(normal)
    for name in ['check.json','check.optimized.json']:
        need(json.loads((ROOT/name).read_text(encoding='utf-8-sig'))==data,'saved receipt '+name)
    rows=[]
    for s in data['seeds']:
        m=s['matrix']; n=len(m); p=s['p']
        for r in range(1,9):
            k=[[m[i%n][j%n] if i//n==j//n else
                (1 if i//n<j//n else -1) for j in range(n*r)] for i in range(n*r)]
            got=valuations(k,p)
            expected=sorted([1]*(2*(r%2))+[2]*(r-r%2)+[0]*(n*r-r-r%2))
            need(got==expected,('Smith layers',p,r))
            rows.append({'p':p,'copies':r,'valuation_multiplicities':
                         {str(e):got.count(e) for e in range(3)}})
    for s in data['two_channel_attachments']:
        m=s['base_matrix']; b=s['cut']; p=s['pfaffian']
        k=[m[i]+b[i] for i in range(4)]+[[-b[j][i] for j in range(4)]+m[i] for i in range(4)]
        need(valuations(k,p)==[0]*6+[1]*2,'two-channel Smith')
    # Independent symbolic polynomial arithmetic for the first-jet claim.
    def times(a,b):
        c=[0]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b):c[i+j]+=x*y
        return c
    for n in range(1,13):
        l=[1]; o=[2]
        for i in range(1,n+1): l=times(l,[1]+[0]*(i-1)+[1])
        for i in range(1,n): o=times(o,[1]+[0]*(i-1)+[1])
        need(sum(l)==sum(o)==2**n,'thin value')
        need(sum(i*x for i,x in enumerate(l))-sum(i*x for i,x in enumerate(o))==n*2**(n-1),'first jet')
    sources={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'check.py',ROOT/'audit.py']}
    return {'status':'PASS','checks':checks,'producer_checks':data['checks'],
            'normal_optimized_identical':True,'prime_power_precision':3,
            'smith_rows':rows,'source_sha256':sources}

if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
