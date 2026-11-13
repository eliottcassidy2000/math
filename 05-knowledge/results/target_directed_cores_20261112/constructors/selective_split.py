"""Exact selective replacement, with literal graphs and symbolic compiler control."""
from pathlib import Path
from itertools import permutations
import json
import subprocess
from paired_scale_probe import evaluate

HERE = Path(__file__).resolve().parent
COMPILER = '05-knowledge/results/session296_multicore_metagraph_20261113/cores/verify_ready.py'

def compiler():
    ns = {'__name__': 'inherited_funnel_compiler'}
    exec(compile(subprocess.check_output(['git', 'show', 'HEAD:' + COMPILER], text=True), COMPILER, 'exec'), ns)
    return ns

def first_funnel(T, word):
    # Root0 is retained as vertex1; fresh root is vertex0, followed by old core, then tail.
    n = len(T); t = len(word); rows = [0] * (n+t+1)
    rows[0] = 2 | sum(1 << (n+1+j) for j,b in enumerate(word) if b == '0')
    tail = ((1 << (n+t+1))-1) ^ ((1 << (n+1))-1)
    for i,row in enumerate(T):
        rows[i+1] = (row << 1) | tail | (1 if i else 0)
    for j,b in enumerate(word):
        rows[n+1+j] = (((1 << (n+t+1))-1) ^ ((1 << (n+j+2))-1)) | (1 if b == '1' else 0)
    return tuple(rows)

def child(T):
    n=len(T)
    return tuple(row | ((1 << n) if i==0 else 0) for i,row in enumerate(T)) + (((1 << n)-1)^1,)

def transitive(m):
    return tuple(((1 << m)-1) ^ ((1 << (i+1))-1) for i in range(m))

def sources(m):
    old = child(first_funnel(first_funnel(first_funnel((0,), '1'*m), '1'), '11'))
    new = child(first_funnel(first_funnel(first_funnel((0,), '1'*(m+1)), '1'), '01'))
    mapping = {i:i for i in range(m+4)}
    mapping.update({m+4:m+5,m+6:m+7,m+7:m+8})
    return old,new,mapping

def cover_bins(T,j):
    # Ordered bins, each carrying one path or empty; exponential only for tiny controls.
    n=len(T); hp=[0]*(1<<n);hp[0]=1
    for mask in range(1,1<<n):
        vs=[i for i in range(n) if mask>>i&1]
        hp[mask]=sum(all(T[a]>>b&1 for a,b in zip(p,p[1:])) for p in permutations(vs))
    dp={0:1};full=(1<<n)-1
    for _ in range(j):
        nd={}
        for mask,count in dp.items():
            rem=full^mask;s=rem
            while True:
                nd[mask|s]=nd.get(mask|s,0)+count*hp[s]
                if s==0:break
                s=(s-1)&rem
        dp=nd
    return dp[full]

def inner(M):
    # root0 -> singleton1 -> M -> root0.
    n=len(M)
    return (2, ((1 << (n+2))-1)^3) + tuple((r<<2)|1 for r in M)

def arbitrary_sources(M):
    old=child(first_funnel(first_funnel(inner(M),'1'),'11'))
    n=len(M);plus=tuple(r|(1<<n) for r in M)+(0,)
    new=child(first_funnel(first_funnel(inner(plus),'1'),'01'))
    return old,new

def main():
    import sympy as z
    ns=compiler();u,v,w=z.symbols('u v w')
    S=(u+v,u,u+v,2*w+4*v-2*u,2*v-u,w+2*v)
    symbolic={}
    for name,base,words in [('old',S,('1','11')),('new',tuple(x.subs({v:2*v,w:3*w},simultaneous=True) for x in S),('1','01'))]:
        for word in words:base=ns['advance'](base,word)
        symbolic[name]=dict(parentH=str(z.expand(base[0])),H=str(z.expand(base[0]+base[4]-u)))
    cases=[];data=[]
    for m in range(1,6):
        old,new,phi=sources(m)
        assert all(((old[i]>>j)&1)==((new[phi[i]]>>phi[j])&1) for i in phi for j in phi)
        for T in (old,new):assert [i for i in range(len(T)-1) if T[i]>>(len(T)-1)&1]==[0]
        cases.extend([(old,len(old)-1),(new,len(new)-1)])
        data.append(dict(m=m,old_rows=old,new_rows=new,retained_map=phi,deleted=[m+5],inserted=[m+4,m+6]))
    vals=iter(evaluate(cases))
    for row in data:
        m=row['m'];a=next(vals);b=next(vals)
        assert a==(105*2**m+52*3**m-7,15+15*2**m,14+15*2**m)
        assert b==(166*2**m+114*3**m+1,17+34*2**m,16+34*2**m)
        row['old_profile_H_q_s']=a;row['new_profile_H_q_s']=b
        row['doubling_defect']=b[0]-(2*a[0]-71)
    modules=[(0,),transitive(2),(2,4,1),(2,4,9,3)]
    hostile=[]
    for M in modules:
        uvw=tuple(cover_bins(M,j) for j in (1,2,3));u,v,w=uvw
        old,new=arbitrary_sources(M)
        a,b=evaluate([(old,len(old)-1),(new,len(new)-1)])
        assert a==(-7*u+105*v+52*w,15*u+15*v,14*u+15*v)
        assert b==(u+166*v+114*w,17*u+34*v,16*u+34*v)
        hostile.append(dict(module=M,bins123=uvw,old_profile=a,new_profile=b))
    out=dict(status='PASS',proof='Universal run/bin expansion plus inherited exact six-coordinate compiler',symbolic=symbolic,literal=data,arbitrary_module_controls=hostile)
    (HERE/'selective_split_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'PASS','symbolic':symbolic,'literal_controls':len(data),'module_controls':len(hostile)}))

if __name__=='__main__':main()
