"""All-height Pareto optimization of the actual Berggren substitution lift.
Transfer matrices and reverse-address convention: ../spinor/PROOF.md §§2–3.
Independent standalone implementation; no external packages or depth cutoff.
"""
from collections import deque
from pathlib import Path
import json
U=((1,0,0,0),(0,2,0,1),(0,0,0,1),(0,0,0,1))
V=((1,1,1,0),(0,2,0,0),(0,0,0,0),(0,0,0,1))
def mul(a,b):return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(4)) for j in range(4)) for i in range(4))
def pack(x,y):return tuple(z for m in (x,y) for row in m for z in row)
def unpack(v):return tuple(tuple(v[i:i+4]) for i in range(0,16,4)),tuple(tuple(v[i:i+4]) for i in range(16,32,4))
def step(v,a):
    x,y=unpack(v)
    if a=='A':xx,yy=mul(y,x),y
    elif a=='B':xx,yy=mul(y,x),mul(mul(x,y),x)
    else:xx,yy=x,mul(mul(y,x),x)
    return pack(xx,yy)
def output(v):
    x,y=unpack(v);m=mul(y,x)
    return m[0][0]+m[0][3]
def dominates(a,b):return all(x<=y for x,y in zip(a,b))
def optimize(modulus,controller=None,certificate_path=None):
    if modulus<1:raise ValueError('modulus must be positive')
    unrestricted=controller is None
    if controller is None:controller=(0,lambda c,a:0,lambda c:True)
    init,tr,accept=controller;seed=pack(U,V)
    def key(c,v):return c,tuple(x%modulus for x in v)
    bank={key(init,seed):{seed:''}};queue=deque([(init,seed)]);added=1
    while queue:
        c,v=queue.popleft();bucket=bank[key(c,v)]
        if v not in bucket:continue
        w=bucket[v]
        for a in 'ABC':
            cc=tr(c,a)
            if cc is None:continue
            z=step(v,a);target=bank.setdefault(key(cc,z),{})
            if any(dominates(old,z) for old in target):continue
            for old in list(target):
                if dominates(z,old):del target[old]
            target[z]=w+a;queue.append((cc,z));added+=1
    retained=[(v,w) for (c,r),b in bank.items() if accept(c) for v,w in b.items()]
    stats={'modulus':modulus,'labels':len(bank),'retained':sum(map(len,bank.values())),
           'accepted':len(retained),'additions':added,'max_reverse_address_length':max((len(w) for v,w in retained),default=0)}
    if certificate_path is not None:
        if not unrestricted:raise ValueError('certificate export currently requires unrestricted controller')
        certificate={'scope':'FINITE-EXACT COMPLETE SEMIGROUP CERTIFICATE',
                     'modulus':modulus,'statistics':stats,'depth_or_size_cap':None,
                     'representatives':[{'vector':list(v),'reverse_address':w} for v,w in retained]}
        Path(certificate_path).write_text(json.dumps(certificate,separators=(',',':'))+'\n')
    return retained,stats
SUB={'A':{'O':'OI','I':'I'},'B':{'O':'OI','I':'OIO'},'C':{'O':'O','I':'OOI'}}
def literal(addr):
    w='OI';m,n=2,1
    for a in addr:
        w=''.join(SUB[a][x] for x in w)
        m,n={'A':(2*m-n,m),'B':(2*m+n,m),'C':(m+2*n,n)}[a]
    h,p,e=1,0,0
    for x in w:
        if x=='O':p,e=2*p+1,1
        else:h,p,e=h+p+e,2*p,0
    rows=[sum(1<<j for j in range(i+1,len(w))) for i in range(len(w))]+[0]
    for i,x in enumerate(w):rows[len(w) if x=='O' else i]|=1<<(i if x=='O' else len(w))
    return {'chronological_address':addr,'cut_word':w,'parameters':[m,n],
            'triple':[m*m-n*n,2*m*n,m*m+n*n],'H':h,'rows':rows}
def summarize(retained,modulus):
    best={}
    for v,w in retained:
        h=output(v);r=h%modulus;candidate=(h,len(w),w)
        if r not in best or candidate<best[r]:best[r]=candidate
    ans=[]
    for r,(h,k,w) in sorted(best.items()):
        rec=literal(w[::-1]);assert rec['H']==h
        rec.update(residue=r,reverse_address=w)
        ans.append(rec)
    return ans

def run():
    experiments=[]
    for q in (2,3,4,8,16):
        records,stats=optimize(q)
        experiments.append({'statistics':stats,'minimum_by_residue':summarize(records,q)})
    # Nontrivial regular restriction: B occurs at least once.
    # This property is reversal-invariant; controller acts on reversed addresses.
    records,stats=optimize(4,(False,lambda seen,a:seen or a=='B',bool))
    constrained={'statistics':stats,'minimum_by_residue':summarize(records,4)}
    inverses=[]
    for target in (1,2,3,5,7,9,13,15):
        # Any modulus greater than the target works; powers of two are
        # substantially smaller reachable quotients for these examples.
        modulus=1<<target.bit_length()
        records,stats=optimize(modulus)
        feasible=[(output(v),len(w),w) for v,w in records if output(v)%modulus==target]
        best=min(feasible,default=None)
        inverses.append({'target':target,'statistics':stats,'exact_member':best is not None and best[0]==target,
                         'fibre_minimum':best[0] if best else None,
                         'minimizer':literal(best[2][::-1]) if best else None})
    from engine import H
    checked={}
    for group in experiments+[constrained]:
        for rec in group['minimum_by_residue']:
            if len(rec['rows'])<=13 and rec['chronological_address'] not in checked:
                assert H(tuple(rec['rows']))==rec['H']
                checked[rec['chronological_address']]=rec['H']
    result={'direct_graph_controls':checked,'experiments':experiments,'contains_B':constrained,'exact_inverse':inverses}
    Path(__file__).with_name('morphic_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='experiments'},indent=2))
if __name__=='__main__':run()
