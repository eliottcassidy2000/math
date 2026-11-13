"""Literal repeated returns, selected escapes, and exact inverse/spectral checks."""
import json,gzip
import sympy as sp
import numpy as np
from typed_glue import HERE,save,need,glue,diff,endpoint_type
from feedback_probe import all_profiles
from port_compiler import mul
from repair23 import strong


def read(name):return json.loads(gzip.decompress((HERE/(name+'.json.gz')).read_bytes()))


def run():
    rules=read('port_search')['rules'];pairs=read('typed_glue')['pairs']
    choice=[('2/3>1/3','1/3','1/2'),('5/7>4/9','4/9','1'),('7/10>6/11','6/11','1/6')]
    escapes=[]
    for label,a,b in choice:
        p=next(x for x in pairs if x['label']==label);c=next(x for x in rules if x['input']==a and x['output']==b)
        rows=[glue(c['rows'],c['entry'],T,p['root']) for T in p['rows']]
        out=[tuple(map(int,all_profiles(np.array(T,dtype=np.int64))[c['exit']])) for T in rows]
        predicted=[mul(c['matrix'],tuple(z[j] for j in (0,1,3,4))) for z in p['profiles']]
        need(predicted==[tuple(z[j] for j in (0,1,3,4)) for z in out],'literal changed-port escape')
        need('E' in endpoint_type(diff(*out)),'actual escaped neutral pair')
        escapes.append(dict(source=label,rule=[a,b],rows=rows,root=c['exit'],profiles=out,orders=[len(T) for T in rows],strong=[strong(T) for T in rows]))
        print('ACTUAL ESCAPE',label,'>',b,'H',out[0][0],'ORDER',len(rows[0]),flush=True)
    p=next(x for x in pairs if x['label']=='2/3>1/3');rows=p['rows'];r=3
    c=next(x for x in rules if x['input']==x['output']=='1/2');M=c['matrix'];matrix=sp.Matrix(M)
    zs=[tuple(int(all_profiles(np.array(T,dtype=np.int64))[r,j]) for j in (0,1,3,4)) for T in rows]
    initial=zs;literal=[];arith=[]
    for k in range(25):
        need(zs[0][0]==zs[1][0] and zs[0][3]-zs[0][1]==zs[1][3]-zs[1][1],'renewed observations')
        need(tuple(y-x for x,y in zip(*zs))==tuple((-1)**k*(-6)*x for x in (0,1,2,1)),'repair eigenmode')
        need(all(z[3]+z[1]==z[0]+z[2] for z in zs),'indegree-one plane')
        h=zs[0][0];E=h+zs[0][3]-zs[0][1];arith.append(dict(k=k,n=11+3*k,H=h,E=E,profiles=zs))
        if k<2:
            rows=[glue(c['rows'],c['entry'],T,r) for T in rows];r=c['exit']
            actual=[tuple(int(all_profiles(np.array(T,dtype=np.int64))[r,j]) for j in (0,1,3,4)) for T in rows]
            need(actual==[mul(M,z) for z in zs],'two literal renewals')
            literal.append(dict(k=k+1,rows=rows,root=r,profiles=actual,strong=[strong(T) for T in rows]))
        zs=[mul(M,z) for z in zs]
    for a,b,c0 in zip(arith,arith[1:],arith[2:]):
        for key in ('H','E'):need(c0[key]==8*b[key]+a[key],'scalar recurrence')
    for a,b in zip(arith,arith[1:]):
        need((4*b['E']-7*b['H'],2*b['H']-b['E'])==(a['H'],a['E']),'exact target inverse')
    lam=sp.Symbol('lambda');poly=lam*(lam+1)*(lam*lam-8*lam-1)
    need(sp.expand(matrix.charpoly(lam).as_expr()-poly)==0,'full characteristic polynomial')
    save('renewal_verify',dict(status='FINITE-EXACT and algebraic certificate',scope='Three selected changed-port escapes literally recounted at orders16,17,18. Two actual repeated returns at14 and17 vertices. Twenty-five symbolic iterations are examples of the proved universal rule, not an all-target claim.',escapes=escapes,initial_profiles=initial,literal_renewals=literal,iterations=arith,characteristic_polynomial=str(poly)))
    print('RETURN HEAD',[(x['H'],x['E']) for x in arith[:5]],flush=True)


if __name__=='__main__':run()
