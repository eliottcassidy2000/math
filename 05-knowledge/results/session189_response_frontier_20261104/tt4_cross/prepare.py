import itertools,json,pathlib
P=pathlib.Path(__file__).resolve().parent
ORDERS=[(0,1,2,3),(1,0,3,2),(1,3,0,2),(2,0,3,1),(2,3,0,1)]
COMBOS=[[(0,1,3,1),(0,1,2,-1)],[(0,1,3,2),(0,1,2,-1)],[(0,1,3,2),(0,1,2,3)],[(0,2,3,1),(0,2,1,-1)]]
CONSTANTS=[2,2,-2,1]
def tt(order):
    n=len(order);r=[0]*n
    for i,a in enumerate(order):
        for b in order[i+1:]:r[a]|=1<<b
    return r
def conv(r):return[((1<<len(r))-1)^(1<<i)^x for i,x in enumerate(r)]
def inflate(Q,B):
    off=[];n=0
    for b in B:off.append(n);n+=len(b)
    r=[0]*n
    for i,b in enumerate(B):
        for u,x in enumerate(b):
            r[off[i]+u]=x<<off[i]
            for j,c in enumerate(B):
                if Q[i]>>j&1:r[off[i]+u]|=((1<<len(c))-1)<<off[j]
    return r
def ambient(r):
    n=len(r);out=[x|(1<<n)for x in r]+[1<<(n+1),(1<<n)-1]
    return out
def kernels(order,p,q,r):
    pos={v:i for i,v in enumerate(order)};L=[0,0,0];R=[0,0,0]
    for perm in itertools.permutations(range(4)):
        for cut in range(1,4):
            a,b=perm[:cut],perm[cut:];z=[1]
            for u,v in list(zip(a,a[1:]))+list(zip(b,b[1:])):
                edge=int(pos[u]<pos[v]);zz=[0]*(len(z)+1)
                for k,w in enumerate(z):zz[k]+=edge*w;zz[k+1]+=w
                z=zz
            if (a[-1]==p and b[0]==q and b[-1]==r)or(a[-1]==q and b[0]==p and b[-1]==r):L=[x+y for x,y in zip(L,z)]
            if (a[0]==r and a[-1]==p and b[0]==q)or(a[0]==r and a[-1]==q and b[0]==p):R=[x+y for x,y in zip(R,z)]
    return L,R
def main():
    keys=[(p,q,r)for p,q in itertools.combinations(range(4),2)for r in range(4)if r not in(p,q)]
    inner=[]
    for order in ORDERS:
        inner.append(dict(order=order,coefficients=[dict(p=p,q=q,r=r,L=kernels(order,p,q,r)[0],R=kernels(order,p,q,r)[1])for p,q,r in keys]))
    (P/'inner_kernels.json').write_text(json.dumps(inner,indent=2)+'\n')
    specs=[dict(Q=[2,4,1],others=[1,1],ambient=False),dict(Q=[2,4,1],others=[2,2],ambient=False),dict(Q=[2,4,1],others=[2,2],ambient=True),dict(Q=[6,12,24,17,3],others=[1,1,1,1],ambient=False)]
    graphs=[];cases=[]
    for spec in specs:
        Q=spec['Q'];others=[tt(tuple(range(n)))for n in spec['others']];n=4+sum(spec['others']);case=dict(spec,core_order=n,source=len(graphs),targets=[],local_replacements=[],group=list(range(4,n)),virtual=inflate(Q,[[0]]+others))
        source=inflate(Q,[tt(ORDERS[0])]+others);graphs.append(ambient(source)if spec['ambient']else source)
        for order in ORDERS[1:]:
            target=inflate(conv(Q),[tt(order)]+others);case['targets'].append(len(graphs));graphs.append(ambient(target)if spec['ambient']else target)
            local=inflate(Q,[tt(tuple(reversed(order)))]+others);case['local_replacements'].append(len(graphs));graphs.append(ambient(local)if spec['ambient']else local)
        cases.append(case)
    (P/'controls.json').write_text(json.dumps(dict(source_order=ORDERS[0],crossed_targets=ORDERS[1:],combinations=COMBOS,constants=CONSTANTS,cases=cases),indent=2)+'\n')
    (P/'graphs.txt').write_text(str(len(graphs))+'\n'+'\n'.join(str(len(r))+' '+' '.join(map(str,r))for r in graphs)+'\n')
if __name__=='__main__':main()
