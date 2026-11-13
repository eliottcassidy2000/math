"""Independent exhaustive small-context audit of the root's module candidate."""
from pathlib import Path
from fractions import Fraction
import json
HERE=Path(__file__).resolve().parent
checks=0

def require(ok,msg):
    global checks
    checks+=1
    if not ok: raise RuntimeError(msg)

def tourney(n,code):
    rows=[0]*n
    k=0
    for i in range(n):
        for j in range(i+1,n):
            if (code>>k)&1: rows[i]|=1<<j
            else: rows[j]|=1<<i
            k+=1
    return rows

def context(ext,mask,module):
    # External vertices precede module labels; bit i means M -> external i.
    m,k=len(ext),len(module)
    rows=ext[:]+[x<<m for x in module]
    for i in range(m):
        for j in range(k):
            if (mask>>i)&1: rows[m+j]|=1<<i
            else: rows[i]|=1<<(m+j)
    return rows

def count(t,start=None,end=None):
    n=len(t)
    states=[{} for _ in range(1<<n)]
    for v in (range(n) if start is None else [start]): states[1<<v][v]=1
    for mask,d in enumerate(states):
        for last,c in d.items():
            choices=t[last]&~mask
            while choices:
                b=choices&-choices
                v=b.bit_length()-1
                target=states[mask|b]
                target[v]=target.get(v,0)+c
                choices-=b
    return sum(states[-1].values()) if end is None else states[-1].get(end,0)

def sub(t,mask):
    vs=[i for i in range(len(t)) if (mask>>i)&1]
    return [sum(((t[u]>>v)&1)<<j for j,v in enumerate(vs)) for u in vs]

def p2(t):
    allmask=(1<<len(t))-1
    h=lambda mask:count(sub(t,mask)) if mask else 1
    return sum(h(mask)*h(allmask^mask) for mask in range(allmask+1))

def tt(n):return [sum(1<<j for j in range(i+1,n)) for i in range(n)]

modules=[tt(4),[2,4,1],[6,12,8,1],[22,4,8,3,14]]
for t in modules:
    require(all(not ((t[i]>>i)&1) for i in range(len(t))),'module loop-free')
    require(all(((t[i]>>j)&1)+((t[j]>>i)&1)==1 for i in range(len(t)) for j in range(i+1,len(t))),'module tournament type')
module_states=[(count(t),p2(t)) for t in modules]
classes={}
max_E=0
for m in range(5):
    num_contexts=0
    all_observers=zero_observers=positive_observers=0
    for code in range(1<<(m*(m-1)//2)):
        ext=tourney(m,code)
        for mask in range(1<<m):
            num_contexts+=1
            observers=[(None,None)]
            if m:observers.append((0,None))
            if m>=2:observers.append((0,1))
            vals=[]
            for k in [1,2,3]:
                t=context(ext,mask,tt(k))
                vals.append([count(t,start,end) for start,end in observers])
            for oi,(start,end) in enumerate(observers):
                all_observers+=1
                f1,f2,f3=[v[oi] for v in vals]
                require((f2-f1)%2==0 and f2>=f1,'A2 nonnegative integral')
                E=f3-3*f2+2*f1
                require(E>=0 and E%6==0,'A3 certificate')
                require(E%12==0,'internal separator permutation divisibility')
                s=int(start is not None)+int(end is not None)
                require(m>=4+s or E==0,'endpoint-refined third-cover cutoff')
                max_E=max(max_E,E)
                if m<4:require(E==0,'two outside vertices per separator')
                if E==0:
                    zero_observers+=1
                    A2=(f2-f1)//2
                    for mod,(h,q) in zip(modules,module_states):
                        actual=count(context(ext,mask,mod),start,end)
                        predicted=(f1-2*A2)*h+A2*q
                        require(actual==predicted,'three-size test arbitrary module')
                else:
                    positive_observers+=1
                    f4=count(context(ext,mask,tt(4)),start,end)
                    fit=f1+(f2-f1)//2*(2**4-2)
                    require(f4-fit>=E//6*(3**4-3*2**4+3),'lower discrepancy')
    classes[m]={'contexts':num_contexts,'observers':all_observers,'zero_E':zero_observers,'positive_E':positive_observers}
receipt={'status':'FINITE-EXACT independent audit; no proof inferred from census','checks':checks,'context_universe':'all labelled tournaments on m=0..4 external vertices, all uniform module attachments','observers':'full H; fixed external start 0; fixed external start/end 0,1 when present','arbitrary_module_controls':[{'masks':t,'H':h,'P2':q} for t,(h,q) in zip(modules,module_states)],'counts':classes,'maximum_E':max_E}
(HERE/'module_audit.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'checks':checks,'status':'PASS','counts':classes,'max_E':max_E}))
