"""Independent review: no producer audit helpers or subset-DP recurrence."""
from itertools import permutations, product
from pathlib import Path
import hashlib, importlib.util, json
S=Path(__file__).resolve().parent.parent/'session26-counter'
spec=importlib.util.spec_from_file_location('target',S/'counter.py'); t=importlib.util.module_from_spec(spec); spec.loader.exec_module(t)
def require(b, msg):
    if not b: raise RuntimeError(msg)
def core(k,z):
    q=[[0]*k for _ in range(k)]
    for v in range(1,k):
        for u in range(v): q[u][v]=z&1; q[v][u]=1-(z&1); z>>=1
    return q
def edge(q,w,u,v):
    k=len(q)
    if u<k and v<k: return q[u][v]
    if u>=k and v>=k: return u<v
    if u<k: return (w[v-k]>>u)&1
    return 1-((w[u-k]>>v)&1)
def oracle(q,w):
    # Enumerate directed vertex words by recursive extension, no memoization.
    n=len(q)+len(w)
    def extend(last,unused):
        if not unused:return 1
        return sum(extend(v,unused-{v}) for v in unused if last is None or edge(q,w,last,v))
    return extend(None,set(range(n)))
def relabel(q,w,p):
    return [[q[u][v] for v in p] for u in p], [sum(((a>>u)&1)<<j for j,u in enumerate(p)) for a in w]
def observed_state(slots,q,w):
    k=len(q)
    ds=[0 if not sl else 2 if w[sl[-1]]>>j&1 else 1 for j,sl in enumerate(slots[:-1])]
    return sum(d*3**j for j,d in enumerate(ds))+(3**k if slots[-1] else 0)
def slot_mass(q,w):
    # Retain explicit partial slot assignments; verify first ports directly.
    k=len(q); mass={}
    for assignment in product(range(k+1),repeat=len(w)):
        slots=[[i for i,s in enumerate(assignment) if s==j] for j in range(k+1)]
        if any(slots[j] and not (w[slots[j][0]]>>(j-1)&1) for j in range(1,k+1)):continue
        st=observed_state(slots,q,w);mass[st]=mass.get(st,0)+1
    return mass
def machine_mass(k,w):
    mass={0:1}
    for a in w:
        nxt={}
        for st,m in mass.items():
            for dst in t.transitions(k)[a][st]:nxt[dst]=nxt.get(dst,0)+m
        mass=nxt
    return mass
out={}; counts=0; covariance=0; duality=0
for k,lim in [(0,6),(1,6),(2,4),(3,3),(4,2)]:
    for z in range(1<<(k*(k-1)//2)):
        q=core(k,z)
        for m in range(lim+1):
            for w in product(range(1<<k),repeat=m):
                h=t.count_hamilton(q,w)
                require(h==oracle(q,w),('oracle',q,w,h));counts+=1
                if k<=3 and m<=2:
                    for p in permutations(range(k)):
                        qq,ww=relabel(q,w,p)
                        require(t.count_hamilton(qq,ww)==h,('covariance',q,w,p));covariance+=1
                dq=[[q[v][u] for v in range(k)] for u in range(k)]
                dw=[a^((1<<k)-1) for a in w[::-1]]
                require(t.count_hamilton(dq,dw)==h,('duality',q,w));duality+=1
out.update(backtracking_oracle_cases=counts,covariance_cases=covariance,duality_cases=duality)
state_cases=0
for k in range(1,4):
    q=core(k,0)
    for m in range(4):
        for w in product(range(1<<k),repeat=m):
            require(slot_mass(q,w)==machine_mass(k,w),('partial_state_mass',k,w));state_cases+=1
out['partial_state_mass_cases']=state_cases
# Direct boundary/multiplicity checks, independent of acceptance totals.
require(t.transitions(1)[0][2]==(1,), 'bad exit must be repairable')
require(t.transitions(1)[0][4]==(4,4),'two distinct slot choices must retain multiplicity')
require(t.transitions(2)[1][0].count(3)==1,'interior singleton first=1 last=0')
require(3 not in t.transitions(2)[0][0], 'reject wrong first port')
require(t.transitions(2)[3][0].count(6)==1,'retain bad exit on valid first port')
require(t.count_hamilton([[0,1],[0,0]],[])==1,'empty internal gap core arc')
require(t.count_hamilton([[0,1,0],[0,0,1],[1,0,0]],[])==3,'cyclic core empty word')
require(t.count_hamilton([[0,1],[0,0]],[1,2,0])==13,'H13')
require([t.count_hamilton(q,[1]) for q in ([[0,1],[0,0]],[[0,0],[1,0]])]==[1,3], 'fixed word hostile')
invalid=[([[0]],[-1]),([[0]],[2]),([[0]],[True]),([[0]],[1.0]),([], [1]),([[0.0]],[]),([[0,1],[1,0]],[]),([[0,0],[0,0]],[]),([[1]],[]),([[0,1]],[]),([[0,'1'],[0,0]],[])]
for q,w in invalid:
    try:t.count_hamilton(q,w)
    except (TypeError,ValueError):pass
    else:raise RuntimeError(('invalid accepted',q,w))
require(t.count_hamilton(((False,True),(False,False)),iter([1,2,0]))==13,'bool core and iterator word')
out['invalid_cases']=len(invalid)
# Existing source manifest and output consistency, without rerunning a producer script.
for line in (S/'SHA256SUMS').read_text().splitlines():
    digest,name=line.split();require(hashlib.sha256((S/name).read_bytes()).hexdigest()==digest,('hash',name))
require((S/'audit.out').read_bytes()==(S/'audit-optimized.out').read_bytes(),'optimized output differs')
out['source_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(S.iterdir()) if p.is_file()}
out['status']='PASS'
print(json.dumps(out,indent=2,sort_keys=True))
