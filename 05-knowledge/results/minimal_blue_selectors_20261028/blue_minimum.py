"""Probe the minimal deterministic blue interface, with exact active checks."""
from pathlib import Path
import importlib.util,json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
spec=importlib.util.spec_from_file_location('prior_boundary',HERE.parent/'singleton_pair_reconstruction_20261028/boundary_profile.py')
import sys
sys.path.insert(0,str(Path(spec.origin).parent))
prior=importlib.util.module_from_spec(spec);spec.loader.exec_module(prior)
r=prior.r;need=r.need;core=r.core

def state(t):
    n=len(t)
    if n<3:return ('tiny',n)
    rows,cols,P,S=prior.raw(t)
    if len(set(rows[:n-2]))==1 and rows[0] in (0,1):return ('flat',rows[0])
    live=[]
    if P[n]:live.append(('center',))
    for d in range(1,n-1):
        if len(set(rows[:d]))==1 and rows[0] in (0,1) and S[n-d]:live.append(('shift',d,rows[0]))
        tail=cols[n-d:]
        if len(set(tail))==1 and tail[0] in (0,1) and P[n-d]:live.append(('shift',-d,tail[0]))
    need(len(live)<=1,('unique live offset',t,live))
    return live[0] if live else ('dead',)

def tau(s):
    if s[0]=='flat':return ('flat',1-s[1])
    if s[0]=='shift':return ('shift',s[1],1-s[2])
    return s

def join(s,left):
    step=1 if left else -1
    if s[0]=='tiny':return ('tiny',2) if s[1]==1 else ('flat',1)
    if s[0]=='dead':return s
    if s==('flat',1):return s
    if s[0] in ('flat','center'):return ('shift',step,1)
    if s[2]==0:return ('dead',)
    d=s[1]+step
    return ('shift',d,1) if d else ('center',)

def observed(s):return s[0] in ('flat','center')

def pure_context(t,d,c):
    for _ in range(abs(d)):
        if not c:t=core.complement(t)
        t=r.adjoin(t,d<0,False)
        if not c:t=core.complement(t)
    return t

def main():
    receipt=[]
    for n in range(1,8):
        reps={};N=1<<((n-1)*(n-2)//2)
        for x in range(N):
            t=core.presentation(n,x);s=state(t);reps.setdefault(s,x)
            need(state(core.complement(t))==tau(s),'tau transition')
            for left in (False,True):need(state(r.adjoin(t,left,False))==join(s,left),'singleton transition')
            need(observed(s)==(n>=3 and r.blue(t)),'blue observation')
        # Every resulting state is separated using these pure end-completions.
        contexts=[(0,0)]+[(d,c) for d in range(-max(1,n-3),max(1,n-3)+1) if d for c in (0,1)]
        profiles={s:tuple(r.blue(pure_context(core.presentation(n,x),d,c)) for d,c in contexts)
                  for s,x in reps.items()}
        need(len(set(profiles.values()))==len(reps),'all blue states distinguished')
        expected=[0,1,1,2,7][n] if n<=4 else 4*n-8
        need(len(reps)==expected,'sharp all-size state-count prediction')
        receipt.append(dict(n=n,presentations=N,blue_states=len(reps),representatives=[dict(state=s,mask=x) for s,x in reps.items()]))
    result=dict(status='FINITE-EXACT',universes=receipt,scope='Deterministic tau/J_L/J_R and blue observation; no reset in this standalone minimality claim.')
    (HERE/'blue-receipt.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps([{k:v for k,v in row.items() if k!='representatives'} for row in receipt]));return result

if __name__=='__main__':main()
