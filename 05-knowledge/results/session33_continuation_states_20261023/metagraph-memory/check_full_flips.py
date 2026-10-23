"""Tiny exact addendum certificate, consuming the independently computed tilings.
Run after check.py. No repository imports or writes.
"""
import json, hashlib
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent

def require(v):
    if not v: raise RuntimeError('Full-flip certificate failed')

def refine(observation, shifts):
    labels=observation[:]; sizes=[len(set(labels))]
    while True:
        sig=[(labels[x],*(labels[x^a] for a in shifts)) for x in range(len(labels))]
        lookup={t:i for i,t in enumerate(sorted(set(sig)))}
        new=[lookup[t] for t in sig]; sizes.append(len(lookup))
        if sizes[-1]==sizes[-2]: return new,sizes
        labels=new

out=[]
for row in json.loads((HERE/'results.json').read_text()):
    n=row['n']
    if n>5: continue
    nodes=row['node']; N=len(nodes); m=(n-1)*(n-2)//2; anchor=N-1
    require(N==1<<m)
    # An explicit directed three-cycle witnesses every nonanchor as nontransitive:
    # choose a shortest backward free arc i<-j. All shorter arcs point forward,
    # so i->i+1->j->i (including the forced path case j=i+2).
    free=[(i,j) for i in range(n) for j in range(i+2,n)]
    witnesses=[]
    for x in range(N-1):
        i,j=min((pair for k,pair in enumerate(free) if not ((x>>k)&1)),key=lambda p:p[1]-p[0])
        k=i+1
        def arc(a,b): return b==a+1 or (b>a+1 and ((x>>free.index((a,b)))&1))
        require(arc(i,k) and arc(k,j) and not arc(i,j))
        witnesses.append([x,i,k,j])
    require([x for x in range(N) if nodes[x]==nodes[anchor]]==[anchor])
    tau,ts=refine(nodes,[anchor]); full,fs=refine(nodes,[1<<k for k in range(m)])
    require(len(set(tau))==row['pair_states']); require(len(set(full))==N)
    tested=0
    for x in range(N):
        shift=x^anchor
        for y in range(x+1,N):
            require(nodes[x^shift]!=nodes[y^shift]); tested+=1
    # Orbit of the singleton anchor indicator is literally a permutation basis.
    translated=s.Matrix([[int((x^a)==anchor) for a in range(N)] for x in range(N)])
    require(translated*translated.T==s.eye(N)); require(translated.rank()==N)
    out.append(dict(n=n,m=m,tilings=N,tau_refinement_sizes=ts,full_flip_refinement_sizes=fs,linear_observer_dimension=N,positive_hidden_state_minimum=N,pairwise_separations=tested,anchor=anchor,nonanchor_cycle_witnesses=witnesses))
receipt=dict(scope='FINITE-EXACT n=3,4,5; full labelled free-coordinate flip continuation',source_results_sha256=hashlib.sha256((HERE/'results.json').read_bytes()).hexdigest(),rows=out)
(HERE/'full_flips_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
for r in out: print({k:v for k,v in r.items() if k!='nonanchor_cycle_witnesses'})
