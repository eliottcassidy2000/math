"""Exact small controls for the contravariant seed guards and Fibonacci charge carry."""
import importlib.util,itertools,json,math
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('seed',HERE.parent/'seed_decision'/'check.py');s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)
def need(x,msg):
    if not x:raise RuntimeError(msg)
def dot(a,b,K):return min(K,sum(x*y for x,y in zip(a,b)))
def pull(f,L,K):
    return tuple(min(K,sum(L[t]*f[s.DEST[i][t]] for t in range(9))) for i in range(9))
FH=(0,1,1,0,0,0,0,0,0);FW=(0,0,0,0,0,0,0,0,1);FU=(0,0,1,0,0,1,0,0,1)
def b(n):
    m=n+1;return (3*m-math.isqrt(5*m*m)-1)//2
def Q(n):return (n-2*b(n),b(n))
def delta(x,y):return b(x)+b(y)-b(x+y)
def digits(n):
    fib=[1,2]
    while fib[-1]<=n:fib.append(sum(fib[-2:]))
    ds=[]
    for f in fib[::-1]:
        d=int(n>=f);ds.append(d);n-=d*f
    while len(ds)>1 and not ds[0]:ds.pop(0)
    return ds
bank=[]
for T in ((0,), (2,4,1),(2,12,9,1)):
    for O in range(1<<len(T)):
        L=s.block_data(T,O)
        if L not in bank:bank.append(L)
start=(1,0,0,0,0,0,0,0,0);checks=0;carrychecks=0;carryvalues=set();neutral=0
prefixes=[start]
for L in bank:prefixes.append(s.advance(start,L))
for A in bank:
    for B in bank:prefixes.append(s.advance(A,B))
q=6;K=q+2
for v in prefixes:
    for L in bank:
        out=s.advance(v,L,K)
        for f in (FH,FW,FU):
            need(dot(out,f,K)==dot(tuple(min(K,x) for x in v),pull(f,L,K),K),'exact contravariant pullback');checks+=1
        h,w,u=s.observables(s.advance(v,L))
        x,y,z=s.advance(v,L)[2::3]
        # Entries 2,5,8 are the actual three U channels.
        k=delta(x,y)+delta(x+y,z)
        need(k==b(x)+b(y)+b(z)-b(u),'carry cocycle telescope')
        raw=tuple(sum(Q(a)[i] for a in (x,y,z)) for i in (0,1))
        need((raw[0]+2*k,raw[1]-k)==Q(u),'corrected three-channel charge');carrychecks+=1;carryvalues.add(k)
        X=Y=prev=0
        for d in digits(u):
            need(not(prev and d),'legal digits');X,Y=Y+d,X+Y+2*d;prev=d
        need(X==u and (2*Y-3*X,2*X-Y)==Q(u),'native two-register charge')
for a in range(20):
    for bb in range(20):
        for c in range(20):
            need(delta(a,bb)+delta(a+bb,c)==delta(bb,c)+delta(a,bb+c),'addition-tree cocycle')
# Historical missing-register hostile is a preserved witness, not an old-result reaudit.
need((2%3,8%3)==(2,2),'same value residue')
need((3%3,13%3)==(0,1),'different shifted residue')
need(Q(6)==(2,2),'neutral charge exists')
out={'status':'PASS','scope':'finite controls for new dual seed guards and inherited charge interface, not universal seed coverage','alphabet_transfers':len(bank),'prefixes':len(prefixes),'dual_guard_checks':checks,'three_channel_charge_checks':carrychecks,'carry_values_seen':sorted(carryvalues),'cocycle_checks':8000,'q':q,'cap':K}
(HERE/'history_bridge_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
