"""Complete target selector below26521 for the guarded23-core prefix language."""
from collections import deque
from math import gcd
import json,gzip
from fibre_probe import HERE,save,need
from repair23 import hp2,join,K9,birth,profile,strong

LIMIT=26520

def run():
    raw=json.loads(gzip.decompress((HERE/'fibre_probe.json.gz').read_bytes()))['records']
    ordinary={}
    for r in raw:ordinary.setdefault((r['n'],r['code']),tuple(r['rows']))
    generators={(1,2):dict(rows=(0,),n=1,code=0)}
    for (n,code),T in sorted(ordinary.items()):
        h,f=hp2(T)
        if 23*h+60*f>LIMIT or gcd(h,30*f)!=1:continue
        generators.setdefault((h,f),dict(rows=T,n=n,code=code))
    gen=[dict(h=h,f=f,**v) for (h,f),v in sorted(generators.items())]
    states={(1,1):dict(rows=(),word=[])};todo=deque([(1,1)])
    while todo:
        h,f=todo.popleft();old=states[h,f]
        for i,g in enumerate(gen):
            H,F=h*g['h'],f*g['f']
            if 23*H+60*F>LIMIT or gcd(H,30*F)!=1 or (H,F) in states:continue
            states[H,F]=dict(rows=join(old['rows'],tuple(g['rows'])),word=old['word']+[i]);todo.append((H,F))
    outputs={}
    for (h,f),s in states.items():
        if not s['rows']:continue
        L=s['rows'];need(hp2(L)==(h,f),'literal prefix recount')
        T=join(L,K9);r=len(T)-1;T=birth(T,r);N=23*h+60*f
        z=profile(T,r);need(z==(N,9*h,9*h,120*f-9*h) and strong(T),'actual target and entry profile')
        need(gcd(50*h,60*f-9*h)==1,'every selected source enters primitively')
        outputs.setdefault(str(N),[]).append(dict(h=h,f=f,prefix_rows=L,word=s['word'],target_rows=T,root=r,profile=z))
    # Exactness for all larger SCCs is proved by the strong-order/P2 bounds.
    eligible=[N for N in range(1,LIMIT+1,2) if gcd(N,15)==1]
    gaps=[N for N in eligible if str(N) not in outputs]
    save('prefix_selector',dict(status='FINITE-EXACT with an all-order exclusion certificate for this language below26521',limit=LIMIT,generators=gen,prefix_pairs=len(states)-1,outputs=outputs,eligible_targets=len(eligible),absent_in_language=gaps,first_target=min(map(int,outputs)),scope='All nonempty finite prefix tournaments L with gcd(H(L),30P_L(2))=1; target23H(L)+60P_L(2). Missing targets are not claimed missing from tournaments.'))
    print('GENERATORS',len(gen),'PREFIX_PAIRS',len(states)-1,'TARGETS',len(outputs),'ELIGIBLE',len(eligible))
    print('FIRST_TARGETS',sorted(map(int,outputs))[:24])
    print('LAST_TARGETS',sorted(map(int,outputs))[-12:])

if __name__=='__main__':run()
