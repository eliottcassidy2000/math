"""Search complete bounded one-clock funnel languages, with exact word coefficients."""
from pathlib import Path
from itertools import product
import sys,json,gzip,time,argparse
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'berggren_target_return_20261111'/'reset'))
import funnel_six_state as f

def need(ok,message):
    if not ok:raise RuntimeError(message)

def run(core_order=8,fixed=12):
    cores=[dict(order=1,rank=0,state=f.SINGLE)]
    for n in range(3,core_order+1):
        cores.extend(dict(order=n,rank=j,state=f.compile_rank(n,j)) for j in range(f.a.class_counts(n)[n]))
    # A target match must hold in the full (A,B) coefficient pair.
    targets=set()
    for a,b in [(45,7),(99,23),(135,23),(297,71),(33,7),(11,1)]:
        for depth in range(8):
            targets.add((a,b));a,b=3*a,3*b+2
    hits={key:None for key in targets};families={};contexts=0;t0=time.time()
    for size in range(fixed+1):
        distinct={}
        for left_len in range(size+1):
            right_len=size-left_len
            for bits in product((0,1),repeat=size):
                left=bits[:left_len];right=bits[left_len:]
                for middle in (0,1):
                    if (right[-1] if right else middle)!=1:continue
                    samples=[]
                    for k in (1,2,3):
                        word=left+(middle,)*k+right
                        l,d,b=f.a.tail_coefficients(word,1);L=f.a.tail_coefficients(word,2)[0]
                        samples.append((L,b,d))
                    A=tuple((v-u)//2 for u,v in zip(samples[0],samples[1]))
                    B=tuple(u-2*a for u,a in zip(samples[0],A))
                    need(tuple(8*a+b for a,b in zip(A,B))==samples[2],'one-clock exact exponential candidate')
                    key=A+B
                    distinct.setdefault(key,(left,middle,right))
        contexts+=sum(1 for _ in distinct)
        for (aL,aB,aD,bL,bB,bD),(left,middle,right) in distinct.items():
            for core in cores:
                h,e,s,p,ee,t=core['state'];G=s+p
                a=aL*G+aB*h-aD*e;b=bL*G+bB*h-bD*e-ee-1
                if a<=0:continue
                key=(a,b)
                cost=core['order']+size+2
                rec=None
                if key in hits and (hits[key] is None or cost<hits[key]['fixed_order']):
                    rec=dict(**core,left=left,middle=middle,right=right,A=a,B=b,fixed_order=cost,word_coefficients=[aL,aB,aD,bL,bB,bD])
                    hits[key]=rec
                # Keep a bounded reusable coefficient grammar, not all huge rays.
                if a<=5000 and abs(b)<=5000 and (key not in families or cost<families[key]['fixed_order']):
                    families[key]=rec or dict(**core,left=left,middle=middle,right=right,A=a,B=b,fixed_order=cost,word_coefficients=[aL,aB,aD,bL,bB,bD])
        print('fixed',size,'contexts',contexts,'families',len(families),'targets',[(k,v['fixed_order']) for k,v in sorted(hits.items()) if v], 'seconds',round(time.time()-t0,1),flush=True)
    edges=[]
    for a,b in families:
        for d in (2,4):
            if (3*a,3*b+d) in families:edges.append([a,b,3*a,3*b+d,d])
    data=dict(scope=f'All unique-start strong rooted cores through{core_order}, plus singleton; all binary contexts with total fixed length <= {fixed}; middle bit0 or1; tail ends1; k>=1',core_order=core_order,fixed=fixed,cores=len(cores),distinct_word_observers=contexts,
              hits=[dict(target=k,recipe=v) for k,v in sorted(hits.items())],families=list(families.values()),edges=edges)
    tag='' if __debug__ else '.optimized'
    (HERE/('family_search'+tag+'.json.gz')).write_bytes(gzip.compress(json.dumps(data,separators=(',',':')).encode(),mtime=0))
    print('PASS',len(families),'families',len(edges),'arithmetic edges; physical compilation audited separately',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--core-order',type=int,default=8);p.add_argument('--fixed',type=int,default=12);a=p.parse_args();run(a.core_order,a.fixed)
