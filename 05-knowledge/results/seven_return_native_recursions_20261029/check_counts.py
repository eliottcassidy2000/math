"""Matched-output audits and benchmarks; no assertion-only checks."""
import argparse, hashlib, json, statistics, time
import baseline87 as b87
import baseline105 as b105
import parked_counts as pc

def main(bench):
    checks=0
    def need(x,msg):
        nonlocal checks
        checks+=1
        if not x:raise ValueError(msg)
    for n in range(31):
        direct=b87.count_direct(n)[0]
        need(b105.count(n)[0]==direct,('unpacked projective',n))
        for park in (False,True):
            for half in (False,True):
                need(pc.projective(n,park,half)[0]==direct,('projective',n,park,half))
        for q in (1,2,3):
            expected=b87.degree_bucket_dp(n,colors=q)[0]
            for park in (False,True):
                need(pc.prefix(n,colors=q,park=park)[0]==expected,('prefix',n,q,park))
    rows=[]
    for n in ((100,200) if bench else (50,100)):
        methods=[('prefix_old',lambda:b87.degree_bucket_dp(n)),
                 ('prefix_parked',lambda:pc.prefix(n)),
                 ('projective_old',lambda:b105.count_packed(n)),
                 ('projective_full',lambda:pc.projective(n,False,False)),
                 ('projective_park',lambda:pc.projective(n,True,False)),
                 ('projective_half',lambda:pc.projective(n,False,True)),
                 ('projective_both',lambda:pc.projective(n,True,True))]
        times={name:[] for name,_ in methods};outputs={};stats={}
        # Reverse alternately to reduce systematic order bias, no other benchmark jobs.
        for rep in range(3 if bench else 1):
            for name,f in (methods if rep%2==0 else list(reversed(methods))):
                start=time.perf_counter();value,stat=f();elapsed=time.perf_counter()-start
                times[name].append(elapsed);outputs[name]=value;stats[name]=stat
        answer=outputs['prefix_old'][-1]
        for name,value in outputs.items():
            need(value==(outputs['prefix_old'] if name.startswith('prefix') else answer),('matched output',n,name))
        compact={name:{k:v for k,v in st.items() if k not in ('rows','stages')} for name,st in stats.items()}
        rows.append(dict(n=n,seconds=times,median_seconds={name:statistics.median(t) for name,t in times.items()},
                         stats=compact,answer_sha256=hashlib.sha256(str(answer).encode()).hexdigest()))
    return dict(status='PASS',explicit_checks=checks,small_universe='all n=0..30; colors 1,2,3; four projective and two prefix modes',benchmarks=rows)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--bench',action='store_true')
    print(json.dumps(main(p.parse_args().bench),indent=2,sort_keys=True))
