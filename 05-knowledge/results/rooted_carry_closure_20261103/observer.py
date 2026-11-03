"""Subset path/partition construction of the eighteen-polynomial named observer."""
from pathlib import Path
from collections import Counter
from math import gcd,factorial
import importlib.util,json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]

def load(name,relative):
    spec=importlib.util.spec_from_file_location(name,ROOT/relative)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

old=load('old_apex','05-knowledge/results/apex_wedge_constructor_20261103/mechanisms.py')
glue=load('old_glue','05-knowledge/results/coverage_machine_induction_20261103/gluing/check.py')
junction=load('old_junction','05-knowledge/results/coverage_machine_induction_20261103/junction.py')

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def path_types(A,ports):
    n=len(A);size=1<<n;types={v:i for i,v in enumerate(ports)}
    table=[Counter() for _ in range(size)]
    for start in range(n):
        dp=[[0]*n for _ in range(size)];dp[1<<start][start]=1
        for mask in range(1,size):
            if not(mask>>start&1):continue
            for last,value in enumerate(dp[mask]):
                if not value:continue
                table[mask][(types.get(start,-1),types.get(last,-1))]+=value
                remaining=A[last]&~mask
                while remaining:
                    bit=remaining&-remaining;remaining-=bit
                    dp[mask|bit][bit.bit_length()-1]+=value
    return table

def observer(A,ports):
    paths=path_types(A,ports);size=1<<len(A);covers=[None]*size;covers[0]=Counter({((),0):1})
    for mask in range(1,size):
        bit=mask&-mask;rest=mask^bit;part=rest;out=Counter()
        while True:
            block=part|bit;remaining=mask^block
            for endpoint,value in paths[block].items():
                for (s,m),c in covers[remaining].items():
                    key=(s,m+1) if endpoint==(-1,-1) else (tuple(sorted(s+(endpoint,))),m)
                    out[key]+=value*c
            if not part:break
            part=(part-1)&rest
        covers[mask]=out
    need(all(s in glue.STATES for s,m in covers[-1]),'eighteen valid skeletons')
    return covers[-1]

def root_profile(V,port):
    h=s=q=F=0
    for (sk,m),v in V.items():
        k=len(sk)+m
        if k==1:
            h+=v
            if any(a==port for a,b in sk):s+=v
        if k==2 and (port,port) in sk:q+=v
        if k in (1,2) and any(b==port for a,b in sk):F+=v
    return h,s,q,F

def delta(z):h,s,q,F=z;return F-h-q+s
def R(z):h,s,q,F=z;return h+F-s,q,q,2*F-q
def N(z):h,s,q,F=z;return h+F-s,h-s,h,h+F
def seal(z,k):
    for _ in range(k):z=R(z)
    return R(N(z))

def records(V):return [{'skeleton':s,'degree':m,'count':v} for (s,m),v in sorted(V.items()) if v]

def main():
    results=[]
    for P in (1,3,5,7):
        w,u=old.chart(P,P+1);plus=old.source(w,u);minus=old.source(w,u,False);ports=(len(plus)-2,len(plus)-1)
        vectors=[observer(A,ports) for A in (plus,minus)];states=[]
        for A,V in zip((plus,minus),vectors):
            z=[root_profile(V,i) for i in range(2)]
            need(z==[junction.profile(A,p) for p in ports],'observer versus independent rooted subset convolutions')
            states.append({'rows':A,'observer':records(V),'root_profiles':z,'defects':list(map(delta,z)),
                           'sealed':[[{'k':k,'state':seal(zz,k),'gcd':gcd(seal(zz,k)[0],seal(zz,k)[1])} for k in range(1,7)] for zz in z]})
        diff=Counter(vectors[1])
        for key,value in vectors[0].items():diff[key]-=value
        results.append({'P':P,'ports':ports,'source':states[0],'target':states[1],'difference':records(diff)})
        print(P,'H',states[0]['root_profiles'][0][0], 'roots',states[0]['root_profiles'],'->',states[1]['root_profiles'],flush=True)
    (HERE/'observer.json').write_text(json.dumps({'status':'FINITE_EXACT','examples':results},indent=2)+'\n',encoding='utf8',newline='\n')
if __name__=='__main__':main()
