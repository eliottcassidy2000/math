"""Independent full rational sums and matched precision benchmarks."""
from fractions import Fraction
from collections import Counter
from math import factorial,gcd
import json,statistics,time
import baseline_precision as old
import pruned_precision as new

def full(n,r,s):
    def partitions(rem,low=1):
        if rem==0:yield ();return
        for a in range(low,rem+1,2):
            for tail in partitions(rem-a,a):yield (a,)+tail
    out=Fraction()
    for parts in partitions(n):
        k=len(parts);E=(n-k)//2+sum(gcd(a,b) for i,a in enumerate(parts) for b in parts[i+1:])
        z=1
        for a,m in Counter(parts).items():z*=a**m*factorial(m)
        out+=Fraction(2**(r*E+s*k),z)
    if out.denominator!=1:raise ArithmeticError('fraction sum not integral')
    v=old.base_valuation(n,r,s);q,rem=divmod(out.numerator,2**v)
    if rem:raise ArithmeticError('nonintegral normalization')
    return q

def verify():
    checks=0
    def need(x,msg):
        nonlocal checks
        checks+=1
        if not x:raise ValueError(msg)
    for n in range(31):
        for r,s in ((1,0),(1,1),(2,1),(3,0)):
            expected=full(n,r,s)
            for bits in (1,2,4,8,12):
                ref=old.normalized_low_bits(n,bits,r,s,shortcut=False)[0]
                need(ref==expected%(1<<bits),('full independent channel',n,r,s,bits))
                for mode in ('collision','completion','sum'):
                    need(new.count(n,bits,r,s,mode)[0]==ref,('pruned channel',n,r,s,bits,mode))
    # A too-low refund gives a real false exclusion: (1,1,1), not just a floor mismatch.
    need(old.term((1,1,1),1,0,True)[0]==1,'repeated-unit hostile')
    need(new.count(3,2)[0]==old.normalized_low_bits(3,2,shortcut=False)[0]==1,'retain multiplicity refund')
    benchmarks=[]
    for n,bits in ((201,8),(401,8),(101,16)):
        times={k:[] for k in ('inherited','collision','completion','sum')};values={};stats={}
        funcs=[('inherited',lambda:old.normalized_low_bits(n,bits,shortcut=False)),
               ('collision',lambda:new.count(n,bits,mode='collision')),
               ('completion',lambda:new.count(n,bits)),
               ('sum',lambda:new.count(n,bits,mode='sum'))]
        for rep in range(3):
            for name,f in (funcs if rep%2==0 else list(reversed(funcs))):
                t=time.perf_counter();v,st=f();times[name].append(time.perf_counter()-t);values[name]=v;stats[name]=st
        need(len(set(values.values()))==1,('benchmark output',n,bits))
        need(stats['completion']['retained']==stats['inherited']['retained'],'no dropped contributing tuple')
        need(stats['sum']['retained']==stats['inherited']['retained'],'sum ceiling retains every contributor')
        benchmarks.append(dict(n=n,bits=bits,residue=values['completion'],seconds=times,
                               median_seconds={k:statistics.median(t) for k,t in times.items()},stats=stats))
    return dict(status='PASS',explicit_checks=checks,channels=[[1,0],[1,1],[2,1],[3,0]],
                full_fraction_orders=[0,30],precisions=[1,2,4,8,12],benchmarks=benchmarks)
if __name__=='__main__':print(json.dumps(verify(),indent=2,sort_keys=True))
