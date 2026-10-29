"""Exact odd Lucas readers, compatible rotation phases, future quotients.

Universe: every cyclic binary no-11 word at the supplied odd lengths.
Only stdlib; explicit checks survive python -O. JSON is regenerated.
"""
from collections import Counter
from math import gcd, lcm
from pathlib import Path
import json

CHECKS = 0
def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok: raise RuntimeError(message)

def fib(n):
    a,b=0,1
    for _ in range(n): a,b=b,a+b
    return a

def factors(n):
    out=[]; p=2
    while p*p<=n:
        if n%p==0:
            q=1
            while n%p==0: n//=p; q*=p
            out.append(q)
        p+=1
    if n>1: out.append(n)
    return out

def words(n):
    def rec(i,w,prev):
        if i==n:
            if not (prev and w&1): yield w
            return
        yield from rec(i+1,w,0)
        if not prev: yield from rec(i+1,w|(1<<i),1)
    yield from rec(0,0,0)

def orbit_lengths(mult, modulus):
    unseen=set(range(modulus)); counts=Counter()
    while unseen:
        x=min(unseen); y=x; size=0
        while True:
            unseen.remove(y); size+=1; y=mult(y)
            if y==x: break
        counts[size]+=1
    return counts

def stable_partition(initial, operations):
    labels=initial; steps=[]
    while True:
        keys={}; new=[]
        for x in range(len(labels)):
            key=(labels[x],)+(tuple(labels[op[x]] for op in operations))
            if key not in keys: keys[key]=len(keys)
            new.append(keys[key])
        steps.append(len(keys))
        if len(keys)==len(set(labels)):
            sizes=Counter(Counter(new).values())
            return {"classes":len(keys),"rounds":steps,"fiber_sizes":dict(sorted(sizes.items()))}
        labels=new

def census(n):
    fn=fib(n); prev=fib(n-1); ln=2*prev+fn
    has4=n%3==0; m=ln//4 if has4 else ln
    need(gcd(fn,m)==1,"reader denominator")
    r=((1-prev)*pow(fn,-1,m))%m
    need((r*r-r-1)%m==0,"golden root")
    ps=factors(m); total=(4 if has4 else 1)*m
    need(total==ln,"norm cardinality")
    def mul4(a): return ((a&1)<<1)^((a>>1)*3)
    pows=[pow(r,i,m) for i in range(n)]
    f4=[1,2,3]
    decoder={}; coords={}
    for w in words(n):
        a=0; b=0
        for i in range(n):
            if w>>i&1:
                b=(b+pows[i])%m
                if has4: a^=f4[i%3]
        x=a*m+b
        need(x not in decoder,"word-reader injection")
        decoder[x]=w; coords[w]=x
    need(len(decoder)==ln,"word-reader surjection")
    rot=[]; ref=[]; add=[]
    for x in range(ln):
        w=decoder[x]; a,b=divmod(x,m)
        rw=((w<<1)&((1<<n)-1))|(w>>(n-1))
        sw=sum(((w>>i)&1)<<((-i)%n) for i in range(n))
        rot.append(coords[rw]); ref.append(coords[sw])
        add.append(((a^1) if has4 else 0)*m+(b+1)%m)
        need(rot[-1]==(mul4(a) if has4 else 0)*m+(r*b)%m,"rotation multiplication")
        if has4: need(ref[-1]//m==((a&1)^(a>>1))+((a>>1)<<1),"reflection Frobenius")
    for x in range(ln):
        need(ref[ref[x]]==x and ref[rot[ref[rot[x]]]]==x,"dihedral identities")
    periods=orbit_lengths(lambda x:rot[x],ln)
    locals_=[]
    if has4: locals_.append(orbit_lengths(mul4,4))
    for q in ps: locals_.append(orbit_lengths(lambda x:(r*x)%q,q))
    # Independent local orbit names omit product(d_i)/lcm(d_i) phase classes.
    phase_hist=Counter(); global_count=0; local_types=0
    def combine(i,prod,period,multiplicity):
        nonlocal global_count,local_types
        if i==len(locals_):
            phase_hist[prod//period]+=multiplicity
            global_count+=multiplicity*prod//period; local_types+=multiplicity
            return
        for d,c in locals_[i].items(): combine(i+1,prod*d,lcm(period,d),multiplicity*c)
    combine(0,1,1,1)
    need(global_count==sum(periods.values()),"compatible CRT orbit formula")
    projections={}
    if has4: projections["4"]=stable_partition([x//m for x in range(ln)],[ref,add])
    for q in ps: projections[str(q)]=stable_partition([(x%m)%q for x in range(ln)],[ref,add])
    return {"n":n,"states":ln,"odd_modulus":m,"root":r,"prime_power_factors":ps,
            "rotation_period_counts":dict(sorted(periods.items())),
            "local_period_counts":[dict(sorted(c.items())) for c in locals_],
            "independent_local_orbit_types":local_types,"global_rotation_orbits":global_count,
            "relative_phase_classes_per_local_type":dict(sorted(phase_hist.items())),
            "future_quotients_for_reflection_and_add_one":projections}

def multiply(x,y,a,m):
    r,s=x; t,v=y
    return ((r*t+s*v)%m,(r*v+s*t+a*s*v)%m)

def power(e,a,m):
    x=(1,0); u=(0,1)
    while e:
        if e&1: x=multiply(x,u,a,m)
        u=multiply(u,u,a,m); e//=2
    return x

def unit_clocks():
    out=[]
    for p in [11,13]:
        for h in range(1,5):
            m=p**h; order=24*11**max(0,h-2) if p==11 else 12*13**(h-1)
            need(power(order,5,m)==(1,0),"unit period")
            for q in [2,3,p]:
                if order%q==0: need(power(order//q,5,m)!=(1,0),"minimal unit period")
            b=power(12,5,m)
            # Trace-matrix determinant is -29 times B; Smith first entry 1.
            fiber=gcd(m,29*71351280)
            out.append({"prime":p,"power":h,"modulus":m,"unit_order":order,
                        "twelfth_power":b,"sampled_fiber":fiber,"sampled_states":m*m//fiber,
                        "primitive_states":m*m})
    return out

if __name__=="__main__":
    result={"status":"FINITE-EXACT","universe":"all cyclic no-11 words at each listed odd length; all residue states",
            "censuses":[census(n) for n in [9,11,13,15,21,23,25]],"unit_clocks":unit_clocks()}
    result["checks"]=CHECKS
    dest=Path(__file__).with_suffix('.json'); dest.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
