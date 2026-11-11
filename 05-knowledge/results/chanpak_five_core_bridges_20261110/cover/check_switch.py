import importlib.util
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parents[1]/'session251_creative_completion_20261111/ternary/compiler.py'
spec=importlib.util.spec_from_file_location('inherited_compiler',SOURCE)
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
core=(2,4,9,19,7)

def require(p,msg):
    if not p:raise RuntimeError(msg)

def literal_count(rows):
    # Independent reverse recurrence: fix the first vertex of each induced path.
    n=len(rows);full=(1<<n)-1
    memo={}
    def f(mask,v):
        key=(mask,v)
        if key in memo:return memo[key]
        if mask==(1<<v):return 1
        rest=mask^(1<<v);choices=rows[v]&rest;total=0
        while choices:
            b=choices&-choices;choices-=b
            total+=f(rest,b.bit_length()-1)
        memo[key]=total;return total
    return sum(f(full,v) for v in range(n))

def main():
    constants=[2,4,-2,0,0,2,-4,-2]
    signature_checks=0
    for R,want in zip(range(0,32,4),constants):
        a=c.signature(core,R|1);b=c.signature(core,R|2)
        delta=tuple(x-y for x,y in zip(a,b))
        require(delta[:7]==(0,0,0,-2,0,2,0),'nonconstant signature delta')
        require(sum(delta[7:])==want,'context constant')
        signature_checks+=1
    graph_checks=0
    for A in range(8):
        for B in range(8):
            low=c.construct(core,2,A,B);high=c.construct(core,1,A,B)
            lo=literal_count(low);hi=literal_count(high)
            require(hi-lo==2*(B-A+1),'two-arc difference')
            require(lo==c.count(c.signature(core,2),A,B),'low compiler')
            require(hi==c.count(c.signature(core,1),A,B),'high compiler')
            graph_checks+=2
    return {'status':'PASS','literal_child_counts':graph_checks,
            'context_signature_checks':signature_checks,
            'scope':'All eight core contexts; actual graph controls A,B in0..7. All-parameter proof by coefficient subtraction.'}

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
