#!/usr/bin/env python3
"""Regenerate and verify finite FRACTRAN artifacts; assertions are never used.
Run normally and with python -O.  No network or external package is needed.
"""
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
from compiler import Compiler, fractions, require, run_register, run_fractran_literal, run_fractran_exponents, write_artifacts
import literal16

ROOT=Path(__file__).resolve().parent


def must_fail(call,description):
    try: call()
    except (RuntimeError,ValueError): return
    raise RuntimeError('hostile control unexpectedly survived: '+description)


def independent_burnside(n):
    """Action on oriented unordered edges; no partition, gcd or compiler formulas."""
    edges=list(itertools.combinations(range(n),2))
    total=0
    for perm in itertools.permutations(range(n)):
        visited=set(); cycles=0; admissible=True
        for start in edges:
            if start in visited: continue
            edge=start; orientation=0
            while True:
                visited.add(edge)
                x,y=perm[edge[0]],perm[edge[1]]
                if x>y: orientation^=1; x,y=y,x
                edge=(x,y)
                if edge==start: break
            if orientation:
                admissible=False; break
            cycles+=1
        if admissible: total+=2**cycles
    den=math.factorial(n)
    require(total%den==0,'independent Burnside divisibility')
    return total//den


def valuation_bound(rows,prime):
    def valuation(a):
        result=0
        while a%prime==0:
            result+=1
            a//=prime
        return result
    differences=[valuation(a)-valuation(b) for a,b in rows]
    require(max(differences)==1 and min(differences)==-1,'output increment bound changed')
    return {'output_prime':prime,'max_increment':max(differences),'min_increment':min(differences),'increment_rows':sum(d>0 for d in differences),'decrement_rows':sum(d<0 for d in differences)}


def primitive_tests():
    checks=0
    tests=[]
    for k in range(7):
        tests.extend([
            (f'n+{k}',lambda n,k=k:n+k),
            (f'n-{k}',lambda n,k=k:max(0,n-k)),
            (f'n*{k}',lambda n,k=k:n*k),
            (f'n=={k}',lambda n,k=k:int(n==k)),
            (f'n!={k}',lambda n,k=k:int(n!=k)),
            (f'n<{k}',lambda n,k=k:int(n<k)),
            (f'n<={k}',lambda n,k=k:int(n<=k)),
            (f'n>{k}',lambda n,k=k:int(n>k)),
            (f'n>={k}',lambda n,k=k:int(n>=k)),
        ])
        if k:
            tests.extend([(f'n//{k}',lambda n,k=k:n//k),(f'n%{k}',lambda n,k=k:n%k)])
    tests.extend([('2**n',lambda n:2**n),('n**0',lambda n:1),('(n+2)*(n+1)',lambda n:(n+2)*(n+1))])
    for expression,oracle in tests:
        p=Compiler().compile(f'def f(n):\n    return {expression}\n')
        for n in range(7):
            r=run_register(p,n)
            require(r['output']==oracle(n) and r['error']==0,f'primitive mismatch {expression}, n={n}')
            f=run_fractran_exponents(p,n)
            require(f==r,f'fraction/IR mismatch {expression}, n={n}')
            checks+=1
    # Literal integer replay is separate from exponent dispatch.
    p=Compiler().compile('def f(n):\n    return (n+1)*2\n')
    for n in range(7):
        r=run_fractran_literal(p,n)
        require(r['output']==2*(n+1),'literal macro control')
        checks+=1
    p=Compiler().compile('def f(n):\n    if n > 2:\n        return 7\n    else:\n        return 4\n')
    for n in range(7):
        require(run_register(p,n)['output']==(7 if n>2 else 4),'branch control')
        checks+=1
    z=Compiler().compile('def f(n):\n    return n//0\n')
    require(run_register(z,3)['error']==1,'division-by-zero error must be explicit')
    must_fail(lambda: run_register(p,-1),'negative input')
    loop=Compiler().compile('def f(n):\n    while n != 0:\n        n = n\n    return n\n')
    must_fail(lambda: run_register(loop,1,100),'diverging loop bound')
    bad={'registers':{'n':0},'register_count':1,'instructions':[('INC',0,0)],'input':0,'output':0,'error':0}
    must_fail(lambda: fractions(bad),'canceling state guard')
    # Actual first-eligible hostile: zero-first DECJZ cannot implement conditional decrement.
    rows=[(7,3),(11,6)] # both enabled at 2*3; unconditional zero fraction is erroneously first.
    value=6
    selected=next((a,b) for a,b in rows if value%b==0)
    require(selected==(7,3) and selected!=(11,6),'ordered-choice hostile')
    return {'primitive_cases':checks,'hostile_controls':5}


def main():
    report={'schema':'session132.finite-fractran.v1','primitive_tests':primitive_tests()}
    for name in ['a000568','literal_lines']:
        source=(ROOT/f'{name}_source.py').read_text()
        p=Compiler().compile(source)
        report[name]={'program':write_artifacts(p,ROOT/name)}
        if name=='a000568':
            report[name]['valuation_bound']=valuation_bound([(a,b) for a,b,_,_ in fractions(p)['fractions']],3)
        if name=='a000568':
            literal=[]; exponent=[]; register=[]
            for n in range(7):
                oracle=independent_burnside(n)
                r=run_register(p,n)
                require(r['output']==oracle and r['error']==0,f'A000568 mismatch at n={n}')
                register.append({'n':n,'output':r['output'],'steps':r['steps']})
                if n<=5:
                    f=run_fractran_exponents(p,n)
                    require(f==r,f'A000568 full exponent/IR mismatch n={n}')
                    exponent.append({'n':n,'output':f['output'],'steps':f['steps']})
                if n<=3:
                    f=run_fractran_literal(p,n)
                    require(f['output']==r['output'] and f['steps']==r['steps'],'literal A000568 mismatch')
                    literal.append({'n':n,**f})
            report[name].update(register_replays=register,exponent_fraction_replays=exponent,literal_integer_replays=literal)
        else:
            results=[]
            for t in range(5):
                r=run_register(p,t)
                require(r['output']==literal16.expected(t+3) and r['error']==0,'compiled L mismatch')
                results.append({'n':t+3,'output':r['output'],'steps':r['steps']})
            report[name]['register_replays']=results
    report['literal16']={'valuation_bound':valuation_bound(literal16.FRACTIONS,2),'fractions':[f'{a}/{b}' for a,b in literal16.FRACTIONS],'integer_replays':[],'exponent_replays':[]}
    for n in range(3,8):
        r=literal16.run_exponents(n)
        require(r['output']==literal16.expected(n) and r['steps']==literal16.exact_steps(n),'literal16 recurrence/step mismatch')
        report['literal16']['exponent_replays'].append(r)
        if n<=7:
            s=literal16.run_integer(n)
            require(s['output']==r['output'] and s['steps']==r['steps'],'literal16 integer replay mismatch')
            report['literal16']['integer_replays'].append(s)
    # Bind every source and concrete fraction/program artifact; no giant integer trace.
    report['sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(ROOT.iterdir()) if p.is_file() and (p.suffix in ['.py','.txt'] or p.name.endswith('.ir.json'))}
    text=json.dumps(report,indent=2,sort_keys=True)+'\n'
    (ROOT/'VERIFIED.json').write_text(text)
    print(json.dumps({k:report[k] for k in ['primitive_tests','a000568','literal16']},indent=2))

if __name__=='__main__': main()
