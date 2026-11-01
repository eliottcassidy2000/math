"""Independent direct permutation Burnside audit of all seven emitted counters."""
from itertools import combinations,permutations
import importlib.util
import json
import math
from pathlib import Path
import sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(Path(__file__).parent))
from audit_compiler import check,ir_run,raw_run
from audit_path_counter import emitted_run

def direct_burnside(n):
    pairs=list(combinations(range(n),2));pindex={p:i for i,p in enumerate(pairs)}
    totals=dict(a000568=0,self_converse=0,pair_tournaments=0,
                colored_tournaments=0,colored_pair_tournaments=0)
    for perm in permutations(range(n)):
        seen=set();vertex_cycles=0
        for v in range(n):
            if v in seen:continue
            vertex_cycles+=1;u=v
            while u not in seen:seen.add(u);u=perm[u]
        seen=set();orbits=0;valid=True;anti=True
        for start in pairs:
            if start in seen:continue
            orbits+=1;pair=start;flip=length=0
            while pair not in seen:
                seen.add(pair);a,b=perm[pair[0]],perm[pair[1]]
                flip^=int(a>b);length+=1;pair=(min(a,b),max(a,b))
            valid=valid and flip==0
            anti=anti and flip==(length%2)
        if valid:
            totals['a000568']+=2**orbits
            totals['pair_tournaments']+=4**orbits
            totals['colored_tournaments']+=2**(orbits+vertex_cycles)
            totals['colored_pair_tournaments']+=2**(2*orbits+vertex_cycles)
        if anti:totals['self_converse']+=2**orbits
    for key,value in totals.items():
        check(value%math.factorial(n)==0,('nonintegral Burnside',n,key))
        totals[key]//=math.factorial(n)
    check((totals['a000568']+totals['self_converse'])%2==0,'converse parity')
    totals['converse_merged']=(totals['a000568']+totals['self_converse'])//2
    return totals

def load_source(name,path):
    spec=importlib.util.spec_from_file_location('counter_'+name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m.count

def load_program(lane,name):
    p=json.loads((lane/(name+'.meta.json')).read_text())
    p['instructions']=sum((json.loads((lane/file).read_text()) for file in p['instruction_parts']),[])
    p['state_primes']=sum((json.loads((lane/file).read_text()) for file in p['state_prime_parts']),[])
    check(len(p['fraction_parts'])==1,'auditor currently expects one fraction part')
    return p,lane/p['fraction_parts'][0]

def main():
    lane=ROOT/'fractran';table=[direct_burnside(n) for n in range(8)]
    for n,row in enumerate(table):
        row['strong']=0 if n==0 else row['a000568']-sum(table[k]['strong']*table[n-k]['a000568'] for k in range(1,n))
        check(row['strong']>=0,'negative strong count')
    cases=[]
    for name in table[0]:
        source=load_source(name,lane/(name+'_source.py'))
        for n,row in enumerate(table):check(source(n)==row[name],('source',name,n,row[name]))
        p,path=load_program(lane,name)
        for n in range(4):
            actual=ir_run(p,n,bound=20_000_000)
            check(actual[0][p['output']]==table[n][name] and not actual[0][p['error']],('IR',name,n))
            if n<=2:
                emitted=emitted_run(p,path,n)
                check(emitted==actual,('emitted',name,n))
            if n==0:check(raw_run(p,n)[:3]==actual,('literal',name,n))
            cases.append({'program':name,'n':n,'output':table[n][name],'steps':actual[2]})
    print(json.dumps({'status':'PASS','independent_Burnside_orders':'0..7',
          'method':'Full symmetric-group permutations; direct signed unordered-pair orbits, including complement-twisted fixed sets.',
          'table':table,'IR_cases':cases,'actual_fraction_exponent_cases':21,
          'literal_integer_cases':7},indent=2))

if __name__=='__main__':main()
