#!/usr/bin/env python3
import itertools,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def require(t,m):
    if not t:raise RuntimeError(m)
def strong(a,S):
    S=set(S);start=next(iter(S));seen={start}
    while True:
        new=seen|{j for i in seen for j in S if a[i]>>j&1}
        if new==seen:break
        seen=new
    if seen!=S:return False
    seen={start}
    while True:
        new=seen|{j for i in seen for j in S if a[j]>>i&1}
        if new==seen:return seen==S
        seen=new

def configurations(a):
    exits=[];entries=[];n=len(a)
    for S in itertools.combinations(range(n),4):
        if not strong(a,S):continue
        for v in range(n):
            if v not in S:
                d=sum(a[x]>>v&1 for x in S)
                if d>=3:exits.append((S,v))
                if d<=1:entries.append((S,v))
    return exits,entries

def canon(a):
    n=len(a)
    return min(sum((a[p[i]]>>p[j]&1)<<k for k,(i,j) in enumerate(itertools.combinations(range(n),2))) for p in itertools.permutations(range(n)))
def extension(core,out):
    n=len(core);a=core[:]+[0]
    for i in range(n):
        if i in out:a[n]|=1<<i
        else:a[i]|=1<<n
    return a
s4=[2,4,9,3]
explicit=[({0,1},(1,2,3,4),0),({0,2},(1,4,2,3),0),({0,3},(1,2,4,3),0),({1,2},(0,4,2,3),1)]
for out,cycle,v in explicit:
    a=extension(s4,out)
    require(all(a[cycle[i]]>>cycle[(i+1)%4]&1 for i in range(4)),'displayed cycle wrong')
    require(sum(a[x]>>v&1 for x in cycle)>=3,'displayed exit wrong')
c5=extension(s4,{2,3});q5=extension(s4,{1,3})
require(all(x.bit_count()==2 for x in c5),'C5 not regular')
# The claimed Q5 labels: a=2, b=1, triangle 0,3,4.
require(all(q5[2]>>v&1 and q5[v]>>1&1 for v in [0,3,4]) and q5[1]>>2&1,'Q5 partition wrong')
require(strong(q5,[0,3,4]),'Q5 triangle wrong')
C={canon(c5):'C5',canon(q5):'Q5'};counts={'strong':0,'exit_or_entry':0,'C5':0,'Q5':0}
for mask in range(1<<10):
    a=[0]*5
    for k,(i,j) in enumerate(itertools.combinations(range(5),2)):
        if mask>>k&1:a[i]|=1<<j
        else:a[j]|=1<<i
    if not strong(a,range(5)):continue
    counts['strong']+=1
    e,b=configurations(a)
    if e or b:counts['exit_or_entry']+=1
    else:
        key=canon(a);require(key in C,'unclassified strong five-core');counts[C[key]]+=1
extensions={}
for name,core in [('C5',c5),('Q5',q5)]:
    ct={'exit':0,'entry':0,'both':0,'neither':0}
    for mask in range(32):
        a=extension(core,{j for j in range(5) if mask>>j&1});e,b=configurations(a)
        ct['exit']+=bool(e);ct['entry']+=bool(b);ct['both']+=bool(e and b);ct['neither']+=not(e or b)
        require(e or b,'core extension lacks entry and exit')
    extensions[name]=ct
result=dict(explicit_balanced_exits_checked=4,five_core_complete_counts=counts,all_external_vertex_extensions=extensions)
(ROOT/'hand_supplement_checks.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
