"""Cross-audit exact receipts, sign guards, small strong words and arithmetic gates."""
from itertools import product
from math import gcd,factorial
import json,gzip
import numpy as np
from typed_glue import HERE,need,save,glue,diff
from port_compiler import mul
from repair23 import profile,literal_profile,strong


def read(name):return json.loads(gzip.decompress((HERE/(name+'.json.gz')).read_bytes()))


def run():
    pairs=read('typed_glue')['pairs'];cone=[]
    for p in pairs:
        d=p['difference'];sign=1 if d[3]>0 else -1;d=[sign*x for x in d]
        need(d[3]>0 and d[6]>0 and d[5]>d[3] and d[5]>d[2],'all inherited pairs in strict cone')
        cone.append(dict(label=p['label'],oriented_difference=d))
    head=read('context_head')['levels'];need([x['neutral'] for x in head]==[2,3,5,13,33,131,994],'exact sequence receipt')
    need(len(read('context_head')['first_rules'])==22,'fixed-root transition count')
    need(len(read('port_search')['rules'])==22,'separate four-input two-port transition count')
    need(len(pairs)==11,'eleven actual input pairs')
    # Literal U/W words, with strongness and independent permutation seed count.
    seed=((2,4,9,3),(2,8,3,5));zs=[profile(T,3) for T in seed]
    for T,z in zip(seed,zs):need(literal_profile(T,3)==z and strong(T),'independent H5 seed')
    cat=read('renewal_catalog')['records'];rules=[next(r for r in cat if r['n']==n) for n in (4,5)]
    words=[]
    for word in [()] + [(0,),(1,)]+list(product(range(2),repeat=2)):
        rows=seed;r=3;states=zs
        for i in word:
            c=rules[i];rows=[glue(c['rows'],c['entry'],T,r) for T in rows];r=c['exit'];states=[mul(c['matrix'],z) for z in states]
        actual=[profile(T,r) for T in rows]
        need(states==actual and all(strong(T) for T in rows),'literal strong small word')
        words.append(dict(word=word,n=len(rows[0]),profiles=actual))
    gates=0
    for h in range(1,101):
        for E in range(1,101):
            if gcd(h,E)!=1:continue
            need(gcd(h+4*E,2*h+7*E)==1,'U primitive gate')
            need(gcd(2*h+9*E,3*h+16*E)==(5 if (h+2*E)%5==0 else 1),'W exact factor5 gate')
            gates+=1
    lang=read('language');need(lang['distinct_H']==1346 and not lang['obstruction_in_language'],'complete small language')
    need(lang['longest_word']==8,'word length dynamic program')
    col=read('collatz_graphs');need(col['greedy']['certified']==100000 and not col['greedy']['unresolved'],'literal greedy census')
    need(col['greedy']['largest_peak']==dict(start=65535,peak=2*3**16-2),'sharp forced family control')
    walk=read('walk_controls');need(set(walk['oriented_gap_witnesses'])=={'7','21'},'one-arc oriented controls')
    need(factorial(18)<2**63,'int64 safety for largest literal DP universe')
    save('audit',dict(status='PASS',scope='Independent permutation seed count; seven literal U/W words; all inherited sign-cone guards; exact arithmetic gcd gates; cross-receipt counts and scope checks. Mathematical proof audit is recorded in AUDIT.md, not delegated.',cone=cone,literal_words=words,primitive_gate_controls=gates,largest_literal_DP_order=18,largest_factorial=factorial(18)))
    print('AUDIT PASS',len(cone),'sign pairs',len(words),'literal words',gates,'gcd gates')


if __name__=='__main__':run()
