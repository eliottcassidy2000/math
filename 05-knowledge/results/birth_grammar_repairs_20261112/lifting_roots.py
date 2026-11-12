"""Exact rational/irrational distinction for the65 two-odd-mode lifting roots."""
from common import *
from fractions import Fraction
from collections import Counter

def power3(q):
    if q<=0:return None
    a=q.numerator;b=q.denominator;k=0
    while a%3==0:a//=3;k+=1
    while b%3==0:b//=3;k-=1
    return k if a==b==1 else None

def run():
    out=[];hist=Counter()
    for r in read('even_valuations.json.gz')['records']:
        if r['type']!='unbounded':continue
        f={int(b):c for b,c in r['normalized_odd_modes'].items()}
        if 5 in f:continue
        ratio=Fraction(-f.get(1,0),f[3]);a=power3(ratio)
        for p in r['unbounded_parities']:
            integer=a is not None and a%2==p
            kind='integer' if integer else 'irrational_in_Z2';hist[kind]+=1
            out.append({'birth_id':r['birth_id'],'port':r['port'],'parity':p,'three_power_ratio':[ratio.numerator,ratio.denominator],'root_type':kind,'a_root':a if integer else None,'t_root':(a-p)//2 if integer else None})
    need(len(out)==65,'all two-mode unbounded branches')
    save('lifting_roots.json.gz',{'scope':'The unique t in Z2 solving C1+C3*3^p*9^t=0 for the65 two-mode unbounded branches. Rational t is equivalent to integer t by prime valuations. The one mixed3/5 branch is excluded from this rationality classification.','root_histogram':dict(hist),'records':out});print(json.dumps({'root_histogram':dict(hist),'integer_roots':[r for r in out if r['root_type']=='integer']}),flush=True)
if __name__=='__main__':run()
