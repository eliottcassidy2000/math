"""Independent recursive path counts for source and sampled hostile partners.
The complete second-moment sum is computed only by hostile.cpp.
"""
from functools import lru_cache
from fractions import Fraction
from pathlib import Path
import json

def need(ok,message):
    if not ok:raise RuntimeError(message)

def count_paths(rows):
    n=len(rows)
    @lru_cache(None)
    def finish(last,remaining):
        if not remaining:return 1
        options=rows[last]&remaining;total=0
        while options:
            bit=options&-options;options-=bit
            total+=finish(bit.bit_length()-1,remaining^bit)
        return total
    return sum(finish(i,((1<<n)-1)^(1<<i)) for i in range(n))

def source(minus):
    rows=[0]*12
    for i in range(12):
        for j in range(i+1,12):
            if i//4!=j//4:forward=(j//4-i//4)%3==1
            elif j%4==3:forward=not(minus and i//4==2)
            else:forward=(j-i)%3==1
            if forward:rows[i]|=1<<j
            else:rows[j]|=1<<i
    return tuple(rows)

def complemented_presentation(code):
    rows=[0]*12;k=0
    for i in range(12):
        for j in range(i+1,12):
            forward=j==i+1 or not(code>>k&1)
            if j>i+1:k+=1
            if forward:rows[i]|=1<<j
            else:rows[j]|=1<<i
    return tuple(rows)

def main():
    data=json.loads((Path(__file__).with_name('hostile.json')).read_text(encoding='utf-8-sig'))
    for i,r in enumerate(data):
        need(count_paths(source(i))==r['H']==462915,'independent source count')
        need(r['H']//r['presentations']==r['Aut'],'free presentation fibre')
        need(r['L2_one']==r['Aut']*r['sum_partner_H_over_presentations'],'orbit multiplicity')
        for code,h in r['sample']:need(count_paths(complemented_presentation(code))==h,'recursive partner control')
    difference=data[0]['L2_one']-data[1]['L2_one']
    need(difference==29097684,'second response separates 81/27')
    print(json.dumps({'independent_source_checks':2,'independent_partner_checks':12,'L2_one_difference':difference,
                      'mean_partner_H_difference':str(Fraction(difference,462915)),
                      'scope':'Full sums: C++ only; sources and twelve partner values: independent Python recursion.'},sort_keys=True,indent=2))
if __name__=='__main__':main()
