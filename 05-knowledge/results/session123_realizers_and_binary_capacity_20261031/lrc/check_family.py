from fractions import Fraction as F
from itertools import combinations
import json
from family_selector import choose_family,attach,R
from witness import norm,need

def full_safe(speeds):
    out=[(F(0),F(1))]
    for v in speeds:out=attach(out,v)
    return out

def main():
    patterns=[(1,1,2,3,13),(14,15,16,17,18),(14,14,14,14,14),(14,27,65,100,101),tuple(2**128+x for x in (1,3,7,11,19))]
    total=literal=0;methods={};maxexpanded=0;maxcut=0;maxprefix=0
    for core in combinations(range(1,14),8):
        for tails in patterns:
            t,st=choose_family(core,tails);need(all(norm(v*t)>=R for v in (*core,*tails)),'actual original row')
            if max(tails)<=18:
                need(any(a<=t<=b for a,b in full_safe((*core,*tails))),'literal full row');literal+=1
            total+=1;methods[st['method']]=methods.get(st['method'],0)+1
            maxexpanded=max(maxexpanded,st['max_expanded_speed']);maxcut=max(maxcut,max(st['cutoffs'],default=0));maxprefix=max(maxprefix,st['finite_tail_intersections'])
    huge=[]
    for bits in (64,256,1024,4096):
        tails=[14,27,2**bits+1,3*2**bits+7,5*2**bits+9]
        t,st=choose_family(tuple(range(1,9)),tails)
        huge.append({'bits':bits,'witness_bits':[t.numerator.bit_length(),t.denominator.bit_length()],'expanded_max':st['max_expanded_speed'],'finite_prefix':st['finite_tail_intersections'],'method':st['method']})
    # All larger bounded-core sizes, including the no-tail base case.
    bigger=0
    for m in range(9,14):
        for core in combinations(range(1,14),m):
            tails=[2**80+17]*(13-m);choose_family(core,tails);bigger+=1
    print(json.dumps({'status':'PASS','eight_core_count':1287,'rows_checked':total,'literal_full_cover_comparisons':literal,'method_counts':methods,'observed_max_expanded_speed':maxexpanded,'observed_max_cutoff':maxcut,'observed_max_finite_prefix':maxprefix,'larger_core_controls':bigger,'large_input_controls':huge,'scope':'finite controls of the implemented selector; universal family theorem remains inherited, and observed maxima are not universal claims'},indent=2,sort_keys=True))
if __name__=='__main__':main()
