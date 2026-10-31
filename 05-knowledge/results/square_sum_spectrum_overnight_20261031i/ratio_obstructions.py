#!/usr/bin/env python3
"""Exact source-collision charts; interval arithmetic is a Farey rectangle."""
from fractions import Fraction
from math import gcd,isqrt
from pathlib import Path
import json
from one_copy_ports import need

def rectangle(n,m):
    return {Fraction(L,R) for L in range(1,isqrt(m)+1) for R in range(1,isqrt(n)+1) if gcd(L,R)==1}

def union_ratios(left,right):
    return {Fraction(b,a)*t for a,n in left for b,m in right for t in rectangle(n,m)}

def literal(left,right):
    result=set()
    for a in left:
        for b in right:
            g=gcd(a,b);u,v=b//g,a//g;U,V=isqrt(u),isqrt(v)
            if U*U==u and V*V==v:result.add(Fraction(U,V))
    return result

def vertices(charts):
    pieces=[{a*a*i for i in range(1,n+1)} for a,n in charts]
    need(sum(map(len,pieces))==len(set().union(*pieces)),'disjoint chart ownership')
    return set().union(*pieces)

def main():
    checks=0
    for n in range(1,65):
        for m in [1,7,16,31,64]:
            need(rectangle(n,m)==literal(range(1,n+1),range(1,m+1)),'complete interval ratio set');checks+=1
    left=[(1,12),(5,7)];right=[(2,8),(11,5)]
    LV,RV=vertices(left),vertices(right)
    need(union_ratios(left,right)==literal(LV,RV),'scaled-chart ratio union')
    native=rectangle(2812,2861)
    phi_sum=sum(sum(gcd(a,k)==1 for a in range(1,k+1)) for k in range(1,54))
    need(len(native)==2*phi_sum-1,'Farey totient identity')
    # Square-class reduction is exact: only equal-class source pairs can collide.
    pair_count=sum(isqrt(a*b)**2==a*b for a in range(1,65) for b in range(1,65))
    need(len(rectangle(64,64))<=pair_count,'distinct ratios versus repeated owner pairs')
    composition=union_ratios([(117,2812),(20,2861)],[(1,2812)])
    print(json.dumps({'status':'passed','complete_interval_controls':checks,'scaled_chart_controls':1,
          'native_interval_pair':[2812,2861],'native_forbidden_ratios':len(native),'native_frame_bound':len(native)+1,
          'raw_pair_bound':2812*2861+1,'totient_sum_through53':phi_sum,
          'second_join_forbidden_ratios':len(composition),'small_equal_class_pairs':pair_count,
          'scope':'Exact finite bad-ratio sets for disjoint square-scaled interval charts; the charts retain source labels and multiplicities without factoring the larger output integers.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
