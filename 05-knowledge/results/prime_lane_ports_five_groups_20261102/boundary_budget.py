#!/usr/bin/env python3
"""Exact controls for the independent-module attachment boundary budget.

No repository imports. Complete one-vertex attachments to k ordered directed
triangles for k=1,2,3. Literal path/stratum controls at k=2; subset DP elsewhere.
Sharp multi-vertex examples use explicit ordered joins of 7-vertex H41 blocks.
"""
from collections import Counter
from itertools import permutations
from math import prod
import json


def need(test, message):
    if not test:
        raise RuntimeError(message)


def triangles(k):
    n=3*k
    rows=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if i//3 != j//3 or (j-i)==1:
                rows[i] |= 1<<j
            else:
                rows[j] |= 1<<i
    return rows


def append_vertex(rows, beats):
    n=len(rows)
    out=list(rows)+[beats]
    for i in range(n):
        if not (beats>>i)&1:
            out[i] |= 1<<n
    return out


def hp(rows):
    n=len(rows)
    if not n:return 1
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    full=(1<<n)-1
    for mask in range(1,full+1):
        for end,value in enumerate(dp[mask]):
            if not value:continue
            choices=rows[end] & (full^mask)
            while choices:
                bit=choices & -choices
                choices-=bit
                dp[mask|bit][bit.bit_length()-1] += value
    return sum(dp[-1])


def strata(rows,k):
    old=3*k
    totals=Counter()
    for word in permutations(range(len(rows))):
        if not all((rows[a]>>b)&1 for a,b in zip(word,word[1:])):continue
        touched=0
        for a,b in zip(word,word[1:]):
            if a<old<=b:touched |= 1<<(a//3)
            if b<old<=a:touched |= 1<<(b//3)
        totals[touched]+=1
    return totals


def join(blocks):
    out=[]
    n=sum(map(len,blocks))
    offset=0
    for rows in blocks:
        end=offset+len(rows)
        later=((1<<n)-1)^((1<<end)-1)
        out.extend((row<<offset)|later for row in rows)
        offset=end
    return out


def val3(n):
    v=0
    while n%3==0:
        v+=1;n//=3
    return v


def main():
    records=[]
    for k in (1,2,3):
        old=triangles(k)
        hist=Counter()
        minimum=100
        stratum_checks=0
        for cut in range(1<<(3*k)):
            rows=append_vertex(old,cut)
            h=hp(rows)
            need(h%3**max(0,k-2)==0,'one-birth divisibility')
            hist[h]+=1
            minimum=min(minimum,val3(h))
            if k==2:
                st=strata(rows,k)
                need(sum(st.values())==h,'literal path vs DP')
                for touched,count in st.items():
                    need(touched.bit_count()<=2,'physical boundary budget')
                    need(count%3**(k-touched.bit_count())==0,'stratum orbit divisor')
                    stratum_checks+=1
        records.append({'k':k,'cuts':1<<(3*k),'min_v3':minimum,
                        'histogram':dict(sorted(hist.items())),
                        'literal_stratum_checks':stratum_checks})
    seed=append_vertex(triangles(2),(1<<0)|(1<<3))
    need(hp(seed)==41,'H41 seed')
    sharp=[]
    for k,r in ((1,0),(2,1),(3,1),(4,2),(5,2)):
        rows=join([seed]*r+[triangles(1)]*(k-2*r))
        expected=41**r*3**(k-2*r)
        actual=hp(rows)
        need(actual==expected,'sharp join count')
        # Delete each new vertex (last in each 7-block), retain actual labels.
        outside={7*j+6 for j in range(r)}
        oldlabels=[i for i in range(len(rows)) if i not in outside]
        induced=[sum(((rows[a]>>b)&1)<<j for j,b in enumerate(oldlabels))
                 for a in oldlabels]
        need(induced==triangles(k),'source graph changed')
        need(val3(actual)==k-2*r,'sharp boundary valuation')
        sharp.append({'k':k,'r':r,'order':len(rows),'H':actual,'v3':val3(actual),
                      'out_rows':rows,'new_vertices':sorted(outside)})
    print(json.dumps({'status':'FINITE-EXACT','complete_one_vertex':records,
                      'sharp_multi_vertex':sharp},indent=2,sort_keys=True))


if __name__=='__main__':main()
