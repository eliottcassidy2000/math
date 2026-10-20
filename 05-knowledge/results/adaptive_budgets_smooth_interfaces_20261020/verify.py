#!/usr/bin/env python3
"""Exact controls; finite universes are printed, no assert-dependent checks."""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import gcd
import argparse
import json

CHECKS = 0
def need(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)

def smooth(n, p):
    for d in range(2, p+1):
        while n % d == 0:
            n //= d
    return n == 1

def phi(n):
    return sum(gcd(n, k) == 1 for k in range(1, n+1))

def determinant(a):
    a = [list(row) for row in a]
    previous = 1
    for k in range(len(a)-1):
        pivot = a[k][k]
        need(pivot != 0, 'positive leading principal pivot')
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                numerator = a[i][j]*pivot-a[i][k]*a[k][j]
                need(numerator % previous == 0, 'Bareiss exact division')
                a[i][j] = numerator//previous
        previous = pivot
    return a[-1][-1] if a else 1

def bank_controls():
    rows = []
    for p, q, z in ((3,5,27),(5,7,25),(7,11,21),(11,13,19)):
        bank = [c for c in range(1,z+1,2) if smooth(c,p)]
        S = sum(bank); horizon = (q+1)*S
        G = [[gcd(a,b) for b in bank] for a in bank]
        for i,a in enumerate(bank):
            for j,b in enumerate(bank):
                need(G[i][j] == sum(phi(d) for d in bank if a%d==b%d==0),
                     'divisor-incidence Gram factorization')
        need(determinant(G) == __import__('functools').reduce(lambda a,c:a*phi(c),bank,1),
             'Smith determinant on divisor-closed smooth bank')
        profiles = set(); cells = defaultdict(list)
        for mask in range(1 << len(bank)):
            A = [bank[i] for i in range(len(bank)) if mask>>i&1]
            degree = q*sum(A)
            response = tuple(sum(gcd(q*a,b) for a in A) for b in bank)
            need(all(b <= horizon-degree for b in bank), 'future probes affordable')
            need(response not in profiles, 'complete response distinguishes every history')
            profiles.add(response); cells[sum(A)].append(A)
        largest = max(cells.values(),key=len)
        # Exact SPD control on one equal-degree cell, plus a full five-bit bank.
        F = largest[:20]
        K = [[2**sum(gcd(a,b) for a in A for b in B) for B in F] for A in F]
        need(determinant(K)>0, 'equal-degree continuation kernel full rank')
        rows.append({'p':p,'q':q,'z':z,'bank':bank,'horizon':horizon,
                     'profiles':len(profiles),'largest_equal_degree_cell':len(largest),
                     'tested_kernel_order':len(F)})
    return rows

def odd_step(n):
    v = 3*n+1
    a = (v & -v).bit_length()-1
    return v >> a, a

def compose(word):
    B=0; K=0
    for a in word:
        B=3*B+(1<<K); K+=a
    return 3**len(word),B,K

def block_table(word):
    m=len(word); table={}
    for i in range(m):
        counts=Counter()
        for j in range(i+1,m+1):
            counts[word[j-1]]+=1
            B=0; K=0
            for a,r in sorted(counts.items(),reverse=True):
                Br=(3**r-(1<<(a*r)))//(3-(1<<a))
                B=3**r*B+(Br<<K); K+=a*r
            table[i,j]=B
    return table

def optimal(word,n):
    m=len(word); P=[0]
    for a in word:P.append(P[-1]+a)
    A,B,K=compose(sorted(word,reverse=True)); target=n<<K
    if A*n+B<target:return [(A*n+B,(m,))]
    maps=block_table(word); powers=[3**j for j in range(m+1)]
    prev={0:(n,())}; rows=[]
    for b in range(1,m+1):
        cur={}
        for j in range(b,m+1):
            cur[j]=min((powers[j-i]*value+(maps[i,j]<<P[i]),cuts+(j,))
                       for i,(value,cuts) in prev.items() if i<j)
        rows.append(cur[m])
        if cur[m][0]<target:return rows
        prev=cur
    raise RuntimeError('actual descent must be certified at singleton refinement')

def source_scan(limit):
    histogram=Counter(); records=[]; hard=[]; max_b=0; max_m=0
    for n in range(3,limit+1,2):
        x=n;word=[]
        while x>=n:
            x,a=odd_step(x);word.append(a)
            if len(word)>1000:raise RuntimeError('declared 1000-step first-exit cap')
        rows=optimal(word,n); b=len(rows); histogram[b]+=1
        A,B,K=compose(word)
        need((A*n+B)==x*(1<<K),'integer-clock endpoint')
        value,cuts=rows[-1]
        need(value<n*(1<<K),'paid block endpoint')
        # Replay the actual selected block envelope, independent of DP numerator transitions.
        v=Q(n);start=0
        for end in cuts:
            a,c,k=compose(sorted(word[start:end],reverse=True))
            v=(a*v+c)/(1<<k);start=end
        need(v==Q(value,1<<K),'chosen cuts direct rational replay')
        if b>1:
            hard.append({'source':n,'minimum_blocks':b,'cuts':cuts,'upper':str(v)})
        if b>max_b or len(word)>max_m:
            records.append({'source':n,'steps':len(word),'clock':K,'endpoint':x,
                            'minimum_blocks':b,'cuts':cuts,'upper':str(v),'word':word})
        max_b=max(max_b,b);max_m=max(max_m,len(word))
    return {'all_odd_sources':[3,limit],'sources':sum(histogram.values()),
            'first_exit_cap':1000,'histogram':dict(sorted(histogram.items())),
            'maximum_blocks':max_b,'records':records,'hard_sources':hard}

def padding_controls():
    records=[]
    for n in (3,7,27,31,63,75,127,703,871,6171,35655,63728127):
        x=n;w=[]
        while x!=1:
            x,a=odd_step(x);w.append(a)
            if len(w)>2000:raise RuntimeError('declared terminal cap')
        A,B,K=compose(sorted(w,reverse=True)); U=Q(A*n+B,1<<K)
        V=U;r=0
        while V>=n:V=(3*V+1)/4;r+=1
        need(V==1+Q(3,4)**r*(U-1),'terminal padding exact debt decay')
        need(r<=2*len(w)+4,'coarse explicit linear clock bound')
        # Retain the actual source cylinder after all padding steps.
        a,b,k=compose(w+[2]*r)
        residue=((1<<k)-b)*pow(a,-1,1<<(k+1)) % (1<<(k+1))
        need(n%(1<<(k+1))==residue,'source remains in the full exact cylinder')
        records.append({'source':n,'steps_to_one':len(w),'clock_to_one':sum(w),
                        'padding_steps':r,'total_integer_clock':sum(w)+2*r,
                        'total_unaccelerated_clock':len(w)+sum(w)+3*r,
                        'coarse_upper':str(U),'paid_upper':str(V)})
    return records

def distance2(a,b):
    # Coordinates (x,y) mean (x/2, sqrt(3)*y/2).
    return Q((a[0]-b[0])**2+3*(a[1]-b[1])**2,4)

def geometry_controls():
    rim=[(2,0),(1,1),(-1,1),(-2,0),(-1,-1),(1,-1)]
    tips=[tuple(rim[2*i][j]+rim[2*i+1][j] for j in range(2)) for i in range(3)]
    pts=[(0,0)]+rim+tips
    edges={(i,j) for i,j in combinations(range(10),2) if distance2(pts[i],pts[j])==1}
    expected={(0,i) for i in range(1,7)}
    expected|={tuple(sorted((1+i,1+(i+1)%6))) for i in range(6)}
    expected|={(1+2*i,7+i) for i in range(3)}|{(2+2*i,7+i) for i in range(3)}
    need(edges==expected,'all 45 squared distances authenticate precisely 18 edges')
    triples=list(permutations(range(6),3)); records=Counter(); extension=Counter()
    projection=set(); all_pairs=set(); palettes={}; hostile=None
    for b in triples:
        for c in triples:
            a=(0,1,2)
            if a[1]==b[0] or b[1]==c[0] or c[1]==a[0]:continue
            ports=a[:2]+b[:2]+c[:2]
            count=6-len(set(ports)); extension[count]+=1
            records['boundary']+=1;records['extendible' if count else 'blocked']+=1
            if count:projection.add(b)
            key=tuple(tuple(sorted(t)) for t in (a,b,c))
            if key in palettes and bool(count)!=palettes[key][0] and hostile is None:
                hostile=[palettes[key][1],(a,b,c)]
            palettes[key]=(bool(count),(a,b,c))
            if a[1]!=b[0]:all_pairs.add(b)
    expected_pairs={b for b in triples if b[0]!=1}
    need(projection==expected_pairs,'two adjacent wings see full proper induced boundary relation')
    need(records=={'boundary':8336,'extendible':7952,'blocked':384},'complete normalized six-colour boundary census')
    need(hostile is not None,'wing palettes alone lose extension information')
    return {'coordinate_convention':'(x/2,sqrt(3)y/2)','points':pts,'edges':sorted(edges),
            'first_wing_fixed':[0,1,2],'boundary_counts':dict(records),
            'extension_count_histogram':dict(sorted(extension.items())),
            'two_wing_projection':len(projection),'same_palette_opposite_extension':hostile}

def main():
    p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=20001)
    args=p.parse_args()
    result={'status':'FINITE-EXACT controls; analytic all-order claims proved separately',
            'smooth_banks':bank_controls(),'geometry':geometry_controls(),
            'padding':padding_controls(),'collatz':source_scan(args.limit)}
    result['checks']=CHECKS
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()
