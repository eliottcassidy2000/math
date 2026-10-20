#!/usr/bin/env python3
"""Separate algorithms, no imported producer code; not an independent agent audit."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import factorial
import json

CHECKS=Counter()
def need(ok,label):
    CHECKS[label]+=1
    if not ok:raise RuntimeError(label)


def triangle_matrix(sign, matching=(0,1,2), complement=False):
    matrix=[[False]*6 for _ in range(6)]
    for i in range(3):
        matrix[i][(i+1)%3]=True
        matrix[3+i][3+(i+sign)%3]=True
        for j in range(3):
            matrix[i][3+j]=(matching[i]==j)^complement
            matrix[3+j][i]=not matrix[i][3+j]
    return matrix


def finite_centers():
    rows=[]
    for sign in (1,-1):
        for matching in permutations(range(3)):
            for complement in (False,True):
                a=triangle_matrix(sign,matching,complement)
                orders=list(permutations(range(6)))
                aut=[p for p in orders if all(a[i][j]==a[p[i]][p[j]]
                                               for i in range(6) for j in range(6))]
                anti=[p for p in orders if all(a[i][j]!=a[p[i]][p[j]]
                                                for i in range(6) for j in range(6) if i!=j)]
                feedback=min(sum(a[p[j]][p[i]] for i in range(6) for j in range(i+1,6))
                             for p in orders)
                paths=[p for p in orders if all(a[p[i]][p[i+1]] for i in range(5))]
                blue=[]
                for p in paths:
                    r=dict(zip(p,p[::-1]))
                    if all(a[i][j]!=a[r[i]][r[j]] for i in range(6) for j in range(6) if i!=j):
                        blue.append(p)
                need(len(aut)==3 and len(anti)==3 and feedback==4,
                     'all 24 matching-oriented centers exact feedback and symmetry')
                need(len(paths)==45,'all matching-oriented centers have 45 Hamilton paths')
                # Matching parity changes the relative cyclic orientation.
                parity=(-1)**sum(matching[i]>matching[j] for i in range(3) for j in range(i+1,3))
                orientation=sign*parity
                need(len(blue)==(27 if orientation==1 else 9),
                     'matching orientation changes blue response by factor three')
                cycle_types=Counter()
                for p in aut:
                    used=set();sizes=[]
                    for v in range(6):
                        if v in used:continue
                        x=v;length=0
                        while x not in used:used.add(x);length+=1;x=p[x]
                        sizes.append(length)
                    cycle_types[tuple(sorted(sizes))]+=1
                need(cycle_types=={(1,1,1,1,1,1):1,(3,3):2},
                     'central cycle index has one shared rotation')
                for ports in combinations(range(6),2):
                    need(sum(all(p[v]==v for v in ports) for p in aut)==1,
                         'any central pair of ports pins the shared phase')
                rows.append({'sign':sign,'matching':matching,'complement':complement,
                             'feedback':feedback,'H':len(paths),'blue_paths':len(blue)})
    return rows


def falling(n,k):
    if n<k:return 0
    return factorial(n)//factorial(n-k)


def transfer_rooks():
    # Inclusion-exclusion over forbidden entries, rather than matrix squaring.
    rows=[]
    for k in range(3,10):
        a=(0,1,2);counts=[]
        for b in permutations(range(k),3):
            forbidden=[{a[i],b[i]} for i in range(3)]
            total=0
            for mask in range(8):
                selected=[i for i in range(3) if mask>>i&1]
                bad_injections=sum(len(set(xs))==len(xs)
                                   for xs in product(*(forbidden[i] for i in selected)))
                total+=(-1)**len(selected)*bad_injections*falling(k-len(selected),3-len(selected))
            direct=sum(all(c[i] not in forbidden[i] for i in range(3))
                       for c in permutations(range(k),3))
            need(total==direct,'boundary transfer rook count equals literal extensions')
            need(k==3 or total>0,'one spare colour destroys ternary orientation obstruction')
            counts.append(total)
        rows.append({'colours':k,'fixed_left_states':len(counts),
                     'minimum':min(counts),'maximum':max(counts),
                     'supported':sum(x>0 for x in counts)})
    return rows


def compose(word):
    A,C=F(1),F(0)
    for a in word:A,C=3*A/2**a,(3*C+1)/2**a
    return A,C


def odd_step(n):
    x=3*n+1;a=0
    while x%2==0:x//=2;a+=1
    return x,a


def cylinder_exhaustion():
    words=0;sources=0
    for m in range(1,5):
        for word in product((1,2,3),repeat=m):
            K=sum(word);A,C=compose(word)
            B=C*2**K
            need(B.denominator==1,'independent composition has integral carry')
            predicted=((2**K-int(B))*pow(3**m,-1,2**(K+1)))%2**(K+1)
            found=[]
            for n in range(1,2**(K+1),2):
                x=n;ok=True
                for a in word:
                    x,b=odd_step(x)
                    if a!=b:ok=False;break
                if ok:found.append(n)
                sources+=1
            need(found==[predicted],'complete residue universe gives exactly the predicted cylinder')
            words+=1
    return {'words':words,'odd_residues_tested':sources}


def envelope(word,n,cuts):
    x=F(n);start=0
    for stop in cuts:
        A,C=compose(sorted(word[start:stop],reverse=True))
        x=A*x+C;start=stop
    return x


def optimal_blocks(word,n,cap=None):
    m=len(word);cap=m if cap is None else min(cap,m)
    maps={(i,j):compose(sorted(word[i:j],reverse=True))
          for j in range(1,m+1) for i in range(j)}
    dp={(0,0):(F(n),())};row=[]
    for b in range(1,cap+1):
        for j in range(b,m+1):
            dp[b,j]=min((maps[i,j][0]*dp[b-1,i][0]+maps[i,j][1],dp[b-1,i][1]+(j,))
                        for i in range(b-1,j) if (b-1,i) in dp)
        row.append(dp[b,m])
        if dp[b,m][0]<n:break
    return row


def refinements():
    # Every contiguous partition of every short word; all immediate splits.
    tested=0
    for m in range(1,7):
        for word in product((1,2,3),repeat=m):
            n=19;values={}
            A,C=compose(word);actual=A*n+C
            for mask in range(1<<(m-1)):
                cuts=tuple(i for i in range(1,m) if mask>>(i-1)&1)+(m,)
                values[mask]=envelope(word,n,cuts)
                need(values[mask]>=actual,'block envelope dominates literal formal endpoint')
                tested+=1
            for mask,value in values.items():
                for i in range(m-1):
                    need(values[mask|(1<<i)]<=value,'every splitting refinement lowers source envelope')
            # Independent all-partition optimum checks the dynamic program.
            rows=optimal_blocks(word,n)
            for b,(value,cuts) in enumerate(rows,1):
                brute=min(v for mask,v in values.items() if mask.bit_count()==b-1)
                need(value==brute,'block optimum DP equals exhaustive partition search')
    examples=[]
    for n in (7,27,31,63,75,127,255,703,871,6171):
        x=n;word=[];first=None
        while x!=1:
            x,a=odd_step(x);word.append(a)
            need(len(word)<=1000,'declared source route cap')
            if first is None and x<n:first=(len(word),x)
        m,endpoint=first
        rows=optimal_blocks(word[:m],n)
        need(rows[-1][0]<n,'listed actual first-exit source has a paid block refinement')
        examples.append({'source':n,'first_descent_steps':m,'first_descent_endpoint':endpoint,
                         'minimum_blocks':len(rows),'cuts':rows[-1][1],
                         'upper':str(rows[-1][0]),'word':word[:m]})
        if n in (27,63):
            # Exhaust all one-fewer-block partitions of these longer words.
            b=len(rows)-1
            maps={(i,j):compose(sorted(word[i:j],reverse=True))
                  for j in range(1,m+1) for i in range(j)}
            best=None;partitions=0
            for interior in combinations(range(1,m),b-1):
                value=F(n);start=0
                for end in interior+(m,):
                    A,C=maps[start,end];value=A*value+C;start=end
                best=value if best is None else min(best,value);partitions+=1
            need(best==rows[-2][0]>=n,'long-source fewer-block optimum by exhaustive cut enumeration')
            examples[-1]['independent_exhausted_partitions']=partitions
            examples[-1]['best_one_fewer_upper']=str(best)
        if n==63:
            need(len(rows)==5 and rows[-1][1]==(5,14,26,30,34),
                 'source 63 refutes a four-block first-exit rule')
            need(rows[-1][0]==F(1095707345239597,17592186044416),
                 'source 63 exact paid five-block endpoint')
        if n==27:
            need(len(rows)==4 and rows[-1][1]==(6,16,32,37),
                 'source 27 requires exactly four blocks on its first descent')
            need(rows[-1][0]==F(7067832544468979611,288230376151711744),
                 'source 27 exact paid four-block endpoint')
            need(all(envelope(word[:j],n,(j,))>=n for j in range(1,len(word)+1)),
                 'source 27 never paid by one block before reaching one')
            u=F(3,2)**word.count(1)
            limit=2*u-1
            need(limit==F(282421147873,8388608)>27,
                 'source 27 one-block all-future hostile limit')
            base=envelope(word,n,(len(word),))
            for r in range(1,21):
                extended=word+[2]*r
                need(envelope(extended,n,(len(extended),))==limit+F(3,4)**r*(base-limit),
                     'all-future one-block hostile exact contraction to unpaid limit')
    return {'all_short_partitions':tested,'source_examples':examples}


def main():
    data={'status':'FINITE-EXACT, separate algorithms by same agent',
          'centers':finite_centers(),'transfer_rooks':transfer_rooks(),
          'cylinders':cylinder_exhaustion(),'refinements':refinements()}
    data['checks']=dict(CHECKS);data['total_checks']=sum(CHECKS.values())
    print(json.dumps(data,indent=2,sort_keys=True))


if __name__=='__main__':main()
