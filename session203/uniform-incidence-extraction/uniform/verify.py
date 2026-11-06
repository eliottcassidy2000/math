#!/usr/bin/env python3
"""Literal controls for uniform binary-incidence extraction obstruction."""
import itertools, json
from fractions import Fraction


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def rows(n, encoding):
    ans=[0]*n
    pos=0
    for i in range(n):
        for j in range(i+1,n):
            if encoding>>pos&1: ans[i]|=1<<j
            else: ans[j]|=1<<i
            pos+=1
    return ans


def hp(a):
    n=len(a)
    if n==0:return 1
    d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for mask in range(1,1<<n):
        for v in range(n):
            q=d[mask][v]
            if q:
                rem=a[v]&~mask
                while rem:
                    lb=rem&-rem;w=lb.bit_length()-1
                    d[mask|lb][w]+=q;rem-=lb
    return sum(d[-1])


def induced(a,verts):
    return [sum(1<<j for j,v in enumerate(verts) if a[u]>>v&1) for u in verts]


def add_apex(a,bits):
    n=len(a);r=a[:]+[0]
    for i in range(n):
        if bits>>i&1:r[n]|=1<<i
        else:r[i]|=1<<n
    return r


def qmatrix(a):
    n=len(a);q=[[0]*n for _ in range(n)]
    for p in itertools.permutations(range(n)):
        bad=[i for i in range(n-1) if not(a[p[i]]>>p[i+1]&1)]
        if len(bad)>1:continue
        for i in range(n-1):
            if not bad or bad==[i]:q[p[i]][p[i+1]]+=1
    return q


def two_apex(m,arc,pa,pb,bits):
    a=[0]*(m+2)
    a[0 if arc else 1]|=1<<(1 if arc else 0)
    for i in range(m):
        x=2+i;b=bits>>i&1
        for v,p in ((0,pa),(1,pb)):
            if b^p:a[v]|=1<<x
            else:a[x]|=1<<v
        for j in range(i+1,m):a[x]|=1<<(2+j)
    return a


def numbins(m):
    return 3**m-(3*2**m+m*2**(m-1))+2*m+3


def main():
    parent_count=cut_count=pair_count=0
    for n in range(1,5):
        for enc in range(1<<(n*(n-1)//2)):
            a=rows(n,enc);q=qmatrix(a);h=hp(a);parent_count+=1
            for i in range(n):
                for j in range(i+1,n):
                    lower=2*hp(induced(a,[v for v in range(n) if v not in (i,j)]))
                    require(q[i][j]+q[j][i]>=lower,'capacity lower bound')
                    pair_count+=1
            for bits in range(1<<n):
                expected=h+sum((q[i][j]+q[j][i]) for i in range(n) for j in range(n)
                               if a[i]>>j&1 and bits>>i&1 and not(bits>>j&1))
                require(hp(add_apex(a,bits))==expected,'cut polynomial')
                cut_count+=1
    mean_results=[];graph_count=0
    for m in range(2,9):
        for arc,pa,pb in itertools.product(range(2),repeat=3):
            total=sum(hp(two_apex(m,arc,pa,pb,bits)) for bits in range(1<<m))
            graph_count+=1<<m
            require(total*16>=numbins(m)*(1<<m),'two-apex mean lower bound')
            mean_results.append({'m':m,'arc':arc,'polarities':[pa,pb],
                                 'mean_numerator':total,'mean_denominator':1<<m,
                                 'bound_numerator':numbins(m),'bound_denominator':16})
    bin_inventories=[];edge_constraint_checks=0
    for m in range(2,8):
        count=0
        for assignment in itertools.product(range(3),repeat=m):
            bins=[[i for i in range(m) if assignment[i]==j] for j in range(3)]
            if not bins[0] or len(bins[1])<2 or not bins[2]:continue
            count+=1
            endpoints=[bins[0][-1],bins[1][0],bins[1][-1],bins[2][0]]
            require(len(set(endpoints))==4,'four distinct constraints')
            # Every polarity pair sets exactly one value at each endpoint.
            for pa,pb in itertools.product(range(2),repeat=2):
                constraints=[pa,1^pa,pb,1^pb]
                good=sum(all(((bits>>i)&1)==v for i,v in zip(endpoints,constraints))
                         for bits in range(1<<m))
                require(good*16==1<<m,'independence probability')
                edge_constraint_checks+=1
        require(count==numbins(m),'three-bin count')
        bin_inventories.append([m,count])
    require(numbins(13)==1516528,'m13 count')
    require(Fraction(numbins(13),16)==94783,'m13 average bound')
    require(6*Fraction(3*2**13+1,2)+1==73732,'m13 plus average')
    require(94783>73732,'m13 contradiction')
    print(json.dumps({'status':'PASS','scope':'Finite literal controls; all-order proof in PROOF.md',
          'parents':parent_count,'cuts':cut_count,'capacity_pairs':pair_count,
          'two_apex_graphs':graph_count,'bin_inventories':bin_inventories,
          'edge_constraint_checks':edge_constraint_checks,'m13':{'N':1516528,
          'lower_mean':94783,'affine_plus_mean':73732,'target_interval':[8193,16384]},
          'mean_checks':mean_results},indent=2,sort_keys=True))

if __name__=='__main__':main()
