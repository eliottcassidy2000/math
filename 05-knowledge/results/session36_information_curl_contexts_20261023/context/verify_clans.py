#!/usr/bin/env python3
"""Exact triangle-detector, clan, SCC and curl-boundary checks."""
import itertools,json,math
from pathlib import Path
from verify_context import tournament,h_subsets,profile,check
from verify_quotient import walk_count,substitute

OUT=Path(__file__).resolve().parent

def triangle_key(a):
    return tuple(all(sum(a[v]>>w&1 for w in tr if w!=v)==1 for v in tr) for tr in itertools.combinations(range(len(a)),3))

def curls(a):
    s=lambda i,j:1 if a[i]>>j&1 else -1
    return tuple(s(i,j)+s(j,k)+s(k,i) for i,j,k in itertools.combinations(range(len(a)),3))

def from_order(order):
    a=[0]*len(order)
    for i,u in enumerate(order):
        a[u]=sum(1<<v for v in order[i+1:])
    return a

def compose(q,blocks):
    ends=[0]
    for a in blocks:ends.append(ends[-1]+len(a))
    out=[]
    for i,a in enumerate(blocks):
        extra=sum(((1<<len(b))-1)<<ends[j] for j,b in enumerate(blocks) if q[i]>>j&1)
        out.extend((row<<ends[i])|extra for row in a)
    return out

def reverse(a):
    n=len(a);return [((1<<n)-1)^(1<<i)^row for i,row in enumerate(a)]

def iso(a,b):
    n=len(a);da=[r.bit_count() for r in a];db=[r.bit_count() for r in b]
    for p in itertools.permutations(range(n)):
        if any(da[i]!=db[p[i]] for i in range(n)):continue
        if all((a[i]>>j&1)==(b[p[i]]>>p[j]&1) for i in range(n) for j in range(n)):return list(p)
    return None

def truncated_walks(a,depth):
    n=len(a);current={((0,)*i+(1,)+(0,)*(n-i-1),i):1 for i in range(n)};total={}
    for _ in range(1,depth+1):
        new={}
        for (c,i),v in current.items():
            total[c]=total.get(c,0)+v
            for j in range(n):
                if a[i]>>j&1:
                    d=list(c);d[j]+=1;key=(tuple(d),j)
                    new[key]=new.get(key,0)+v
        current=new
    return total

def main():
    stats={};detectorchecks=0
    for n in range(1,6):
        bank={};groups={}
        for mask in range(1<<(n*(n-1)//2)):
            a=tournament(n,mask);key=triangle_key(a);w=truncated_walks(a,7)
            if key in bank:check(w==bank[key],('same C3 different W truncation',n,mask))
            else:bank[key]=w
            groups[key]=groups.get(key,0)+1
            if n>=3:
                for tr,cyclic in zip(itertools.combinations(range(n),3),key):
                    c=tuple(n-1 if i in tr else 1 for i in range(n))
                    check((walk_count(a,c)>0)==cyclic,('triangle detector',n,mask,tr));detectorchecks+=1
        stats[n]={'tournaments':1<<(n*(n-1)//2),'C3_classes':len(bank),'walk_degree':7}
        print('C3 and W checks',n,stats[n],flush=True)
    stats['detector_checks']=detectorchecks
    t=from_order([0,1,2,3]);u=from_order([0,1,3,2])
    check(curls(t)!=curls(u) and curls(t)!=tuple(-v for v in curls(u)),'minimal transitive curl hostile')
    stats['transitive_hostile']={'rows':[t,u],'curls':[curls(t),curls(u)]}
    c3=[2,4,1]
    t=compose(c3,[[2,0],[0],[0]]);u=compose(c3,[[0,1],[0],[0]])
    check(curls(t)!=curls(u) and curls(t)!=tuple(-v for v in curls(u)),'strong labelled hostile')
    check(truncated_walks(t,10)==truncated_walks(u,10),'strong labelled W')
    stats['strong_labelled_hostile']={'rows':[t,u],'curls':[curls(t),curls(u)]}
    a4=compose([2,0],[[0],c3]);b2=[2,0]
    t=compose(c3,[a4,b2,[0]]);u=compose(c3,[reverse(a4),b2,[0]])
    check(iso(t,u) is None and iso(reverse(t),u) is None,'strong unlabelled nonconverse hostile')
    check(profile(t)==profile(u),'clan profiles')
    values=[]
    for sizes in [(1,)*7,(2,)*7,(1,2,1,1,2,1,1)]:
        pair=[h_subsets(substitute(v,sizes))[-1] for v in (t,u)]
        check(pair[0]==pair[1],('clan blowup',sizes));values.append({'sizes':sizes,'H':pair[0]})
    stats['strong_nonisomorphic_nonconverse_hostile']={'rows':[t,u],'curls':[curls(t),curls(u)],'profile':profile(t),'blowups':values}
    # Minimal labelled n<=3: same unsigned C3 forces curl equality up to global sign.
    for n in (1,2,3):
        groups={}
        for mask in range(1<<(n*(n-1)//2)):
            a=tournament(n,mask);groups.setdefault(triangle_key(a),[]).append(curls(a))
        for group in groups.values():
            for x in group:check(x==group[0] or x==tuple(-v for v in group[0]),'n<=3 minimality')
    (OUT/'clan_verification.json').write_text(json.dumps(stats,indent=2)+'\n')
    print('detector checks',detectorchecks)
    print('strong nonisomorphic nonconverse rows',t,u)
    print('ALL CLAN CHECKS PASSED')

if __name__=='__main__':main()
