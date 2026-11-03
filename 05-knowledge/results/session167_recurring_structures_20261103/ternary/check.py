#!/usr/bin/env python3
"""Self-contained exact checks; no external census or package dependencies."""
from itertools import combinations, permutations, product
from collections import Counter, deque
from math import comb, factorial
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
SEEDS=(153,177,181)
EXPECTED=(0,29,131,162,77,15,1)

def need(test, message):
    if not test:
        raise RuntimeError(message)

def matrix(n,mask):
    a=[[0]*n for _ in range(n)]
    for b,(i,j) in enumerate(combinations(range(n),2)):
        a[i][j]=(mask>>b)&1
        a[j][i]=1-a[i][j]
    return a

def encode(a):
    return sum(a[i][j]<<b for b,(i,j) in enumerate(combinations(range(len(a)),2)))

def relabel(n,mask,p):
    a=matrix(n,mask);b=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):b[p[i]][p[j]]=a[i][j]
    return encode(b)

def sign(p):
    return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))

def strong(a):
    n=len(a)
    for start in range(n):
        seen={start};todo=[start]
        while todo:
            i=todo.pop()
            for j in range(n):
                if a[i][j] and j not in seen:seen.add(j);todo.append(j)
        if len(seen)!=n:return False
    return True

def h_subsets(a):
    n=len(a);dp=[[0]*n for _ in range(1<<n)];h=[0]*(1<<n)
    for v in range(n):dp[1<<v][v]=1
    for s in range(1,1<<n):
        for v in range(n):
            if not (s>>v&1):continue
            r=s^(1<<v)
            if r:dp[s][v]=sum(dp[r][u] for u in range(n) if r>>u&1 and a[u][v])
        h[s]=sum(dp[s])
    return h

def profile_dp(a):
    n=len(a);h=h_subsets(a);dp=[[0]*(n+1) for _ in range(1<<n)];dp[0][0]=1
    for s in range(1,1<<n):
        bit=s&-s;rest=s^bit;t=rest
        while True:
            block=t|bit;r=s^block
            for j in range(r.bit_count()+1):dp[s][j+1]+=h[block]*dp[r][j]
            if not t:break
            t=(t-1)&rest
    return tuple(dp[-1])

def profile_permutations(a):
    n=len(a);p=[0]*(n+1)
    for word in permutations(range(n)):
        b=sum(not a[word[i]][word[i+1]] for i in range(n-1))
        for j in range(b+1,n+1):p[j]+=comb(n-1-b,j-1-b)
    return tuple(x//factorial(j) for j,x in enumerate(p))

def seed_checks():
    result={'masks':SEEDS,'profile':EXPECTED,'pairs':[],'triangle_lift':[]}
    need(all(f_at(EXPECTED,r)==r**6+12*r**4+16*r**2 for r in range(7)), 'expanded F polynomial mismatch')
    for mask in SEEDS:
        a=matrix(6,mask);need(strong(a),'seed not strong')
        need(profile_dp(a)==EXPECTED,'wrong subset profile')
        need(profile_permutations(a)==EXPECTED,'independent profile mismatch')
    for x,y in combinations(SEEDS,2):
        distances=[((x^relabel(6,y,p)).bit_count(),p) for p in permutations(range(6))]
        d,p=min(distances)
        need(d==1,'not distinct adjacent classes')
        other=relabel(6,y,p);flip=(x^other).bit_length()-1
        result['pairs'].append({'from':x,'to':y,'minimum_distance':d,'old_to_new':p,'relabelled_target':other,'flipped_pair':list(combinations(range(6),2))[flip]})
    current=SEEDS[0]
    for target in (SEEDS[1],SEEDS[2],SEEDS[0]):
        p=next(p for p in permutations(range(6)) if (current^relabel(6,target,p)).bit_count()==1)
        other=relabel(6,target,p)
        result['triangle_lift'].append({'from_labelled_mask':current,'to_labelled_mask':other,'target_class_seed':target,'old_to_new':p})
        current=other
    need(sign(p)==-1,'triangle holonomy should be odd')
    result['closing_permutation_sign']=sign(p)
    return result

def minimality_checks():
    result=[]
    for n in range(1,6):
        m=comb(n,2);which=[-1]*(1<<m);reps=[];perms=list(permutations(range(n)))
        for mask in range(1<<m):
            if which[mask]>=0:continue
            i=len(reps);reps.append(mask)
            for p in perms:which[relabel(n,mask,p)]=i
        adj=[set() for _ in reps]
        for i,x in enumerate(reps):
            for b in range(m):
                j=which[x^(1<<b)]
                if j!=i:adj[i].add(j)
        profiles=[profile_dp(matrix(n,x)) for x in reps]
        equal_edges=[(i,j) for i in range(len(reps)) for j in adj[i] if i<j and profiles[i]==profiles[j]]
        triangles=[(i,j,k) for i,j,k in combinations(range(len(reps)),3) if j in adj[i] and k in adj[i] and k in adj[j] and profiles[i]==profiles[j]==profiles[k]]
        need(not triangles,'smaller profile triangle found')
        result.append({'n':n,'labelled_tournaments':1<<m,'classes':len(reps),'profile_equal_edges':[(reps[i],reps[j]) for i,j in equal_edges],'profile_equal_triangles':0})
    need([r['classes'] for r in result]==[1,1,2,4,12],'class census mismatch')
    return result

def ordinal(a,b):
    n=len(a);m=len(b);c=[[0]*(n+m) for _ in range(n+m)]
    for i in range(n):
        for j in range(n):c[i][j]=a[i][j]
        for j in range(m):c[i][n+j]=1
    for i in range(m):
        for j in range(m):c[n+i][n+j]=b[i][j]
    return c

def f_at(p,r):
    fall=1;total=p[0]
    for j in range(1,len(p)):
        fall*=r-j+1;total+=p[j]*fall
    return total

def ordinal_checks():
    profiles=[]
    for x,y in product(SEEDS,repeat=2):
        p=profile_dp(ordinal(matrix(6,x),matrix(6,y)))
        need(p[1]==29**2,'wrong ordinal Hamilton count')
        need(all(f_at(p,r)==f_at(EXPECTED,r)**2 for r in range(13)),'ordinal interface product fails')
        profiles.append(p)
    need(len(set(profiles))==1,'nine ordinal sums not profile equivalent')
    return {'states':9,'order':12,'profile':profiles[0]}

def budget_checks():
    result=[]
    for k in range(1,9):
        states=[w for w in product(range(3),repeat=k) if all(not(w[i] and w[i+1]) for i in range(k-1))]
        ids={w:i for i,w in enumerate(states)};adj=[set() for _ in states]
        for i,w in enumerate(states):
            for pos in range(k):
                for c in range(3):
                    if c==w[pos]:continue
                    v=w[:pos]+(c,)+w[pos+1:]
                    if v in ids:adj[i].add(ids[v])
        vertices=sum(comb(k-j+1,j)*2**j for j in range((k+1)//2+1))
        edges=3*sum(j*comb(k-j+1,j)*2**(j-1) for j in range(1,(k+1)//2+1))
        need(len(states)==vertices==(2**(k+2)-(-1)**k)//3,'budget vertex count')
        need(sum(map(len,adj))==2*edges,'budget edge count')
        for start,w in enumerate(states):
            dist={start:0};todo=deque([start])
            while todo:
                x=todo.popleft()
                for y in adj[x]:
                    if y not in dist:dist[y]=dist[x]+1;todo.append(y)
            need(len(dist)==len(states),'budget graph disconnected')
            need(all(dist[i]==sum(a!=b for a,b in zip(w,v)) for i,v in enumerate(states)),'budget graph not Hamming-isometric')
        result.append({'slots':k,'states':vertices,'edges':edges,'all_pairs_Hamming_isometry':True})
    # Triple interaction which cannot be represented by pairwise conflicts.
    weights=(2,2,2);capacity=5
    faces=[s for s in range(8) if sum(weights[i] for i in range(3) if s>>i&1)<=capacity]
    need(all((1<<i)|(1<<j) in faces for i,j in combinations(range(3),2)) and 7 not in faces,'higher budget interaction control')
    # All nonempty downward-closed families on a three-element ground set.
    tested=0
    for bits in range(1,256):
        family={s for s in range(8) if bits>>s&1}
        if any((s^(1<<i)) not in family for s in family for i in range(3) if s>>i&1):continue
        tested+=1
        states=[w for w in product(range(3),repeat=3) if sum((w[i]!=0)<<i for i in range(3)) in family]
        ids={w:i for i,w in enumerate(states)};adj=[set() for _ in states]
        for i,w in enumerate(states):
            for pos in range(3):
                for c in range(3):
                    v=w[:pos]+(c,)+w[pos+1:]
                    if v!=w and v in ids:adj[i].add(ids[v])
        need(len(states)==sum(2**s.bit_count() for s in family),'general family vertex formula')
        need(sum(map(len,adj))==6*sum(s.bit_count()*2**(s.bit_count()-1) for s in family if s),'general family edge formula')
        for start,w in enumerate(states):
            dist={start:0};todo=deque([start])
            while todo:
                x=todo.popleft()
                for y in adj[x]:
                    if y not in dist:dist[y]=dist[x]+1;todo.append(y)
            need(len(dist)==len(states) and all(dist[i]==sum(a!=b for a,b in zip(w,v)) for i,v in enumerate(states)),'general family Hamming isometry')
    need(tested==19,'three-variable downward family census')
    return {'all_nonempty_three_variable_downward_families_checked':tested,'path_exclusion':result,'three_two_fifths_pieces':{'faces':faces,'states':sum(2**s.bit_count() for s in faces),'edges':3*sum(s.bit_count()*2**(s.bit_count()-1) for s in faces if s)}}

def ca_checks():
    result=[]
    for s in range(5):
        k=3**s;x=(1,)+(0,)*(k-1);initial=x;seen=[]
        while x not in seen:
            seen.append(x);x=tuple((x[i]+x[(i+1)%k])%3 for i in range(k))
        need(x==initial and len(seen)==2*k,'CA delta clock order')
        need(seen[k]==tuple(-a%3 for a in initial),'CA half-turn negation')
        result.append({'slots':k,'delta_period':len(seen),'unwrapped_tournament_order':6*k,'constant_H':str(29**k)})
    # Exhaustive matrix identity on all 3^9 states, not needed for proof.
    k=9
    for x in product(range(3),repeat=k):
        y=x
        for _ in range(k):y=tuple((y[i]+y[(i+1)%k])%3 for i in range(k))
        need(y==tuple(-a%3 for a in x),'CA all-state half-turn')
    return result

def carrier_checks():
    # A singleton has weight 1; a larger terminal block of size j has C(j,2).
    a=[1]
    for n in range(1,6):
        a.append(sum(comb(n-1,j-1)*(1 if j==1 else comb(j,2))*a[n-j] for j in range(1,n+1)))
    need(a[3:]==[7,28,126],'terminal-pair carrier counts')
    return a

def main():
    out={'seed':seed_checks(),'minimality':minimality_checks(),'ordinal':ordinal_checks(),'budget':budget_checks(),'cellular_automaton':ca_checks(),'unweighted_carriers':carrier_checks()}
    (HERE/'certificate.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
    print('ALL CHECKS PASSED')

if __name__=='__main__':main()
