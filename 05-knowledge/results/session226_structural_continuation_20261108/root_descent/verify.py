#!/usr/bin/env python3
import itertools, json, math

def require(x, msg):
    if not x: raise RuntimeError(msg)

def profile(rows, root):
    n=len(rows); D=E=U=0
    def path(p): return all(rows[a]>>b&1 for a,b in zip(p,p[1:]))
    for p in itertools.permutations(range(n)):
        good=path(p)
        D+=good; E+=good and p[-1]==root
        j=p.index(root)
        U+=path(p[:j]) and path(p[j:])
    d=sum(path(p) for p in itertools.permutations([i for i in range(n) if i!=root]))
    return D,E,U,d

def move(s, k):
    D,E,U,d=s
    return {'O':(D,0,2*U,d),'Z':(D-E+U,d,2*U-d,d),'J':(D-E+U,D-E,D+U,D)}[k]

def cone(s):
    D,E,U,d=s
    return 0<=E<=d and 0<=D-E<=U-d and 1<=D<=U

profiles=0; unique=set(); seeds=[]
for n in range(1,6):
    pairs=list(itertools.combinations(range(n),2))
    for mask in range(1<<len(pairs)):
        rows=[0]*n
        for j,(a,b) in enumerate(pairs):
            if mask>>j&1: rows[a]|=1<<b
            else: rows[b]|=1<<a
        for root in range(n):
            s=profile(rows,root); require(cone(s),('physical cone',rows,root,s))
            profiles+=1; unique.add(s)
            if n<=3 and s not in seeds: seeds.append(s)
require((1,1,1,1) in seeds, 'singleton/sink hostile absent')
require((1,0,2,1) in seeds, 'original seed absent')
word_checks=0; C=15
for seed in seeds:
    stack=[(seed,'')]
    while stack:
        s,w=stack.pop(); D,E,U,d=s
        require(cone(s),('descendant cone',seed,w,s))
        h,g=D+2*U,2*U
        require(3*g>=2*h and g%2==0, 'endpoint constraint')
        end=move(move(s,'O'),'J')
        require(end==(h,D,h,D), 'endpoint formula')
        if D<=C:
            v=U
            while v%2==0: v//=2
            sv=seed[2]
            while sv%2==0: sv//=2
            require(v<=2*C or (D==seed[0] and v==sv),('dyadic envelope',seed,w,s))
            if any(k!='O' for k in w):
                last=max(i for i,k in enumerate(w) if k!='O')
                p=seed
                for k in w[:last]: p=move(p,k)
                q=move(p,w[last])
                require(p[2]<=C and q[2]<=2*C,('last non-O',seed,w,p,q))
        word_checks+=1
        if len(w)<7:
            for k in 'OZJ': stack.append((move(s,k),w+k))
congruence_checks=0
for b in range(2,13):
    for a in range(1,b):
        if math.gcd(a,b)!=1: continue
        for c in range(-20,21):
            for N in range(-30,61):
                den=a+b; num=a*N-c
                integer_even=(num%den==0 and (num//den)%2==0)
                require(integer_even==((a*N-c)%(2*den)==0),'residue equivalence')
                if abs(c)<a: require((a-c)%(2*den)!=0,'small intercept')
                congruence_checks+=1
require(all(any((2*r-c)%10==0 for c in (0,-2,-4,-6,-8)) for r in (1,3,5,7,9)), 'full-cover hostile')
print(json.dumps({'status':'PASS','labelled_rooted_tournament_profiles':profiles,'distinct_profiles':len(unique),'physical_seed_profiles':seeds,'word_checks':word_checks,'congruence_checks':congruence_checks,'max_tournament_order':5,'max_word_length':7,'boundary_cap':C,'hostile_residue_cover':[0,-2,-4,-6,-8]},sort_keys=True,indent=2))
