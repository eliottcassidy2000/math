#!/usr/bin/env python3
"""Exact consequence checks. No asserts; runs identically under python -O."""
from itertools import product, permutations, combinations
from fractions import Fraction as F
from math import isqrt, gcd, comb
from collections import Counter
import json

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def pairs(n):return [(i,j) for i in range(1,n+1) for j in range(1,i-1)]
def tile(n,mask):return dict(zip(pairs(n),[(mask>>i)&1 for i in range(len(pairs(n)))]))
def corner(t,n,h):return tuple(t[n-b,1+a] for a in range(h+1) for b in range(h+1-a))
def cd(t,n,h):return {(a,b):t[n-b,1+a] for a in range(h+1) for b in range(h+1-a)}
def delete(t,n,low):
    return {(i-(1 if low else 0),j-(1 if low else 0)):v for (i,j),v in t.items() if (j>1 if low else i<n)}
def attach(t,n,row,low):
    if low:
        out={(i+1,j+1):v for (i,j),v in t.items()}
        out.update({(i+1,1):row[i-2] for i in range(2,n+1)})
    else:
        out=dict(t);out.update({(n+1,j):row[j-1] for j in range(1,n)})
    return out

def legal(t,n,w):
    def edge(a,b):
        i,j=max(a,b),min(a,b);bit=0 if i-j==1 else t[i,j]
        return (bit==0)==(a>b)
    return all(edge(a,b) for a,b in zip(w,w[1:]))

def corners():
    counts=[];updates=0;attachment_updates=0
    for n in range(3,8):
        rows=[set() for _ in range(n-2)]
        for mask in range(1<<len(pairs(n))):
            t=tile(n,mask)
            for h in range(n-2):
                c=cd(t,n,h);rows[h].add(corner(t,n,h))
                # Independent decoding by actual endpoint deletion sequences.
                for a,b in c:
                    u=t;size=n
                    for _ in range(a):u=delete(u,size,True);size-=1
                    for _ in range(b):u=delete(u,size,False);size-=1
                    need(u[size,1]==c[a,b],'deletion query decoder')
                refl={(n+1-j,n+1-i):v for (i,j),v in t.items()}
                need(cd(refl,n,h)=={(a,b):c[b,a] for a,b in c},'reflection')
                need(cd({k:1-v for k,v in t.items()},n,h)=={k:1-v for k,v in c.items()},'complement')
                if h:
                    for low in [False,True]:
                        target=cd(delete(t,n,low),n-1,h-1)
                        pred={(a,b):c[a+int(low),b+int(not low)] for a,b in target}
                        need(target==pred,'deletion compiler');updates+=1
                    if n<=6:
                        for row in product([0,1],repeat=n-1):
                            for low in [False,True]:
                                target=cd(attach(t,n,row,low),n+1,h-1)
                                pred={}
                                for a,b in target:
                                    if low:pred[a,b]=row[n-b-2] if a==0 else c[a-1,b]
                                    else:pred[a,b]=row[a] if b==0 else c[a,b-1]
                                need(target==pred,'attachment compiler');attachment_updates+=1
        for h,seen in enumerate(rows):
            count=1<<((h+1)*(h+2)//2)
            need(len(seen)==count,'minimal corner carrier cardinality')
            counts.append([n,h,count])
    a=tile(5,0);a[5,2]=1;b=dict(a);b[3,1]=1
    need(corner(a,5,1)==corner(b,5,1) and legal(a,5,[3,2,5,4,1]) and legal(b,5,[3,2,5,4,1]),'common legal same-clock remark')
    need(a[3,1]!=b[3,1],'remark separates endpoint')
    languages=[]
    for n in range(3,7):
        words=list(permutations(range(1,n+1)));profiles=set()
        for mask in range(1<<len(pairs(n))):
            t=tile(n,mask);profiles.add(tuple(i for i,w in enumerate(words) if legal(t,n,w)))
        need(len(profiles)==1<<len(pairs(n)),'labelled HP legality must separate')
        languages.append([n,len(profiles)])
    return {'corner_counts':counts,'deletion_updates':updates,'attachment_updates':attachment_updates,'legality_profiles':languages,'remark_hostile':{'path':[3,2,5,4,1],'initial_corner':list(corner(a,5,1)),'after_remark':[0,1]}}

def cc(n):return n//2 if n%2==0 else 3*n+1
def branches(n):
    return [(k,((1<<k)*n-1)//3) for k in range(1,7) if ((1<<k)*n-1)%3==0 and (((1<<k)*n-1)//3)%3]
def collatz():
    residue={}
    for r in [1,2,4,5,7,8]:
        n=r+18 if r%2 else r+9
        z=branches(n);need(len(z)==2,'two unit inverse branches');residue[str(r)]=[k for k,_ in z]
    level={85};seen={85};counts=[1];edge_checks=0
    for depth in range(1,17):
        nxt=set()
        for n in level:
            z=branches(n);need(len(z)==2,'inverse arity')
            for k,y in z:
                need(y>1 and y%2==1 and y%3!=0,'odd unit child')
                cur=y
                for step in range(k+1):
                    cur=cc(cur)
                    if step<k:need(cur%2==0,'no intermediate odd node')
                need(cur==n,'native clock')
                need(y not in nxt and y not in seen,'source labels do not collide')
                nxt.add(y);edge_checks+=1
        level=nxt;seen|=level;need(len(level)==1<<depth,'binary level cardinality');counts.append(len(level))
    # An exact finite witness to the L=1 polynomial-growth bound failing.
    r=1
    while (1<<r)<=2*(14*r+1)**2:r+=1
    need(cc(5)==cc(32)==16,'branching hostile')
    need(1+(F(1,4)**2+F(1,2)**2)/(2*F(1,2))==F(21,16),'hyperbolic constant-distance identity')
    return {'unit_residue_exponents':residue,'depth_counts':counts,'native_macro_edges_checked':edge_checks,'unit_dilation_contradiction_depth':r,'largest_depth16_integer_bits':max(level).bit_length()}

def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def genus_prime(p):
    need(prime(p),'prime-level input')
    e2=1 if p==2 else 2*int(p%4==1)
    e3=1 if p==3 else 2*int(p%3==1)
    g=F(p+1,12)-F(e2,4)-F(e3,3)
    need(g.denominator==1 and g>=0,'integral nonnegative genus');return int(g)
def power2(n):return n>0 and n&(n-1)==0

def genus():
    O=[2,3,5,7,11,13,17,19,23,29,31,41,47,59,71]
    need([p for p in O if power2(p-1)]==[2,3,5,17],'power-two p-1 includes2; Fermat definition excludes2')
    need([p for p in O if p>2 and power2(p-1)]==[3,5,17],'Ogg Fermat intersection')
    zero=[p for p in range(2,200) if prime(p) and genus_prime(p)==0]
    need(zero==[2,3,5,7,13],'genus zero prime levels')
    rows=[]
    for m in range(2,5):
        p=(1<<(1<<m))+1;g=genus_prime(p);k=(1<<(m-1))-1
        need(g==((1<<(2*k))-1)//3 and 3*g+1==1<<(2*k),'Fermat genus Collatz identity')
        cur=g
        for _ in range(2*k+1):cur=cc(cur)
        need(cur==1,'certified standard return')
        rows.append({'m':m,'p':p,'genus':g,'trunk_index':k,'power':3*g+1})
    need(5461==43*127 and prime(43) and prime(127),'nonconstructible factor example')
    F5=(1<<32)+1;need(F5==641*6700417 and prime(641) and prime(6700417),'composite Fermat control')
    actual=F(642*6700418,12)-2;formal=F(F5-5,12)
    need(actual==358472361 and formal==357913941 and not power2(3*int(actual)+1),'composite genus hostile')
    channels=[]
    for k in range(1,9):
        g=((1<<(2*k))-1)//3;ps=[(1<<(2*k+2))+d for d in [9,1,3,-5]]
        for p in ps:
            if prime(p):need(genus_prime(p)==g,'four-channel genus')
        channels.append({'k':k,'genus':g,'prime_levels':[p for p in ps if prime(p)]})
    return {'genus_zero_prime_levels_through199':zero,'Fermat_genus_rows':rows,'Ogg_Fermat_intersection':[3,5,17],'channels':channels,'composite_F5':{'actual_genus':int(actual),'formal_trunk':int(formal),'actual_3g_plus1':3*int(actual)+1}}

def graph(n,mask):
    E=list(combinations(range(n),2));return {e for i,e in enumerate(E) if mask>>i&1}
def ep(n,E):
    pair=Counter();start=Counter()
    for w in permutations(range(n)):
        if all(tuple(sorted((a,b))) in E for a,b in zip(w,w[1:])):pair[w[0],w[-1]]+=1;start[w[0]]+=1
    return [start[i] for i in range(n)],dict(pair)
def ports():
    controls=0
    for n in range(2,5):
        for mask in range(1<<(n*(n-1)//2)):
            E=graph(n,mask);es,ps=ep(n,E)
            for p in range(n):
                aug=E|{(p,n)};a,_=ep(n+1,aug)
                need(sum(a)==2*es[p],'single-leaf exact count');controls+=1
                for q in range(p+1,n):
                    aug=E|{(p,n),(q,n+1)};a,_=ep(n+2,aug)
                    need(sum(a)==2*ps.get((p,q),0),'two-leaf joint count');controls+=1
    G={(0,1),(1,2),(2,3),(0,3)};H={(0,1),(1,3),(2,3),(0,2)}
    ge,gp=ep(4,G);he,hp=ep(4,H)
    need(ge==he==[2,2,2,2] and gp.get((0,3),0)==1 and hp.get((0,3),0)==0,'endpoint marginal hostile')
    Q=lambda n:{(a,b) for a,b in combinations(range(1,n+1),2) if isqrt(a+b)**2==a+b}
    path=[8,1,15,10,6,3,13,12,4,5,11,14,2,7,9]
    need(set(path)==set(range(1,16)) and all(tuple(sorted(e)) in Q(15) for e in zip(path,path[1:])),'square-sum witness')
    for n,target in [(16,9),(17,8),(18,7)]:need([a for a in range(1,n) if (a,n) in Q(n)]==[target],'native leaf port')
    leaves=[v for v in range(1,19) if sum(v in e for e in Q(18))==1]
    need(leaves==[16,17,18],'three-leaf obstruction')
    # Exact independent subset DP for Q15 endpoint multiplicities.
    adj=[[] for _ in range(15)]
    for a,b in Q(15):adj[a-1].append(b-1);adj[b-1].append(a-1)
    dp=[{} for _ in range(1<<15)]
    for i in range(15):dp[1<<i][i]=1
    for mask,row in enumerate(dp):
        for v,count in row.items():
            for w in adj[v]:
                if not mask>>w&1:
                    row2=dp[mask|1<<w];row2[w]=row2.get(w,0)+count
    endpoint={str(v+1):c for v,c in dp[-1].items()}
    need(endpoint=={'8':1,'9':1},'Q15 exact endpoint counts')
    return {'leaf_attachment_controls':controls,'two_port_hostile':{'marginals':ge,'joint_03':[1,0],'two_leaf_H':[2,0]},'Q15_endpoint_counts':endpoint,'Q18_leaves':leaves}

def fermat_neighbours():
    rows=[]
    for n in [3,5,7,9,11,13,15,17]:
        units=[a for a in range(1,n) if gcd(a,n)==1]
        reps=set();stabs=Counter()
        for bits in product([0,1],repeat=(n-1)//2):
            S=frozenset(i+1 if b==0 else n-i-1 for i,b in enumerate(bits))
            images=[tuple(sorted(a*x%n for x in S)) for a in units]
            reps.add(min(images));stabs[images.count(tuple(sorted(S)))]+=1
        rows.append({'n':n,'connection_sets':1<<((n-1)//2),'unit_group_order':len(units),'multiplier_orbits':len(reps),'stabilizer_histogram':dict(stabs)})
    need(rows[-2]['multiplier_orbits']==rows[-1]['multiplier_orbits']==16,'15/17 orbit equality')
    need(all(list(r['stabilizer_histogram'])==[1] for r in rows if r['n'] in [3,5,15,17]),'free constructible unit actions')
    controls=[]
    for a,b in [(3,5),(5,3)]:
        n=a*b;S={x for x in range(1,n) if x%a in range(1,(a+1)//2) or (x%a==0 and x//a in range(1,(b+1)//2))}
        group=set()
        for shift in range(a):
            for rots in product(range(b),repeat=a):
                perm=tuple((v%a+shift)%a+a*((v//a+rots[v%a])%b) for v in range(n))
                need(len(set(perm))==n,'block bijection')
                need(all(((y-x)%n in S)==((perm[y]-perm[x])%n in S) for x in range(n) for y in range(n) if x!=y),'block automorphism')
                group.add(perm)
        need(len(group)==a*b**a,'explicit block subgroup size')
        need([u for u in range(1,n) if gcd(u,n)==1 and {u*x%n for x in S}==S]==[1],'no multiplier symmetry despite block symmetry')
        controls.append({'block_count':a,'block_size':b,'S':sorted(S),'subgroup_order':len(group)})
    ladder=[]
    for m in range(1,5):
        top=(1<<(1<<m))+1;bottom=top-2
        exponent=(1<<((1<<m)-1))-(1<<m)
        need((bottom-1)//2-((1<<m)-1)==exponent==(top-1)//2-(1<<m),'orbit exponent identity')
        ladder.append([bottom,top,exponent])
    def fermat_product(s):
        out=1
        for i in range(s.bit_length()):
            if s>>i&1:out*=1+(1<<(1<<i))
        return out
    carry_zeros=[]
    for s in range(129):
        pascal_row=sum((comb(s,t)%2)<<t for t in range(s+1))
        need(fermat_product(s)==pascal_row,'independent Pascal parity row')
    for s in range(4096):
        m=0
        while s>>m&1:m+=1
        A=fermat_product(s & ~((1<<(m+1))-1))
        left,right=fermat_product(s),fermat_product(s+1)
        need(left==A*((1<<(1<<m))-1) and right==A*((1<<(1<<m))+1),'Fermat binary carry product')
        delta=(right-left)//2-1
        need(delta==A-1 and delta>=0,'monotone orbit exponent')
        need((delta==0)==(s==(1<<m)-1),'exact carry equality boundary')
        if delta==0:carry_zeros.append(s)
    collisions={}
    for s in range(1,32):
        n=fermat_product(s);e=(n-1)//2-s
        collisions.setdefault(e,[]).append(n)
    pairs_known=sorted(v for v in collisions.values() if len(v)>1)
    need(pairs_known==[[3,5],[15,17],[255,257],[65535,65537]],'all known-prime subsets collision census')
    # Missing F5 blocks an initial segment, not the auxiliary carry identity.
    need(fermat_product(32)==641*6700417,'first inadmissible carry partner')
    return {'small_orbit_census':rows,'block_symmetry_controls':controls,'ladder_log2_counts':ladder,'binary_carry_controls':4096,'Pascal_parity_rows':129,'zero_increments_including_auxiliary_zero':carry_zeros,'all_known_Fermat_subset_collisions':pairs_known}

def monodromy():
    def cycles(perm):
        seen=set();counts=Counter()
        for x in range(len(perm)):
            if x in seen:continue
            y=x;length=0
            while y not in seen:seen.add(y);length+=1;y=perm[y]
            need(y==x,'permutation closes on source');counts[length]+=1
        return counts
    out=[]
    for p in [5,7,11,13,17,19,59,67,73,251,257,65537]:
        need(prime(p),'prime permutation field')
        S=[p if x==0 else -pow(x,-1,p)%p for x in range(p)]+[0]
        R=[p if x==0 else (1-pow(x,-1,p))%p for x in range(p)]+[1]
        T=[R[S[x]] for x in range(p+1)]
        need(T==[(x+1)%p for x in range(p)]+[p],'actual modular generator product')
        cs,cr,ct=cycles(S),cycles(R),cycles(T)
        chi=sum(cs.values())+sum(cr.values())+sum(ct.values())-(p+1)
        need(chi%2==0 and 1-chi//2==genus_prime(p),'independent permutation genus')
        out.append({'p':p,'S_cycles':dict(cs),'R_cycles':dict(cr),'cusp_cycles':dict(ct),'genus':1-chi//2})
    return out

def main():
    print(json.dumps({'status':'passed','corner':corners(),'Collatz':collatz(),'genus':genus(),'ports':ports(),'Fermat_neighbours':fermat_neighbours(),'monodromy':monodromy(),'scope':'Exact finite controls supplement the written all-size proofs; no independent reviewer, all-order enumeration or universal Collatz claim.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
