#!/usr/bin/env python3
"""Session73 independent exact audit; all checks survive python -O."""
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import factorial, gcd
import json

def check(ok, why):
    if not ok:
        raise RuntimeError(why)

def w(L,K):
    return Q(2*factorial(K)*factorial(L+1),factorial(L+K+2))

def odd_step(n):
    t=3*n+1; a=0
    while t%2==0:
        a+=1; t//=2
    return t,a

def root_controls():
    rows=[]
    for j in range(1,81):
        n=(4**(j+1)-1)//3
        nxt,a=odd_step(n)
        check(nxt==1 and a==2*j+2,'root ray actual word')
        W=w(0,j)
        check(W==Q(2,(j+1)*(j+2)),'root ray first-hit counters')
        check((n%3==0)==(j%3==2),'native lambda support')
        b=Q(3*n,3*n+1); d=1-b*b
        check(1/d>=Q(3*n,2),'passive linear cost lower bound')
        f=1-b*W
        check(f>=Q(2,3),'root score forcing lower bound')
        killed=-b*W
        if n%3==0:
            bill=W*f*f/d
            killed_bill=W*killed*killed/d
            check(bill>=Q(2*n,3)*W,'live first-step paid lower bound')
            check(killed_bill>=Q(3*n,2)*b*b*W**3,'killed first-step paid lower bound')
            rows.append((j,n,bill,killed_bill))
    # On this subsequence j advances by3: n grows at least64-fold;
    # (j+1)(j+2) grows by less than4-fold for j>=2. Hence n*W^3
    # grows by more than one and in fact eventually at least8-fold.
    for x,y in zip(rows,rows[1:]):
        check(y[1]>64*x[1],'exponential sibling height')
        check(y[2]>x[2],'nonkilled bills increase')
        check(y[3]>x[3],'killed bills increase')
    sums=0
    for N in range(1,101):
        sums+=w(0,N)
        check(sums==1-Q(2,N+2),'exceptional root incoming row')
    tails=0
    for L,K,J,N in product(range(5),range(1,6),range(5),range(1,9)):
        total=sum((w(L+1,K+j) for j in range(J,J+N)),Q())
        check(total==w(L,K+J)-w(L,K+J+N),'complete incoming row finite tail')
        tails+=1
    identities=0
    for n,e,f in product(range(1,256,2),[Q(-2),Q(0),Q(1,3),Q(2)],[Q(-1),Q(0),Q(2,5)]):
        b=Q(3*n,3*n+1); d=1-b*b
        out=b*e+f; r=d*e-b*f
        check(out*out+r*r/d==e*e+f*f/d,'passive identity')
        identities+=1
    # Hostile: un-killed root self-loop would add W(1)=1 to the root
    # incoming row and violate the contraction at f=1_{1}.
    check(1+sum((w(0,j) for j in range(1,10)),Q())>1,'root-loop hostile')
    return dict(ray_terms=80,lambda_subsequence_terms=len(rows),first_terms=[dict(j=j,n=n,live_bill=str(a),killed_bill=str(b)) for j,n,a,b in rows[:4]],incoming_tail_checks=tails,passive_identities=identities,hostile_root_loop='rejected')

def rooted_data(n):
    word=[]; x=n
    while x!=1:
        check(len(word)<10000,'bounded test failed to reach root')
        x,a=odd_step(x); word.append(a)
    return word,(w(len(word)-1,sum((a-1)//2 for a in word)) if word else Q(1))

def signed_and_phase_controls():
    samples=[]
    for k in range(32):
        n=(10*64**k-1)//3
        word,W=rooted_data(n)
        check(word==[6*k+1,4] and n%3==0,'signed native ray guards')
        check(W==Q(4,(3*k+2)*(3*k+3)*(3*k+4)),'signed native lambda')
        sigma=-1;b=Q(3*n,3*n+1);d=1-b*b
        e=Q(-1);ep=Q(1)
        F=ep-b*e
        check(F*F/d==6*n+1,'unsigned-gauge signed price')
        Fp=ep-sigma*b*e
        check(Fp*Fp/d==Q(1,6*n+1),'phase-aware signed price')
        if k<4:samples.append(dict(k=k,n=n,lambda_mass=str(W),unsigned_bill=str(W*F*F/d),phase_bill=str(W*Fp*Fp/d)))
    controls=0
    occupancy={}; pathbill=Q(); literalclock=shortcutclock=Q()
    for src in range(3,768,6):
        word,W=rooted_data(src)
        n=src
        for a in word:
            n1,actual_a=odd_step(n)
            rem,_=rooted_data(n);rem1,_=rooted_data(n1)
            e=Q((-1)**sum(rem));ep=Q((-1)**sum(rem1));sigma=(-1)**a
            check(actual_a==a and e==sigma*ep,'signed shortcut recurrence')
            b=Q(3*n,3*n+1);d=1-b*b;c=sigma*b;F=ep-c*e;r=d*e-c*F
            check(F*F/d==Q(1,6*n+1),'phase gain exact local price')
            check(ep*ep+r*r/d==e*e+F*F/d,'signed passive identity')
            occupancy[n]=occupancy.get(n,Q())+W
            pathbill+=W*F*F/d
            literalclock+=W*(a+1);shortcutclock+=W*a
            n=n1;controls+=1
    statebill=sum((mass/Q(6*n+1) for n,mass in occupancy.items()),Q())
    check(statebill==pathbill,'source-weight occupation payment')
    for n,mass in occupancy.items():
        _,W=rooted_data(n)
        check(mass<=W,'finite source occupation below actual W')
    check(pathbill<=Q(13,57),'finite positive-source paid bound')
    check(literalclock>shortcutclock,'literal versus shortcut clocks')
    # Finite authenticated horizon A_T uses zero outside. At its boundary
    # a known target has sign +/-1 but the source is still zero, so every
    # gain choice gives forcing magnitude1 and passive price1/d.
    target=1;frontiers=[]
    for T in range(6):
        tw,_=rooted_data(target)
        check(len(tw)==T and target%3,'unit target at each odd-step depth')
        base=None
        for a in [1,2]:
            if (2**a*target-1)%3==0:
                base=(2**a*target-1)//3;break
        check(base is not None,'inverse base')
        tested=0
        for j in range(1,25):
            n=(4**j*(3*base+1)-1)//3
            if n%3:continue
            word,W=rooted_data(n)
            check(word[1:]==tw and len(word)==T+1,'horizon frontier actual word')
            b=Q(3*n,3*n+1);d=1-b*b
            check(W/d>=W*Q(3*n,2),'phase truncation frontier lower bound')
            tested+=1
        check(tested>0,'frontier lambda siblings')
        frontiers.append(dict(T=T,target=target,lambda_siblings=tested))
        # Pick first strictly nonroot inverse unit, retaining a concrete certificate.
        for a in range(1,20):
            num=2**a*target-1
            if num%3==0 and num//3>1 and (num//3)%3:
                target=num//3;break
    n=21;word,_=rooted_data(n)
    check(sum(word)==6 and sum(word)+len(word)==7,'parity clock hostile')
    return dict(ray_cases=32,ray_first_terms=samples,local_signed_identities=controls,original_lambda_sources=128,finite_source_occupation_bill=str(pathbill),horizon_frontiers=frontiers,clock_parity_hostile=dict(n=21,shortcut=6,literal=7))

def psd(a):
    a=[list(r) for r in a]; N=len(a)
    for k in range(N):
        check(a[k][k]>=0,'negative PSD pivot')
        if not a[k][k]:
            check(all(a[k][i]==0 for i in range(k+1,N)),'zero PSD pivot with nonzero row')
            continue
        for i in range(k+1,N):
            for j in range(i,N):
                a[i][j]-=a[i][k]*a[k][j]/a[k][k]
                a[j][i]=a[i][j]

def hardy_controls():
    tests=0
    for L in range(1,5):
        for u,v in product(product([Q(1),Q(2)],repeat=L),repeat=2):
            for gains in product([Q(-1,2),Q(2)],repeat=L-1):
                p=[Q(1)]
                for a in gains:p.append(p[-1]*a)
                cuts=[sum((u[j]*p[j]**2 for j in range(k,L)),Q())*sum((1/(v[i]*p[i]**2) for i in range(k+1)),Q()) for k in range(L)]
                B=max(cuts)
                R=[[Q() if i>j else p[j]/p[i] for i in range(L)] for j in range(L)]
                a=[[4*B*v[i]*(i==j)-sum((u[k]*R[k][i]*R[k][j] for k in range(L)),Q()) for j in range(L)] for i in range(L)]
                psd(a)
                for k in range(L):
                    f=[1/(v[i]*p[i]) if i<=k else Q() for i in range(L)]
                    out=[sum((R[j][i]*f[i] for i in range(L)),Q()) for j in range(L)]
                    ratio=sum((u[j]*out[j]**2 for j in range(L)),Q())/sum((v[i]*f[i]**2 for i in range(L)),Q())
                    check(ratio>=cuts[k],'Hardy necessary bound test')
                tests+=1
    # Coarse atom mass1 copied to children1/3,2/3 isometric, child contrast
    # vector(2,-1) has weighted conditional mean0 and can carry arbitrary norm.
    child=[Q(1,3),Q(2,3)]
    check(sum(child)==1 and sum((a*x for a,x in zip(child,[2,-1])),Q())==0,'zero mean fine mode')
    check(Q(1,2)*(1**2+2**2)==Q(5,2),'false microclock hostile')
    return dict(hardy_cases=tests,refinement_hostiles=['two-step copied energy5/2','zero-mean fine mode unbounded'])

def clock_controls():
    words=cells=0; non_surjective=0
    for h in range(4):
        for word in product(range(1,5),repeat=h):
            S=sum(word); P=3**h; Qd=2**S
            B=0; pref=0
            for i,a in enumerate(word):
                B+=3**(h-1-i)*2**pref; pref+=a
            n0=((Qd-B)*pow(P,-1,2*Qd))%(2*Qd)
            z0=(n0-1)//2
            d=((P*n0+B)//Qd-1)//2
            for q in range(1,19):
                g=gcd(P,q); period=Qd*q//g
                image=set()
                for t in range(q//g):
                    z=z0+Qd*t; n=2*z+1; raw=n; steps=0
                    for a in word:
                        n1,k=odd_step(n)
                        check(k==a,'native branch cylinder')
                        raw=3*raw+1; steps+=1
                        for _ in range(a):raw//=2; steps+=1
                        check(raw==n1,'raw clock replay')
                        n=n1
                    target=(n-1)//2
                    check(target==d+P*t,'affine endpoint index')
                    check(steps==h+S,'integer clock')
                    image.add(target%q); cells+=1
                check(len(image)==q//g,'image size')
                check(image=={r for r in range(q) if (r-d)%g==0},'image guard')
                check(P*(period//Qd)%q==0,'period sufficiency')
                for t in range(1,q//g):check(P*t%q!=0,'minimal source period')
                non_surjective+=g>1
            words+=1
    # An odd-modulus representative alone cannot supply parity-selected division.
    check(1%5==6%5 and (1//2)%5!=(6//2)%5,'division precision hostile')
    # Two positive lifts with the same k-bit residue can have different later bits.
    for k in range(1,10):
        check(0%(2**k)==2**k%(2**k),'finite tape collision')
        check(((2**k)>>k)&1,'delayed tape guard')
    return dict(words=words,endpoint_cells=cells,non_surjective_word_moduli=non_surjective,precision_hostiles='rejected')

def paley(p):
    residues={x*x%p for x in range(1,p)}
    return [[False if x==y else (y-x)%p in residues for y in range(p)] for x in range(p)]

def cyclic(T,a,b,c):
    return (T[a][b] and T[b][c] and T[c][a]) or (T[b][a] and T[c][b] and T[a][c])

def fas(T,order):
    return sum(T[order[j]][order[i]] for i in range(len(T)) for j in range(i+1,len(T)))

def components(n,edges):
    par=list(range(n))
    def rt(v):
        while par[v]!=v:v=par[v]
        return v
    for a,b in edges:par[rt(a)]=rt(b)
    return len({rt(i) for i in range(n)})

def isforest(n,edges):return len(edges)==n-components(n,edges)

def tournament_controls():
    result={}
    for p in [7,11]:
        T=paley(p); triangles=[e for e in combinations(range(p),3) if cyclic(T,*e)]
        incidences={e:sum(set(e)<=set(t) for t in triangles) for e in combinations(range(p),2)}
        check(set(incidences.values())=={(p+1)//4},'Paley triangle incidence')
        check(len(triangles)==p*(p*p-1)//24,'Paley triangle count')
        result[str(p)]=dict(triangles=len(triangles),triangles_per_edge=(p+1)//4,triangle_FAS_lower=(len(triangles)+(p+1)//4-1)//((p+1)//4))
    T=paley(7); fs=[fas(T,o) for o in permutations(range(7))]
    check(min(fs)==7,'Paley7 all-rank FAS')
    order=(0,3,1,6,4,2,5); rank={v:i for i,v in enumerate(order)}
    S=[[False if a==b else T[a][b] if abs(a-b)==1 else not T[a][b] for b in range(7)] for a in range(7)]
    D={e for e in combinations(range(7),2) if S[e[0]][e[1]]!=(rank[e[0]]<rank[e[1]])}
    check(D=={(0,4),(0,2),(2,3),(3,5),(1,5),(5,6)} and isforest(7,D),'Paley7 forest predecessor')
    checks=0
    # All fixed-rank forest frames and all Hamilton guards at orders4,5.
    for n in [4,5]:
        pairs=list(combinations(range(n),2))
        for mask in range(1<<len(pairs)):
            D={e for i,e in enumerate(pairs) if mask>>i&1}
            if not isforest(n,D):continue
            S=[[False if a==b else (a<b)!=((min(a,b),max(a,b)) in D) for b in range(n)] for a in range(n)]
            for P in permutations(range(n)):
                if not all(S[a][b] for a,b in zip(P,P[1:])):continue
                E={tuple(sorted(e)) for e in zip(P,P[1:])}; Dp=D^E
                T=[[False if a==b else S[a][b] if tuple(sorted((a,b))) in E else not S[a][b] for b in range(n)] for a in range(n)]
                order=tuple(range(n-1,-1,-1))
                if not D&E:
                    check(order==P[::-1],'zero-overlap reverse rank')
                    new=list(order)
                    for i in range(0,n-1,2):new[i],new[i+1]=new[i+1],new[i]
                    check(fas(T,new)==len(Dp)-n//2,'disjoint adjacent repair')
                    order=new
                check(fas(T,order)<=2*n-4,'forest one-flip upper bound')
                checks+=1
    # Ternary lexicographic C3 tower, independent subset-DP FAS minimum.
    for k in [1,2]:
        n=3**k
        def edge(a,b):
            if a==b:return False
            for j in range(k-1,-1,-1):
                x=(a//3**j)%3;y=(b//3**j)%3
                if x!=y:return (y-x)%3==1
        T=[[edge(a,b) for b in range(n)] for a in range(n)]
        dp=[0]+[10**9]*((1<<n)-1)
        for mask in range(1,1<<n):
            for v in range(n):
                if mask>>v&1:
                    rest=mask^(1<<v)
                    dp[mask]=min(dp[mask],dp[rest]+sum(T[v][u] for u in range(n) if rest>>u&1))
        check(dp[-1]==n*(n-1)//6,'C3 tower FAS')
        result['tower'+str(n)]=dict(FAS=dp[-1],beta_lower=max(0,dp[-1]-n+1))
    result['forest_one_flip_cases']=checks
    result['Paley11_one_flip_from_forest']='impossible:19>18'
    result['Paley7_orders']=len(fs)
    return result

def weighted_fas(T,cost):
    n=len(T); dp=[Q()]+[None]*((1<<n)-1)
    for mask in range(1,1<<n):
        vals=[]
        for v in range(n):
            if mask>>v&1:
                rest=mask^(1<<v)
                vals.append(dp[rest]+sum((cost[v]*cost[u] for u in range(n) if rest>>u&1 and T[v][u]),Q()))
        dp[mask]=min(vals)
    return dp[-1]

def substitution_controls():
    cases=0
    blocks=[[0,1],[2],[3,4]]
    for qbits in product([False,True],repeat=3):
        Qt=[[False]*3 for _ in range(3)]
        for (i,j),bit in zip(combinations(range(3),2),qbits):Qt[i][j]=bit;Qt[j][i]=not bit
        for internal in product([False,True],repeat=2):
            T=[[False]*5 for _ in range(5)]
            for i,j in combinations(range(3),2):
                for u in blocks[i]:
                    for v in blocks[j]:T[u][v]=Qt[i][j];T[v][u]=Qt[j][i]
            for (u,v),bit in zip([(0,1),(3,4)],internal):T[u][v]=bit;T[v][u]=not bit
            for costs in product([Q(1),Q(2)],repeat=5):
                actual=weighted_fas(T,costs)
                Acost=[sum(costs[v] for v in b) for b in blocks]
                child=sum((weighted_fas([[T[u][v] for v in b] for u in b],[costs[v] for v in b]) for b in blocks),Q())
                check(actual==child+weighted_fas(Qt,Acost),'weighted substitution exact law')
                cases+=1
    def ctower(k):
        if not k:return [[False]]
        C=ctower(k-1);m=len(C)
        return [[False if i==j else C[i%m][j%m] if i//m==j//m else (j//m-i//m)%3==1 for j in range(3*m)] for i in range(3*m)]
    connected=[]
    for k in range(2,5):
        T=ctower(k);n=len(T);m=n//3
        order=[0]+list(range(m,2*m))+list(range(2*m,3*m))+list(range(1,m))
        D={tuple(sorted((order[i],order[j]))) for i in range(n) for j in range(i+1,n) if T[order[j]][order[i]]}
        check(len(D)==n*(n-1)//6 and components(n,D)==1,'connected optimum C3 tower')
        connected.append(dict(n=n,FAS=len(D),exact_beta=len(D)-n+1))
    paley_towers=[]
    P=paley(7)
    for k in range(4):
        C=ctower(k);m=len(C);n=7*m
        T=[[False if i==j else C[i%m][j%m] if i//m==j//m else P[i//m][j//m] for j in range(n)] for i in range(n)]
        D={(i,j) for i,j in combinations(range(n),2) if T[j][i]}
        check(len(D)==n*(n-1)//6 and components(n,D)==1,'connected P7 C3 tower optimum certificate')
        paley_towers.append(dict(n=n,FAS=len(D),exact_beta=len(D)-n+1))
    return dict(weighted_cases=cases,connected_C3_towers=connected,connected_P7_C3_towers=paley_towers,scope='minimum proof uses audited substitution law; exact edge-count/connectivity certificates checked')

def main():
    print(json.dumps(dict(status='PASS',arithmetic='exact integers and Fraction',root=root_controls(),signed_phase=signed_and_phase_controls(),hardy=hardy_controls(),clocks=clock_controls(),tournaments=tournament_controls(),substitutions=substitution_controls()),indent=2,sort_keys=True))
if __name__=='__main__':main()
