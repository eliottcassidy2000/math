"""Standalone exact controls: arithmetic signs, exclusion responses and eleven."""
from fractions import Fraction as Q
from itertools import combinations
from math import gcd,isqrt,prod
import json

def need(b,s):
    if not b:raise RuntimeError(s)

def factors(n):
    out={};p=2
    while p*p<=n:
        while n%p==0:out[p]=out.get(p,0)+1;n//=p
        p+=1
    if n>1:out[n]=out.get(n,0)+1
    return out

def pair(n):
    e=list(factors(n).values());l=(-1)**sum(e)
    return (l if max(e,default=0)<2 else 0),l

def sieve(N):
    mu=[1]*(N+1);lam=mu.copy();prime=[True]*(N+1)
    for p in range(2,N+1):
        if prime[p]:
            for k in range(p,N+1,p):prime[k]=False;mu[k]*=-1
            for k in range(p*p,N+1,p*p):mu[k]=0
            d=p
            while d<=N:
                for k in range(d,N+1,d):lam[k]*=-1
                d*=p
    mu[0]=lam[0]=0
    return mu,lam

def conv(a,b):
    N=min(len(a),len(b))-1;c=[0]*(N+1)
    for d in range(1,N+1):
        for k in range(1,N//d+1):c[d*k]+=a[d]*b[k]
    return c

def arithmetic():
    N=10000;mu,lam=sieve(N);one=[0]+[1]*N;eps=[0,1]+[0]*(N-1)
    sq=[0]*(N+1);msq=sq.copy()
    for d in range(1,isqrt(N)+1):sq[d*d]=1;msq[d*d]=mu[d]
    need(conv(mu,one)==eps,'Mobius inverse')
    need(conv(lam,list(map(abs,mu)))==eps,'Liouville squarefree inverse')
    need(conv(mu,sq)==lam and conv(lam,msq)==mu,'square dilation inverses')
    for n in range(1,N+1):need(pair(n)==(mu[n],lam[n]),'independent trial factors')
    M=[0];L=[0]
    for n in range(1,N+1):M.append(M[-1]+mu[n]);L.append(L[-1]+lam[n])
    for x in range(1,2001):
        need(L[x]==sum(M[x//(d*d)] for d in range(1,isqrt(x)+1)),'summatory square lift')
        need(M[x]==sum(mu[d]*L[x//(d*d)] for d in range(1,isqrt(x)+1)),'summatory inverse')
    for a in range(1,101):
        for b in range(1,101):
            need(mu[a*b]==mu[a]*mu[b]*(gcd(a,b)==1),'coprimality guard')
            need(lam[a*b]==lam[a]*lam[b],'sign character')
    for m in range(1,100):
        need(lam[m]==-lam[2*m],'no fixed Liouville clock')
        need(pair((m+1)**2)[0]==0 and (m+1)**2%m==1%m,'no fixed Mobius clock')
    return {'factor_cross_checks':N,'convolution_range':N,'summatory_range':2000,'multiplication_pairs':10000,
            'eleven_powers':{str(11**j):pair(11**j) for j in range(1,4)},'hostile_H_factorization':factors(462915)}

def matrix_product(A,B):return [[sum(a*b for a,b in zip(r,c)) for c in zip(*B)] for r in A]

def exclusion():
    results=[]
    for k in range(1,7):
        N=1<<k
        A=[[0 if a&b else (-1)**(a.bit_count()+b.bit_count()) for b in range(N)] for a in range(N)]
        inv=[[prod((0 if ((a|b)>>j)&1==0 else -1) for j in range(k)) for b in range(N)] for a in range(N)]
        need(matrix_product(A,inv)==[[int(i==j) for j in range(N)] for i in range(N)],'unimodular tensor inverse')
        signatures={tuple(r) for r in A}|{(0,)*N}
        need(len(signatures)==N+1,'minimal support-plus-failure states')
        signs=[(-1)**a.bit_count() for a in range(N)]
        joint={tuple(A[a])+tuple(signs[a]*s for s in signs) for a in range(N)}
        joint|={(0,)*N+tuple(signs),(0,)*N+tuple(-s for s in signs)}
        need(len(joint)==N+2,'joint carrier retains phase after squarefree failure')
        results.append({'primes':k,'live_states':N,'mu_total_states':N+1,'joint_mu_lambda_states':N+2,
                        'mu_rank':N,'lambda_rank':1,'joint_output_linear_rank':N+1})
    bank=[2,3,11];numbers=[prod(p for j,p in enumerate(bank) if mask>>j&1) for mask in range(8)]
    need(all(pair(n)[0]==pair(n)[1] for n in numbers),'squarefree bank agreement')
    need(pair(6)==pair(22) and pair(3*6)[0]==0 and pair(3*22)[0]==-1,'same observer opposite continuation')
    correlations=[]
    for k in range(3,8):
        C=(1<<(k-1))-1;S=(1<<k)-1
        E=[a for a in range(1<<k) if (a&C).bit_count()%2==0]
        O=[a for a in range(1<<k) if (a&C).bit_count()%2==1]
        for t in range(1<<k):
            if t.bit_count()>k-2:continue
            need(sorted(a&t for a in E)==sorted(a&t for a in O),'all bounded-order joint marginals')
        need(sorted((-1)**a.bit_count() for a in E)==sorted((-1)**a.bit_count() for a in O),'same full sign distributions')
        e=sum(0 if a&S else (-1)**(a.bit_count()+k) for a in E)
        o=sum(0 if a&S else (-1)**(a.bit_count()+k) for a in O)
        need(e==(-1)**k and o==0,'full guard distinguishes parity laws')
        correlations.append({'primes':k,'equal_marginal_order':k-2,'after_full_multiplier':[str(Q(e,len(E))),str(Q(o,len(O)))]})
    return {'tensor_controls':results,'correlation_controls':correlations,'bank_numbers':numbers,'bank_mu':[pair(n)[0] for n in numbers],
            'hostile':{'states':[6,22],'same_mu_lambda':[1,1],'multiplier':3,'mu_outputs':[0,-1]}}

def mul_series(a,b,N):
    out=[0]*(N+1)
    for i,x in enumerate(a):
        if x:
            for j in range(min(len(b)-1,N-i)+1):out[i+j]+=x*b[j]
    return out

def inverse_series(a,N):
    need(a[0]==1,'unit series');out=[1]+[0]*N
    for n in range(1,N+1):out[n]=-sum(a[k]*out[n-k] for k in range(1,min(n,len(a)-1)+1))
    return out

def euler(b,N):
    out=[1]+[0]*N
    for d in range(1,N+1):
        for _ in range(abs(b[d])):
            if b[d]>0:
                for n in range(N,d-1,-1):out[n]-=out[n-d]
            else:
                for n in range(d,N+1):out[n]+=out[n-d]
    return out

def series_controls():
    N=160;mu,lam=sieve(N);P=euler(mu,N);R=euler(lam,N)
    forward=[1]+[0]*N;back=forward.copy()
    for d in range(1,isqrt(N)+1):
        dd=d*d;term=[P[n//dd] if n%dd==0 else 0 for n in range(N+1)]
        forward=mul_series(forward,term,N)
        if mu[d]:
            term=[R[n//dd] if n%dd==0 else 0 for n in range(N+1)]
            back=mul_series(back,term if mu[d]==1 else inverse_series(term,N),N)
    need(forward==R and back==P,'ordinary-product square-dilation duality')
    b=[0]+[2+2*(n%11==0) for n in range(1,N+1)]
    ghost=[0]*(N+1)
    for d in range(1,N+1):
        for n in range(d,N+1,d):ghost[n]+=d*b[d]
    inv=conv(mu,ghost);wrong=conv(lam,ghost)
    need(all(inv[n]==n*b[n] for n in range(1,N+1)),'eta exponent recovery')
    F=euler(b,N);rec=[1]+[0]*N
    for n in range(1,N+1):
        s=-sum(ghost[k]*rec[n-k] for k in range(1,n+1));need(s%n==0,'Newton integrality');rec[n]=s//n
    need(F==rec,'independent coefficient recurrence')
    need(Q(wrong[4],4)==Q(5,2) and Q(wrong[121],121)==Q(486,121),'wrong inverse hostile')
    return {'degree':N,'mu_product_prefix':P[:17],'lambda_product_prefix':R[:17],'eta_divided_by_q_prefix':F[:12],
            'wrong_lambda_inverse':{'b4':str(Q(wrong[4],4)),'b121':str(Q(wrong[121],121))}}

def binary_mod(a,b):
    while a.bit_length()>=b.bit_length():a^=b<<(a.bit_length()-b.bit_length())
    return a

def irreducible(f):
    d=f.bit_length()-1
    return all(binary_mod(f,g)!=0 for k in range(1,d//2+1) for g in range((1<<k)|1,1<<(k+1),2))

def binary_mul(a,b):
    c=0
    while b:
        if b&1:c^=a
        a<<=1;b>>=1
    return c

def polynomial_quotient(a,b):
    q=0
    while a.bit_length()>=b.bit_length():
        k=a.bit_length()-b.bit_length();q^=1<<k;a^=b<<k
    need(a==0,'exact polynomial division');return q

def field_controls():
    mu,lam=sieve(121);neck=[]
    for n in range(1,13):
        primitive=0
        for w in range(1<<n):
            mask=(1<<n)-1
            period=next(k for k in range(1,n+1) if (((w<<k)|(w>>(n-k)))&mask)==w)
            primitive+=(period==n)
        formula=sum(mu[d]*(1<<(n//d)) for d in range(1,n+1) if n%d==0)
        need(formula==primitive and primitive%n==0,'necklace primitive period')
        neck.append({'n':n,'primitive_orbits':primitive//n})
    exact=Q((1<<121)-(1<<11),121);false=exact+Q(2,121)
    need(exact.denominator==1 and false.denominator==121,'prime-square necklace hostile')
    phi11=(1<<11)-1;phi23=(1<<23)-1
    need(irreducible(phi11),'Phi11 irreducible')
    f=next(g for g in range((1<<11)|1,1<<12,2) if binary_mod(phi23,g)==0)
    h=polynomial_quotient(phi23,f);need(irreducible(f) and irreducible(h) and f!=h,'Phi23 split')
    need(binary_mul(f,h)==phi23,'binary factor multiplication')
    orders=[]
    for p in [11,23,89]:
        x=2%p;d=1
        while x!=1:x=x*2%p;d+=1
        seen=set();orbits=[]
        for a in range(1,p):
            if a in seen:continue
            orbit=set();x=a
            while x not in orbit:orbit.add(x);x=2*x%p
            seen|=orbit;orbits.append(orbit)
        pairs=[];used=set();impossible=False
        for i,o in enumerate(orbits):
            if i in used:continue
            neg={(-a)%p for a in o};j=orbits.index(neg)
            if i==j:impossible=True;break
            pairs.append((o,neg));used|={i,j}
        valid=0
        if not impossible:
            for mask in range(1<<len(pairs)):
                D=set().union(*(ab[(mask>>i)&1] for i,ab in enumerate(pairs)))
                need(not D&{(-a)%p for a in D} and len(D)==(p-1)//2 and {2*a%p for a in D}==D,'native circulant orientation')
                valid+=1
        need(valid==(0 if d%2==0 else 1<<((p-1)//(2*d))),'orbit count criterion')
        orders.append({'prime':p,'ord_2':d,'doubling_invariant_tournaments':valid})
    return {'primitive_necklaces':neck,'degree121_primitive_orbits':str(exact),'lambda_wrong_degree121':str(false),
            'Phi11_binary_mask':phi11,'Phi23_binary_factors':[f,h],'doubling_clocks':orders}

def H(rows):
    n=len(rows);dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1<<n):
        for v,a in enumerate(dp[mask]):
            if a:
                left=rows[v]&(((1<<n)-1)^mask)
                while left:b=left&-left;left-=b;dp[mask|b][b.bit_length()-1]+=a
    return sum(dp[-1])

def tournament(mask,n):
    rows=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        u,v=(j,i) if mask>>k&1 else (i,j);rows[u]|=1<<v
    return rows

def join(a,b):
    return [r|(((1<<len(b))-1)<<len(a)) for r in a]+[r<<len(a) for r in b]

def geometric_controls():
    escapes=[]
    for p in range(2,100001):
        if factors(p)!={p:1}:continue
        j=(isqrt(p)+1)//2;hi=(2*j+1)**2
        for k in range(2,hi//p+1):escapes.append([k,p])
    need(escapes==[[2,2],[3,2],[4,2],[2,3],[3,3],[2,11]],'odd-square bracket escapes')
    c=tournament(2,3);need(H(c)==3,'C3 control')
    mask=next(m for m in range(1<<10) if H(tournament(m,5))==11);a=tournament(mask,5)
    examples=[]
    for label,t in [('C3+C3',join(c,c)),('C3+H11',join(c,a)),('H11+H11',join(a,a))]:
        h=H(t);examples.append({'join':label,'n':len(t),'H':h,'mu_lambda':pair(h)})
    need([x['H'] for x in examples]==[9,33,121],'literal join multiplication')
    return {'prime_search_bound':100000,'bracket_escape_pairs':escapes,'H11_mask_downward':mask,'H11_out_rows':a,'joins':examples}

def projected_interval(a,b,d):
    # Push forward the indicator under x -> d*x (mod 1), with weight 1/d.
    pieces=[]
    for r in range(d):
        lo=max(Q(0),d*a-r);hi=min(Q(1),d*b-r)
        if lo<hi:pieces.append((lo,hi,Q(1,d)))
    return pieces

def extra_controls():
    projectors=[]
    for q in [6,30,66,1001]:
        scores=[]
        for d in range(1,q+1):
            if q%d:continue
            A=projected_interval(Q(0),Q(1,2*q),d)
            B=projected_interval(Q(1,q),Q(3,2*q),d)
            I=sum((max(Q(0),min(b,y)-max(a,x))*u*v for a,b,u in A for x,y,v in B),Q(0))
            need(I==(Q(1,2*q*q) if d==q else 0),'one-shot interval corner')
            scores.append((d,I))
        z=sum(pair(d)[0]*I for d,I in scores)
        need(z==Q(pair(q)[0],2*q*q),'Mobius corner amplifier')
        projectors.append({'Q':q,'nonzero_overlap':str(scores[-1][1]),'mobius_sum':str(z)})
    cycles=[]
    for start in [1,-5,-17,-1]:
        n=start;seen=set();odd=even=0;word=[]
        while n not in seen:
            seen.add(n);word.append(n)
            if n%2:n=3*n+1;odd+=1
            else:n//=2;even+=1
        need(n==start,'signed Collatz exact return')
        cycles.append({'root':start,'odd_steps':odd,'halvings':even,'period':len(word),'word':word})
    need([(c['odd_steps'],c['halvings']) for c in cycles]==[(1,2),(2,3),(7,11),(1,1)],'known cycle counts')
    need(pair(3)==pair(17)==(-1,-1) and pair(10)==(1,1) and pair(52)==(0,-1),'affine operation hostile')
    return {'one_shot_projectors':projectors,'signed_collatz_cycles':cycles,
            'affine_hostile':{'odd_inputs':[3,17],'same_pair':[-1,-1],'outputs':[10,52],'output_pairs':[[1,1],[0,-1]]}}

if __name__=='__main__':
    print(json.dumps({'arithmetic':arithmetic(),'exclusion':exclusion(),'series':series_controls(),
                      'fields':field_controls(),'geometry':geometric_controls(),'extra':extra_controls()},sort_keys=True,indent=2))
