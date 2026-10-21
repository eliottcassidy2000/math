#!/usr/bin/env python3
"""Independent exact renewal-law controls. Actual tournaments and formal words
are deliberately separate universes. No assertions; no repository imports."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import permutations, product
from math import comb, factorial
from pathlib import Path
import hashlib, json
import mpmath as mp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT/'session23-root/RENEWAL_LIMIT_LAWS.md'
checks=0

def check(cond,msg):
    global checks
    checks += 1
    if not cond: raise RuntimeError(msg)

def triangle_dist(N,l):
    if N<0:return {}
    return {m:comb(N-2*m,m)*comb(N-2*m+l,l) for m in range(N//3+1)}

def conv(a,b,N):
    out=[0]*(N+1)
    for i,x in enumerate(a):
        for j,y in enumerate(b[:N+1-i]):out[i+j]+=x*y
    return out

def neutral_counts(N):
    a=[0]*(N+1);a[0]=1
    for n in range(1,N+1):a[n]=a[n-1]+(a[n-3] if n>=3 else 0)
    return a

# ACTUAL universe: every labelled tournament through order 6, quotiented by
# explicit vertex-permutation orbits; H and b are literal Hamilton words.
actual_rows=[]; actual1={};actual3={};actual_sc=defaultdict(Counter)
known_counts=[1,1,2,4,12,56]
for n in range(1,7):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    index={p:k for k,p in enumerate(pairs)}
    ps=list(permutations(range(n)))
    def arc(mask,i,j):
        bit=(mask>>index[min(i,j),max(i,j)])&1
        return bool(bit) if i<j else not bool(bit)
    maps=[]
    for p in ps:
        maps.append([(index[min(p[i],p[j]),max(p[i],p[j])], int(p[i]>p[j])) for i,j in pairs])
    seen=set();reps=[];h1=Counter();h3=Counter();degrees=Counter();labels_h1=Counter()
    for mask in range(1<<len(pairs)):
        if mask in seen:continue
        orbit=set()
        for mapping in maps:
            out=0
            for k,(dest,flip) in enumerate(mapping):out|=(((mask>>k)&1)^flip)<<dest
            orbit.add(out)
        seen.update(orbit)
        aut=factorial(n)//len(orbit)
        hwords=[p for p in ps if all(arc(mask,p[i],p[i+1]) for i in range(n-1))]
        H=len(hwords);check(H%aut==0,'actual H/Aut integrality')
        h=H//aut
        blue=0
        for p in hwords:
            rev={p[i]:p[n-1-i] for i in range(n)}
            if all(arc(mask,rev[i],rev[j]) != arc(mask,i,j) for i,j in pairs):blue+=1
        check(blue%aut==0,'actual blue orbit integrality')
        b=blue//aut
        converse=((1<<len(pairs))-1)^mask
        selfconverse=converse in orbit
        check(selfconverse or b==0,'nonself blue count')
        if selfconverse:check((h-b)%2==0,'self converse degree parity')
        degree=(h-b)//2 if selfconverse else h
        degrees[degree]+=Fraction(1,1 if selfconverse else 2)
        reach=[[i==j or (i!=j and arc(mask,i,j)) for j in range(n)] for i in range(n)]
        for k in range(n):
            for i in range(n):
                if reach[i][k]:
                    for j in range(n):reach[i][j] = reach[i][j] or reach[k][j]
        unused=set(range(n));sccs=[]
        while unused:
            i=min(unused);c={j for j in unused if reach[i][j] and reach[j][i]}
            unused-=c;sccs.append(c)
        sizes=sorted(map(len,sccs));m=sizes.count(3)
        if selfconverse:actual_sc[n,h-b][m]+=1
        if h==1:
            check(all(s in (1,3) for s in sizes),'actual h1 neutral SCC classification')
            check(aut==3**m,'actual h1 Aut=3^M')
            h1[m]+=1;labels_h1[m]+=len(orbit)
        if h==3:
            check(sizes.count(5)==1 and all(s in (1,3,5) for s in sizes),'actual h3 atom classification')
            check(aut==5*3**m,'actual h3 Aut=5*3^M')
            h3[m]+=1
        reps.append(mask)
    check(len(seen)==1<<len(pairs),'complete labelled census')
    check(len(reps)==known_counts[n-1],'actual ordinary type count')
    check(all(v.denominator==1 for v in degrees.values()),'actual converse-merged degree integrality')
    check(dict(h1)==triangle_dist(n,0),'actual h1 exact M distribution')
    check(dict(h3)==triangle_dist(n-5,1),'actual h3 exact M distribution')
    actual1[n]=dict(h1);actual3[n]=dict(h3)
    actual_rows.append({'n':n,'labelled':len(seen),'types':len(reps),'h1_M':dict(h1),'h3_M':dict(h3),'degree_nodes':{k:int(v) for k,v in degrees.items()},'h1_labelled_M':dict(labels_h1)})
check(actual_rows[3]['h1_labelled_M']=={0:24,1:16},'actual labelled weighting hostile n4')
check(actual_rows[3]['h1_M']=={0:1,1:2},'actual ordinary weighting n4')

# FORMAL universe only: no assertion that any extra atom is realizable.
# A,B have identical size/h but distinct identities and swap under converse.
# C is larger and fixed; D,E are shorter-length skeleton alternatives.
alphabet={'x':(1,1),'t':(3,1),'A':(5,3),'B':(5,3),'C':(7,3),'D':(9,9),'E':(11,9)}
opp={'x':'x','t':'t','A':'B','B':'A','C':'C','D':'D','E':'E'}
formal=defaultdict(Counter);fgaps=defaultdict(Counter);wordcount=0

def visit(word,n,h):
    global wordcount
    wordcount+=1
    skeleton=tuple(c for c in word if alphabet[c][1]>1)
    m=word.count('t');formal[n,h][skeleton,m]+=1
    gaps=[0]
    for c in word:
        if alphabet[c][1]>1:gaps.append(0)
        else:gaps[-1]+=alphabet[c][0]
    fgaps[n,h,skeleton][tuple(gaps)]+=1
    for c,(size,weight) in alphabet.items():
        if n+size<=18 and h*weight<=27:visit(word+(c,),n+size,h*weight)
visit((),0,1)
for (n,h),dist in formal.items():
    for skeleton in {s for s,m in dist}:
        N=n-sum(alphabet[c][0] for c in skeleton);l=len(skeleton)
        got={m:dist.get((skeleton,m),0) for m in range(N//3+1)}
        check(got==triangle_dist(N,l),'formal literal words vs multinomial exact M')
        for g,count in fgaps[n,h,skeleton].items():
            a=neutral_counts(max(g));expected=1
            for x in g:expected*=a[x]
            check(count==expected,'formal independent gap coefficient product')
        revskel=tuple(opp[c] for c in reversed(skeleton))
        check(all(formal[n,h].get((revskel,m),0)==v for (s,m),v in dist.items() if s==skeleton),'formal converse multiplicities')
# Shorter skeletons are genuinely present for fixed h=9 in the formal model.
check(formal[11,9][(('E',),0)]==1,'formal shorter skeleton retained')
check(sum(v for (s,m),v in formal[10,9].items() if len(s)==2)==4,'formal identical-size atom multiplicity retained')

# Independent exact coefficient recurrence from powers of denominator:
# A(z)^q satisfies p(z)^q times series = 1, not a binomial expansion.
NMAX=1500
a=neutral_counts(NMAX)
coeffs={1:a}
for q in range(2,5):
    p=[1,-1,0,-1];den=[1]
    for _ in range(q):
        out=[0]*(len(den)+3)
        for i,x in enumerate(den):
            for j,y in enumerate(p):out[i+j]+=x*y
        den=out
    v=[1]+[0]*NMAX
    for n in range(1,NMAX+1):v[n]=-sum(den[j]*v[n-j] for j in range(1,min(n,len(den)-1)+1))
    coeffs[q]=v
for q in range(1,5):
    for n in list(range(61))+[100,200,400,800,1500]:
        check(coeffs[q][n]==sum(triangle_dist(n,q-1).values()),'denominator recurrence vs multinomial coefficients')

mp.mp.dps=70
rho=mp.findroot(lambda x:x+x**3-1,mp.mpf('.68'))
c=1/(rho*(1+3*rho*rho));a0=rho*rho
mu=a0/(1+3*a0);sig=a0*(1+a0)/(1+3*a0)**3
check(abs(mu-mp.mpf('0.1942540040'))<mp.mpf('1e-10'),'displayed mu')
check(abs(sig-mp.mpf('0.04956151754'))<mp.mpf('1e-11'),'displayed sigma squared')
# Coefficients / moments / characteristic function finite controls, not proofs.
asymptotic=[]
for q in range(1,5):
    for n in (200,800,1500):
        d=triangle_dist(n,q-1);Z=sum(d.values())
        mean=Fraction(sum(m*v for m,v in d.items()),Z)
        var=Fraction(sum(m*m*v for m,v in d.items()),Z)-mean*mean
        mm=mp.mpf(mean.numerator)/mean.denominator
        vv=mp.mpf(var.numerator)/var.denominator
        ratio=mp.mpf(Z)*rho**n*factorial(q-1)/(c**q*n**(q-1))
        characteristic=sum(mp.mpf(v)*mp.exp(1j*(m-mu*n)/mp.sqrt(sig*n)) for m,v in d.items())/Z
        row={'q':q,'n':n,'leading_ratio':str(ratio),'mean_offset':str(mm-mu*n),'variance_offset':str(vv-sig*n),'cf_error_s1':str(abs(characteristic-mp.exp(mp.mpf('-.5'))))}
        asymptotic.append(row)
        check(abs(mm-mu*n)<q+1,'mean bounded-offset numerical control')
        check(abs(vv-sig*n)<q+1,'variance bounded-offset numerical control')
        if n==1500:check(abs(ratio-1)<mp.mpf('.02'),'fixed-q leading constant convergence control')

# Formal h9 full law: each ordered maximal skeleton retains its identity.
skels=[s for s in product('ABC',repeat=2)]+[('D',),('E',)]
weights={s:rho**sum(alphabet[t][0] for t in s) for s in skels if len(s)==2}
W=sum(weights.values());synthetic=[]
for n in (100,400,1500):
    totals={s:coeffs[len(s)+1][n-sum(alphabet[t][0] for t in s)] for s in skels}
    Z=sum(totals.values());short=sum(totals[s] for s in skels if len(s)==1)
    per={''.join(s):str(mp.mpf(totals[s])/Z) for s in weights}
    maxerr=max(abs(mp.mpf(totals[s])/Z-weights[s]/W) for s in weights)
    synthetic.append({'n':n,'short_probability_times_n':str(mp.mpf(short)*n/Z),'maximum_skeleton_weight_error':str(maxerr),'maximal_skeleton_probabilities':per})
check(mp.mpf(synthetic[-1]['maximum_skeleton_weight_error'])<mp.mpf('.002'),'formal Boltzmann skeleton weight control')
check(weights[('A','A')]==weights[('B','B')],'same-size identities have same individual weights')
check(weights[('A','A')]>weights[('C','C')],'longer skeleton not uniformly weighted')

# Exact simplex moment coefficients: one gap marked by derivatives of A.
# Weighted vector enumeration is unnecessary: sum_g g^k a_g * [z^(N-g)]A^(q-1).
simplex=[]
for q in (2,3,4):
    for n in (100,400,1500):
        Z=coeffs[q][n]
        v1=Fraction(sum(g*a[g]*coeffs[q-1][n-g] for g in range(n+1)),Z*n)
        v2=Fraction(sum(g*g*a[g]*coeffs[q-1][n-g] for g in range(n+1)),Z*n*n)
        expected=Fraction(2,q*(q+1))
        check(v1==Fraction(1,q),'exact exchangeability normalized first moment')
        error=abs(mp.mpf(v2.numerator)/v2.denominator-mp.mpf(expected.numerator)/expected.denominator)
        simplex.append({'q':q,'n':n,'second_moment':str(float(v2)),'Dirichlet_second_moment':str(float(expected)),'error':str(error)})
        if n==1500:check(error<mp.mpf('.002'),'simplex second-moment numerical control')
# r=0 at h=1 is a point mass; for a positive skeleton the finite-n
# sum of normalized gap sizes is (n-s)/n, not exactly one.
check(Fraction(10,15)!=1,'finite-size atom normalization hostile')
# At finite N=3,q=2 gap vectors do not have uniform weight: (0,3) has2, (1,2)1.
check(a[0]*a[3]==2 and a[1]*a[2]==1,'finite-gap uniformity hostile')
# Uniform among all compositions of component COUNTS would yield wrong M rate.
check(abs(mu-mp.mpf('.25'))>mp.mpf('.01'),'renewal vertex clock hostile')

# MIRRORED actual pure-blue family: the class is determined by its first
# half and one of the four inherited centers. Keep parity and triangle mark.
pureblue_centers=[(0,0),(1,0),(3,1),(5,0)]
def mirror_dist(n,centers):
    out=Counter()
    for c,e in centers:
        if n>=c and (n-c)%2==0:
            for m,v in triangle_dist((n-c)//2,0).items():out[2*m+e]+=v
    return out
for n in range(1,7):
    check(actual_sc[n,0]==mirror_dist(n,pureblue_centers),'actual d0 mirror histogram')
    check(actual_sc[n,2]==mirror_dist(n,[(4,0)]),'actual d2 mirror histogram')
# Exact finite reflected-word construction, separate from rational coefficients.
mirror_literal=defaultdict(Counter)
def mirror_visit(half,n):
    for center,c,e in [((),0,0),(('x',),1,0),(('t',),3,1),(('R5',),5,0)]:
        total=2*n+c
        if total<=24:
            literal=half+center+tuple(reversed(half))
            check(sum(3 if z=='t' else 5 if z=='R5' else 1 for z in literal)==total,'literal mirror size')
            mirror_literal[total][literal.count('t')]+=1
    for letter,size in [('x',1),('t',3)]:
        if 2*(n+size)<=24:mirror_visit(half+(letter,),n+size)
mirror_visit((),0)
for n in range(25):check(mirror_literal[n]==mirror_dist(n,pureblue_centers),'literal pureblue mirrors vs exact half formula')
mirror_controls=[]
for n in (100,101,400,401,1500,1501):
    d=mirror_dist(n,pureblue_centers);Z=sum(d.values())
    mean=Fraction(sum(m*v for m,v in d.items()),Z)
    var=Fraction(sum(m*m*v for m,v in d.items()),Z)-mean*mean
    mm=mp.mpf(mean.numerator)/mean.denominator
    vv=mp.mpf(var.numerator)/var.denominator
    check(abs(mm-mu*n)<2,'mirror mean bounded offset')
    check(abs(vv-2*sig*n)<2,'mirror doubled variance bounded offset')
    if n%2==0:
        check(all(m%2==0 for m in d),'empty-center even parity support')
        half=triangle_dist(n//2,0);halfZ=sum(half.values())
        halfmean=Fraction(sum(m*v for m,v in half.items()),halfZ)
        halfvar=Fraction(sum(m*m*v for m,v in half.items()),halfZ)-halfmean*halfmean
        check(mean==2*halfmean and var==4*halfvar,'exact mirror vs independent half variance hostile')
    mirror_controls.append({'n':n,'mean_offset':str(mm-mu*n),'variance_offset_from_2_sigma2_n':str(vv-2*sig*n),'variance_ratio_to_sigma2_n':str(vv/(sig*n))})

# EXACT SAMPLER controls on actual K1/C3/R5 component words. These are
# classified tournament words, unlike the deliberately synthetic alphabet above.
sampler_bins=Counter();sampler_pal=Counter();sampler_orbits=defaultdict(set)
sampler_word_count=0
sizes={'x':1,'t':3,'R':5}
def sampler_visit(word,n,r):
    global sampler_word_count
    sampler_word_count+=1
    a=word.count('x');m=word.count('t');key=(n,r,m)
    sampler_bins[key]+=1
    if word==tuple(reversed(word)):sampler_pal[key]+=1
    else:sampler_orbits[n,r].add(min(word,tuple(reversed(word))))
    left=Counter(word);remaining=len(word);prob=Fraction(1)
    for letter in word:
        prob*=Fraction(left[letter],remaining)
        left[letter]-=1;remaining-=1
    w=factorial(a+m+r)//(factorial(a)*factorial(m)*factorial(r))
    check(prob==Fraction(1,w),'exact sequential multiset probability')
    for letter,size in sizes.items():
        if n+size<=18 and r+(letter=='R')<=3:
            sampler_visit(word+(letter,),n+size,r+(letter=='R'))
sampler_visit((),0,0)
for n in range(19):
    for r in range(4):
        if n<5*r:continue
        N=n-5*r;allwords=pal=0
        for m in range(N//3+1):
            a=N-3*m;w=factorial(a+m+r)//(factorial(a)*factorial(m)*factorial(r))
            odd=a%2+m%2+r%2
            p=0 if odd>1 else factorial((a+m+r)//2)//(factorial(a//2)*factorial(m//2)*factorial(r//2))
            check(sampler_bins[n,r,m]==w,'sampler exact multinomial bin')
            check(sampler_pal[n,r,m]==p,'sampler exact palindrome bin')
            allwords+=w;pal+=p
        check(len(sampler_orbits[n,r])==(allwords-pal)//2 and (allwords-pal)%2==0,'sampler exact reversal quotient')
        if r<2 and n>=max(4,5*r+1):check(allwords>pal,'sampler nonempty surviving subfamily')
# The finite uniform-subset TV identity with explicitly counted full/subsets.
for B in range(1,20):
    for A in range(1,B+1):
        tv=(A*abs(Fraction(1,A)-Fraction(1,B))+(B-A)*Fraction(1,B))/2
        check(tv==1-Fraction(A,B),'exact uniform-subset TV identity')

result={'verdict'  :'PASS exact finite controls; asymptotic observations are not proofs',
'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
'checks':checks,'actual_universe':'All 33,867 labelled tournaments on 1 through 6 vertices; every vertex-permutation orbit checked',
'actual':actual_rows,'formal_universe':{'alphabet':alphabet,'converse':opp,'word_count':wordcount,'max_order':18,'max_h':27,'realizability_claim':False},
'constants':{'rho':str(rho),'mu':str(mu),'sigma2':str(sig),'pole_c':str(c)},
'sampler_controls':{'literal_word_count':sampler_word_count,'max_order':18,'r_range':[0,3],'actual_alphabet':['K1','C3','R5']},
'mirror_controls':mirror_controls,
'asymptotic_controls':asymptotic,'synthetic_skeleton_controls':synthetic,'simplex_moment_controls':simplex}
print(json.dumps(result,indent=2,sort_keys=True))
