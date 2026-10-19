#!/usr/bin/env python3
"""Exact controls; universal all-m proof is in the manuscript."""
from fractions import Fraction
import hashlib,json
checks=0

def need(ok,label):
    global checks
    checks+=1
    if not ok:raise RuntimeError(label)

def branch(w):
    a,b,M=1,0,1
    for z in w:
        if z=='1':a,b=3*a,3*b+M
        M*=2
    return a,b,M

def comp(u,v):
    a,b,M=u;c,d,Q=v
    return c*a,c*b+d*M,Q*M

def least_pair(w,k):
    seen={};best=None
    for j in range(len(w)-k+1):
        f=w[j:j+k]
        if f in seen:
            pair=(seen[f],j)
            if best is None or pair<best:best=pair
        else:seen[f]=j
    return best,len(seen)

# Complete finite-word universe, all binary words through length ten and
# every valid factor length. This audits the minimal-left-occurrence lemma,
# not a bounded numerical substitute for the all-m theorem.
finite_words=0;repeated_cases=0;overlap_cases=0;interior_cases=0
for m in range(1,11):
    for z in range(1<<m):
        w=format(z,'0%db'%m);finite_words+=1
        for k in range(1,m+1):
            pair,P=least_pair(w,k)
            need((pair is not None)==(P<m-k+1),'finite pigeonhole iff repeated factor')
            if pair is None:continue
            i,j=pair;repeated_cases+=1
            need(w[i:i+k]==w[j:j+k],'actual equal finite factors')
            if j-i<k:overlap_cases+=1
            if i>0:
                interior_cases+=1
                need(w[i-1]!=w[j-1],'minimal earlier occurrence has different predecessors')
            for u in range(i):
                need(all(w[u:u+k]!=w[v:v+k] for v in range(u+1,m-k+1)),
                     'chosen earlier occurrence is globally minimal')
K='110110110110110101011011101010';V=4*K+'10';F='10'
data={'V':branch(V),'F':branch(F)}
aV,bV,mV=data['V']
need(bV*pow(mV,-1,243)%243==82,'V image class')
need((3*bV+mV)*pow(4*mV,-1,243)%243==1,'VF image class')
need(Fraction(3,4)*Fraction(aV,mV)**9<1<Fraction(3,4)*Fraction(aV,mV)**10,
     'greedy construction nine-ten bound')
records=[]
for m in (1,2,3,10,101,901,1024):
    u=(1,0,1);symbols=[];macro=[]
    for _ in range(m):
        u=comp(u,data['V']);macro.append('V')
        a,b,M=u;reset=3*a>=4*M
        symbols.append('1' if reset else '0')
        if reset:u=comp(u,data['F']);macro.append('F')
        a,b,M=u
        need(M<=a and 3*a<4*M,'all complete iteration multipliers in invariant interval')
    u=comp(u,data['F']);macro.append('F')
    a,b,M=u;s=''.join(symbols)
    literal=''.join(V+('10' if z=='1' else '') for z in s)+'10'
    N=len(literal);L=122*m;k=(m+99)//100
    need(a<M and M==2**N,'final first crossing and length')
    pair,P=least_pair(s,k)
    need(P<=k+1,'mechanical factor complexity control')
    x=-b*pow(a,-1,M)%M
    need(x>0,'positive least lawful source')
    y=x;states=[x]
    for z in s:
        need((aV*y+bV)%mV==0,'V branch guard')
        y=(aV*y+bV)//mV
        if z=='1':
            need((3*y+1)%4==0,'conditional F branch guard')
            y=(3*y+1)//4
        need(y%243==(82 if z=='0' else 1),'typed output for every iteration')
        states.append(y)
    need(len(set(states))==m+1,'reverse-injective checkpoint distinctness')
    need((3*y+1)%4==0,'final F guard');endpoint=(3*y+1)//4
    need(M*endpoint==a*x+b and endpoint<x,'exact endpoint descent in control')
    need((M-a)*x>b,'all positive lifts descend in control cylinder')
    if m<=10:
        z=x
        for bit in literal:
            need(z%2==int(bit),'literal Terras replay guard')
            z=(3*z+1)//2 if z%2 else z//2
        need(z==endpoint,'literal Terras replay endpoint')
    selected=None
    if m>=3:
        need(m>2*k and m-k+1>k+1,'strict finite pigeonhole margin')
        need(pair is not None,'repeated reset factor exists')
        i,j=pair
        future=''.join(V+('10' if z=='1' else '') for z in s[i:i+k])
        future2=''.join(V+('10' if z=='1' else '') for z in s[j:j+k])
        r=len(future)
        need(future==future2 and r>=122*k and 100*r>=L,'variable-length future match reaches height scale')
        need(j+k<=m,'matched futures before final extra F')
        need(states[i]!=states[j] and (states[i]-states[j])%(2**r)==0,
             'actual unequal states share future cylinder')
        if i==0:need(states[j]>x,'anchored later checkpoint strictly above source')
        else:
            need(s[i-1]!=s[j-1],'selected pair has distinct preceding symbols')
            need(states[i]%243!=states[j]%243,'past-image guards remove return alternative')
        selected={'i':i,'j':j,'k':k,'matched_bits':r,'case':'source_anchor' if i==0 else 'different_image_types'}
    else:need(N<2**52,'small words in proved short range')
    records.append({'m':m,'N':N,'L':L,'factor_count':P,'selected_pair':selected,
        'source_bits':x.bit_length(),'literal_replay':m<=10,
        'word_sha256':hashlib.sha256(literal.encode()).hexdigest()})
# Fixed-template controls include templates needing multiple F resets.
def valuation3(z):
    if z==0:return None
    z=abs(z);v=0
    while z%3==0:v+=1;z//=3
    return v

templates={'K':K,'V':V,'KV':K+V,'VK':V+K,'K32':K*32,'V20':V*20,'KV7':(K+V)*7}
fixed=[]
for name,U in templates.items():
    aa,bb,MM=branch(U);ell=len(U);base=3 if name in ('K','VK','K32') else 4
    ca,cM=aa,MM;floor_theta=0
    while 3*ca>=4*cM:ca*=3;cM*=4;floor_theta+=1
    for m in (2,3,17):
        carrier=(1,0,1);resets=[]
        for _ in range(m):
            carrier=comp(carrier,(aa,bb,MM));j=0
            while 3*carrier[0]>=4*carrier[2]:carrier=comp(carrier,data['F']);j+=1
            need(carrier[2]<=carrier[0] and 3*carrier[0]<4*carrier[2],'fixed-template interval invariant')
            need(j in (floor_theta,floor_theta+1),'two adjacent reset counts')
            resets.append(j)
        a,b,M=comp(carrier,data['F']);x=-b*pow(a,-1,M)%M
        need(a<M and x>0,'fixed-template first crossing and positive source')
        values=[x];y=x
        for j in resets:
            need((aa*y+bb)%MM==0,'fixed-template U guard');y=(aa*y+bb)//MM
            for _ in range(j):
                need((3*y+1)%4==0,'fixed-template conditional reset guard');y=(3*y+1)//4
            need(valuation3(y-1)==base+j,'fixed-template disjoint terminal types')
            values.append(y)
        need(len(set(values))==m+1,'fixed-template reverse-injective checkpoints')
        need((3*y+1)%4==0,'fixed-template final reset guard');endpoint=(3*y+1)//4
        need(M*endpoint==a*x+b and endpoint<x,'fixed-template direct endpoint control')
        if m==2:
            need(values[1]>x and (values[1]-x)%(2**ell)==0,'two-copy anchored U match')
        else:
            symbols=''.join(str(j-floor_theta) for j in resets);k=(m+99)//100
            pair,P=least_pair(symbols,k)
            need(P<=k+1 and pair is not None,'fixed-template factor complexity and repetition')
            i,j=pair;r=k*ell+2*sum(resets[i:i+k])
            need(100*r>=m*ell and values[i]!=values[j] and (values[i]-values[j])%(2**r)==0,
                 'fixed-template unequal states with long matched future')
        fixed.append({'template':name,'m':m,'U_bits':ell,'base_type':base,
            'reset_choices':sorted(set(resets)),'floor_theta':floor_theta,'N':M.bit_length()-1})
print(json.dumps({'status':'PASS exact controls; all-m and reverse-injective proofs in manuscript',
 'checks':checks,'complete_binary_words':finite_words,'repeated_factor_cases':repeated_cases,
 'overlapping_cases':overlap_cases,'interior_minimal_cases':interior_cases,'greedy_controls':records,
 'fixed_template_controls':fixed},
 indent=2,sort_keys=True))
