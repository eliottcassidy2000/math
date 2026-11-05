#!/usr/bin/env python3
"""Exact response decoding, affine-code obstructions, and quartic surgery.
Standalone standard-library certificate. All checks survive python -O.
"""
import hashlib,itertools,json,random
from fractions import Fraction as Q
from pathlib import Path

def require(ok,message):
    if not ok: raise ValueError(message)

TRANS=[[[] for _ in range(36)] for _ in range(4)]
ACCEPT=[]
for o,x,y,z in itertools.product(range(2),range(3),range(3),range(2)):
    s=o*18+x*6+y*2+z
    if x!=2 and (y==1 or (y==0 and o==0)): ACCEPT.append(s)
    for c in range(4):
        p=(c>>o)&1; q=(c>>(1-o))&1
        TRANS[c][s].append(o*18+(2 if p else 1)*6+y*2+z)
        if y or p: TRANS[c][s].append(o*18+x*6+(2 if q else 1)*2+z)
        if z or q: TRANS[c][s].append(o*18+x*6+y*2+1)
INIT=(1,)+(0,)*17+(1,)+(0,)*17

def step(v,c):
    out=[0]*36
    for i,a in enumerate(v):
        if a:
            for j in TRANS[c][i]: out[j]+=a
    return tuple(out)
def state(w):
    v=INIT
    for c in w: v=step(v,c)
    return v
def value(w): return sum(state(w)[i] for i in ACCEPT)
def graph(w):
    n=len(w)+2; rows=[0]*n; rows[0]=2
    for j,c in enumerate(w,2):
        for a in (0,1):
            if c>>a&1: rows[a]|=1<<j
            else: rows[j]|=1<<a
        for i in range(2,j): rows[i]|=1<<j
    return rows
def dp(rows):
    n=len(rows); d=[[0]*n for _ in range(1<<n)]
    for v in range(n): d[1<<v][v]=1
    for s in range(1,1<<n):
        for v in range(n):
            if s>>v&1:
                old=s^(1<<v)
                if old: d[s][v]=sum(d[old][u] for u in range(n) if rows[u]>>v&1)
    return sum(d[-1])
def strong(w):
    return bool(w) and w[0]!=0 and w[-1]!=3 and (2 in w or any(a==3 and b==0 for i,a in enumerate(w) for b in w[i+1:]))
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def row(s):
    out=[int(i in ACCEPT) for i in range(36)]
    for c in s[::-1]: out=[sum(out[j] for j in TRANS[c][i]) for i in range(36)]
    return tuple(out)

# Integer coefficients; D_A/36=A, D_B/12=B, D_C/72=1.
DECODERS={
 'A':(36,{'0':-5,'1':-3,'2':6,'00':15,'10':-9,'20':-18,'30':6,'01':9,'000':-4,'100':6,'200':6,'010':-9}),
 'B':(12,{'0':-1,'1':3,'00':3,'10':-3,'01':-3,'000':-2,'010':3}),
 'C':(72,{'':18,'0':13,'1':-21,'2':-12,'00':-39,'10':27,'20':36,'30':-30,'01':45,'11':-18,'21':-36,'31':18,'000':14,'100':6,'200':-12,'010':-27,'210':18})}
SUFFIXES=tuple(sorted(set().union(*(set(c) for den,c in DECODERS.values())),key=lambda x:(len(x),x)))
ROWS={s:row(tuple(map(int,s))) for s in SUFFIXES}
READOUTS={name:tuple(sum(c*ROWS[s][i] for s,c in co.items()) for i in range(36)) for name,(den,co) in DECODERS.items()}
A_POS=(2,4,8,10,14,16)
B_POS=(20,22,26,28,32,34)
C_POS=(0,6,12)

def bits(w,b):
    x=0
    for c in w: x=2*x+((c>>b)&1)
    return x
def canonical(a,b):
    m=max(a.bit_length(),b.bit_length())
    return tuple(((a>>j)&1)+2*((b>>j)&1) for j in range(m-1,-1,-1))
def decode(v):
    vals=[]
    for name in ('A','B'):
        den=DECODERS[name][0]; q,r=divmod(dot(READOUTS[name],v),den)
        require(r==0,'integral decoder'); vals.append(q)
    return canonical(*vals)
def profile(w):
    st=state(w)
    return tuple(dot(ROWS[s],st) for s in SUFFIXES)

class Span:
    def __init__(self): self.rows={}
    def add(self,v):
        a=list(map(Q,v))
        for p,b in sorted(self.rows.items()):
            if a[p]:
                q=a[p]; a=[x-q*y for x,y in zip(a,b)]
        p=next((i for i,x in enumerate(a) if x),None)
        if p is None: return False
        q=a[p]; self.rows[p]=[x/q for x in a]; return True

def reach():
    out=[INIT]; words=[()]; s=Span(); s.add(INIT); i=0
    while i<len(out):
        for c in range(4):
            v=step(out[i],c)
            if s.add(v): out.append(v); words.append(words[i]+(c,))
        i+=1
    return out,words

def physical_assignments(w,order,alloc):
    rr=graph(w); hw=[]
    for j in range(3):
        hw.extend(i+2 for i,g in enumerate(alloc) if g==j)
        if j<2: hw.append(order[j])
    return all(rr[u]>>v&1 for u,v in zip(hw,hw[1:]))

def bulk_check(prefix,kind,k):
    """Independently check every free assignment in the lower-bound injection."""
    tail=tuple(map(int,kind))*k; w=prefix+tail; m=len(prefix)
    alloc=[0]*len(w)
    if kind=='12':
        alloc[m]=1; alloc[m+1]=2; alloc[-2]=1; alloc[-1]=0
        free=range(m+2,len(w)-2)
    else:
        targetbit=1 if kind=='10' else 0
        special=next(i for i,c in enumerate(prefix) if c>>targetbit&1)
        alloc[special]=2 if kind=='10' else 1
        alloc[m]=1 if kind=='10' else 2
        alloc[-3]=0; alloc[-2]=2; alloc[-1]=1
        free=range(m+1,len(w)-3)
    count=0
    for choice in itertools.product(range(3),repeat=2*k-4):
        a=alloc[:]
        for j,g in zip(free,choice): a[j]=g
        require(physical_assignments(w,(0,1),a),'bulk injection creates valid path')
        count+=1
    require(count==3**(2*k-4),'bulk injection cardinality')
    return count

def bank_word(options,mask):
    out=(2,)
    for i,pair in enumerate(options): out+=pair[(mask>>i)&1]
    return out+(0,)
def mobius(values):
    out=[list(v) for v in values]; m=(len(out)-1).bit_length()
    for b in range(m):
        for mask in range(len(out)):
            if mask>>b&1: out[mask]=[a-b for a,b in zip(out[mask],out[mask^(1<<b)])]
    return out

def main():
    receipt={}
    # Literal Hamilton path DP independent of run-state transfer.
    checks=0
    for n in range(6):
        for w in itertools.product(range(4),repeat=n):
            require(value(w)==dp(graph(w)),'independent physical path count'); checks+=1
    receipt['independent_physical_presentations']=checks
    R,W=reach(); require(len(R)==28,'invariant reachable dimension')
    # All 28 basis vectors have all four children checked in reach(), so this
    # is an invariant-span certificate, not a bounded-depth conjecture.
    for st in R:
        require(dot(READOUTS['A'],st)==36*sum(st[i] for i in A_POS),'A decoder identity')
        require(dot(READOUTS['B'],st)==12*sum(st[i] for i in B_POS),'B decoder identity')
        require(dot(READOUTS['C'],st)==72*sum(st[i] for i in C_POS),'constant decoder identity')
        for c in range(4):
            nxt=step(st,c)
            require(dot(READOUTS['A'],nxt)==2*dot(READOUTS['A'],st)+36*(c&1)*sum(st[i] for i in C_POS),'A recurrence')
            require(dot(READOUTS['B'],nxt)==2*dot(READOUTS['B'],st)+12*(c>>1)*sum(st[i] for i in C_POS),'B recurrence')
            require(dot(READOUTS['C'],nxt)==dot(READOUTS['C'],st),'constant recurrence')
    require(dot(READOUTS['C'],INIT)==72,'constant initial value')
    for den,co in DECODERS.values(): require(sum(co.values())==0,'zero-sum decoder weights')
    recon=0
    for n in range(8):
        for w in itertools.product(range(4),repeat=n):
            st=state(w); stripped=tuple(w)
            while stripped and stripped[0]==0: stripped=stripped[1:]
            require(decode(st)==stripped,'canonical word reconstruction')
            require(sum(st[i] for i in A_POS)==bits(w,0),'first physical integer')
            require(sum(st[i] for i in B_POS)==bits(w,1),'second physical integer')
            recon+=1
    receipt.update(reachable_dimension=len(R),reachable_basis_words=[''.join(map(str,w)) for w in W],
        decoder_coefficients=DECODERS,decoder_suffixes=SUFFIXES,canonical_reconstruction_checks=recon)
    # Constant letter readout annihilators after one step. They give no 3^k
    # channel for c=1,2; c=3 is exactly neutral for every physical prefix.
    for st in R:
        for c in range(4):
            v1=step(st,c); v2=step(v1,c); v3=step(v2,c); v4=step(v3,c)
            require(all(a-6*b+11*d-6*e==0 for a,b,d,e in zip(v4,v3,v2,v1)),'full constant-clock annihilator')
        for c in (1,2):
            v1=step(st,c); v2=step(v1,c); v3=step(v2,c)
            require(sum(v3[i]-3*v2[i]+2*v1[i] for i in ACCEPT)==0,'no ternary channel')
        require(sum(step(st,3)[i] for i in ACCEPT)==sum(st[i] for i in ACCEPT),'sink neutrality')
    receipt['bulk_injection_checks']=[{'prefix':''.join(map(str,p)),'tail':kind,'k':k,'count':bulk_check(p,kind,k)}
        for p,kind in [((), '12'),((1,2,0),'12'),((2,),'10'),((1,),'20'),((3,0,2),'10'),((3,0,1),'20')]
        for k in (2,3,4)]
    # Exact frozen-bit formula for the remaining binary pair {0,2}.
    frozen=0
    for n in range(8):
        for w in itertools.product((0,2),repeat=n):
            A=bits(w,1); p=int(bool(w) and w[-1]==2); h=value(w)
            require(value(w+(0,))==h+2*A,'frozen bit append zero')
            require(value(w+(2,))==h+A+2-p,'frozen bit append two')
            frozen+=1
    receipt['frozen_bit_formula_checks']=frozen
    # Three initials plus the full constant-clock annihilator above prove the
    # following distinct-root formulas for every k>=1; later terms are controls.
    clock=[]
    for k in range(1,41):
        u=(2,1,3)+(0,)*k; w=(1,3,2)+(0,)*k
        require(value(u)==-11+3*2**k+12*3**k,'first exact correction family')
        require(value(w)==-9+3*2**k+12*3**k,'second exact correction family')
        require(value(w+(1,))-value(u+(1,))==10-3*2**k,'hidden correction response')
        require(strong(u) and strong(w),'strong correction families')
        if k<=5:
            require(value(u)==dp(graph(u)) and value(w)==dp(graph(w)),'independent correction counts')
            clock.append({'k':k,'H_low':value(u),'H_high':value(w),'difference_after_1':10-3*2**k})
    receipt['exact_plus_two_clock_families']=clock
    # Quartic response, not just scalar H: every retained suffix coordinate.
    rng=random.Random(193); banks=[]
    for m,length in ((8,1),(7,2),(6,3)):
        opts=[]
        for i in range(m):
            pair=tuple(tuple(rng.randrange(4) for _ in range(length)) for _ in range(2)); opts.append(pair)
        vals=[profile(bank_word(opts,e)) for e in range(1<<m)]
        coeff=mobius(vals)
        for mask,c in enumerate(coeff):
            if mask.bit_count()>4: require(not any(c),'all fifth-and-higher interactions vanish')
        for e,v in enumerate(vals):
            got=[0]*len(SUFFIXES)
            for sub,c in enumerate(coeff):
                if sub.bit_count()<=4 and sub&e==sub: got=[a+b for a,b in zip(got,c)]
            require(tuple(got)==v,'quartic reconstruction of full response')
        banks.append({'switches':m,'block_length':length,'all_choices':len(vals),
            'samples_sufficient':sum(__import__('math').comb(m,j) for j in range(5))})
    receipt['quartic_banks']=banks
    sharp=[value((2,)+tuple(c*((e>>i)&1) for i,c in enumerate((1,1,2,2)))+(0,)) for e in range(16)]
    delta=sum((-1)**(4-e.bit_count())*h for e,h in enumerate(sharp))
    require(delta==4,'sharp fourth derivative')
    for e in range(16):
        w=(2,)+tuple(c*((e>>i)&1) for i,c in enumerate((1,1,2,2)))+(0,)
        require(strong(w) and value(w)==dp(graph(w)),'strong physical sharpness bank')
    receipt['sharp_strong_quartic']={'values_mask_order':sharp,'fourth_difference':delta}
    raw=(json.dumps(receipt,sort_keys=True,indent=2)+'\n').encode()
    Path(__file__).with_name('receipt.json').write_bytes(raw)
    print(json.dumps({'status':'PASS','physical_dp_checks':checks+26,'reachable_dimension':len(R),
        'suffix_probes':len(SUFFIXES),'canonical_reconstructions':recon,'quartic_sharpness':delta,
        'receipt_sha256':hashlib.sha256(raw).hexdigest()},sort_keys=True))
if __name__=='__main__': main()
