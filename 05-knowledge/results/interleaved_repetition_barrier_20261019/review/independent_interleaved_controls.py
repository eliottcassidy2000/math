#!/usr/bin/env python3
"""Independent exact audit controls. No imports from the submitted verifiers.

Uses the closed sum over odd-bit positions for affine coefficients, and a
parity-lifting construction (no modular inverse) for lawful sources. Finite
controls are explicitly separated from proof inputs and universal claims.
"""
from fractions import Fraction as Q
from itertools import product
from math import gcd
import hashlib
import json

checks = 0

def check(ok, label):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(label)

K = '110110110110110101011011101010'
bank = {'K':K, 'V':K*4+'10', 'F':'10', 'H':'110101000'}

def coeff(bits):
    # b = sum over odd positions i of 2^i * 3^(number of later odd bits).
    a, b = 1, 0
    for i in range(len(bits)-1, -1, -1):
        if bits[i] == '1':
            b += a << i
            a *= 3
    return a, b, 1 << len(bits)

def parity_lift(bits):
    # At stage i, r is the least source for the first i parity bits and z
    # its exact endpoint. Adding 2^i to r adds odd 3^h to that endpoint.
    r, z, oddpower = 0, 0, 1
    for i, bit in enumerate(bits):
        bit = int(bit)
        if z % 2 != bit:
            r += 1 << i
            z += oddpower
        if bit:
            z = (3*z+1)//2
            oddpower *= 3
        else:
            z //= 2
    return r, z

def guarded_replay(x, bits):
    for bit in bits:
        if x % 2 != int(bit):
            return None
        x = (3*x+1)//2 if bit == '1' else x//2
    return x

def literal(macro):
    return ''.join(bank[c] for c in macro)

def slopes(macro):
    out = [Q(1)]
    for name in macro:
        bits = bank[name]
        out.append(out[-1]*Q(3**bits.count('1'), 1<<len(bits)))
    return out

def digest_int(x):
    return hashlib.sha256(str(x).encode()).hexdigest() if x.bit_length()<14000 else hashlib.sha256(x.to_bytes((x.bit_length()+7)//8,'big')).hexdigest()

data = {k:coeff(v) for k,v in bank.items()}
aK,bK,mK = data['K']
aV,bV,mV = data['V']
lam = Q(aK,mK)
C = Q(bV,aV-mV)
check(data['K'] == (1162261467,4912087879,1073741824), 'closed-sum K coefficients')
check(data['H'] == (81,101,512), 'closed-sum H coefficients')
check(C == Q(83829391903560839836419813590128166431,157489106280555890461927692448373459), 'V center')
check(1<lam<Q(13,12), 'K slope interval')
check(Q(aV,mV)>1 and Q(aV,mV)**30 < lam**122, 'V expansion and bit rate')
check(C<533, 'center ceiling')
for name,bits in bank.items():
    a,b,m = data[name]
    r,z = parity_lift(bits)
    check(a*r+b == m*z, 'parity-lift endpoint')
    check(guarded_replay(r,bits)==z, 'parity-lift direct replay')
    check(Q(b,a)<=Q(len(bits),6), 'macro normalized carry tariff')
    if name != 'F':
        check(Q(b,a)<Q(len(bits),6), 'strict tariff')
    if name in 'KV':
        check(0<Q(b,a-m)<=C, 'individual center bound')
        for i in range(1,len(bits)+1):
            aa,bb,mm=coeff(bits[:i])
            check(aa>mm, 'expansion internal prefix slope')
    else:
        check(Q(a,m)<=Q(3,4) and b<=m-a, 'reset ceiling and value at 1')

# Independent complete finite test of the fuel iff, including all-even and
# contracting words outside the positive-fuel regime of the main theorem.
fuel_words = 0
fuel_residues = 0
for n in range(1,5):
    for chars in product('01', repeat=n):
        bits=''.join(chars)
        a,b,m=coeff(bits)
        for k in (2,3,4):
            repeated=bits*k
            rho,endpoint=parity_lift(repeated)
            check(((a-m)*rho+b)%(m**k)==0, 'fuel at lawful source')
            for y in range(m**k):
                lawful=guarded_replay(y,repeated) is not None
                fuel=((a-m)*y+b)%(m**k)==0
                check(lawful==fuel, 'fuel iff at every residue')
                fuel_residues+=1
            fuel_words+=1

r1,_=parity_lift(K)
r4,_=parity_lift(K*4)
check(r1==815464443 and r4==588943162011782273059357097228697595, 'source floors without modular inverse')
short_rows=[]
for k in (1,2,3):
    for reset in 'FH':
        a,b,m=coeff(K*k+bank[reset])
        gap=(m-a)*r1-b
        check(a<m and gap>0, 'short initial expansion run descends')
        short_rows.append({'word':'K'*k+reset,'descent_gap_at_r1':gap})

# A different logarithm decomposition from the submitted verifier:
# log 3 = log 2 + log(3/2), using z=1/3 and z=1/5, with only 64 terms.
def log_bracket(z, count=64):
    power=z
    total=Q(0)
    for j in range(count):
        total += 2*power/(2*j+1)
        power *= z*z
    return total, total+2*power/((2*count+1)*(1-z*z))
l2,u2=log_bracket(Q(1,3))
l32,u32=log_bracket(Q(1,5))
l3,u3=l2+l32,u2+u32
lo,hi=l2/u3,u2/l3
p,q=5750934602875680,9115015689657667
H0,L0,B=2**52,2**48,10**10
check(Q(1,2)<l2<u2<1 and 1<l3<u3<2, 'elementary log bounds')
check(gcd(p,q)==1 and q>2*H0, 'spacing denominator')
check(Q(p,q)-Q(1,q*q)<lo<hi<Q(p,q)+Q(1,q*q), 'independent rational log certificate')
check(2*H0*q<3*r4, 'short-range exclusion margin')
check(2136*H0*q<r4, 'coarse margin')
check(Q(H0,6*r4)<Q(1,2), 'log conversion smallness')
check(2*30**5*32*2<B, 'cited coefficient majorant')
check(15*L0<H0, 'interval overlap')
check(Q(L0,450)>53*B+52, 'long-range base')
check(L0>450*(B+1), 'long-range derivative')
check(Q(1,200)-Q(1,360)==Q(1,450), 'height lower rate')
check(Q(3,2)**2<Q(4,3)**3, 'reset clock ratio')

# Complete expansion-first bank forest through nine macroblocks. Every leaf
# is constructed from its actual parity cylinder; no repeated-word power or
# affine-composition helper from the submitted script is used.
counts={'first_crossing_leaves':0,'open_depth_9':0,'earlier_context_repeated_factors':0,'repeated_factors':0}
hist={}
leafhash=hashlib.sha256()

def check_leaf(macro, inspect_factors=True):
    bits=literal(macro)
    a,b,m=coeff(bits)
    prefslopes=slopes(macro)
    check(all(s>=1 for s in prefslopes[:-1]) and prefslopes[-1]<1, 'actual first crossing')
    N=len(bits)
    L=sum(len(bank[c]) for c in macro if c in 'KV')
    check(N<=15*L, 'episode length clock')
    check(Q(b,a)<Q(N,6), 'normalized carry ceiling')
    rho,end=parity_lift(bits)
    check(rho>0 and a*rho+b==m*end and end<rho, 'least lawful source descends')
    check(guarded_replay(rho,bits)==end, 'leaf literal replay')
    for lift in (1,3):
        xx=rho+lift*m
        yy=end+lift*a
        check(yy<xx and guarded_replay(xx,bits)==yy, 'positive lift replay')
    Z=Q(rho)+Q(N,6)+C
    prefix_bits=''
    for i in range(len(macro)+1):
        ap,bp,mp=coeff(prefix_bits)
        check(Q(bp,ap)<=Q(len(prefix_bits),6), 'prefix normalized carry ceiling')
        if inspect_factors:
            for width in range(1,(len(macro)-i)//2+1):
                U=macro[i:i+width]
                if any(c not in 'KV' for c in U):
                    continue
                k=2
                while macro[i:i+k*width]==U*k:
                    Ubits=literal(U)
                    au,bu,mu=coeff(Ubits)
                    y=guarded_replay(rho,prefix_bits)
                    sP=Q(ap,mp)
                    fuel=(au-mu)*y+bu
                    check(fuel>0 and fuel%mu**k==0, 'interior integer fuel')
                    check(Q(y)+C<=sP*(Q(rho)+Q(len(prefix_bits),6)+C), 'transport through context')
                    check(Q(y)+C>Q(mu**k,au), 'local fuel height')
                    r=(k-1)*len(Ubits)
                    check(Z>Q(1<<r)/(sP*Q(au,mu)), 'source-preserving height')
                    # Raise the fractional exponent comparison to power 30.
                    check((sP*Q(au,mu))**30<=lam**L, 'full expansion slope ceiling')
                    contextL=sum(len(bank[c]) for c in macro[:i+width] if c in 'KV')
                    check((sP*Q(au,mu))**30<=lam**contextL, 'position-sensitive slope ceiling')
                    counts['repeated_factors']+=1
                    if any(c in 'FH' for c in macro[:i]):
                        counts['earlier_context_repeated_factors']+=1
                    k+=1
        if i<len(macro):
            prefix_bits+=bank[macro[i]]
    episodes=sum(c in 'KV' and (i==0 or macro[i-1] in 'FH') for i,c in enumerate(macro))
    return {'word':macro,'N':N,'L':L,'episodes':episodes,'source_bits':rho.bit_length()}, (rho,end)

stack=['K','V']
while stack:
    w=stack.pop()
    s=slopes(w)[-1]
    if s<1:
        rec,_=check_leaf(w)
        counts['first_crossing_leaves']+=1
        hist[rec['episodes']]=hist.get(rec['episodes'],0)+1
        leafhash.update(json.dumps(rec,sort_keys=True).encode()+b'\n')
    elif len(w)==9:
        counts['open_depth_9']+=1
    else:
        stack.extend(w+c for c in 'KVFH')

# Large genuinely interior repetitions following 10 or 100 explicit resets.
interior_rows=[]
for initial,prior_resets,U,k in ((80,10,'KVK',3),(600,100,'VK',4)):
    context='K'*initial+'FK'*prior_resets
    repeated=U*k
    w=context+repeated
    ss=slopes(w)
    check(all(s>=1 for s in ss), 'interior example context uncrossed')
    while ss[-1]>=1:
        w+='F'
        ss.append(ss[-1]*Q(3,4))
    rec,(rho,end)=check_leaf(w,inspect_factors=False)
    ubits=literal(U);au,bu,mu=coeff(ubits)
    pbits=literal(context);ap,bp,mp=coeff(pbits)
    y=guarded_replay(rho,pbits)
    r=(k-1)*len(ubits)
    Z=Q(rho)+Q(rec['N'],6)+C
    fuel=(au-mu)*y+bu
    check(100*r>=rec['L'], 'interior example density')
    check(fuel>0 and fuel%mu**k==0, 'large interior repetition fuel')
    check(Z>Q(1<<r)/(Q(ap,mp)*Q(au,mu)), 'large context transported height')
    rec.pop('word')
    rec.update({'prior_explicit_resets':prior_resets,'factor':U,'power':k,'duplicated_bits':r})
    interior_rows.append(rec)

# More than a thousand episodes in the stated nonvacuity family, using
# closed-sum coefficients plus actual parity lifting and a full replay.
ell=4096
j=1
s=lam**ell*Q(3,4)
while s>=1:
    j+=1
    s*=lam*Q(3,4)
family='K'*ell+'FK'*(j-1)+'F'
family_bits=literal(family)
a,b,m=coeff(family_bits)
rho,end=parity_lift(family_bits)
check(a<m and (a*rho+b)==m*end and end<rho, 'large family endpoint')
check(guarded_replay(rho,family_bits)==end, 'large family full chronological replay')
check(j>1000 and j<=ell, 'large family episodes')
check(lam*Q(3,4)<1 and lam**2*Q(3,4)<1, 'unbounded family symbolic facts')
check(100*30*(ell-1)>=30*(ell+j-1), 'family density')

# Explicit scope probes. F itself is an expansion-free first crossing and
# returns 1, so the broad short-range lemma needs expansion-first scope.
check(slopes('F')[-1]<1 and guarded_replay(1,'10')==1, 'reset-only counterexample to unqualified short lemma')
check(guarded_replay(1,'1010')==1, 'contracting-square equality hostile')
check(guarded_replay(1,K+'10') is None, 'wrong guard rejected')
vals=[4,1,1,1,1,2,2,1,2,1,1,2,1,1,1,2,3]
hostile=''.join('1'+'0'*(v-1) for v in vals)
aa,bb,mm=coeff(hostile)
check(aa<mm and guarded_replay(165,hostile)==167 and guarded_replay(165,'1000')==31, 'earlier-crossing overshoot hostile')

report={
  'verdict':'All independent controls pass; primary theorem proof audit passes subject to the explicitly cited Matveev bound. Qualify the standalone short-range lemma as expansion-first.',
  'method':'Independent closed-sum affine coefficients, parity-lift source construction, literal replay, and log(3)=log(2)+log(3/2) rational certificate; no submitted code imports.',
  'fuel_equivalence':{'words_and_powers':fuel_words,'exhaustive_residues':fuel_residues,'word_length_max':4,'powers':[2,3,4]},
  'short_initial_controls':short_rows,
  'source_floor':r4,
  'short_certificate':{'terms_per_log':64,'log_arguments':['2','3/2'],'p':p,'q':q,'q_minus_2H0':q-2*H0,'short_margin':3*r4-2*H0*q},
  'long_certificate':{'B':B,'L0':L0,'base_margin':str(Q(L0,450)-(53*B+52)),'derivative_margin':L0-450*(B+1),'safe_Matveev_coefficient':2*30**5*32*2},
  'finite_forest':dict(counts,maximum_macro_depth=9,episode_histogram=hist,leaf_digest_sha256=leafhash.hexdigest()),
  'interior_repetition_controls':interior_rows,
  'large_family':{'initial_K_power':ell,'episodes':j,'N':len(family_bits),'source_bits':rho.bit_length(),'full_parity_replay':True,'source_sha256_bytes':hashlib.sha256(rho.to_bytes((rho.bit_length()+7)//8,'big')).hexdigest()},
  'scope_probe':'W=F, x=1 gives S=3/4 and endpoint 1: outside the main theorem; exclude it explicitly from standalone short-range wording.',
  'checks':checks
}
print(json.dumps(report,sort_keys=True,indent=2))
