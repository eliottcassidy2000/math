#!/usr/bin/env python3
"""Independent finite controls for the all-greedy theorem, no submitted imports."""
from fractions import Fraction as Q
from itertools import product,combinations
from hashlib import sha256
from pathlib import Path
import json

checks=0
def check(ok,label):
    global checks
    checks+=1
    if not ok:
        raise RuntimeError(label)

def coeff(w):
    a,b=1,0
    for i in range(len(w)-1,-1,-1):
        if w[i]=='1':
            b+=(1<<i)*a
            a*=3
    return a,b,1<<len(w)

def valuation(z):
    if z==0:
        return None
    v=0
    while z%3==0:
        z//=3
        v+=1
    return v

def replay(x,w):
    for z in w:
        check(x%2==int(z),'actual literal parity guard')
        x=(3*x+1)//2 if x%2 else x//2
    return x

def logs(z,n=192):
    lo=Q(0);power=z
    for i in range(n):
        lo+=2*power/(2*i+1)
        power*=z*z
    return lo,lo+2*power/((2*n+1)*(1-z*z))

l2,u2=logs(Q(1,3));l15,u15=logs(Q(1,5))

# Every repeated pair, not merely the algorithm's selected minimal pair,
# can be pulled left until anchored or preceded by unequal symbols.
word_count=0;pair_count=0;overlap_count=0;interior_count=0;shift_count=0
for alphabet,limit in [('01',10),('012',6)]:
    for length in range(1,limit+1):
        for letters in product(alphabet,repeat=length):
            w=''.join(letters);word_count+=1
            for k in range(1,length+1):
                groups={}
                for i in range(length-k+1):
                    groups.setdefault(w[i:i+k],[]).append(i)
                repeat_pairs=[(i,j) for positions in groups.values() for i,j in combinations(positions,2)]
                check(bool(repeat_pairs)==(len(groups)<length-k+1),'exact finite repetition count')
                if repeat_pairs:
                    i,j=min(repeat_pairs)
                    check(i==0 or w[i-1]!=w[j-1],'globally minimal left occurrence lemma')
                for i,j in repeat_pairs:
                    pair_count+=1
                    overlap_count+=j-i<k
                    while i>0 and w[i-1]==w[j-1]:
                        i-=1;j-=1;shift_count+=1
                        check(w[i:i+k]==w[j:j+k],'constructive one-place left shift')
                    check(i==0 or w[i-1]!=w[j-1],'terminal anchored-or-disjoint-predecessor case')
                    interior_count+=i>0

# Hostiles showing why neither overlapping images nor a lost source anchor
# can be silently admitted to the abstract affine-code lemma.
E=lambda x:2*x+1
C=lambda x:Q(x+3,2)
check(E(1)==C(E(1))==3,'overlapping-image maps can repeat at pre-crossing checkpoints')
check(2*Q(1,2)==1,'overlap hostile preserves coefficient prefix invariant')
D=lambda x:Q(x+2,2)
check(D(2)==2 and Q(1,2)<1,'disjoint-image code can return without pre-crossing anchor')
# D on x=2 mod 4 has even outputs, E has odd outputs, so its fixed point
# does not rely on overlapping images; it violates the multiplier anchor.
check(all(D(4*j+2)%2==0 and E(j)%2==1 for j in range(10)),
      'anchor hostile retains disjoint output parity classes')

K='110110110110110101011011101010'
V=K*4+'10'
F='10'
bank={'K':K,'V':V}
for name,w in bank.items():
    h=0
    for i,bit in enumerate(w,1):
        h+=int(bit)
        check(3**h>2**i,'every nonempty literal expansion prefix has multiplier above one')
    a,b,M=coeff(w)
    check(valuation(b-M)==(3 if name=='K' else 4),'last-letter finite image type')

templates=[]
for size in range(1,5):
    templates.extend(''.join(t) for t in product('KV',repeat=size))
templates+=['K'*32,'V'*20,'K'*128,'KV'*7]

def make_control(template,m,force_replay=False):
    U=''.join(bank[z] for z in template)
    ell=len(U);h=U.count('1')
    aU,bU,MU=coeff(U)
    slope=Q(aU,MU)
    t=3 if template[-1]=='K' else 4
    check(valuation(bU-MU)==t<h,'arbitrary-template universal output type')
    lo=(h*l15+(h-ell)*u2)/(u2-l15)
    hi=(h*u15+(h-ell)*l2)/(l2-u15)
    base=lo.__floor__()
    check(base==hi.__floor__(),'exact floor of template logarithmic rotation')
    counts=[];S=Q(1);total=0
    for n in range(1,m+1):
        S*=slope;j=0
        while S>=Q(4,3):
            S*=Q(3,4);j+=1
            check(S>=1,'each conditional reset retains noncrossing')
        check(1<=S<Q(4,3),'completed iteration invariant')
        check(j in (base,base+1),'two adjacent reset counts')
        total+=j
        check((n*lo).__floor__()==(n*hi).__floor__()==total,'all mechanical floors certified independently')
        counts.append(j)
    check(S*Q(3,4)<1,'final forced reset is crossing')
    codes={j:U+F*j for j in (base,base+1)}
    code_coeff={j:coeff(w) for j,w in codes.items()}
    for j,(a,b,M) in code_coeff.items():
        check(b-M==3**j*(bU-MU),'exact affine image type increment across F powers')
        check(valuation(b-M)==t+j,'two coded images have distinct fixed types')
        check(a==3**j*aU and M==4**j*MU,'code affine slope and denominator')
    literal=''.join(codes[j] for j in counts)+F
    a,b,M=coeff(literal)
    check(a<M and len(literal)==m*ell+2*(total+1),'full first-crossing coefficients and length')
    x=(-b*pow(a,-1,M))%M
    check(x>0,'least positive cylinder representative')
    y=x;states=[x];prefix_a,prefix_b,prefix_M=1,0,1
    for j in counts:
        aa,bb,MM=code_coeff[j]
        check((aa*y+bb)%MM==0,'actual complete code endpoint integrality')
        y=(aa*y+bb)//MM
        check(valuation(y-1)==t+j,'actual coded checkpoint image type')
        prefix_a,prefix_b,prefix_M=aa*prefix_a,aa*prefix_b+bb*prefix_M,MM*prefix_M
        check(prefix_a>=prefix_M and prefix_b>0,'strict source anchor keeps complete actual past')
        check(y>x,'every nonempty completed prefix above original source')
        states.append(y)
    check(len(set(states))==m+1,'all actual iteration checkpoints are pairwise distinct')
    check((3*y+1)%4==0,'last F guard')
    endpoint=(3*y+1)//4
    check(M*endpoint==a*x+b and endpoint<x,'direct exact endpoint descent')
    check((M-a)*x>b,'least representative exclusion implies all positive lifts descend')
    check((a*(x+M)+b)//M==endpoint+a<x+M,'second cylinder lift control')
    if force_replay:
        check(replay(x,literal)==endpoint,'full literal replay agrees with branch endpoint')
    symbols=''.join(str(j-base) for j in counts)
    k=(m+99)//100
    factors={symbols[i:i+k] for i in range(m-k+1)}
    check(len(factors)<=k+1,'finite mechanical factor complexity')
    selected=None
    if m>=3:
        check(m>2*k and len(factors)<m-k+1,'strict finite pigeonhole')
        allpairs=[(i,j) for i in range(m-k+1) for j in range(i+1,m-k+1)
                  if symbols[i:i+k]==symbols[j:j+k]]
        i,j=min(allpairs)
        check(i==0 or symbols[i-1]!=symbols[j-1],'minimal-left choice on actual reset sequence')
        future_i=''.join(codes[c] for c in counts[i:i+k])
        future_j=''.join(codes[c] for c in counts[j:j+k])
        r=len(future_i)
        check(future_i==future_j and r>=k*ell and 100*r>=m*ell,'variable-length encoding supplies long literal match')
        check(j+k<=m,'both futures inside finite pre-final-reset stream')
        check(states[i]!=states[j] and (states[j]-states[i])%(1<<r)==0,
              'actual unequal states in same future cylinder')
        if i==0:
            check(states[j]>states[0],'selected source-anchored strict inequality')
        else:
            check(valuation(states[i]-1)!=valuation(states[j]-1),'selected incompatible incoming finite types')
        selected={'i':i,'j':j,'k':k,'r':r}
    elif m==2:
        first=len(codes[counts[0]])
        check(literal[:ell]==literal[first:first+ell]==U,'two-copy match uses start of iteration two')
        check(states[1]>states[0] and (states[1]-states[0])%(1<<ell)==0,
              'two-copy source-anchored unequal checkpoints')
        check(100*ell>=m*ell,'two-copy height threshold')
    else:
        check(template=='V' and len(literal)<=250<2**52,'only small V one-copy case invoked')
    return {'template':template,'m':m,'U_bits':ell,'N':len(literal),'floor_theta':base,
            'reset_counts':sorted(set(counts)),'selected_pair':selected,'literal_replay':force_replay}

records=[]
for template in templates:
    for m in (2,3,17):
        records.append(make_control(template,m,len(template)<=2 and m<=3))
for template,m in [('V',1),('V',101),('V',205),('K',101),('KV',205)]:
    records.append(make_control(template,m,template=='V' and m==1))

# Literal code prefix-freeness is deliberately not assumed.
check(V== (V+F)[:len(V)],'the valid V/VF code is not prefix-free')

basepath=Path(__file__).resolve().parent.parent
manifest=json.loads((basepath/'MANIFEST.json').read_text())
hashes={}
for entry in manifest['files']:
    data=(basepath/entry['path']).read_bytes()
    digest=sha256(data).hexdigest()
    check(digest==entry['sha256'] and len(data)==entry['bytes'],'frozen submitted manifest')
    hashes[entry['path']]=digest
check(hashes['ALL_GREEDY_ENDPOINT_DESCENT.md']=='c7ae20aa091b0d8e56bc50bbaa7ae21f87b25146509678256bc23e4e7cd7cd83',
      'assigned frozen proof identity')
print(json.dumps({'status':'PASS independent all-greedy controls','checks':checks,
    'finite_words':word_count,'all_repeated_pairs':pair_count,'overlapping_pairs':overlap_count,
    'left_shifts_checked':shift_count,'pairs_ending_at_different_predecessors':interior_count,
    'fixed_templates':len(templates),'whole_word_controls':len(records),
    'records':records,'source_hashes':hashes,
    'universal_statements_proved_in_review':True,'arbitrary_template_m1_claimed':False},
    indent=2,sort_keys=True))
