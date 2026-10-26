"""81/27, marked Frobenius, Vitali witnesses, and distinct ABC normalizations."""
from collections import Counter
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path
import json
from core import require, edge, converse, relabel, sigma, complement, hamilton, substitute

def triples(t):
    return tuple(c for c in combinations(range(len(t)),3)
                 if all(sum(edge(t,u,v) for v in c if u!=v)==1 for u in c))

def lambdas(t):
    cs=triples(t)
    return tuple(sum(u in c and v in c for c in cs) for u,v in combinations(range(len(t)),2))

def ordinary_mask(n,m):
    t=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if m>>k&1:t[i]|=1<<j
        else:t[j]|=1<<i
    return tuple(t)

def endpoint_matrix(t):
    n=len(t);full=(1<<n)-1;out=[]
    for s in range(n):
        dp=[[0]*n for _ in range(1<<n)];dp[1<<s][s]=1
        for mask in range(1<<n):
            if not mask>>s&1:continue
            for v,val in enumerate(dp[mask]):
                if not val:continue
                allowed=t[v] & (full^mask)
                while allowed:
                    bit=allowed&-allowed;allowed^=bit
                    dp[mask|bit][bit.bit_length()-1]+=val
        out.append(dp[full])
    return out

def gf_mul(a,b):
    r=0
    while b:
        if b&1:r^=a
        b>>=1;a<<=1
        if a&8:a^=11
    return r

def gf_f(x):return gf_mul(x,x)

def amplification():
    receipt=json.loads((Path(__file__).resolve().parent/'atlas-normal.json').read_text(encoding='utf-8-sig'))
    w=receipt['first_paired_port_guard_failure'];x,y=map(tuple,w['rows']);probe=tuple(w['permutation'])
    targets=[complement(x),complement(y)]
    def rooted_iso(a,b):
        return next((0,)+p+(6,) for p in permutations(range(1,6)) if relabel(b,(0,)+p+(6,))==a)
    src=rooted_iso(x,y);dst=rooted_iso(*targets)
    out=[]
    for k in range(1,7):
        def wrapped(bits):
            n=7*k+2;u=n-2;v=n-1;rows=[0]*n
            for i,b in enumerate(bits):
                a=[x,y][b]
                for j in range(7):
                    rows[7*i+j]=(a[j]<<(7*i))|sum(1<<z for z in range(7*(i+1),7*k))|(1<<u)
            rows[u]=1<<v;rows[v]=(1<<(7*k))-1
            return tuple(rows)
        base=wrapped([0]*k);target_base=complement(base);fingerprints=set()
        for bits in product([0,1],repeat=k):
            t=wrapped(bits);s=complement(t);n=len(t)
            fs=tuple(7*i+(src[j] if bits[i] else j) for i in range(k) for j in range(7))+(n-2,n-1)
            ft=tuple(7*i+(dst[j] if bits[i] else j) for i in range(k) for j in range(7))+(n-2,n-1)
            require(relabel(t,fs)==base and relabel(s,ft)==target_base,'amplified paired rooted types')
            require(sigma(t)!=t,'amplification is black')
            fingerprints.add(t)
            for i,b in enumerate(bits):
                p=tuple(range(7*i))+tuple(7*i+j for j in probe)+tuple(range(7*(i+1),n))
                require(all(edge(t,a,z) for a,z in zip(p,p[1:]))==bool(b),'one local guard reads each bit')
        require(len(fingerprints)==2**k,'independent alignment bits')
        out.append({'copies':k,'n':7*k+2,'states':len(fingerprints),'one_bit_guards':k})
    return {'source_rooted_isomorphism':src,'target_rooted_isomorphism':dst,'controls':out}

def two_outer_ports_fail():
    # Same source child, quotient run word, number of runs, and complemented
    # child with its two outer roots. Inner ports/coincidences still matter.
    from core import canonical
    a=(2,4,9,3);t,owner=substitute((2,4,1),[a,(0,),(0,)])
    ps=[(0,1,2,4,5,3),(1,2,4,5,3,0)];records=[];keys=[]
    for p in ps:
        runs=[]
        for v in p:
            if not runs or owner[runs[-1][-1]]!=owner[v]:runs.append([v])
            else:runs[-1].append(v)
        child_runs=[r for r in runs if owner[r[0]]==0];b=list(converse(a))
        for r in child_runs:
            for u,v in zip(r,r[1:]):b[u]^=1<<v;b[v]^=1<<u
        outer=(child_runs[0][0],child_runs[-1][-1]);ports=tuple(v for r in child_runs for v in (r[0],r[-1]))
        keys.append(canonical(tuple(b),outer));u=complement(t,p)
        records.append({'path':p,'runs':runs,'child_target':b,'outer_roots':outer,'all_ports':ports,
                        'target':u,'H':hamilton(u),'blue':sigma(relabel(t,p))==relabel(t,p)})
    require(keys[0]==keys[1] and [r['H'] for r in records]==[15,23],'two outer ports lose run attachment')
    return {'source_child':a,'source':t,'quotient_run_word':[0,1,2,0],'records':records}

def main():
    ap=(10,12,9,0);am=converse(ap);q=(2,4,1)
    ts=[substitute(q,[ap]*3)[0],substitute(q,[am,ap,ap])[0]]
    require(triples(ts[0])==triples(ts[1]),'same labelled C3 field')
    control=[]
    for t,word,aut in zip(ts,[(0,0,0),(1,0,0)],[81,27]):
        # Canonical four-block recovery in the proof makes this the full group.
        autos=[]
        for shift in range(3):
            if any(word[i]!=word[(i+shift)%3] for i in range(3)):continue
            for internal in product(range(3),repeat=3):
                p=tuple(4*((i+shift)%3)+(v+internal[i])%3 if v<3 else 4*((i+shift)%3)+3
                        for i in range(3) for v in range(4))
                require(relabel(t,p)==t,'structured automorphism');autos.append(p)
        require(len(autos)==aut,'81/27 group orders')
        em=endpoint_matrix(t);H=sum(map(sum,em));require(H==hamilton(t)==462915,'independent path counts')
        flux=Counter()
        for s in range(12):
            for e in range(12):flux[t[s].bit_count()-t[e].bit_count()-1]+=em[s][e]
        flux={k:v for k,v in sorted(flux.items()) if v}
        require(all(v%aut==0 for v in flux.values()),'flux orbit integrality')
        root_orbits=[];unseen=set(range(12))
        while unseen:
            v=min(unseen);orb={p[v] for p in autos};unseen-=orb;root_orbits.append(sorted(orb))
        control.append({'word':word,'H':H,'Aut':aut,'presentations':H//aut,
                        'root_orbits':root_orbits,'flux':flux,
                        'root_stabilizer_orders':[sum(p[o[0]]==o[0] for p in autos) for o in root_orbits]})
    # All three normal-basis vectors in F8, with chosen normal element 3.
    basis=[3,gf_f(3),gf_f(gf_f(3))]
    encode=lambda w:(basis[0]*w[0])^(basis[1]*w[1])^(basis[2]*w[2])
    def rotate(w):return (w[2],w[0],w[1])
    def reflect(w):return (w[0],w[2],w[1])
    def S(w):return tuple(1-b for b in reflect(w))
    frames=[]
    for eps in [1,-1]:
        for i in range(3):
            w=tuple(int(j==i) if eps==1 else int(j!=i) for j in range(3))
            f=(i,(i+eps)%3,(i+2*eps)%3)
            fw=rotate(w);sw=S(w)
            want_f=tuple(int(j==(i+1)%3) if eps==1 else int(j!=(i+1)%3) for j in range(3))
            want_s=tuple(int(j==(-i)%3) if eps==-1 else int(j!=(-i)%3) for j in range(3))
            require(fw==want_f and sw==want_s,'S3 equivariant bijection')
            require(encode(fw)==gf_f(encode(w)),'normal basis Frobenius')
            require(rotate(rotate(rotate(w)))==w and S(S(w))==w and S(rotate(S(w)))==rotate(rotate(w)),'dihedral action')
            require(len({w,rotate(w),rotate(rotate(w)),S(w),rotate(S(w)),rotate(rotate(S(w)))})==6,'regular S3')
            frames.append({'sign':eps,'path':f,'word':w,'field':encode(w)})
    marked=[]
    for bits in product([0,1],repeat=3):
        w=tuple(bits);x=encode(w);d=1 if gf_f(x)==x else 3
        for i in range(3):
            pair=(x,basis[i]);orbit=[]
            while pair not in orbit:orbit.append(pair);pair=tuple(gf_f(a) for a in pair)
            require(len(orbit)==3,'marked phase always full orbit')
        marked.append({'word':w,'field':x,'degree':d})
    # Independently documented labelled-lambda twins; do not replace simple
    # cycle counts by traces or treat finite observer loss as nonmeasurability.
    vt=[ordinary_mask(7,m) for m in [1354195,1351956]]
    require(lambdas(vt[0])==lambdas(vt[1]),'Vitali lambda equality')
    vh=[hamilton(t) for t in vt];require(vh==[109,111],'Vitali H witness')
    def divisor_sum(n):return sum(d for d in range(1,n+1) if n%d==0)
    def norm(a,b):
        g=gcd(a,b);return [a//g,b//g,(a+b)//g]
    abc=[]
    for n in [6,12,24,28,72,120,496,8128]:
        power=n&-n;odd=n//power;s=divisor_sum(n)-n
        abc.append({'n':n,'odd_core_split':norm(power,n-power),
                    'inherited_square_later_map':[odd-power,power,odd] if power>=2 and odd>power else None,
                    'aliquot_split':norm(n,s),'abundancy_num_den':[divisor_sum(n)//gcd(n,divisor_sum(n)),n//gcd(n,divisor_sum(n))],
                    'perfect':s==n,'odd_core':odd})
    # Same odd core does not descend through aliquot continuation.
    require((6//(6&-6))==(12//(12&-12))==3,'same doubling quotient')
    require((divisor_sum(6)-6)==6 and (divisor_sum(12)-12)==16,'different odd cores after aliquot')
    require(norm(12,divisor_sum(12)-12)==[3,4,7] and norm(56,divisor_sum(56)-56)==[7,8,15],
            'same perfect normalization fails after doubling')
    words=[(1,0,0)*3,(1,1,1,0,0,0,0,0,0)]
    periods=[next(d for d in range(1,10) if all(w[i]==w[(i+d)%9] for i in range(9))) for w in words]
    require([sum(w) for w in words]==[3,3] and periods==[3,9],'weight moments miss Frobenius degree')
    return {'carousel_control':control,'labelled_cyclic_triples':len(triples(ts[0])),
            'F8_basis':basis,'S3_frame_bridge':frames,'marked_Frobenius':marked,
            'Vitali_control':{'masks':[1354195,1351956],'H':vh,'same_labelled_lambda':True,
                              'same_labelled_C3':triples(vt[0])==triples(vt[1])},
            'alignment_amplification':amplification(),
            'two_outer_ports_failure':two_outer_ports_fail(),
            'equal_weight_different_degree':{'words':words,'periods':periods},
            'ABC_maps':abc,'ABC_scope':'square-later map inherited from incoming session69; odd-core split is a distinct comparison',
            'ABC_continuation_hostiles':{'same_odd_core':[6,12],'aliquot_odd_cores':[3,1],
                                        'perfect_sources':[6,28],'after_doubling_aliquot_triples':[[3,4,7],[7,8,15]]}}

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
