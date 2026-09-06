"""Connected-decoder gcd descent and sharp actual-entry controls.

No mathematical producer imports. Incoming hereditary caps are explicit
inherited hypotheses; their large classifier is not rerun here.
"""
from itertools import combinations
from math import gcd
from functools import reduce
import hashlib,json,sys
from pathlib import Path
sys.stdout.reconfigure(newline='\n')
GATES=0
Q=91**6
CAPS={7:90,8:30,9:9,10:4,11:2,12:1}
SIX_TAIL_GCDS={1,2,3,4,5,6,8,9,10,11,12,15,16,17,18,20,22,23,24,25,
              27,29,30,32,33,34,36,40,44,45,46,48,50,51,54,58,60,64,
              66,72,88,90}

def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise RuntimeError(label)

def cgcd(values):return reduce(gcd,values,0)

def allowed_sum(n):
    if n>356:return False
    d=2
    while d*d<=n:
        if n%d==0:
            power=0
            while n%d==0:n//=d;power+=1
            if d%3!=2 or power>2:return False
        d+=1
    return n==1 or n%3==2

def atlas_edge(a,b):
    d=gcd(a,b)
    return allowed_sum((a+b)//d)

def components(values):
    unseen=set(range(len(values)));out=[]
    while unseen:
        todo=[min(unseen)];seen=set()
        while todo:
            i=todo.pop()
            if i in seen:continue
            seen.add(i)
            todo.extend(j for j in unseen-seen if atlas_edge(values[i],values[j]))
        unseen-=seen;out.append(sorted(seen))
    return out

def main():
    inherited_path=Path(__file__).with_name('overnight12_20260906_lrc_decoder_descent_inherited_profiles.json')
    inherited_raw=inherited_path.read_bytes()
    inherited_hash=hashlib.sha256(inherited_raw).hexdigest()
    need(inherited_hash=='935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f',
         'frozen incoming joint-shadow classifier data')
    inherited=json.loads(inherited_raw)
    profiles={int(d):{(c,tuple(gs)) for c,gs in level['profiles']}
              for d,level in inherited['levels'].items()}
    need(set(inherited['levels']['6']['gcds'])==SIX_TAIL_GCDS,'exact inherited seven-body gcd set')
    atlas=[(a,b) for a in range(1,178) for b in range(a+1,357-a)
           if gcd(a,b)==1 and allowed_sum(a+b)]
    need(len(atlas)==5855,'inherited all-scale primitive atlas')
    coefficients={a for pair in atlas for a in pair}
    divisors={d for a in coefficients for d in range(1,a+1) if a%d==0}
    clocks={a*g for a in divisors for g in SIX_TAIL_GCDS}
    need(max(coefficients)==355 and (1,355) in atlas,'sharp oriented edge multiplier')
    need(coefficients==set(range(1,356)) and divisors==coefficients,'complete oriented multiplier universe')
    need(len(SIX_TAIL_GCDS)==42 and max(clocks)==31950,'finite six-subset clock restriction')
    need(len(clocks)==6121,'exact necessary six-subset clock universe size')
    for u in range(1,181):
        for d in range(1,u+1):
            if u%d:continue
            for v in range(1,181):
                ell=d//gcd(d,v);a=u//gcd(u,v)
                need(a%ell==0,'exact outgoing-edge gcd-loss divisibility')
    reports=[]
    patterns={4:(1,3,4,9),5:(1,3,4,9,10),6:(1,3,4,9,10,16)}
    for k,U in patterns.items():
        D=90*355**(7-k)
        top=[D*u for u in U]
        chain=[90*355**j for j in range(7-k)]
        V=top+chain+[2,4,6,12]
        K=max(V);g=2*Q*K+1
        row=V+[g,3*g]
        need(cgcd(V)==2 and cgcd(v//2 for v in V)==1,'physical core t=2 and primitive normalization')
        need(len(V)==11 and len(set(row))==13 and cgcd(row)==1,'primitive distinct physical entry')
        need(sum(row)<=Q*Q,'actual THM3818 physical box')
        need(all(v%7 for v in row),"literal strict safe phase x=1/7")
        need(g>2*Q*K,'every support-three crossing excluded by dominance')
        need(all(gcd(g,v)==1 for v in V),'outside coordinates coprime to core')
        comps=components(row)
        need(sorted(map(len,comps))==[2,11], 'actual decoder components')
        need(any(set(c)==set(range(11)) for c in comps),'claimed core is the actual component')
        need(cgcd(top)==D,'attainment of exact descent bound')
        observed={}
        for size,cap in CAPS.items():
            maximum=0
            for indices in combinations(range(13),size):
                selected=set(indices)
                value=cgcd(row[i] for i in indices);maximum=max(maximum,value)
                need(value<=cap,'complete hereditary subset cap control')
                signature=tuple(sorted(gcd(value,row[j]) for j in range(13) if j not in selected))
                need((value,signature) in profiles[13-size],
                     'complete incoming joint-shadow profile control, including complement gcd word')
            observed[size]=maximum
        for size in (4,5,6):
            for indices in combinations(range(11),size):
                S=set(indices);d=cgcd(V[i] for i in S)
                need(d<=90*355**(7-size),'connected descent cap for every core subset')
                # Retain an explicit adaptive growth path, not just a final scalar.
                chain_losses=[]
                while len(S)<7:
                    candidates=[]
                    for i in S:
                        for j in set(range(11))-S:
                            if atlas_edge(V[i],V[j]):
                                nextd=gcd(d,V[j]);loss=d//nextd
                                candidates.append((loss,i,j,nextd))
                    need(bool(candidates),'connected component supplies outgoing edge')
                    loss,i,j,nextd=min(candidates)
                    a=V[i]//gcd(V[i],V[j])
                    need(a%loss==0 and loss<=355,'actual oriented gcd-loss multiplier')
                    chain_losses.append(loss);S.add(j);d=nextd
                need(d in SIX_TAIL_GCDS,'grown seven-subset lies in inherited exact gcd set')
        reports.append(dict(k=k,top_gcd=D,physical_core=V,t=2,primitive_core=[v//2 for v in V],outside_scale=g,sum=sum(row),
                            maximum=K,hereditary_maxima=observed))
    need(90*355**2==11342250<Q//(42*177),'five-support gcd below pair-closure threshold')
    need(90*355**3==4026498750<Q,'four-support gcd below coefficient height')
    print('STATUS: PASS; connected-decoder gcd descent and sharp safe actual-entry controls')
    print('INHERITED ATLAS:',len(atlas),'pairs;',len(coefficients),'oriented coefficients')
    print('NECESSARY SIX-SUBSET CLOCK UNIVERSE:',len(clocks),'values; maximum',max(clocks))
    print('INHERITED PROFILE SHA256:',inherited_hash)
    print('PHYSICAL CORE GCD BOUNDS: six=31950, five=11342250, four=4026498750')
    for record in reports:print(json.dumps(record,sort_keys=True))
    print('SCOPE: no support-five Bezout to support-three crossing inference; no actual failure claim')
    print('SEMANTIC SHA256:',hashlib.sha256(json.dumps([reports,sorted(clocks)],sort_keys=True).encode()).hexdigest())
    print('ACTIVE GATES:',GATES)

if __name__=='__main__':main()
