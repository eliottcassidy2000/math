#!/usr/bin/env python3
"""Second-round exact certificates: graph rerooting, seed grammar, slot polynomials."""
import importlib.util
import itertools
import json
import math
from collections import deque
from pathlib import Path
spec=importlib.util.spec_from_file_location('first_round',Path(__file__).with_name('verify.py'))
v=importlib.util.module_from_spec(spec);spec.loader.exec_module(v)
require=v.require


def graph_step(a,r,c):
    a=list(a);n=len(a);a.append(0 if c=='O' else 1<<r)
    for i in range(n):
        if c=='O' or i!=r:a[i]|=1<<n
    return tuple(a), (n if c=='J' else r)


def jstep(s):
    D,E,U,d=s
    return D+U-E,D-E,D+U,D


def seed_witnesses(nmax):
    out={}
    for n in range(2,nmax+1):
        pairs=list(itertools.combinations(range(n),2));full=(1<<n)-1
        for bits in range(1<<len(pairs)):
            a=[0]*n
            for k,(i,j) in enumerate(pairs):
                if bits>>k&1:a[i]|=1<<j
                else:a[j]|=1<<i
            H,dp=v.counts(a)
            rev=tuple(sum(1<<j for j in range(n) if a[j]>>i&1) for i in range(n))
            _,back=v.counts(rev)
            for r in range(n):
                if not v.reachable(a,r):continue
                U=sum(back[m][r]*H[full^m] for m in range(1,full+1) if m>>r&1)
                s=H[full],dp[full][r],U,H[full^(1<<r)]
                out.setdefault(s,(tuple(a),r))
    return out


def grammar(bank,bound):
    seen={s:None for s in sorted(bank) if s[2]<=bound};todo=deque(seen)
    image={s[0] for s in bank if s[0]<=bound};cert={}
    while todo:
        s=todo.popleft()
        for c,t in [('O',v.step(s,'O')),('Z',v.step(s,'Z')),('J',jstep(s))]:
            if t[0]<=bound:
                image.add(t[0]);cert.setdefault(t[0],(s,c))
            if t[2]<=bound and t not in seen:
                seen[t]=(s,c);todo.append(t)
    return image,seen,cert


def factors_closure(image,bound):
    result=set(image)|{1};der={}
    for h in range(3,bound+1,2):
        if h in result:continue
        for a in range(3,math.isqrt(h)+1,2):
            if h%a==0 and a in result and h//a in result:
                result.add(h);der[h]=(a,h//a);break
    return result,der


# Core = directed triangle 0,1,2 followed by 3 beating only 0.
CORE=v.graph_word('ZZ')

def slot_data(path, distinguished=False, fixed_end=False):
    if not path:return 1,0
    s=q=0
    if not distinguished and path[0]==3:s+=1
    if path[-1]!=3 and not fixed_end:s+=1
    for a,b in zip(path,path[1:]):
        allowed=a!=3 and b==3
        direct=bool(CORE[a]>>b&1)
        if not allowed and not direct:return None
        if allowed:s+=1;q+=not direct
    return s,q


def addpoly(out,s,q):
    for j in range(q+1):out[s-j]=out.get(s-j,0)+(-1)**j*math.comb(q,j)


def slot_polynomials():
    result={}
    for name,vertices,end in [('D',range(4),None),('E',range(4),2),('d',[0,1,3],None)]:
        out={}
        for p in itertools.permutations(vertices):
            if end is not None and p[-1]!=end:continue
            sq=slot_data(p,fixed_end=end is not None)
            if sq:addpoly(out,*sq)
        result[name]={b:c for b,c in sorted(out.items()) if c}
    out={}
    for mask in range(16):
        if not mask>>2&1:continue
        P=[i for i in range(4) if mask>>i&1 and i!=2]
        Q=[i for i in range(4) if not mask>>i&1]
        for p in itertools.permutations(P):
            sp=slot_data((2,)+p,True)
            if sp is None:continue
            for q in itertools.permutations(Q):
                sq=slot_data(q)
                if sq is not None:addpoly(out,sp[0]+sq[0],sp[1]+sq[1])
    result['U']={b:c for b,c in sorted(out.items()) if c}
    return result


def main():
    polys=slot_polynomials()
    require(polys=={'D':{0:-1,1:3,2:3},'E':{1:1},'d':{1:1,2:2},'U':{1:-1,2:6,3:4}},('slot polynomials',polys))
    for k in range(8):
        a=v.graph_word('Z'+'O'*k+'Z');s=v.profile(a,2)
        predicted=tuple(sum(c*b**k for b,c in polys[name].items()) for name in ['D','E','U','d'])
        require(s==predicted,('reroot family',k,s,predicted))
        if 1<=k<=5:
            for t in range(1,4):
                aa=a;rr=2
                for _ in range(t):aa,rr=graph_step(aa,rr,'Z')
                expected=4*(2**t-1)*3**k+(4*2**t+1)*2**k+5-2**(t+1)
                require(v.counts(aa)[0][-1]==expected,('closed count',k,t))
    bank5=seed_witnesses(5);bank6=seed_witnesses(6)
    require(set(bank5)==v.seeds(5)[0],'independent seed profile implementation')
    bound=2001
    im5,seen5,cert5=grammar(bank5,bound)
    im6,seen6,cert6=grammar(bank6,bound)
    misses5=[h for h in range(1,bound+1,2) if h not in im5]
    misses6=[h for h in range(1,bound+1,2) if h not in im6]
    require(misses5==[7,21,173],('five-seed grammar',misses5))
    require(misses6==[7,21],('six-seed grammar',misses6))
    baseoz,_=v.bank_image(set(bank5),bound);baseoz|={s[0] for s in bank5}
    mult,der=factors_closure(baseoz,bound)
    require(105 in mult and 173 not in mult,'multiplication control')
    # Explicit all-edge-generated witness with one reroot, ending strongly connected.
    a=v.graph_word('ZOOZ');r=2;sequence=[v.profile(a,r)]
    for _ in range(2):a,r=graph_step(a,r,'Z');sequence.append(v.profile(a,r))
    require(sequence==[(15,1,59,9),(73,9,109,9),(173,9,209,9)],'173 sequence')
    require(all(v.reachable(a,r) for r in range(len(a))),'173 strongness')
    a13=v.graph_word('ZZ');a13,r13=graph_step(a13,2,'Z')
    require(v.counts(a13)[0][-1]==13,'13 reroot repair')
    # Exact finite check of the full-graph minimum-response potential, independent DP.
    controls=0
    for s,(aa,r) in bank5.items():
        mu=min(v.profile(aa,q)[2] for q in range(len(aa)) if v.reachable(aa,q))
        for c in 'OZ':
            bb,rr=graph_step(aa,r,c)
            hh=v.counts(bb)[0][-1]
            mm=min(v.profile(bb,q)[2] for q in range(len(bb)) if v.reachable(bb,q))
            require(mm==2*mu if c=='O' else mm>=mu and hh>s[0],('potential',s,c,mu,mm,hh))
            controls+=1
    print(json.dumps({'status':'PROVED formulas and decreasing potential; finite exact grammar census',
      'slot_polynomials':polys,'literal_family_profiles_checked':8,'literal_suffix_counts_checked':15,
      'seed_profiles':{'5':len(bank5),'6':len(bank6)},
      'grammar_states_through_response_2001':{'5':len(seen5),'6':len(seen6)},
      'all_OZJ_tail_misses_through_2001':{'5':misses5,'6':misses6},
      'OZ_bank5_product_closure_misses_through_2001':[h for h in range(1,bound+1,2) if h not in mult],
      '105_product_factorization':der.get(105),
      'minimum_response_potential_checks':controls,
      '173_witness':{'starting_word':'ZOOZ','reroot_to':2,'suffix':'ZZ','profiles':sequence,'adjacency_bitmasks':a},
      '13_witness':{'starting_word':'ZZ','reroot_to':2,'suffix':'Z','adjacency_bitmasks':a13}},indent=2,sort_keys=True))


if __name__=='__main__':main()
