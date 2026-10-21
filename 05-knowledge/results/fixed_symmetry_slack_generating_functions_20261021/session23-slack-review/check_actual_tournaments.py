#!/usr/bin/env python3
"""Independent focused audit: direct actual tournaments, no imports from candidate.
Universe: every labelled tournament <=4; seeded arbitrary tournaments 5..8;
structured actual tournament controls <=15; all 27 C3 weight profiles 1..3.
No large general tournament census. All claims checked exactly.
"""
from collections import Counter
from itertools import product, permutations
from pathlib import Path
from random import Random
import hashlib, json, time

CHECKS=Counter()
RECORDS=[]
def need(x,label):
    CHECKS[label]+=1
    if not x: raise AssertionError(label)
def bits(m):
    while m:
        b=m&-m; yield b.bit_length()-1; m-=b

def tournament(n,code):
    a=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if code&1:a[i]|=1<<j
            else:a[j]|=1<<i
            code>>=1
    return tuple(a)
def converse(a):
    return tuple(((1<<len(a))-1)^r^(1<<i) for i,r in enumerate(a))
def sub(a,vs):
    return tuple(sum(1<<j for j,v in enumerate(vs) if a[u]>>v&1) for u in vs)
def inflate(q,parts):
    offsets=[]; n=0
    for a in parts: offsets.append(n); n+=len(a)
    rows=[0]*n
    for i,a in enumerate(parts):
        for u,r in enumerate(a):
            rows[offsets[i]+u]=r<<offsets[i]
            for j,b in enumerate(parts):
                if q[i]>>j&1:rows[offsets[i]+u]|=((1<<len(b))-1)<<offsets[j]
    return tuple(rows)
def transitive(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
C3=(2,4,1)
def word(w):return inflate(transitive(len(w)),[C3 if x==3 else (0,) for x in w])
def join(parts):return inflate(transitive(len(parts)),parts)

def feedback(a,weights=None,long=frozenset()):
    """Subset DP over the last vertex, not representative sampling or orders."""
    n=len(a); weights=weights or (1,)*n
    dp=[0]+[10**9]*((1<<n)-1); last=[None]*(1<<n)
    for m in range(1,1<<n):
        for v in bits(m):
            p=m^(1<<v); backward=a[v]&p
            if backward and (v in long or any(j in long for j in bits(backward))):continue
            val=dp[p]+weights[v]*sum(weights[j] for j in bits(backward))
            if val<dp[m]:dp[m]=val;last[m]=v
    order=[];m=(1<<n)-1
    while m and last[m] is not None:
        v=last[m];order.append(v);m^=1<<v
    return dp[-1],tuple(reversed(order))

def maps(a,b=None,colors_a=None,colors_b=None):
    """All actual graph isomorphisms by adjacency-consistent vertex extension."""
    b=a if b is None else b
    n=len(a)
    if len(b)!=n:return []
    ca=colors_a or (None,)*n; cb=colors_b or (None,)*n
    da=[r.bit_count() for r in a];db=[r.bit_count() for r in b]
    candidates=[tuple(j for j in range(n) if da[i]==db[j] and ca[i]==cb[j]) for i in range(n)]
    order=sorted(range(n),key=lambda i:(len(candidates[i]),-da[i]))
    out=[];p=[-1]*n
    def rec(k,used):
        if k==n:out.append(tuple(p));return
        i=order[k]
        for j in candidates[i]:
            if used>>j&1:continue
            if all(((a[i]>>h)&1)==((b[j]>>p[h])&1) for h in order[:k]):
                p[i]=j;rec(k+1,used|(1<<j));p[i]=-1
    rec(0,0);return out

def module(a,m):
    return all((a[v]&m) in (0,m) for v in range(len(a)) if not(m>>v&1))
def components(a,m):
    reach={}
    for v in bits(m):
        r=1<<v;todo=r
        while todo:
            u=(todo&-todo).bit_length()-1;todo^=1<<u
            new=a[u]&m&~r;r|=new;todo|=new
        reach[v]=r
    left=m;out=[]
    while left:
        v=(left&-left).bit_length()-1
        c=sum(1<<u for u in bits(left) if reach[v]>>u&1 and reach[u]>>v&1)
        out.append(c);left^=c
    out.sort(key=lambda c:-(a[next(bits(c))]&(m^c)).bit_count())
    return out

def w_type(a,m):
    cs=components(a,m)
    if any(c.bit_count() not in (1,3) for c in cs):return None
    return tuple(c.bit_count() for c in cs)
def canonical(a):
    wm=[m for m in range(1,1<<len(a)) if module(a,m) and w_type(a,m) is not None]
    maximal=[m for m in wm if not any(m!=b and m&b==m for b in wm)]
    need(sum(m.bit_count() for m in maximal)==len(a),'maximal W modules partition by size')
    need(all(not(x&y) for i,x in enumerate(maximal) for y in maximal[i+1:]),'maximal W modules disjoint')
    for i,x in enumerate(wm):
        for y in wm[i+1:]:
            if x&y:
                need(module(a,x|y) and w_type(a,x|y) is not None,'all overlapping W modules have W union')
    q=sub(a,[next(bits(m)) for m in maximal]);ws=tuple(w_type(a,m) for m in maximal)
    return maximal,q,ws

def inflation_is_word_by_flags(q,ws,mask):
    for c in components(q,mask):
        k=c.bit_count()
        if k==1:continue
        if k!=3 or any(ws[i]!=(1,) for i in bits(c)):return False
    return True

def analyze(a,name):
    n=len(a);aut=maps(a);s=sum(any(p[i]!=i for p in aut) for i in range(n));f,order=feedback(a);d=3*f-s
    need(d>=0,'nonnegative slack')
    ms,q,ws=canonical(a);weights=tuple(sum(w) for w in ws);ts=tuple(w.count(3) for w in ws)
    Uwords=tuple(w.count(1) for w in ws)
    kg=maps(q,colors_a=ws,colors_b=ws)
    J={i for i in range(len(q)) if any(p[i]!=i for p in kg)}
    C,qorder=feedback(q,weights)
    U=sum(Uwords[i] for i in J);V=sum(weights[i] for i in J)
    need(f==sum(ts)+C,'actual feedback equals straightened weighted cost')
    need(s==3*sum(ts)+U,'actual support equals internal triangle support plus moved singleton mass')
    need(d==3*C-U,'actual slack equals 3C-U')
    need(len(aut)==len(kg)*3**sum(ts),'actual full automorphism group order factors')
    need(V<=2*d,'moved canonical block whole mass is bounded by 2delta')
    need(C<=d,'canonical quotient weighted cost at most delta')
    need(len(q)<=4*d+1,'canonical quotient at most 4delta+1')
    for mask in range(1,1<<len(q)):
        if mask.bit_count()<2 or not module(q,mask):continue
        actual=w_type(a,sum(ms[i] for i in bits(mask))) is not None
        need(actual==inflation_is_word_by_flags(q,ws,mask),'maximality flags match actual induced inflation')
        need(not actual,'canonical blocks pass exact maximality test')
    actual_sc=bool(maps(a,converse(a)))
    if d:
        L={i for i,w in enumerate(weights) if w>d};colors=tuple(None if i in L else ws[i] for i in range(len(q)))
        C0,_=feedback(q,weights,L)
        K0=maps(q,colors_a=colors,colors_b=colors)
        need(C0==C,'admissible-order compiler cost equals actual weighted cost')
        need(set(K0)==set(kg),'short-word long-flag stabilizer equals actual stabilizer')
        need(all(all(p[i]==i for i in L) for p in K0),'every long slot fixed by configuration automorphisms')
        need(len(L)<=2*d+1,'long slot bound')
        reverse_colors=tuple(None if i in L else tuple(reversed(ws[i])) for i in range(len(q)))
        antis=maps(q,converse(q),colors_a=colors,colors_b=reverse_colors)
        predicted_sc=False
        if antis:
            tau=antis[0]
            need(all(tuple(p[i] for i in sorted(L))==tuple(tau[i] for i in sorted(L)) for p in antis),'configuration anti-isomorphism has unique action on long slots')
            need(all(tau[tau[i]]==i for i in L),'anti-isomorphism action squares to identity on long slots')
            need(sum(tau[i]==i for i in L)==len(L)%2,'anti-isomorphism long slots have at most one middle slot')
            predicted_sc=all(ws[tau[i]]==tuple(reversed(ws[i])) for i in L)
        need(predicted_sc==actual_sc,'palindrome and paired-reversal test equals actual self-converse')
    elif n:
        need(len(q)==1,'zero slack canonical quotient is one word')
        need(actual_sc==(ws[0]==tuple(reversed(ws[0]))),'zero slack self-converse iff word palindrome')
    r={'name':name,'n':n,'f':f,'s':s,'delta':d,'aut_order':len(aut),'canonical_block_words':[list(w) for w in ws], 'canonical_quotient_size':len(q),'weighted_feedback':C,'moved_singleton_mass_U':U,'moved_block_mass_V':V,'self_converse':actual_sc}
    RECORDS.append(r);return r

def main():
    start=time.time()
    for n in range(1,5):
        for code in range(1<<(n*(n-1)//2)):analyze(tournament(n,code),f'all-labelled-{n}-{code}')
    rng=Random(23016271)
    for n in range(5,9):
        for j in range(12):analyze(tournament(n,rng.getrandbits(n*(n-1)//2)),f'seeded-arbitrary-{n}-{j}')
    # Mandatory synchronized triangles, with row matching rather than homogeneous cross arcs.
    sync=list(join([C3,C3]))
    for i in range(3):
        for j in range(3):
            if i!=j:sync[i]^=1<<(3+j);sync[3+j]^=1<<i
    r=analyze(tuple(sync),'diagonal-C3-synchronization')
    need((r['f'],r['s'],r['delta'],r['aut_order'])==(4,6,6,3),'mandatory synchronization hostile exact invariants')
    need(r['canonical_quotient_size']==3 and r['canonical_block_words']==[[1,1]]*3,'synchronization canonical quotient is C3 on three transitive pairs')
    for j,w in enumerate([(1,3),(3,1),(1,3,1),(3,3),(1,1,1,1)]):
        core=inflate(C3,[(0,),word(w),(0,)])
        r=analyze(core,f'long-central-delta3-{j}')
        need(r['delta']==3,'C3[1,longword,1] positive delta3 family')
    for j,(left,central,right) in enumerate([
        ((1,1,1,1),(1,1),(1,1,1,1)),
        ((1,3),(1,1),(3,1)),
        ((1,3),(1,1),(1,3)),
        ((1,1,1,1),(1,3,1),()),
        ((),(1,3,1),(1,1,1,1)),
        ((1,1,1,1),(1,3,1),(1,1,1,1)),
    ]):
        a=join([word(left),inflate(C3,[(0,),word(central),(0,)]),word(right)])
        r=analyze(a,f'long-palindrome-channel-{j}')
        need(r['delta']==3,'outer words preserve delta3')
    for j,w in enumerate([(3,),(1,3)]):
        r=analyze(inflate(C3,[word(w)]*3),f'moved-whole-block-not-singleton-mass-{j}')
        need(r['moved_block_mass_V']>r['moved_singleton_mass_U'],'hostile distinguishes moved whole block mass from moved singleton mass')
    r5=tuple(sum(1<<((i+k)%5) for k in (1,2)) for i in range(5))
    need(analyze(r5,'regular-five')['delta']==4,'regular-five delta4')
    analyze(join([word((1,3)),r5,word((3,1))]),'regular-five-with-word-tails')
    # Test exact maximality on proposed partitions, including deliberately nonmaximal ones.
    # This is distinct from checking only the already-canonical partition above.
    proposed_count=0
    short_words=[(1,),(1,1),(1,1,1),(3,)]
    for qn in range(1,4):
        for code in range(1<<(qn*(qn-1)//2)):
            q=tournament(qn,code)
            for ws in product(short_words,repeat=qn):
                a=inflate(q,[word(w) for w in ws]); offset=0; bs=[]
                for w in ws:
                    bs.append(((1<<sum(w))-1)<<offset);offset+=sum(w)
                predicted=not any(module(q,m) and inflation_is_word_by_flags(q,ws,m)
                    for m in range(1,1<<qn) if m.bit_count()>=2)
                actual=not any(module(a,m) and w_type(a,m) is not None and
                    any(m!=b and m&b==b for b in bs)
                    for m in range(1,1<<len(a)))
                need(predicted==actual,'proposed-block maximality iff flag test including rejected configurations')
                proposed_count+=1
    # Direct feedback oracle for arbitrary homogeneous modules, including non-W ones.
    profile_count=0
    for weights in product(range(1,4),repeat=3):
        parts=[transitive(w) for w in weights]
        a=inflate(C3,parts)
        need(feedback(a)[0]==feedback(C3,weights)[0],'all 27 cyclic weight profiles straighten exactly')
        profile_count+=1
    for q in [C3,tournament(4,19),r5]:
        parts=[(C3 if i%2==0 else (0,)) for i in range(len(q))]
        a=inflate(q,parts)
        need(feedback(a)[0]==sum(feedback(p)[0] for p in parts)+feedback(q,tuple(map(len,parts)))[0],'cyclic nontransitive homogeneous modules straighten exactly')
    # Test finite-configuration isomorphism counting on different long fillings.
    choices=[(1,1,1,1),(1,3),(3,1)]
    family=[inflate(C3,[(0,),word(w),(0,)]) for w in choices]
    for i,a in enumerate(family):
        for j,b in enumerate(family):
            need(bool(maps(a,b))==(i==j),'distinct long words give distinct actual tournament classes')
            need(bool(maps(a,converse(b)))==(choices[i]==tuple(reversed(choices[j]))),'actual converse pairing matches long-word reversal')
    receipt={'status':'PASS','scope':__doc__,'actual_tournaments':len(RECORDS),'proposed_partitions_tested':proposed_count,'checks':dict(sorted(CHECKS.items())),'total_checks':sum(CHECKS.values()),'records':RECORDS,'seconds':round(time.time()-start,3),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(__file__).with_name('actual-tournament-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='records'},indent=2))
    print('NAMED HOSTILES:')
    for r in RECORDS:
        if not r['name'].startswith(('all-labelled','seeded')):print(json.dumps(r))
if __name__=='__main__':main()
