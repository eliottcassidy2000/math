#!/usr/bin/env python3
"""Literal rooted gate, continuant compiler, exact sidecar decoder, and controls.
Python 3 standard library only. Run python3 verify_gate.py > VERIFY.json.
"""
from itertools import product
from collections import defaultdict
import json

def require(condition, context):
    """Optimization-safe verification: checks remain active under python -O."""
    if not condition:
        raise AssertionError(context)

def tournament(n,code):
    a=[0]*n; k=0
    for j in range(n):
        for i in range(j):
            u,v=(j,i) if code>>k&1 else (i,j)
            a[u]|=1<<v; k+=1
    return a

def all_path_data(a,root):
    n=len(a); full=(1<<n)-1
    dp=[[0]*n for _ in range(1<<n)]
    rp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    rp[1<<root][root]=1
    for mask in range(1,1<<n):
        for v in range(n):
            c=dp[mask][v]; r=rp[mask][v]
            if not c and not r:continue
            todo=a[v]&~mask
            while todo:
                bit=todo&-todo;todo-=bit;u=bit.bit_length()-1
                dp[mask|bit][u]+=c;rp[mask|bit][u]+=r
    hs=[sum(row) for row in dp]; hs[0]=1
    R=sum(sum(rp[mask])*hs[full^mask] for mask in range(1,1<<n))
    return hs[-1],R,dp[-1][root]

def hp(a):
    n=len(a);dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for mask in range(1,1<<n):
        for v,c in enumerate(dp[mask]):
            if not c:continue
            todo=a[v]&~mask
            while todo:
                bit=todo&-todo;todo-=bit;dp[mask|bit][bit.bit_length()-1]+=c
    return sum(dp[-1])

def sink(a):
    v=len(a)
    return [row|(1<<v) for row in a]+[0]

def gate(a,root,m):
    require(m>=0, 'Verification check failed on source line 51')
    a=a[:]
    for _ in range(m+1):a=sink(a)
    b=len(a)-1;c=len(a)
    a=[row|(1<<c) for row in a]+[1<<root]
    a[root]^=1<<c
    return a,b

def jmove(a,root):
    v=len(a); b=[row|(1<<v) for row in a]+[1<<root]
    b[root]^=1<<v
    return b,v

def compile_program(ms):
    a=[0];root=0
    for m in ms:a,root=gate(a,root,m)
    return a,root

def recurrence(ms):
    old=0;h=1
    for m in ms:
        require(m>=1, 'Verification check failed on source line 72')
        q=1<<(m+1);old,h=h,(q+1)*h+q*old
    return old,h

def decode(old,h):
    """Return unique positive-m program for the paired continuation state, or None."""
    rev=[]
    while old:
        if old<0 or h<=old:return None
        v=(h-old)//old
        if v<4:return None
        q=1<<(v.bit_length()-1)
        numerator=h-(q+1)*old
        if numerator<0 or numerator%q:return None
        earlier=numerator//q
        if earlier>=old:return None
        rev.append(q.bit_length()-2)
        h,old=old,earlier
    return tuple(reversed(rev)) if h==1 else None

def decode_twoadic(h,R):
    """Independent inverse using v2(R_old)=1 after the first gate."""
    rev=[]
    while (h,R)!=(1,1):
        old=R-h; delta=2*h-R
        if old<1 or delta<=0:return None
        if old==1:
            q=delta
            if q<4 or q&(q-1):return None
            rev.append(q.bit_length()-2)
            h,R=1,1
            continue
        v2=(delta&-delta).bit_length()-1
        if v2<3:return None
        q=1<<(v2-1)
        RR=delta//q
        if not old<RR<2*old:return None
        rev.append(q.bit_length()-2)
        h,R=old,RR
    return tuple(reversed(rev))

def strong(a):
    n=len(a); full=(1<<n)-1
    def reach(rows):
        seen=1;todo=1
        while todo:
            bit=todo&-todo;todo-=bit;new=rows[bit.bit_length()-1]&~seen
            seen|=new;todo|=new
        return seen==full
    rev=[sum(1<<j for j in range(n) if a[j]>>i&1) for i in range(n)]
    return reach(a) and reach(rev)

def reverse_arcs(a):
    return [(j,i) for j in range(len(a)) for i in range(j) if a[j]>>i&1]

def main():
    checks=0; zero_m=0; root_presentations=0
    for n in range(1,5):
        for code in range(1<<(n*(n-1)//2)):
            a=tournament(n,code)
            for root in range(n):
                root_presentations+=1
                h,R,e=all_path_data(a,root)
                for m in range(4):
                    b,r=gate(a,root,m); hp_,Rp,_=all_path_data(b,r)
                    q=1<<(m+1)
                    expected=(h+q*R,2*h+q*R-(e if m==0 else 0))
                    require((hp_,Rp)==expected, (a,root,m,(hp_,Rp),expected))
                    # Physically equal to O^m J^2 after swapping the last two labels.
                    c=a[:];rr=root
                    for _ in range(m):c=sink(c)
                    c,rr=jmove(c,rr);c,rr=jmove(c,rr)
                    perm=list(range(len(c)));perm[-1],perm[-2]=perm[-2],perm[-1]
                    relabeled=[0]*len(c)
                    for u,row in enumerate(c):
                        for v in range(len(c)):
                            if row>>v&1:relabeled[perm[u]]|=1<<perm[v]
                    require(relabeled==b and perm[rr]==r, 'Verification check failed on source line 149')
                    checks+=1;zero_m+=int(m==0)
    program_checks=0;max_n=0
    for k in range(5):
        for ms in product(range(1,5),repeat=k):
            n=1+sum(m+2 for m in ms)
            if n>16:continue
            a,r=compile_program(ms);old,h=recurrence(ms)
            require(hp(a)==h, (ms,a,h))
            require(decode(old,h)==ms, 'Verification check failed on source line 158')
            require(decode_twoadic(h,h+old if ms else 1)==ms, 'Verification check failed on source line 159')
            require(h%8 in (1,5), 'Verification check failed on source line 160')
            require(strong(a), 'Verification check failed on source line 161')
            edges=reverse_arcs(a)
            require(len(edges)==k, 'Verification check failed on source line 163')
            require(len({v for e in edges for v in e})==2*k, 'Verification check failed on source line 164')
            require(all(sum(i<=cut<j for j,i in edges)<=2 for cut in range(n-1)), 'Verification check failed on source line 165')
            program_checks+=1;max_n=max(max_n,n)
    # All words with output at most cap. Exact monotone pruning.
    cap=1000000; seen={}; collisions=[]; nodes=0; bydepth=defaultdict(int)
    def visit(old,h,ms):
        nonlocal nodes
        nodes+=1;bydepth[len(ms)]+=1
        require(decode(old,h)==ms, (old,h,ms))
        require(decode_twoadic(h,h+old if ms else 1)==ms, 'Verification check failed on source line 173')
        R=h+old if ms else 1
        if h in seen:
            om,op=seen[h]
            if op!=old:
                collisions.append({'h':h,'first':list(om),'first_previous':op,
                                   'second':list(ms),'second_previous':old,
                                   'q4_first_continuation':5*h+4*op,
                                   'q4_second_continuation':5*h+4*old})
        else:seen[h]=(ms,old)
        q=4;m=1
        while h+q*R<=cap:
            hh=h+q*R
            visit(h,hh,ms+(m,));q*=2;m+=1
    visit(0,1,())
    require(min(c['h'] for c in collisions)==89, 'Verification check failed on source line 188')
    for ms in ((2,2),(3,1)):
        a,r=compile_program(ms);old,h=recurrence(ms)
        require(all_path_data(a,r)[:2]==(89,89+old), 'Verification check failed on source line 191')
        b,rr=gate(a,r,1)
        require(hp(b)==5*89+4*old, 'Verification check failed on source line 193')
    fixed=[]
    for k in range(9):fixed.append(recurrence((1,)*k)[1])
    print(json.dumps({'status':'all assertions passed','rooted_presentations':root_presentations,
        'literal_gate_checks':checks,'m_zero_corrected_controls':zero_m,
        'literal_program_checks':program_checks,'max_program_vertices':max_n,
        'program_census_cap':cap,'program_words_including_empty':nodes,
        'distinct_scalar_values':len(seen),'by_depth':dict(sorted(bydepth.items())),
        'scalar_collision_controls':sorted(collisions,key=lambda c:c['h'])[:5],'constant_m1_counts':fixed,
        'm0_false_reset_control':{'input':'singleton','actual_h_R':[3,3],'false_h_R':[3,4]},
        'noncommutation_control':{'m_word_1_2':recurrence((1,2))[1],
                                  'm_word_2_1':recurrence((2,1))[1]}},indent=2))
if __name__=='__main__':main()
