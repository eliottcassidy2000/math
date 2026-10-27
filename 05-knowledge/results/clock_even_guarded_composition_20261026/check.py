"""Exact, standalone controls for signed graph responses, clocks and tape stacks."""
from itertools import combinations, permutations, product
from collections import Counter
from fractions import Fraction as Q
from math import factorial
import json

def need(ok, message):
    if not ok: raise RuntimeError(message)

def matmul(a,b):
    return [[sum(x*y for x,y in zip(r,c)) for c in zip(*b)] for r in a]

def inverse(a):
    n=len(a); z=[[Q(v) for v in r]+[Q(i==j) for j in range(n)] for i,r in enumerate(a)]
    for j in range(n):
        k=next(i for i in range(j,n) if z[i][j]);z[j],z[k]=z[k],z[j]
        d=z[j][j];z[j]=[v/d for v in z[j]]
        for i in range(n):
            if i!=j:
                d=z[i][j];z[i]=[v-d*w for v,w in zip(z[i],z[j])]
    return [r[n:] for r in z]

def rank(rows):
    basis={}
    for r in rows:
        v=list(map(Q,r))
        for j,w in sorted(basis.items()):
            if v[j]:
                d=v[j];v=[a-d*b for a,b in zip(v,w)]
        if any(v):
            j=next(j for j,x in enumerate(v) if x);d=v[j];basis[j]=[x/d for x in v]
    return len(basis)

def fwht(v):
    a=list(v);s=1
    while s<len(a):
        for k in range(0,len(a),2*s):
            for j in range(s):
                x,y=a[k+j],a[k+j+s];a[k+j],a[k+j+s]=x+y,x-y
        s*=2
    return a

def parity(x):return -1 if x.bit_count()%2 else 1

class Atlas:
    def __init__(self,n,ports=()):
        self.n=n;self.edges=list(combinations(range(n),2));self.m=len(self.edges);self.N=1<<self.m
        self.index={e:i for i,e in enumerate(self.edges)}
        self.perms=[p for p in permutations(range(n)) if all(p[i]==i for i in ports)]
        self.maps=[]
        for p in self.perms:
            inds=[];rev=0
            for k,(i,j) in enumerate(self.edges):
                inds.append(self.index[tuple(sorted((p[i],p[j])))])
                if p[i]>p[j]:rev|=1<<k
            self.maps.append((inds,rev))
        self.tc={};self.orbits={}
        for t in range(self.N):
            if t in self.tc:continue
            orb={self.move(t,g,True) for g in range(len(self.perms))};c=min(orb)
            self.orbits[c]=orb
            for x in orb:self.tc[x]=c
        self.ts=sorted(self.orbits);self.even=[];seen=set()
        for f in range(self.N):
            if f in seen:continue
            signs={};odd=False
            for g in range(len(self.perms)):
                y=self.move(f,g);s=parity(f&self.maps[g][1])
                if y in signs and signs[y]!=s:odd=True
                signs[y]=s
            seen.update(signs)
            if not odd:self.even.append((f,signs))
        need(len(self.even)==len(self.ts),'twisted orbit dimension')
        self.M=[[sum(s*parity(t&f) for f,s in orb.items()) for _,orb in self.even] for t in self.ts]
        for a,(_,oa) in enumerate(self.even):
            for b in range(len(self.even)):
                dot=sum(len(self.orbits[t])*self.M[i][a]*self.M[i][b] for i,t in enumerate(self.ts))
                need(dot==(self.N*len(oa) if a==b else 0),'signed orbit Parseval')

    def move(self,x,g,tournament=False):
        inds,rev=self.maps[g];x=x^rev if tournament else x
        return sum(((x>>i)&1)<<j for i,j in enumerate(inds))

    def forward(self,t,i,j):
        k=self.index[tuple(sorted((i,j)))];return bool(((t>>k)&1)==(i>j))

    def path_data(self,p):
        mask=sum(1<<self.index[tuple(sorted(e))] for e in zip(p,p[1:]))
        negative=sum(1<<self.index[tuple(sorted((a,b)))] for a,b in zip(p,p[1:]) if a>b)
        return mask,negative

def tournament_controls():
    results=[]
    for n in range(2,7):
        a=Atlas(n);ti={t:i for i,t in enumerate(a.ts)};L=[[0]*len(a.ts) for _ in a.ts]
        for i,t in enumerate(a.ts):
            for p in permutations(range(n)):
                if all(a.forward(t,u,v) for u,v in zip(p,p[1:])):
                    mask,_=a.path_data(p);u=t^((a.N-1)^mask);L[i][ti[a.tc[u]]]+=1
        for i,t in enumerate(a.ts):
            for j,u in enumerate(a.ts):
                need(len(a.orbits[t])*L[i][j]==len(a.orbits[u])*L[j][i],'path-complement detailed balance')
        B=matmul(matmul(inverse(a.M),L),a.M)
        need(matmul(a.M,B)==matmul(L,a.M),'operator intertwiner')
        one=[1]*len(a.ts);space=[one];v=one
        for _ in a.ts:
            v=[sum(x*y for x,y in zip(r,v)) for r in L];space.append(v)
        conv=[ti[a.tc[t^(a.N-1)]] for t in a.ts]
        need(all(L[conv[i]][conv[j]]==L[i][j] for i in range(len(a.ts)) for j in range(len(a.ts))),'converse intertwining')
        row={'n':n,'tournament_types':len(a.ts),'even_graph_types':len(a.even),
             'converse_orbits':len({min(i,conv[i]) for i in range(len(a.ts))}),
             'negative_operator_entries':sum(x<0 for r in B for x in r),'counted_response_rank_from_constant':rank(space)}
        if n<=5:row.update({'graph_masks':[x for x,_ in a.even],'tournament_masks':a.ts,'counted_kernel':L,'walsh_basis':a.M,'signed_operator':[[str(x) for x in r] for r in B]})
        if n<=4:
            rooted=Atlas(n,(0,n-1));row['joint_port_types']=len(rooted.ts);row['relative_even_graph_types']=len(rooted.even)
        results.append(row)
    # Character-level path-guard formula by independent truth-table transforms.
    a=Atlas(4);cases=0
    for p in permutations(range(4)):
        mask,negative=a.path_data(p);flip=(a.N-1)^mask
        for f in range(a.N):
            values=[int(all(a.forward(t,u,v) for u,v in zip(p,p[1:])))*parity(f&(t^flip)) for t in range(a.N)]
            got=fwht(values);expected=[0]*a.N;sub=mask
            while True:
                expected[f^sub]=(a.N//8)*parity(f&flip)*parity(sub&negative)
                if sub==0:break
                sub=(sub-1)&mask
            need(got==expected,'guarded complement Fourier convolution');cases+=1
    # Ordered-join restriction of every raw character, both split sizes at n=5.
    joins=0;a=Atlas(5)
    for cut in [2,3]:
        left=Atlas(cut);right=Atlas(5-cut)
        for u in range(left.N):
            for v in range(right.N):
                t=0
                for k,(i,j) in enumerate(a.edges):
                    if j<cut:t|=((u>>left.index[i,j])&1)<<k
                    elif i>=cut:t|=((v>>right.index[i-cut,j-cut])&1)<<k
                for f in range(a.N):
                    fl=fr=0
                    for k,(i,j) in enumerate(a.edges):
                        if f>>k&1:
                            if j<cut:fl|=1<<left.index[i,j]
                            elif i>=cut:fr|=1<<right.index[i-cut,j-cut]
                    need(parity(f&t)==parity(fl&u)*parity(fr&v),'ordered join character restriction');joins+=1
    # Products of full labelled Hamilton-word guards, independent Boolean algebra.
    guard_ranks=[]
    for n in range(2,5):
        a=Atlas(n)
        guards=[sum(int(all(a.forward(t,u,v) for u,v in zip(p,p[1:])))<<t for t in range(a.N)) for p in permutations(range(n))]
        basis={};queue=[(1<<a.N)-1]
        while queue:
            v=queue.pop()
            while v:
                j=v.bit_length()-1
                if j in basis:v^=basis[j]
                else:
                    basis[j]=v;queue.extend(v&g for g in guards);break
        need(len(basis)==a.N,'named full-cube guard algebra')
        guard_ranks.append({'n':n,'labelled_states':a.N,'guard_algebra_rank_mod2':len(basis)})
    return {'atlases':results,'guard_character_truth_tables':cases,'ordered_join_character_checks':joins,'named_full_cube_guard_algebra':guard_ranks}

def clock_controls():
    affine=0;division=0;parities=[]
    for m in range(2,13):
        for a,b,c,d,r in product(range(m),repeat=5):
            need((a*(c*r+d)+b)%m==((a*c)*r+a*d+b)%m,'affine composition');affine+=1
        for a,b,q in product(range(1,5),range(4),range(2,5)):
            for r in range(q*m):
                need(((a*(r+q*m)+b)//q)%m==((a*r+b)//q)%m,'division precision lift');division+=1
    for k in range(1,11):
        words={}
        for n in range(1<<k):
            x=n;w=0
            for j in range(k):w|=(x%2)<<j;x=(3*x+1)//2 if x%2 else x//2
            need(w not in words,'shortcut parity coding bijection');words[w]=n
        parities.append({'k':k,'distinct_parity_words':len(words)})
    m=5
    need((1%m)==(6%m) and (1//2)%m!=(6//2)%m,'fixed-clock division hostile')
    fib=Counter();x,y=0,1
    while (x,y) not in fib:fib[x,y]+=1;x,y=y,(x+y)%5
    arrows={}
    for x,y in fib:arrows.setdefault(x,set()).add(y)
    need(len(fib)==20 and len(arrows[1])>1,'Fibonacci scalar arrow loses state')
    # Two-level refinement: forgetting a hidden bit is not compatible with a guard.
    even_guard=[int(n%2==0) for n in range(10)]
    need(any(even_guard[r]!=even_guard[r+5] for r in range(5)),'clock needs branch guard')
    return {'affine_composition_checks':affine,'division_lift_checks':division,'shortcut_parity_bijections':parities,
            'division_mod5_hostile':[1,6,0,3],'Fibonacci_mod5_period':len(fib),'Fibonacci_scalar_successors':{str(k):sorted(v) for k,v in arrows.items()}}

def stacks(tape,h):
    return (sum(1<<(h-i-1) for i in tape if i<h),sum(1<<(i-h) for i in tape if i>=h))

def stack_step(L,R,w,d):
    if d==1:return 2*L+w,R//2
    return L//2,(L%2)+2*w+4*(R//2)

def tape_controls():
    cases=0
    for word in range(1<<7):
        tape={i-3 for i in range(7) if word>>i&1}
        for h,w,d in product(range(-3,4),[0,1],[-1,1]):
            L,R=stacks(tape,h);out=set(tape)
            if w:out.add(h)
            else:out.discard(h)
            need(stack_step(L,R,w,d)==stacks(out,h+d),'two-stack/native tape equality');cases+=1
    finite=[]
    for k in range(1,11):
        # A right scanner outputs each bit, leaving the tape unchanged.
        words={}
        for r in range(1<<k):
            x=r;out=[]
            for _ in range(k):out.append(x%2);x//=2
            words[r]=tuple(out)
        need(len(set(words.values()))==1<<k,'scanner distinguishes all tape words')
        finite.append({'window_bits':k,'distinct_guard_traces':len(words)})
    return {'native_stack_checks':cases,'scanner':finite,'hidden_bit_hostile':{'k':6,'same_residue_mod64':[0,64],'different_read_at_step':6}}

def main():
    a=Atlas(3)
    scores={tuple(sum(a.forward(t,i,j) for j in range(3) if i!=j) for i in range(3)) for t in range(a.N)}
    need(len(scores)==7>4,'score is not the cut bit vector')
    a4=Atlas(4);cut=sum(1<<k for k,(i,j) in enumerate(a4.edges) if (i<2)!=(j<2))
    need(cut!=0 and all(sum(bool(cut>>k&1) for k,e in enumerate(a4.edges) if i in e)%2==0 for i in range(4)),'even-order cut/cycle intersection')
    print(json.dumps({'tournaments':tournament_controls(),'clocks':clock_controls(),'tapes':tape_controls(),
                      'historical_score_cut_hostile':{'n3_labelled_score_vectors':len(scores),'n3_cut_vectors':4,'n4_nonzero_cut_and_cycle_mask':cut}},sort_keys=True,indent=2))
if __name__=='__main__':main()
