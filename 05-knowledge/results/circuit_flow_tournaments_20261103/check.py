"""Literal graph checks versus the exterior-word compiler and exact algebra."""
import hashlib
import importlib.util
import itertools as it
import json
from pathlib import Path
from math import factorial
import sys
import sympy as sp
from flow_operator import apply, compile_operator, need, symbolic

HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent
LO=(16,17,3,7,12); HI=(16,25,3,5,12)


def rows(n, mask):
    a=[0]*n
    for b,(i,j) in enumerate(it.combinations(range(n),2)):
        a[i if mask>>b&1 else j] |= 1 << (j if mask>>b&1 else i)
    return tuple(a)


def substitute(seed, module, root=0):
    outside=[i for i in range(len(seed)) if i!=root]
    n=len(module); a=list(module)+[0]*len(outside)
    for j,y in enumerate(outside,n):
        for x in range(n):
            if seed[root]>>y&1: a[x]|=1<<j
            else:a[j]|=1<<x
    for j,y in enumerate(outside,n):
        for k,z in enumerate(outside,n):
            if seed[y]>>z&1:a[j]|=1<<k
    return tuple(a)


def subset_hp(a):
    n=len(a); d=[[0]*n for _ in range(1<<n)]; h=[1]+[0]*((1<<n)-1)
    for v in range(n):d[1<<v][v]=1
    for mask in range(1,1<<n):
        rem=((1<<n)-1)^mask
        for v in range(n):
            if not d[mask][v]:continue
            dest=rem & a[v]
            while dest:
                b=dest&-dest;dest-=b;w=b.bit_length()-1
                d[mask|b][w]+=d[mask][v]
        h[mask]=sum(d[mask])
    return h


def slots(a, upto):
    h=subset_hp(a); full=len(h)-1; prev=[1]+[0]*full; out=[prev[full]]
    for _ in range(upto):
        cur=[0]*(full+1)
        for mask in range(full+1):
            sub=mask
            while True:
                cur[mask]+=prev[mask^sub]*h[sub]
                if not sub:break
                sub=(sub-1)&mask
        out.append(cur[full]);prev=cur
    return out


def literal_slots(a, r):
    # Ordered permutations and weak cuts: exactly one history per labelled slot cover.
    n=len(a); total=0
    if r==0:return int(n==0)
    for word in it.permutations(range(n)):
        for cuts in it.combinations_with_replacement(range(n+1),r-1):
            stops=(0,)+cuts+(n,)
            if all(all(a[word[i]]>>word[i+1]&1 for i in range(stops[j],stops[j+1]-1)) for j in range(r)):
                total+=1
    return total


def strong(a):
    n=len(a)
    for v in range(n):
        seen=1<<v
        while True:
            new=seen
            for u in range(n):
                if seen>>u&1:new |= a[u]
            if new==seen:break
            seen=new
        if seen!=(1<<n)-1:return False
    return True


def program(word, module=(0,)):
    a=tuple(module);names=['source:'+str(i) for i in range(len(a))];cuts=[]
    for level,letter in enumerate(word):
        old=len(a);a=substitute(LO if letter=='L' else HI,a)
        for v in range(old,len(a)):
            cuts.append(sum(1<<u for u in range(v) if a[u]>>v&1))
            names.append('gate:'+str(level)+':'+str(v-old+1))
    # Replay the actual birth cuts independently of substitution.
    b=list(module)
    for incoming in cuts:
        n=len(b);b.append(((1<<n)-1)^incoming)
        for i in range(n):
            if incoming>>i&1:b[i]|=1<<n
    need(tuple(b)==a,'actual source preserving births')
    return {'word_inner_to_outer':word,'rows':a,'names':names,'birth_cuts':cuts,
            'birth_clock':len(cuts),'native_complement_clock':0,
            'source_names_retained':names[:len(module)]}


def polynomial(op, p, r):
    return sp.expand(apply(op,r,lambda x:p.subs(r,x)))


def main():
    L=compile_operator(LO);U=compile_operator(HI)
    r=sp.Symbol('r'); checks=0; literal=0
    # Universal compiler on every labelled 4-vertex seed and every root;
    # modules: empty, singleton, transitive pair, both 3-vertex profile types.
    modules=[(),(0,),rows(2,1),rows(3,7),(2,4,1)]
    for mask in range(64):
        seed=rows(4,mask)
        for root in range(4):
            op=compile_operator(seed,root)
            for a in modules:
                source=slots(a,6)
                f=sp.interpolate([(i,z) for i,z in enumerate(source[:len(a)+1])],r)
                b=substitute(seed,a,root);actual=slots(b,3)
                for j in range(1,4):
                    need(apply(op,j,lambda x:f.subs(r,x))==actual[j],('compiler',mask,root,a,j))
                checks+=1
    # Gate full profile for all 75 nonempty labelled tournaments through order 4.
    gate_checks=0
    for n in range(1,5):
        for mask in range(1<<(n*(n-1)//2)):
            a=rows(n,mask);fv=slots(a,4)
            for seed,op in ((LO,L),(HI,U)):
                b=substitute(seed,a);actual=slots(b,3)
                for j in range(1,4):
                    need(apply(op,j,lambda x:fv[x])==actual[j],('gate',n,mask,j))
                need(strong(b),'strong actual wrapper');gate_checks+=1
    # A completely independent tiny permutation-and-cut check.
    for seed in (LO,HI):
        for a in ((0,),rows(2,1)):
            b=substitute(seed,a);fv=slots(b,2)
            for j in (1,2):
                need(literal_slots(b,j)==fv[j],'literal slot enumeration');literal+=1
    # Polynomial identities, never interpolation of the context operator.
    algebra=[]
    for m in range(1,13):
        f=r**m; lu=polynomial(L,polynomial(U,f,r),r);ul=polynomial(U,polynomial(L,f,r),r)
        need((lu-ul).subs(r,1)==-12*(3**m-3*2**m+3),'three-run commutator')
        G=compile_operator((2,4,1))
        jump=polynomial(G,polynomial(U,f,r)-polynomial(L,f,r),r).subs(r,1)
        need(jump==8*2**m-2*3**m-4,'contextual sign change')
        algebra.append({'m':m,'LU_H':int(lu.subs(r,1)),'UL_H':int(ul.subs(r,1)),
                        'commutator':int((lu-ul).subs(r,1)),'triangle_read_jump':int(jump)})
    # Actual two-gate sources, all old labels retained.
    sources=[]
    for module in ((0,),rows(2,1),rows(3,7),(2,4,1),rows(4,63)):
        p=program('UL',module);q=program('LU',module)
        h=subset_hp(p['rows'])[-1];g=subset_hp(q['rows'])[-1]
        f=slots(module,3);q3=f[3]-3*f[2]+3*f[1]
        need(h-g==-12*q3 and q3%6==0,'literal commutator')
        need(strong(p['rows']) and strong(q['rows']),'literal nested strong')
        sources.append({'module':module,'first':p,'second':q,'counts':[h,g],'pc3':q3//6})
    neutral=[]
    for word in ('LL','UL','LU','UU','ULL','LUL'):
        d=program(word);d['H']=subset_hp(d['rows'])[-1]
        d['first_slots']=slots(d['rows'],2) if len(d['rows'])<=9 else None
        neutral.append(d)
    need([z['H'] for z in neutral[:4]]==[589,611,611,637],'actual shared source square')
    need(neutral[4]['H']-neutral[5]['H']==2376,'actual one-context separator')
    # Triangular observation matrix gives a rational linear decoder.
    kernels=[[0,1]]
    for depth in range(1,13):
        prev=kernels[-1];nxt=[0]*(len(prev)+1)
        for j,c in enumerate(prev):
            if c:
                for (shift,degree),v in L.items():
                    from flow_operator import falling
                    if j+shift>=0:nxt[j+shift]+=c*v*falling(j,degree)
        top=factorial(depth+1)*factorial(2*depth+1)//2**depth
        need(nxt[-1]==top,'leading observer coefficient');kernels.append(nxt)
    for depth in range(8):
        M=sp.Matrix([[sum(c*j**m for j,c in enumerate(row)) for m in range(1,depth+2)] for row in kernels[:depth+1]])
        need(M.det()!=0,'actual transitive sources witness sharp rank')
    out={'status':'PROVED identities; FINITE_EXACT graph comparisons','operators':{'L':symbolic(L),'U':symbolic(U)},
         'universal_compiler_cases':checks,'gate_profile_cases':gate_checks,'literal_slot_cases':literal,
         'transitive_algebra':algebra,'commutator_sources':sources,'neutral_square_and_readers':neutral,
         'observer_kernels':kernels,'source_scope':'uniform module substitutions preserve old labels and arcs; arbitrary named-port continuations require full source'}
    (HERE/'certificate.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',encoding='utf8',newline='\n')
    print('PASS compiler',checks,'gate profiles',gate_checks,'literal slots',literal,'commutator sources',len(sources),'611 reader',neutral[4]['H'],neutral[5]['H'])


if __name__=='__main__':
    main()
