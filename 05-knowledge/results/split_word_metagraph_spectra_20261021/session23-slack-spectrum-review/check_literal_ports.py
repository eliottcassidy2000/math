#!/usr/bin/env python3
"""Independent exact audit using literal Hamilton paths in their own vertex coordinates.
No imports from candidate code. Path-coordinate tournament masks quotient the paths
without computing or assuming the automorphism group. All valid middle words of
vertex weight <=9, plus selected longer controls; padded cases have total order <=9.
"""
from collections import Counter
from itertools import product, permutations
from pathlib import Path
import argparse, hashlib, json, math

CHECKS=Counter()
def need(value,label):
    CHECKS[label]+=1
    if not value: raise RuntimeError(label)
def bits(x):
    while x:
        b=x&-x;x-=b;yield b.bit_length()-1

def words(cap):
    yield ()
    for s in (1,3):
        if cap>=s:
            for w in words(cap-s):yield (s,)+w

def word_graph(w):
    blocks=[];off=0
    for s in w:blocks.append(tuple(range(off,off+s)));off+=s
    a=[0]*off
    for k,B in enumerate(blocks):
        for v in B:
            a[v]=sum(1<<u for C in blocks[k+1:] for u in C)
            if len(B)==3:a[v]|=1<<B[(B.index(v)+1)%3]
    return tuple(a)
def join(graphs):
    n=sum(map(len,graphs));a=[];off=0
    for g in graphs:
        tail=((1<<n)-1)^((1<<(off+len(g)))-1)
        a.extend((r<<off)|tail for r in g);off+=len(g)
    return tuple(a)
def core(w):
    a=list(join([(0,),word_graph(w),(0,)]));n=len(a)
    a[0]^=1<<(n-1);a[-1]^=1
    return tuple(a)
def ham_paths(a):
    full=(1<<len(a))-1
    def rec(P,remain):
        if not remain:yield P;return
        for v in bits(remain&a[P[-1]]):yield from rec(P+(v,),remain^(1<<v))
    for v in range(len(a)):yield from rec((v,),full^(1<<v))
def path_mask(a,P):
    # All path edges are one and omitted. This mask is the COMPLETE source graph
    # in Hamilton-path coordinates, so equality is exactly an automorphism orbit.
    out=0;k=0
    for i in range(len(P)):
        for j in range(i+2,len(P)):
            out|=((a[P[i]]>>P[j])&1)<<k;k+=1
    return out
def mask_graph(mask,n):
    a=[0]*n;k=0
    for i in range(n):
        if i+1<n:a[i]|=1<<(i+1)
        for j in range(i+2,n):
            if mask>>k&1:a[i]|=1<<j
            else:a[j]|=1<<i
            k+=1
    return tuple(a)
def reflect(mask,n):
    pairs=[(i,j) for i in range(n) for j in range(i+2,n)]
    ix={p:k for k,p in enumerate(pairs)}
    return sum(((mask>>k)&1)<<ix[(n-1-j,n-1-i)] for k,(i,j) in enumerate(pairs))
def full_maps(a,anti=False):
    n=len(a);scores=[r.bit_count() for r in a]
    classes={d:tuple(i for i,s in enumerate(scores) if s==d) for d in set(scores)}
    target=[classes.get(n-1-s if anti else s,()) for s in scores]
    order=sorted(range(n),key=lambda v:len(target[v]));p=[None]*n;out=[]
    def rec(k,used):
        if k==n:out.append(tuple(p));return
        v=order[k]
        for u in target[v]:
            if used>>u&1:continue
            if all((((a[v]>>x)&1)^anti)==((a[u]>>p[x])&1) for x in order[:k]):
                p[v]=u;rec(k+1,used|1<<u);p[v]=None
    rec(0,0);return out

def direct_target(a,P):
    n=len(a)
    t=[((1<<n)-1)^r^(1<<i) for i,r in enumerate(a)]
    for u,v in zip(P,P[1:]):t[u]|=1<<v;t[v]&=~(1<<u)
    return tuple(sum(((t[P[i]]>>P[j])&1)<<j for j in range(n)) for i in range(n))

def graph_key(a,merged=False):
    # Canonical vertex ordering sorted by outdegree, then all permutations within
    # tied degrees. Independent of Hamilton-path mask encoding.
    def ordinary(g):
        groups={d:[] for d in sorted(set(r.bit_count() for r in g))}
        for v,r in enumerate(g):groups[r.bit_count()].append(v)
        best=None
        for ps in product(*(permutations(vs) for vs in groups.values())):
            p=tuple(v for q in ps for v in q)
            code=tuple((g[p[i]]>>p[j])&1 for i in range(len(g)) for j in range(i+1,len(g)))
            if best is None or code<best:best=code
        return best
    key=ordinary(a)
    if merged:
        n=len(a);op=tuple(((1<<n)-1)^r^(1<<i) for i,r in enumerate(a))
        key=min(key,ordinary(op))
    return key

def split_coordinate(P,w):
    n=len(P)
    if P.index(n-1)+1!=P.index(0):return None
    left=set(P[:P.index(n-1)]);off=1;js=[]
    for s in w:js.append(sum(v in left for v in range(off,off+s)));off+=s
    return tuple(js)

def analyze(w,targets=False):
    a=core(w);n=len(a);t=w.count(3);e=w.count(1)+2*t
    masks={};coords={};H=0
    for P in ham_paths(a):
        H+=1;m=path_mask(a,P);coord=split_coordinate(P,w)
        if m in masks:
            need(coords[m]==coord,'path-coordinate orbit has one split coordinate')
        else:masks[m]=P;coords[m]=coord
    autos=full_maps(a);anti=full_maps(a,True)
    need(len(autos)==3**t,'full actual automorphism order')
    need(H==len(masks)*len(autos),'literal free-orbit quotient')
    need(H==3**t*(1+2**e),'literal raw H formula')
    need(len(masks)==1+2**e,'literal path-mask h formula')
    need(Counter(coords.values())==Counter([None]+list(product(*(range(s+1) for s in w)))),'complete split-mask bijection')
    is_pal=w==w[::-1]
    need(bool(anti)==is_pal,'actual anti-isomorphism palindrome criterion')
    b=sum(reflect(m,n)==m for m in masks)
    need(b==(0 if not is_pal else 1 if n%2 else 1+2**(e//2)),'blue masks equal parity formula')
    if is_pal:
        need(all(reflect(m,n) in masks for m in masks),'reflection stays in actual path fibre')
        for m,j in coords.items():
            reflected=coords[reflect(m,n)]
            need(reflected==(None if j is None else tuple(s-k for s,k in zip(w,j[::-1]))),'literal reflection is complement-reversal')
        if n%2==0:
            r=len(w)//2;X=list(product(*(range(s+1) for s in w[:r])));pairs={};fixed=set()
            for m,j in coords.items():
                if j is None:continue
                x=j[:r];y=tuple(s-k for s,k in zip(w[:r],j[r:][::-1]))
                jm=coords[reflect(m,n)]
                xx=jm[:r];yy=tuple(s-k for s,k in zip(w[:r],jm[r:][::-1]))
                need((xx,yy)==(y,x),'X times X coordinates convert reflection to swap')
                if x==y:fixed.add(x)
                else:pairs[frozenset((x,y))]=frozenset((m,reflect(m,n)))
            need(len(fixed)==len(X),'blue diagonal X complete')
            need(len(pairs)==math.comb(len(X),2),'black local ports are two-subsets X')
            need(len(set(pairs.values()))==len(pairs),'two-subset indexing is injective at local incidences')
    row={'word':list(w),'n':n,'e':e,'H':H,'aut':len(autos),'h':len(masks),'b':b,'self_converse':is_pal}
    if targets:
        full=(1<<((n-1)*(n-2)//2))-1
        source=graph_key(a,True);ports={};dest_counts=Counter();loop_ports=0;port_receipts=[]
        for m in masks:
            if reflect(m,n)==m:continue
            pair=min(m,reflect(m,n))
            if pair in ports:continue
            need(reflect(m^full,n)==(reflect(m,n)^full),'nonpath complement commutes with reflection')
            target=mask_graph(m^full,n);target_key=graph_key(target,True)
            need(graph_key(mask_graph(reflect(m,n)^full,n),True)==target_key,'twinned local ports have same actual merged destination')
            need(path_mask(target,tuple(range(n)))==(m^full),'target adjacency is exact nonpath complement')
            need(direct_target(a,masks[m])==target,'actual-vertex nonpath flip agrees with path-coordinate target')
            ports[pair]=target_key;dest_counts[target_key]+=1;loop_ports+=target_key==source
            port_receipts.append({'split':coords[m],'twin_split':coords[reflect(m,n)],'literal_path':masks[m],
                'source_path_mask':m,'target_path_mask':m^full,'target_adjacency_rows':target,
                'target_merged_canonical_key':target_key})
        row['actual_black_incidence_ports']=len(ports)
        row['distinct_actual_merged_destinations']=len(dest_counts)
        row['destination_incidence_multiplicities']=sorted(dest_counts.values())
        row['black_loop_incidences']=loop_ports
        row['source_adjacency_rows']=a
        row['literal_port_receipts']=port_receipts
    return row

def fixed_coefficients():
    rows=[]
    for e in range(2,21):
        allc=Counter();even=Counter();odd=Counter()
        def ew(rem,p=()):
            if rem==0:yield p
            if rem>=1:yield from ew(rem-1,p+(1,))
            if rem>=2:yield from ew(rem-2,p+(3,))
        for w in ew(e):
            n=2+sum(w);allc[n]+=1
            if w==w[::-1]:(odd if n%2 else even)[n]+=1
        expected={e+t+2:math.comb(e-t,t) for t in range(e//2+1)}
        need(dict(allc)==expected,'fixed-e binomial polynomial literal-word enumeration')
        P=lambda r:Counter({r+t:math.comb(r-t,t) for t in range(r//2+1)})
        ev=Counter();od=Counter()
        if e%2==0:
            r=e//2;ev=Counter({2+2*k:v for k,v in P(r).items()});od=Counter({5+2*k:v for k,v in P(r-1).items()})
        else:
            r=(e-1)//2;od=Counter({3+2*k:v for k,v in P(r).items()})
        need(even==ev,'fixed-e even-palindrome coefficient formula')
        need(odd==od,'fixed-e odd-palindrome coefficient formula')
        rows.append({'e':e,'ordinary':dict(sorted(allc.items())),'SC_even':dict(sorted(even.items())),'SC_odd':dict(sorted(odd.items()))})
    return rows

def padding_controls():
    result=[]
    for w in [(1,1),(3,),(1,3),(3,1),(1,3,1)]:
        for left,right in product([(),(1,),(3,),(1,1),(1,3),(3,1)],repeat=2):
            if sum(left)+sum(w)+sum(right)+2>9:continue
            a=join([word_graph(left),core(w),word_graph(right)])
            masks={path_mask(a,P) for P in ham_paths(a)};e=w.count(1)+2*w.count(3)
            need(len(masks)==1+2**e,'arbitrary neutral padding preserves h')
            actual=bool(full_maps(a,True));pal=w==w[::-1] and right==left[::-1]
            need(actual==pal,'padded actual self-converse iff reversed tails and palindromic core')
            b=sum(reflect(m,len(a))==m for m in masks)
            need(b==(0 if not pal else 1 if (sum(w)+2)%2 else 1+2**(e//2)),'mirrored padding preserves b including parity channels')
            result.append({'left':left,'core':w,'right':right,'n':len(a),'h':len(masks),'b':b})
    return result

def degree_series_controls(N=14):
    # Enumerate complete classified (left,core,right) signatures and quotient
    # their converse involution directly; then compare coefficientwise with
    # the rational degree channels, using integer series convolution only.
    counted=Counter();seen=set();classes=0
    for w in words(N-2):
        if sum(w)<2:continue
        for left in words(N-2-sum(w)):
            for right in words(N-2-sum(w)-sum(left)):
                classes+=1
                sig=(left,w,right);rev=(right[::-1],w[::-1],left[::-1]);rep=min(sig,rev)
                if rep in seen:continue
                seen.add(rep);n=2+sum(w)+sum(left)+sum(right);e=w.count(1)+2*w.count(3)
                if sig!=rev:k=1+2**e
                elif (2+sum(w))%2:k=2**(e-1)
                else:k=math.comb(2**(e//2),2)
                counted[k,n]+=1
    neutral=[0]*(N+1);neutral[0]=1
    for n in range(1,N+1):neutral[n]=(neutral[n-1] if n>=1 else 0)+(neutral[n-3] if n>=3 else 0)
    p2=[sum(neutral[i]*neutral[n-i] for i in range(n+1)) for n in range(N+1)]
    q=[neutral[n//2] if n%2==0 else 0 for n in range(N+1)]
    expected=Counter()
    for e in range(2,N-1):
        P=Counter({e+t+2:math.comb(e-t,t) for t in range(e//2+1)});E=Counter();O=Counter()
        for w in words(N-2):
            if w.count(1)+2*w.count(3)==e and w==w[::-1]:
                ((O if (2+sum(w))%2 else E))[2+sum(w)]+=1
        for n in range(N+1):
            F=sum(c*p2[n-i] for i,c in P.items() if i<=n)
            Ev=sum(c*q[n-i] for i,c in E.items() if i<=n)
            Od=sum(c*q[n-i] for i,c in O.items() if i<=n)
            need((F-Ev-Od)%2==0,'ordinary minus self-converse coefficients are even')
            expected[2**(e-1),n]+=Od
            if e%2==0:expected[math.comb(2**(e//2),2),n]+=Ev
            expected[1+2**e,n]+=(F-Ev-Od)//2
    expected=+expected
    need(expected==counted,'complete degree-channel rational coefficients match direct converse quotient')
    return {'maximum_order':N,'ordinary_signatures':classes,'merged_signatures':len(seen),
        'nonzero_degree_order_coefficients':len(counted),
        'coefficients':[{'degree':k,'n':n,'count':c} for (k,n),c in sorted(counted.items())]}

def matching_bridge_controls(R=10):
    def exponent(w):return sum(1 if s==1 else 2 for s in w)
    def ew(r):
        if r==0:yield ()
        if r>=1:
            for w in ew(r-1):yield (1,)+w
        if r>=2:
            for w in ew(r-2):yield (3,)+w
    wordsets={r:tuple(ew(r)) for r in range(2*R+2)}
    def cut(w,r):
        acc=0
        for j,s in enumerate(w):
            if acc==r:return w[:j],w[j:]
            acc+=1 if s==1 else 2
            if acc>r:return None
        return (w,()) if acc==r else None
    def decode(A,B):
        e=exponent(A)
        while e and B[0]==3:A,B=B[1:],A;e-=1
        need(B and B[0]==1,'tiling decoder stops at a monomer or unit')
        if e==0:
            need(A==() and B==(1,),'actual empty unit at decoder boundary')
            return None
        return e,A+tuple(reversed(B[1:]))
    def encode(obj,r):
        if obj is None:e=0;A=();B=(1,)
        else:
            e,w=obj;left,right=cut(w,e);A=left;B=(1,)+right[::-1]
        for _ in range(r-e):A,B=B,(3,)+A
        return A,B
    rows=[]
    for r in range(R+1):
        image=set()
        for A,B in product(wordsets[r],wordsets[r+1]):
            obj=decode(A,B)
            need(encode(obj,r)==(A,B),'matching bridge decoder and inverse recover literal pair')
            image.add(obj)
        expected={None}
        for e in range(1,r+1):
            expected.update((e,U+V[::-1]) for U,V in product(wordsets[e],repeat=2))
        need(image==expected,'matching bridge image is exactly cumulative balanced-cut cores plus unit')
        need(len(image)==len(wordsets[r])*len(wordsets[r+1]),'matching bridge injective and onto')
        if r:
            lhs=Counter(map(lambda w:2+sum(w),wordsets[2*r]));rhs=Counter();balanced=0;straddling=0
            for w in wordsets[2*r]:
                split=cut(w,r)
                if split is not None:balanced+=1
                else:
                    acc=0
                    for j,s in enumerate(w):
                        nxt=acc+(1 if s==1 else 2)
                        if acc<r<nxt:
                            need(s==3 and acc==r-1,'unique straddler is a triangle at half exponent')
                            need(exponent(w[j+1:])==r-1,'straddler leaves equal lower-weight sides')
                            straddling+=1;break
                        acc=nxt
            for U,V in product(wordsets[r],repeat=2):rhs[2+sum(U)+sum(V)]+=1
            for U,V in product(wordsets[r-1],repeat=2):rhs[5+sum(U)+sum(V)]+=1
            need(lhs==rhs,'weighted fixed-even-exponent polynomial splits into two Fibonacci squares')
            need(balanced==len(wordsets[r])**2 and straddling==len(wordsets[r-1])**2,'balanced and straddling channels have correct cardinalities')
        else:balanced=straddling=0
        rows.append({'r':r,'unequal_pair_count':len(image),'balanced_layer':balanced,'straddling_layer':straddling})
    need([r['unequal_pair_count'] for r in rows[:7]]==[1,2,6,15,40,104,273],'inherited Fibonacci-product sequence with exact index')
    return {'max_half_exponent':R,'rows':rows}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--receipt',default='normal.json');args=parser.parse_args()
    rows=[]
    for w in words(9):
        if sum(w)<2:continue
        rows.append(analyze(w,targets=w in [(1,1),(1,1,1,1),(3,3),(1,3,3,1),(3,1,1,3)]))
    for w in [(1,3,3,1),(3,1,1,3),(3,3,3,3)]:
        if sum(w)>9:rows.append(analyze(w))
    gf=fixed_coefficients();padding=padding_controls();degree=degree_series_controls();matching=matching_bridge_controls()
    # Explicit exclusions/hostiles: at W=1 the endpoints cease to be intrinsic.
    boundary=core((1,));bp=list(ham_paths(boundary));bm={path_mask(boundary,p) for p in bp}
    need((len(bp),len(full_maps(boundary)),len(bm))==(3,3,1),'excluded W=1 is triangle not h=3 delta-three atom')
    need(len(bm)!=1+2**1,'unrestricted formula fails at excluded singleton boundary')
    # Same h need not imply same b or order: 11 versus 3 have e=2.
    lookup={tuple(r['word']):r for r in rows}
    need(lookup[(1,1)]['h']==lookup[(3,)]['h']==5 and lookup[(1,1)]['b']==3 and lookup[(3,)]['b']==1,'same h different blue counts hostile')
    need(lookup[(1,3)]['h']==lookup[(3,1)]['h']==9 and lookup[(1,3)]['b']==0,'same e does not ensure palindrome or nonzero b')
    same_order=[lookup[w] for w in [(1,3,3,1),(3,1,1,3)]]
    need(all(tuple(r[k] for k in ('n','H','aut','h','b'))==(10,585,9,65,9) for r in same_order),'same-order wiring hostile has identical full scalar signature')
    targetsets=[{tuple(p['target_merged_canonical_key']) for p in r['literal_port_receipts']} for r in same_order]
    need(all(len(s)==28 for s in targetsets),'same-order wiring hostile has 28 distinct destinations on each side')
    need(not targetsets[0].intersection(targetsets[1]),'same-order equal-scalar cores have disjoint actual merged destination sets')
    need(lookup[(1,1,1,1)]['destination_incidence_multiplicities']==[1,1,1,1,2] and lookup[(3,3)]['destination_incidence_multiplicities']==[1]*6,'six-port wiring hostile has unequal neighbor multiplicity profiles')
    data={'status':'PASS','method':__doc__,'literal_cores':len(rows),'literal_hamilton_paths':sum(r['H'] for r in rows),'literal_path_orbits':sum(r['h'] for r in rows),'maximum_core_order':max(r['n'] for r in rows),'padding_cases':len(padding),'fixed_e_range':[2,20],'degree_series_summary':{k:v for k,v in degree.items() if k!='coefficients'},'matching_bridge':matching,'checks':dict(sorted(CHECKS.items())),'total_checks':sum(CHECKS.values()),'rows':rows,'fixed_e':gf,'padding':padding,'degree_series':degree,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(args.receipt).write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({k:v for k,v in data.items() if k not in ('rows','fixed_e','padding','degree_series')},indent=2))
    print('LOCAL INCIDENCE TARGET CONTROLS:')
    for r in rows:
        if 'actual_black_incidence_ports' in r:print(json.dumps({k:v for k,v in r.items() if k!='literal_port_receipts'}))
if __name__=='__main__':main()
