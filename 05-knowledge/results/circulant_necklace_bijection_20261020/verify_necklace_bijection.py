#!/usr/bin/env python3
"""Exact binary-necklace path bijection and prime circulant controls.
No dependencies; no randomization; no catalogue ranking or Hall matching.
"""
from itertools import product
from collections import Counter
import argparse
import json

CHECKS=0
def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(message)

def canonical(word):
    return min(word[i:]+word[:i] for i in range(len(word)))

def period(word):
    n=len(word)
    return next(d for d in range(1,n+1) if n%d==0 and word==word[:d]*(n//d))

def oddpart(n):
    while n%2==0:n//=2
    return n

def dyadicpart(n):return n//oddpart(n)

def even_root(word):return sum(word[:period(word)])%2==0

def derivative(word):return tuple(word[i]^word[(i+1)%len(word)] for i in range(len(word)))

def stable(word):
    e=dyadicpart(len(word))
    return all(sum(word[r::e])%2==0 for r in range(e))

def preimages(word):
    """All predecessor necklaces, represented by least rotations."""
    if sum(word)%2:return ()
    ans=[0]
    for bit in word[:-1]:ans.append(ans[-1]^bit)
    ans=tuple(ans)
    return tuple(sorted({canonical(ans),canonical(tuple(1-bit for bit in ans))}))

def cut_predecessor(word, rule='stable'):
    pred=preimages(word)
    if len(pred)!=2:return None
    if rule=='stable' and stable(word):
        return next(w for w in pred if not stable(w))
    return pred[-1]

def forward_word(word, rule='stable'):
    """Odd-total word -> even-root word. Retains the input rotation phase."""
    if sum(word)%2!=1:raise ValueError('The source must have odd total weight.')
    bound=dyadicpart(len(word)) if rule=='stable' else 2**len(word)
    for steps in range(1,bound+1):
        target=derivative(word)
        if canonical(word)==cut_predecessor(target,rule):return target,steps
        word=target
    raise AssertionError('The proved termination bound failed.')

def inverse_necklace(word, rule='stable'):
    """Even-root word -> canonical odd-total source necklace."""
    if not even_root(word):raise ValueError('The target must have even primitive-root weight.')
    current=cut_predecessor(word,rule)
    bound=dyadicpart(len(word)) if rule=='stable' else 2**len(word)
    for steps in range(1,bound+1):
        pred=preimages(current)
        if not pred:return current,steps
        cut=cut_predecessor(current,rule)
        kept=[w for w in pred if w!=cut]
        if len(kept)!=1:raise AssertionError('Invalid retained predecessor count.')
        current=kept[0]
    raise AssertionError('The proved inverse termination bound failed.')

def tournament_word(bits):return bits+tuple(1-bit for bit in bits)

def tournament_forward(full):
    m=len(full)//2
    if len(full)%2 or full[m:]!=tuple(1-bit for bit in full[:m]):
        raise ValueError('The tournament word must be antiperiodic.')
    source=derivative(full)[:m]
    return forward_word(source)

def tournament_inverse(graph_word):
    source,steps=inverse_necklace(graph_word)
    first=[0]
    for bit in source[:-1]:first.append(first[-1]^bit)
    return canonical(tournament_word(tuple(first))),steps

def primes_to(n):
    return [p for p in range(3,n+1,2) if all(p%d for d in range(2,int(p**.5)+1))]

def primitive_root(p):
    return next(g for g in range(2,p) if len({pow(g,k,p) for k in range(1,p)})==p-1)

def graph_mask(p,g,word):
    S={pow(g,i,p) for i in range(p-1) if word[i%len(word)]}
    return tuple(tuple(int((j-i)%p in S) if i!=j else 0 for j in range(p)) for i in range(p))

def graph_even_by_automorphisms(adj):
    """Independent automorphism enumeration, no multiplier assumption.
    Empty and complete graphs have direct certificates; otherwise backtrack.
    """
    n=len(adj); edges=[(i,j) for i in range(n) for j in range(i+1,n) if adj[i][j]]
    if not edges:return True,1
    if len(edges)==n*(n-1)//2:return False,1
    image=[];used=set();count=0;odd=False
    def rec():
        nonlocal count,odd
        i=len(image)
        if i==n:
            count+=1
            if sum(image[u]>image[v] for u,v in edges)%2:odd=True
            return
        for j in range(n):
            if j not in used and all(adj[i][a]==adj[j][b] for a,b in enumerate(image)):
                image.append(j);used.add(j);rec();used.remove(j);image.pop()
                if odd:return
    rec()
    return not odd,count

def run(max_length,prime_max,automorphism_max):
    report={'necklaces':[],'prime_circulants':[],'hostiles':{}}
    for n in range(1,max_length+1):
        necklaces={canonical(w) for w in product((0,1),repeat=n)}
        incoming={w:[] for w in necklaces}
        for w in necklaces:incoming[canonical(derivative(w))].append(w)
        source={w for w in necklaces if sum(w)%2}
        target={w for w in necklaces if even_root(w)}
        images={};hist=Counter();lex_hist=Counter()
        for w in necklaces:
            pred=preimages(w)
            check(set(pred)==set(incoming[w]),('independent indegree',n,w))
            check(len(pred)==(2 if even_root(w) else 0 if sum(w)%2 else 1),('indegree trichotomy',n,w))
            d=period(w);e=period(derivative(w))
            check(d in (e,2*e),('period loss',n,w))
            z=w
            for _ in range(dyadicpart(n)):z=derivative(z)
            check(stable(z),('stable entry',n,w))
            if stable(w):
                check(len(pred)==2 and sum(stable(v) for v in pred)==1,('stable predecessor',n,w))
        for w in source:
            out,steps=forward_word(w);out=canonical(out)
            check(out in target and out not in images,('forward injection',n,w,out))
            back,backsteps=inverse_necklace(out)
            check((back,backsteps)==(w,steps),('inverse',n,w,out,back))
            check(oddpart(period(w))==oddpart(period(out)),('odd period',n,w,out))
            check(steps<=dyadicpart(n),('step bound',n,w,steps))
            for i in range(n):
                shifted=w[i:]+w[:i]
                phase,sk=forward_word(shifted)
                original,k=forward_word(w)
                check((phase,sk)==(original[i:]+original[:i],k),('rotation equivariance',n,w,i))
            lexout,lexsteps=forward_word(w,'lex')
            check(inverse_necklace(lexout,'lex')==(w,lexsteps),('generic lex rule',n,w))
            images[out]=w;hist[steps]+=1;lex_hist[lexsteps]+=1
        check(set(images)==target,('surjectivity',n))
        report['necklaces'].append({'m':n,'necklaces':len(necklaces),'sources':len(source),'targets':len(target),'max_steps':max(hist),'step_histogram':dict(sorted(hist.items())),'generic_lex_max_steps':max(lex_hist)})
    for p in primes_to(prime_max):
        m=(p-1)//2;g=primitive_root(p)
        tourn_words=[tournament_word(w) for w in product((0,1),repeat=m)]
        tourn_types={canonical(w) for w in tourn_words}
        even_words={w for w in product((0,1),repeat=m) if even_root(w)}
        even_types={canonical(w) for w in even_words}
        type_map={};label_images=set();fibers=Counter()
        for word in tourn_words:
            out,k=tournament_forward(word);label_images.add(out);fibers[out]+=1
            a,b=canonical(word),canonical(out)
            check(a not in type_map or type_map[a]==b,('well-defined prime map',p,a))
            type_map[a]=b
            inverse,ik=tournament_inverse(out)
            check((inverse,ik)==(a,k),('prime inverse',p,word,out))
            # Every multiplier is a cyclic shift; translations fix both objects.
            for i in range(2*m):
                out2,k2=tournament_forward(word[i:]+word[:i])
                j=i%m
                check((out2,k2)==(out[j:]+out[:j],k),('AGL equivariance',p,word,i))
        check(set(type_map)==tourn_types and set(type_map.values())==even_types and len(type_map)==len(even_types),('prime type bijection',p))
        check(label_images==even_words,('labelled surjection',p))
        input_period_by_image={image:period(source) for source,image in type_map.items()}
        for w,sz in fibers.items():
            d=input_period_by_image[canonical(w)]
            check(sz==d//period(w) and sz&(sz-1)==0,('dyadic labelled fibers',p,w,sz,d))
        aut_checked=0;aut_perms=0
        if p<=automorphism_max:
            for w in product((0,1),repeat=m):
                exact,count=graph_even_by_automorphisms(graph_mask(p,g,w))
                check(exact==even_root(w),('independent graph evenness',p,w))
                aut_checked+=1;aut_perms+=count
        report['prime_circulants'].append({'p':p,'primitive_root':g,'tournament_labelled':len(tourn_words),'even_graph_labelled':len(even_words),'types':len(tourn_types),'fibers':dict(sorted(Counter(fibers.values()).items())),'graph_words_independently_checked':aut_checked,'automorphisms_examined':aut_perms})
    # Fixed-rule failures: D alone never lands in target; squaring fails at p=5.
    f=tournament_word((0,0));u=derivative(f)[:2];d2=derivative(u)
    check(not even_root(u) and not even_root(d2),('fixed derivative hostile',f,u,d2))
    report['hostiles']['p5_fixed_square']={'tournament':f,'D':u,'D2':d2,'D3':derivative(d2)}
    # Odd root repeated twice has even total weight but is not an even graph.
    check(sum((0,1,0,1))%2==0 and not even_root((0,1,0,1)),'primitive parity hostile')
    report['hostiles']['primitive_parity']='0101 has even total weight and odd primitive-root weight.'
    report['checks']=CHECKS
    return report

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--max-length',type=int,default=16);ap.add_argument('--prime-max',type=int,default=31);ap.add_argument('--automorphism-max',type=int,default=13);ap.add_argument('--output',default=None)
    args=ap.parse_args()
    result=run(args.max_length,args.prime_max,args.automorphism_max)
    text=json.dumps(result,indent=2,sort_keys=True)
    if args.output:
        with open(args.output,'w') as f:f.write(text+'\n')
    print(text)
