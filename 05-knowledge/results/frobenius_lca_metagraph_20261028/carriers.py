"""Affine path-presentation operations, LCA query certificates and 81/27."""
import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
import json
from pathlib import Path
from fields import need, rank, span


def basis(values):
    result=[]
    for v in sorted(set(values)):
        if rank(result+[v])>len(result): result.append(v)
    return result


def affine_model(points):
    points=set(points)
    if not points: return None
    root=min(points); directions=basis(x^root for x in points)
    need({root^x for x in span(directions)}==points,'affine closure')
    return root,tuple(directions)


def all_affine(m):
    spaces={frozenset({0})}
    frontier=list(spaces)
    while frontier:
        s=frontier.pop()
        for x in range(1<<m):
            if x not in s:
                new=s|{v^x for v in s}
                if new not in spaces: spaces.add(new); frontier.append(new)
    cells=set()
    for s in spaces:
        for x in range(1<<m): cells.add(frozenset(v^x for v in s))
    return sorted(cells,key=lambda s:(len(s),tuple(sorted(s))))


def pairs(n): return list(combinations(range(n),2))
def free(n): return [(i,j) for i,j in pairs(n) if j>i+1]


def from_free(n,word):
    edges=pairs(n); ids={e:i for i,e in enumerate(free(n))}
    return sum((int(j==i+1) if j==i+1 else ((word>>ids[i,j])&1))<<k for k,(i,j) in enumerate(edges))


def edge(n,word,u,v):
    ids={p:i for i,p in enumerate(pairs(n))}
    return ((word>>ids[min(u,v),max(u,v)])&1) ^ int(u>v)


def relabel(n,word,p):
    return sum(edge(n,word,p[i],p[j])<<k for k,(i,j) in enumerate(pairs(n)))


def legal(n,word,p): return all(edge(n,word,p[i],p[i+1]) for i in range(n-1))


def join(n,x,m,y):
    return sum((edge(n,x,u,v) if v<n else edge(m,y,u-n,v-n) if u>=n else 1)<<k
               for k,(u,v) in enumerate(pairs(n+m)))


def affine_controls():
    summary=[]; families={}; operation_tests=0
    for n in range(1,5):
        cells=[frozenset(from_free(n,x) for x in c) for c in all_affine(len(free(n)))]
        families[n]=cells
        tau=sum(1<<k for k,(i,j) in enumerate(pairs(n)) if j>i+1)
        pp=list(permutations(range(n))); full=(1<<len(pairs(n)))-1
        for cell in cells:
            model=affine_model(cell)
            need(len(cell)==1<<len(model[1]),'state information dimension')
            affine_model(x^tau for x in cell); operation_tests+=1
            # Reflection-fixed blue cells are affine slices.
            affine_model(x for x in cell if relabel(n,x,tuple(reversed(range(n))))^full==x)
            operation_tests+=1
            for p in pp:
                image={relabel(n,x,p) for x in cell if legal(n,x,p)}
                affine_model(image)
                need(all(legal(n,y,tuple(range(n))) for y in image),'output path guard')
                operation_tests+=1
        summary.append(dict(n=n,affine_cells=len(cells),presentations=1<<len(free(n))))
    joins=0
    for n in range(1,5):
        for m in range(1,5):
            if n+m>6: continue
            tau=sum(1<<k for k,(i,j) in enumerate(pairs(n+m)) if j>i+1)
            for a in families[n]:
                for b in families[m]:
                    image={join(n,x,m,y) for x in a for y in b}
                    model=affine_model(image)
                    need(len(image)==len(a)*len(b),'join injectivity')
                    affine_model(x^tau for x in image)
                    joins+=1
    # A direct example of an unrecorded cross-child parity.
    code=[x for x in range(8) if x.bit_count()%2==0]
    left={x&3 for x in code}; right={(x>>1)&3 for x in code}
    invented=[x for x in range(8) if x&3 in left and (x>>1)&3 in right]
    need(len(code)==4 and len(invented)==8,'separator joint loss')
    # Once the parent parity is owned at the LCA, the third bit is reconstructed.
    need(all(((x>>2)&1)==((x&1)^((x>>1)&1)) for x in code),'LCA recovery')
    # Literal blue half-edge control inherited from the order-four packet.
    x,y=from_free(4,2),from_free(4,5)
    tau5=sum(1<<k for k,(i,j) in enumerate(pairs(5)) if j>i+1)
    glued=[join(1,0,4,z)^tau5 for z in (x,y)]
    hh=[hamilton(rows(5,z)) for z in glued]
    need(sorted(hh)==[11,13],'half-edge hostile')
    return dict(universes=summary,operation_tests=operation_tests,ordered_join_cells=joins,
                missing_separator_bits=1,blue_half_edge_hamilton_counts=hh)


def rows(n,word):
    return tuple(sum(edge(n,word,u,v)<<v for v in range(n) if u!=v) for u in range(n))


def hamilton(t):
    n=len(t); dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for mask in range(1,1<<n):
        for v,count in enumerate(dp[mask]):
            if count:
                available=t[v]&~mask
                while available:
                    bit=available&-available; available-=bit
                    dp[mask|bit][bit.bit_length()-1]+=count
    return sum(dp[-1])


def aut_count(t,keep_elements=False):
    n=len(t); degrees=[x.bit_count() for x in t]
    candidates=[[j for j in range(n) if degrees[j]==degrees[i]] for i in range(n)]
    order=sorted(range(n),key=lambda i:len(candidates[i])); assignment={}; total=0; elements=[]
    def rec(k,used):
        nonlocal total
        if k==n:
            total+=1
            if keep_elements: elements.append(tuple(assignment[i] for i in range(n)))
            return
        u=order[k]
        for v in candidates[u]:
            if used>>v&1: continue
            if all(((t[u]>>a)&1)==((t[v]>>b)&1) for a,b in assignment.items()):
                assignment[u]=v; rec(k+1,used|(1<<v)); del assignment[u]
    rec(0,0); return elements if keep_elements else total


def hostile():
    c3=(2,4,1)
    child=tuple(c3[i]|8 for i in range(3))+(0,)
    converse=tuple(((1<<4)-1)^child[i]^(1<<i) for i in range(4))
    def assemble(bits):
        t=[]
        for i in range(3):
            b=converse if bits>>i&1 else child
            for j in range(4):
                row=b[j]<<(4*i)
                row |= 15<<(4*((i+1)%3))
                t.append(row)
        return tuple(t)
    def triples(t):
        return tuple(c for c in combinations(range(12),3) if all(sum((t[u]>>v)&1 for v in c if v!=u)==1 for u in c))
    a,b=assemble(0),assemble(1)
    need(triples(a)==triples(b),'labelled C3 hostile')
    result=[]; groups={}
    for word in range(8):
        t=assemble(word); h=hamilton(t)
        if word<2:
            groups[word]=set(aut_count(t,True)); aut=len(groups[word])
        else: aut=aut_count(t)
        stabilizer=sum(all(((word>>i)&1)==((word>>((i+s)%3))&1) for i in range(3)) for s in range(3))
        need(aut==27*stabilizer and h==462915,'decorated LCA stabilizer')
        result.append(dict(decoration=word,H=h,Aut=aut,presentations=h//aut,
                           reset_diagonal=str(Fraction(aut,h)),cyclic_colour_stabilizer=stabilizer))
    need([result[i]['Aut'] for i in (0,1)]==[81,27],'81/27 control')
    G,H=groups[0],groups[1]
    need(H<G and all(all(g[v]!=v for v in range(12)) for g in G-H),'nonempty pointwise stabilizer theorem')
    need(all({g for g in G if g[v]==v}=={h for h in H if h[v]==v} for v in range(12)),
         'all named vertex stabilizers identical')
    transport=[sum(g[3]==7 for g in group) for group in (G,H)]
    need(transport==[27,0],'same local stabilizers lose inter-port transport')
    return dict(all_decorations=result,identical_labelled_cyclic_triples=len(triples(a)),
                transporter_3_to_7=transport,identical_nonempty_pointwise_stabilizers=True,
                all_elements_in_difference_fixed_point_free=len(G-H),
                LCA_invariant_constraints=[3,6],constraint_rank=rank([3,6]),
                invariant_decorations=[z for z in range(8) if all((z&r).bit_count()%2==0 for r in (3,6))])


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output'); args=parser.parse_args()
    result=dict(affine_compiler=affine_controls(),hostile=hostile())
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output: Path(args.output).write_text(text,encoding='utf-8',newline='\n')
    else: print(text,end='')


if __name__=='__main__': main()
