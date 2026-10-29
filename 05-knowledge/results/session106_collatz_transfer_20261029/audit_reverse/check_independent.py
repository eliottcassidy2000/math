#!/usr/bin/env python3
"""Independent native dictionary and literal atlas audit; no producer imports."""
import json
from pathlib import Path

def need(p,m):
    if not p:raise RuntimeError(m)

def step(n):
    t=3*n+1;a=(t&-t).bit_length()-1
    return t>>a,a
WORDS=((1,),(3,1,1,1),(3,3,1,1,1,1,1),(3,1,3,1,1,1,1),(3,1,1,3,1,1,1))
def affine(w):
    a,b,c=1,1,0
    for v in w:a,b,c=3*a,b*2**v,3*c+b
    return a,b,c

def inverse(y):
    ans=[]
    for i,w in enumerate(WORDS):
        a,b,c=affine(w)
        numerator=b*(2*y-1)-c
        if numerator<=0 or numerator%a:continue
        n=numerator//a
        if n%2==0:continue
        x=(n+1)//2
        if x%81:continue
        z=n;got=[]
        for _ in w:z,v=step(z);got.append(v)
        if tuple(got)==w and z==2*y-1:ans.append((x,i))
    return sorted(ans)

def main():
    # Arbitrary guarded macro lift: compare composed native words from a nonroot cylinder.
    prefix=(1,3); aa,bb,cc=affine(prefix)
    M=2*bb;r=((bb-cc)*pow(aa,-1,M))%M;t=(aa*r+cc)//bb
    macro_checks=0
    for w in WORDS:
        alpha,P,c=affine(w)
        g=((P-c)*pow(alpha,-1,2*P))%(2*P)
        s=((g-t)//2*pow(aa,-1,P))%P
        rr,MM,AA,tt=r+M*s,M*P,alpha*aa,(alpha*(t+2*aa*s)+c)//P
        for q in range(12):
            source=rr+MM*q; y=source
            for want in prefix+w:
                y,actual=step(y);need(actual==want,'macro native guard')
            need(y==tt+2*AA*q,'macro joint source endpoint')
            # Mixed source/endpoint inequalities reduce without changing truth.
            need((2*source-y>=7)==((2*MM-2*AA)*q>=7-2*rr+tt),'mixed affine guard')
            need((y+source<=100000)==((2*AA+MM)*q<=100000-tt-rr),'upper affine guard')
            macro_checks+=1
    limit=300000; incoming={}; edges=0
    for x in range(81,limit+1,81):
        n=2*x-1;q=1;w=[];z=n
        for _ in range(10000):
            z,a=step(z);w.append(a);q+=2 if a%3==0 else -1
            if q==0:break
        need(q==0,'bounded local return audit')
        if z>n:
            need(tuple(w) in WORDS,'five-growing language')
            i=WORDS.index(tuple(w));y=(z+1)//2
            if y<=limit:incoming.setdefault(y,[]).append((x,i));edges+=1
    for y in range(81,limit+1,81):
        need(inverse(y)==sorted(incoming.get(y,[])),'reverse vs native source dictionary')
    nodes=0
    for y in range(81,limit+1,81):
        row=[(y,0,0)]
        for k in range(1,16):
            row=[(x,j+len(WORDS[i]),d+sum(WORDS[i])) for t,j,d in row for x,i in inverse(t)]
            need(len(row)<=3,'reverse width')
            need(len({d-j for _,j,d in row})==len(row),'clock uniqueness')
            need(all(5*j-3*d==2*k for _,j,d in row),'clock recovery')
            nodes+=len(row)
    # Full literal coverage audit of delivered depth16 atlas.
    atlas=json.loads((Path(__file__).resolve().parents[1]/'descent'/'atlas16.json').read_text())
    bits=atlas['bits'];mod=1<<bits;covered={}
    for c in atlas['certificates']:
        w=c['word'];a,b,cc=affine(w)
        need((a,b,cc)==(c['A'],c['B'],c['C']),'atlas affine recomputation')
        need(c['modulus']==2*b,'atlas modulus')
        for residue in range(c['residue'],mod,c['modulus']):
            need(residue not in covered,'atlas disjoint')
            covered[residue]=tuple(w)
        # independent endpoint coefficient proves infinite cylinder descent
        n=c['smallest_nonterminal'];y=n
        for v in w:
            y,actual=step(y);need(v==actual,'atlas native guard')
        need(y<n and a<b,'atlas all-ray payment')
    unresolved=[]
    for r in range(1,mod,2):
        n=r if r>1 else r+mod;x=n;word=[];clock=0;paid=None
        while clock<bits:
            x,a=step(x)
            if clock+a+1>bits:break
            clock+=a;word.append(a)
            if x<n:
                paid=tuple(word);break
        if paid is None:unresolved.append(r)
        else:need(covered.get(r)==paid,'literal first-descent atlas')
    need(unresolved==atlas['unresolved_residues'],'literal unresolved atlas')
    print(json.dumps({'status':'PASS','native_return_sources':limit//81,'native_incoming_edges':edges,
                      'reverse_endpoints':limit//81,'reverse_depth':15,'reverse_nodes':nodes,
                      'macro_checks':macro_checks,'atlas_bits':bits,'atlas_certificates':len(atlas['certificates']),
                      'atlas_unresolved':len(unresolved)},indent=2,sort_keys=True))
if __name__=='__main__':main()
