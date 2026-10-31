#!/usr/bin/env python3
"""All signed offsets independent, within a fixed individual-label phase fibre."""
import hashlib,json,sys
from math import isqrt,prod
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'session124_mixed_radix_square_rays_20261031'))
from core import Constructor,OPS,prefix,vertices
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031d'))
from leaf_switches import apply
BANK=json.loads((HERE/'phase_switches.json').read_text())

def arity(word):
    if any(s not in 'AC' for s in word):raise ValueError('only audited A/C')
    return prod(OPS[s][0] for s in word)
def threshold(word):return (225*arity(word)-1)//2
def digest(p):return hashlib.sha256(''.join(f'{x}\n' for x in p).encode()).hexdigest()
def phase(p):return tuple(j%2 for v,j in sorted((v,j) for j,v in enumerate(p)))

class Trees:
    def __init__(self):
        self.c=Constructor();self.supply_cache={};self.decodes=0

    def supply(self,n,bit):
        if n<112 or bit not in '01':raise ValueError('same-phase supply domain')
        key=(n,bit)
        if key in self.supply_cache:return self.supply_cache[key]
        if n<2812:
            p=self.c.pair(n)
            out=p if bit=='0' else [apply(v,m) for v,m in zip(p,BANK['moves'][str(n)])]
        else:
            m=(n-12)//25;d=n-25*m
            p=self.supply(m,bit);r=self.c.certificate('A',m%2,d)
            out=[list(vertices(25,p,w)) for w in r['words']]
        self.c.bases.check_pair(out)
        self.supply_cache[key]=out
        return out

    def pair(self,n,word,bits):
        size=arity(word)
        if n<threshold(word) or len(bits)!=size or any(b not in '01' for b in bits):
            raise ValueError('tree domain')
        if not word:return self.supply(n,bits)
        q=OPS[word[0]][0];half=(q-1)//2;m=(n-half)//q;d=n-q*m;chunk=size//q
        children={c:self.pair(m,word[1:],bits[(c+half)*chunk:(c+half+1)*chunk]) for c in range(-half,half+1)}
        record=self.c.certificate(word[0],m%2,d)
        out=[self.expand(w,q,children) for w in record['words']]
        self.c.bases.check_pair(out,prefix(word))
        return out

    @staticmethod
    def expand(word,q,children):
        out=[]
        for a in word:
            kind,c=a[:2]
            if kind=='single':out.append(c)
            elif kind=='head':out.append(q-c)
            else:
                _,c,i,o=a;p=children[c][i];start=int(kind=='tail')
                indices=range(start,len(p)) if o==1 else range(len(p)-1,start-1,-1)
                out.extend(q*p[j]+(-1)**(j+1)*c for j in indices)
        return out

    def decode(self,path,n,word,index):
        if len(path)!=n+index:raise RuntimeError('member length')
        if not word:
            matches=[b for b in '01' if path==self.supply(n,b)[index]]
            if len(matches)!=1:raise RuntimeError('same-phase leaf signature')
            return matches[0]
        q=OPS[word[0]][0];half=(q-1)//2;m=(n-half)//q;d=n-q*m
        r=self.c.certificate(word[0],m%2,d)
        return self.decode_word(path,r['words'][index],q,m,word[1:])

    def decode_word(self,path,atoms,q,m,suffix):
        half=(q-1)//2;at=0;children={}
        for a in atoms:
            kind,c=a[:2]
            if kind in ('single','head'):
                if path[at]!=(c if kind=='single' else q-c):raise RuntimeError('fixed atom')
                at+=1
            else:
                _,c,i,o=a;length=m+i-int(kind=='tail');raw=path[at:at+length]
                if o==-1:raw=raw[::-1]
                child=[(v+half)//q for v in raw]
                if kind=='tail':child.insert(0,1)
                if c in children:raise RuntimeError('duplicate readable channel')
                children[c]=self.decode(child,m,suffix,i);self.decodes+=1;at+=length
        if at!=len(path) or set(children)!=set(range(-half,half+1)):raise RuntimeError('signed coverage')
        return ''.join(children[c] for c in range(-half,half+1))

    def cycle(self,n,word,bits):
        size=arity(word)
        if not word or n<threshold(word) or len(bits)!=size:raise ValueError('cycle domain')
        q=OPS[word[0]][0];half=(q-1)//2;m=(n-half)//q;d=n-q*m;chunk=size//q
        children={c:self.pair(m,word[1:],bits[(c+half)*chunk:(c+half+1)*chunk]) for c in range(-half,half+1)}
        r=self.c.certificate(word[0],m%2,d,True)
        atoms=r['word'] if word[0]=='C' else r['words'][0]
        out=self.expand(atoms,q,children)
        if len(out)!=n or sorted(out)!=list(range(1,n+1)) or out[:len(prefix(word))]!=prefix(word):
            raise RuntimeError('cycle ownership or prefix')
        if any(isqrt(v+out[j-1])**2!=v+out[j-1] for j,v in enumerate(out)):
            raise RuntimeError('cycle edge')
        if self.decode_word(out,atoms,q,m,word[1:])!=bits:raise RuntimeError('cycle decoder')
        return out

def code(word,salt):
    return ''.join(str(hashlib.sha256(f'{salt}:{a}'.encode()).digest()[0]%2) for a in range(arity(word)))

def main():
    t=Trees();supply=[]
    for n in [112,1111,2811,2812,2813,70312,70313]:
        p=t.supply(n,'0');z=t.supply(n,'1');standard=t.c.pair(n)
        if p!=standard or any(a==b or phase(a)!=phase(b) for a,b in zip(p,z)):
            raise RuntimeError('uniform same canonical phase')
        supply.append({'n':n,'member_sha256':[digest(x) for x in z]})
    cubes=[]
    for w in ['A','C']:
        q=arity(w);half=(q-1)//2;channels=[-half,-1,0,1,half];hashes=[set(),set()];colors=None
        for k in range(32):
            bits=list(code(w,'fixed'))
            for j,c in enumerate(channels):bits[c+half]=str((k>>j)&1)
            bits=''.join(bits);p=t.pair(threshold(w),w,bits)
            for i,x in enumerate(p):
                if t.decode(x,threshold(w),w,i)!=bits:raise RuntimeError('cube decode')
                hashes[i].add(digest(x))
            ph=[phase(x) for x in p]
            if colors is not None and colors!=ph:raise RuntimeError('phase across independent signed choices')
            colors=ph
        if any(len(s)!=32 for s in hashes):raise RuntimeError('cube injectivity')
        cubes.append({'word':w,'n':threshold(w),'channels':channels,'choices':32,'distinct_members':[len(s) for s in hashes]})
    deep=[]
    for w,extra in [('A',1),('C',1),('AA',0),('AC',0),('CA',1),('CC',0)]:
        n=threshold(w)+extra;bits=code(w,'deep');p=t.pair(n,w,bits)
        for i,x in enumerate(p):
            if t.decode(x,n,w,i)!=bits:raise RuntimeError('deep decode')
        cyc=t.cycle(n,w,bits)
        deep.append({'word':w,'n':n,'bits':len(bits),'pair_sha256':[digest(x) for x in p],'cycle_sha256':digest(cyc)})
    print(json.dumps({'status':'passed','uniform_supply':supply,'signed_cubes':cubes,'deep':deep,'decoded_channels':t.decodes},indent=2,sort_keys=True))

if __name__=='__main__':main()
