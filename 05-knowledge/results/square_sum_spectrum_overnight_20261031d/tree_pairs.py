#!/usr/bin/env python3
"""Twenty-five independently chosen weak-B children, with a member decoder."""
import hashlib
import json
from pathlib import Path
import sys
from math import isqrt

HERE=Path(__file__).resolve().parent
ROOTS=HERE.parent/'session124_mixed_radix_square_rays_20261031'
sys.path.insert(0,str(ROOTS))
from core import Constructor,prefix
from leaf_switches import apply
SWITCHES=json.loads((HERE/'leaf_switches.json').read_text())
BASE=SWITCHES['first_uniform_base_index']


def threshold(depth):return ((2*BASE+1)*49**depth-1)//2


class Trees:
    def __init__(self):
        self.constructor=Constructor()
        self.leaves={}
        self.decoded=0

    def pair(self,n,depth,bits):
        if n<threshold(depth) or len(bits)!=25**depth or any(b not in 'AB' for b in bits):
            raise ValueError('tree domain')
        if depth==0:
            key=(n,bits)
            if key not in self.leaves:
                if n>=1062:self.leaves[key]=self.constructor.pair(n,bits)
                else:
                    pair=self.constructor.bases.base(n)
                    self.leaves[key]=pair if bits=='A' else [apply(p,m) for p,m in zip(pair,SWITCHES['moves'][str(n)])]
            return self.leaves[key]
        m=(n-24)//49;d=n-49*m;size=25**(depth-1)
        children=[self.pair(m,depth-1,bits[a*size:(a+1)*size]) for a in range(25)]
        record=self.constructor.certificate('C',m%2,d)
        result=[]
        for word in record['words']:
            output=[]
            for atom in word:
                kind,c=atom[:2]
                if kind=='single':output.append(c)
                elif kind=='head':output.append(49-c)
                else:
                    _,c,i,o=atom;p=children[abs(c)][i]
                    indices=range(int(kind=='tail'),len(p)) if o==1 else range(len(p)-1,int(kind=='tail')-1,-1)
                    output.extend(49*p[j]+(-1)**(j+1)*c for j in indices)
            result.append(output)
        self.constructor.bases.check_pair(result,prefix('C'*depth))
        return result

    def decode(self,path,n,depth,index):
        if len(path)!=n+index:raise RuntimeError('decoder member length')
        if depth==0:
            if path[:2]!=[1,3]:raise RuntimeError('leaf start')
            if n<1062:
                original=self.constructor.bases.base(n)[index]
                if path==original:return 'A'
                if path==apply(original,SWITCHES['moves'][str(n)][index]):return 'B'
                raise RuntimeError('switched leaf signature')
            if path[2]==13:return 'A'
            if path[2]==6:return 'B'
            raise RuntimeError('leaf signature')
        m=(n-24)//49;d=n-49*m
        record=self.constructor.certificate('C',m%2,d)
        return self.decode_word(path,record['words'][index],m,depth)

    def decode_word(self,path,word,m,depth):
        at=0;children={}
        for atom in word:
            kind,c=atom[:2]
            if kind in ('single','head'):
                expected=c if kind=='single' else 49-c
                if path[at]!=expected:raise RuntimeError('fixed singleton/head')
                at+=1
            else:
                _,c,i,o=atom;length=m+i-int(kind=='tail')
                if c>=0:
                    raw=path[at:at+length]
                    if o==-1:raw=raw[::-1]
                    child=[(v+24)//49 for v in raw]
                    if kind=='tail':child.insert(0,1)
                    children[c]=self.decode(child,m,depth-1,i)
                    self.decoded+=1
                at+=length
        if at!=len(path) or set(children)!=set(range(25)):raise RuntimeError('all child channels')
        return ''.join(children[a] for a in range(25))

    def cycle(self,n,depth,bits):
        if depth<1 or n<threshold(depth) or len(bits)!=25**depth:
            raise ValueError('terminal tree domain')
        m=(n-24)//49;d=n-49*m;size=25**(depth-1)
        children=[self.pair(m,depth-1,bits[a*size:(a+1)*size]) for a in range(25)]
        word=self.constructor.certificate('C',m%2,d,True)['word']
        output=[]
        for atom in word:
            kind,c=atom[:2]
            if kind=='single':output.append(c)
            elif kind=='head':output.append(49-c)
            else:
                _,c,i,o=atom;p=children[abs(c)][i]
                indices=range(int(kind=='tail'),len(p)) if o==1 else range(len(p)-1,int(kind=='tail')-1,-1)
                output.extend(49*p[j]+(-1)**(j+1)*c for j in indices)
        seen=bytearray(n+1)
        P=prefix('C'*depth)
        if len(output)!=n or output[:len(P)]!=P:raise RuntimeError('terminal size/root prefix')
        for j,v in enumerate(output):
            if not 1<=v<=n or seen[v]:raise RuntimeError('terminal ownership')
            seen[v]=1;s=v+output[j-1]
            if isqrt(s)**2!=s:raise RuntimeError('terminal square edge')
        if self.decode_word(output,word,m,depth)!=bits:raise RuntimeError('terminal decoder')
        digest=hashlib.sha256()
        for v in output:digest.update((str(v)+'\n').encode())
        return digest.hexdigest()


def code(depth,salt):
    return ''.join('AB'[hashlib.sha256(f'{salt}:{a}'.encode()).digest()[0]%2] for a in range(25**depth))


def main():
    trees=Trees();records=[]
    # Complete32-choice cube on five channels; the other20 remain marked/fixed.
    hashes=[set(),set()]
    for choice in range(32):
        bits=''.join('AB'[(choice>>j)&1] for j in range(5))+code(1,'fixed')[5:]
        pair=trees.pair(threshold(1),1,bits)
        for i,path in enumerate(pair):
            if trees.decode(path,threshold(1),1,i)!=bits:raise RuntimeError('cube decoder')
            hashes[i].add(hashlib.sha256(','.join(map(str,path)).encode()).hexdigest())
    if any(len(h)!=32 for h in hashes):raise RuntimeError('member injectivity')
    records.append({'depth':1,'N':threshold(1),'complete_tested_choice_cube':32,'distinct_each_member':[len(h) for h in hashes]})
    for depth,extra,salt in [(1,1,'opposite parity'),(2,0,'deep')]:
        n=threshold(depth)+extra;bits=code(depth,salt);pair=trees.pair(n,depth,bits)
        digests=[]
        for i,path in enumerate(pair):
            if trees.decode(path,n,depth,i)!=bits:raise RuntimeError('deep decoder')
            digest=hashlib.sha256()
            for v in path:digest.update((str(v)+'\n').encode())
            digests.append(digest.hexdigest())
        records.append({'depth':depth,'N':n,'independent_leaf_choices':len(bits),'code_sha256':hashlib.sha256(bits.encode()).hexdigest(),'member_sha256':digests})
    cycles=[]
    for depth,extra in [(1,0),(1,1),(2,0)]:
        n=threshold(depth)+extra;bits=code(depth,f'cycle:{extra}')
        cycles.append({'depth':depth,'N':n,'leaf_choices':len(bits),'cycle_sha256':trees.cycle(n,depth,bits)})
    print(json.dumps({'status':'passed','base':BASE,'records':records,'decoded_child_channels':trees.decoded,
                      'cycles':cycles,'pair_checks':trees.constructor.bases.checks,'certificate_checks':trees.constructor.checks},indent=2,sort_keys=True))


if __name__=='__main__':main()
