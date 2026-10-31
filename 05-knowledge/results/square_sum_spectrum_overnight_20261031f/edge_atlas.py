#!/usr/bin/env python3
"""A product of disjoint edge switches, with an exact boundary-query compiler."""
import json,sys,hashlib
from dataclasses import dataclass
from math import isqrt
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031e'))
from phase_trees import Trees,OPS,arity,threshold,code,digest

def need(x,msg):
    if not x:raise RuntimeError(msg)
def edge(a,b):return (min(a,b),max(a,b))
def edges(p,cycle=False):
    out={edge(a,b) for a,b in zip(p,p[1:])}
    if cycle:out.add(edge(p[-1],p[0]))
    return out
@dataclass
class Carrier:
    path:list
    moves:dict

class Atlas:
    def __init__(self):self.t=Trees();self.leaves={}

    def pair(self,n,word,offset=0):
        need(n>=threshold(word),'atlas threshold')
        if not word:
            if n not in self.leaves:
                p=self.t.supply(n,'0');z=self.t.supply(n,'1');raw=[]
                for a,b in zip(p,z):
                    old,new=edges(a),edges(b);r=old-new;s=new-old
                    need(r and s and edge(1,3) not in r|s,'readable leaf switch away from head cut')
                    raw.append((a,(r,s)))
                self.leaves[n]=raw
            return [Carrier(p,{offset:move}) for p,move in self.leaves[n]]
        q=OPS[word[0]][0];half=(q-1)//2;m=(n-half)//q;d=n-q*m;chunk=arity(word[1:])
        children={c:self.pair(m,word[1:],offset+(c+half)*chunk) for c in range(-half,half+1)}
        r=self.t.c.certificate(word[0],m%2,d)
        return [self.expand(w,q,children) for w in r['words']]

    def expand(self,atoms,q,children):
        p=self.t.expand(atoms,q,{c:[v.path for v in pair] for c,pair in children.items()})
        moves={}
        for atom in atoms:
            if atom[0] in ('single','head'):continue
            kind,c,i,o=atom;child=children[c][i]
            ph={v:j%2 for j,v in enumerate(child.path)}
            def lift(v):return q*v+(-1 if ph[v]==0 else 1)*c
            for b,(removed,added) in child.moves.items():
                need(b not in moves,'one selected source member per bit')
                if kind=='tail':need(edge(1,3) not in removed|added,'cut is invariant')
                moves[b]=tuple({edge(lift(u),lift(v)) for u,v in part} for part in (removed,added))
        return Carrier(p,moves)

    def cycle(self,n,word):
        need(bool(word) and n>=threshold(word),'terminal atlas domain')
        q=OPS[word[0]][0];half=(q-1)//2;m=(n-half)//q;d=n-q*m;chunk=arity(word[1:])
        children={c:self.pair(m,word[1:],(c+half)*chunk) for c in range(-half,half+1)}
        r=self.t.c.certificate(word[0],m%2,d,True);atoms=r['word'] if word[0]=='C' else r['words'][0]
        out=self.expand(atoms,q,children)
        return EdgeCube(out)

class EdgeCube:
    def __init__(self,carrier,anchor=3):
        self.reference=carrier.path;self.moves=carrier.moves;self.base=edges(self.reference,True)
        self.anchor=anchor
        need(edge(1,anchor) in self.base,'oriented cycle anchor')
        self.controls={};self.owners={};self.bits=len(self.moves)
        need(set(self.moves)==set(range(self.bits)),'continuous bit addresses')
        for b,(removed,added) in self.moves.items():
            need(removed<=self.base and not added&self.base,'genuine edge replacement')
            need(len(removed)==len(added)>=2,'nontrivial balanced move')
            for value,part in [(0,removed),(1,added)]:
                for e in part:
                    need(e not in self.controls,'each changed edge has one bit')
                    self.controls[e]=(b,value)
                    for v in e:
                        need(v not in self.owners or self.owners[v]==b,'disjoint changed-vertex supports')
                        self.owners[v]=b

    def edge_set(self,bits):
        need(len(bits)==self.bits and set(bits)<={'0','1'},'exact bit string')
        out=set(self.base)
        for b,bit in enumerate(bits):
            if bit=='1':
                removed,added=self.moves[b];out.difference_update(removed);out.update(added)
        return out

    def compile(self,queries):
        pins={}
        for e,wanted in queries:
            e=edge(*e)
            if e in self.controls:
                b,present=self.controls[e];value=present if wanted else 1-present
                if b in pins and pins[b]!=value:return None
                pins[b]=value
            elif (e in self.base)!=wanted:return None
        return {'pins':pins,'free_bits':self.bits-len(pins),'one_code':''.join(str(pins.get(b,0)) for b in range(self.bits))}

    def walk(self,bits):
        es=self.edge_set(bits);n=len(self.reference);adj={v:[] for v in range(1,n+1)}
        for a,b in es:
            need(1<=a<=n and 1<=b<=n and isqrt(a+b)**2==a+b,'native square edge')
            adj[a].append(b);adj[b].append(a)
        need(len(es)==n and all(len(v)==2 for v in adj.values()),'2-regular exact interval')
        out=[1];last=1;v=self.anchor
        while v!=1:
            need(v not in (1,) and len(out)<n,'no premature duplicate/overflow')
            need(last in adj[v],'walk follows a retained edge')
            out.append(v);nxt=adj[v][0] if adj[v][1]==last else adj[v][1];last,v=v,nxt
        need(len(out)==n and len(set(out))==n,'one Hamiltonian cycle')
        return out

def main():
    a=Atlas();records=[]
    bank=json.loads((HERE.parent/'square_sum_spectrum_overnight_20261031e/phase_switches.json').read_text())
    exchange_n=min(int(n) for n,ms in bank['moves'].items() if any(m[0]=='exchange' for m in ms))
    for word,extra in [('A',0),('C',1),('A',25*exchange_n+12-threshold('A')),('A',25*exchange_n+36-threshold('A')),('A',67500),('AA',0),('AC',1),('CC',0)]:
        n=threshold(word)+extra;cube=a.cycle(n,word);bits=code(word,'independent switches')
        actual=a.t.cycle(n,word,bits);recovered=cube.walk(bits)
        need(edges(actual,True)==edges(recovered,True),'independent edge-union realization')
        # Any marked vertex's complete reference neighbourhood costs <=one bit.
        pinned=[1,3]+[min(v for e in cube.moves[b][0]|cube.moves[b][1] for v in e) for b in range(4)]
        constraints=[(e,True) for e in cube.base if any(v in e for v in pinned)]
        compiled=cube.compile(constraints)
        need(compiled is not None and len(compiled['pins'])<=len(pinned),'vertex-port budget')
        need(len(compiled['pins'])==4,'changed-vertex controls genuinely pay four bits')
        need(all((e in cube.edge_set(compiled['one_code']))==wanted for e,wanted in constraints),'compiled witness')
        b=0;removed,added=cube.moves[b]
        hostile=cube.compile([(min(removed),True),(min(added),True)])
        need(hostile is None,'two incompatible literals on one bit')
        weights={}
        for r,s in cube.moves.values():weights[len(r)+len(s)]=weights.get(len(r)+len(s),0)+1
        records.append({'word':word,'N':n,'bits':cube.bits,'changed_vertex_support':len(cube.owners),
                        'edge_hamming_weights':weights,'pinned_vertices':pinned,'paid_bits':len(compiled['pins']),
                        'cycle_sha256':digest(recovered)})
    # Complete64-choice fibre and complete single/two-query checks on its first two bits.
    cube=a.cycle(threshold('A'),'A');universe=[]
    for k in range(64):
        bits=format(k,'06b')+'0'*(cube.bits-6);es=cube.edge_set(bits)
        need(edges(cube.walk(bits),True)==edges(a.t.cycle(threshold('A'),'A',bits),True),'complete cube actual paths')
        universe.append((bits,es))
    need(cube.compile([(edge(1,3),False)]) is None,'cannot forbid invariant anchor')
    need(cube.compile([(edge(1,2),True)]) is None,'cannot require edge outside family union')
    need(cube.compile([(edge(1,3),True)])['free_bits']==cube.bits,'fixed edge costs no bit')
    need(cube.compile([([3,1],True)])['free_bits']==cube.bits,'unordered JSON edge query')
    queries=[(e,v) for b in range(2) for part in cube.moves[b] for e in part for v in (False,True)]
    tests=0
    for qa in queries:
        for qb in queries:
            found=cube.compile([qa,qb]);matching=[bits for bits,es in universe if all((e in es)==want for e,want in (qa,qb))]
            need((found is None)==(not matching),'complete conflicting query result')
            if found is not None:need(len(matching)==1<<(6-len(found['pins'])),'exact free-face count')
            tests+=1
    print(json.dumps({'status':'passed','records':records,'complete_choice_cube':64,'complete_two_query_tests':tests,
                      'scope':'Exact extension relation inside the authenticated cycle family; arbitrary external forests are not asserted extendable.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
