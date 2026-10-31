#!/usr/bin/env python3
"""Exact two-seam gluing of closed square-sum carriers with square ports."""
from dataclasses import dataclass
from itertools import product
from math import isqrt,gcd
from pathlib import Path
import hashlib,json,sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031f'))
from edge_atlas import Atlas,edge,need

def square(n):return isqrt(n)**2==n
def scale(es,f):return {edge(f*a,f*b) for a,b in es}

@dataclass
class Closed:
    vertices:set
    base:set
    moves:list
    port:tuple
    square_class:int=1
    norm_seed:tuple=None

    def edges(self,bits):
        need(len(bits)==len(self.moves) and set(bits)<={'0','1'},'branch clock')
        es=set(self.base)
        for z,(r,a) in zip(bits,self.moves):
            if z=='1':es.difference_update(r);es.update(a)
        return es

    def walk(self,bits):
        es=self.edges(bits);adj={v:[] for v in self.vertices}
        for u,v in es:
            need(u in adj and v in adj and square(u+v),'native square edge')
            adj[u].append(v);adj[v].append(u)
        need(all(len(ns)==2 for ns in adj.values()),'exact degree two')
        u,v=self.port;out=[u];last=u
        while v!=u:
            need(len(out)<len(adj),'no overflow')
            out.append(v);a,b=adj[v];last,v=v,b if a==last else a
        need(len(out)==len(adj) and len(set(out))==len(adj),'connected full traversal')
        return out

    def check_port(self):
        u,v=self.port
        d=self.square_class
        need(d>=1 and all(d%(p*p) for p in range(2,isqrt(d)+1)),'normalized square class')
        need(self.port in self.base and u%d==v%d==0 and all(square(x) for x in (u//d,v//d,u+v)),'same-class square-sum edge')
        need(all(self.port not in r|a for r,a in self.moves),'invariant marked port')
        if self.norm_seed is None:
            a,b=isqrt(u//d),isqrt(v//d);g=gcd(a,b);a//=g;b//=g
            w=isqrt((a*a+b*b)//d)
            need(a*a+b*b==d*w*w,'rational norm seed')
            self.norm_seed=(a,b,w)
        a,b,w=self.norm_seed
        need(gcd(a,b)==1 and a*a+b*b==d*w*w,'primitive norm seed')

def from_cube(c,port=(9,16),square_class=1):
    out=Closed(set(c.reference),set(c.base),[c.moves[b] for b in range(c.bits)],edge(*port),square_class)
    out.check_port();return out

def condition_cube(c,port,square_class):
    compiled=c.compile([(port,True)]);need(compiled is not None,'feasible seam edge')
    free=[b for b in range(c.bits) if b not in compiled['pins']]
    out=Closed(set(c.reference),c.edge_set(compiled['one_code']),[c.moves[b] for b in free],edge(*port),square_class)
    out.check_port()
    return out,{'pins':compiled['pins'],'free_source_bits':free}

def glue(A,B,common_frame=None):
    A.check_port();B.check_port();need(A.port==B.port,'same normalized port type')
    need(A.square_class==B.square_class,'same square class')
    u,v=A.port
    M=max(max(A.vertices),max(B.vertices))
    d=A.square_class;m=None;t=None;budget=M;seed=A.norm_seed
    if common_frame is not None:
        need(common_frame['square_class']==d,'common norm class')
        m,t=common_frame['m'],common_frame['t']
        r,k,l=(common_frame[z] for z in ('r','k','l'))
        constant=common_frame['label_constant'];seed=tuple(common_frame['norm_seed'])
        budget=common_frame['frame_maximum_label']
        need(budget>=M,'common frame covers both source sizes')
    elif d==1:
        m=isqrt(isqrt(M))+2;m+=m%2
        r=m*m-1;k=2*m;l=m*m+1;constant=256
    else:
        a0,b0,w=A.norm_seed;L=2*d*w
        T=max(4*b0+4,2*(isqrt(isqrt(M))+1))
        t=L*((T+L-1)//L)
        r=(t*t-1)*a0-2*t*b0;k=(t*t-1)*b0+2*t*a0;l=w*(t*t+1)
        K=4*b0+8+L
        constant=max(a0*K*K,b0*K*K+2*a0*K)**2
    need(r>0 and k>0 and r*r+k*k==d*l*l and gcd(r,k)==1 and r*r>M,'coprime disjoint scale guard')
    a,b=r*r,k*k
    va={a*x for x in A.vertices};vb={b*x for x in B.vertices}
    need(not va&vb,'disjoint source copies')
    left,right=scale(A.base,a),scale(B.base,b)
    left.remove(edge(a*u,a*v));right.remove(edge(b*u,b*v))
    seams={edge(a*u,b*u),edge(a*v,b*v)}
    need(all(square(x+y) for x,y in seams),'two exact arithmetic seams')
    moves=[(scale(R,a),scale(S,a)) for R,S in A.moves]+[(scale(R,b),scale(S,b)) for R,S in B.moves]
    out=Closed(va|vb,left|right|seams,moves,edge(a*u,b*u),d,seed);out.check_port()
    union=out.base|set().union(*(S for R,S in moves));degrees=CounterDegrees(union)
    need(max(degrees.values())<=3,'subcubic carrier preserved')
    need({e for e in union if bool(e[0] in va)!=bool(e[1] in va)}==seams,'exact two-edge boundary')
    need(max(out.vertices)<=constant*budget*budget,'quadratic maximum-label bound')
    record={'m':m,'t':t,'square_class':d,'norm_seed':seed,'label_constant':constant,'r':r,'k':k,'l':l,'gcd_legs':gcd(r,k),'input_maximum_label':M,'frame_maximum_label':budget,'source_scales':[a,b],'source_vertices':[len(va),len(vb)],
            'source_bits':[len(A.moves),len(B.moves)],'seams':sorted(seams),'new_port':out.port,
            'output_vertices':len(out.vertices),'output_bits':len(moves),'maximum_label':max(out.vertices)}
    return out,record

def CounterDegrees(es):
    out={}
    for e in es:
        for v in e:out[v]=out.get(v,0)+1
    return out

def digest(p):return hashlib.sha256(','.join(map(str,p)).encode()).hexdigest()

def main():
    atlas=Atlas();A=from_cube(atlas.cycle(2812,'A'));B=from_cube(atlas.cycle(2862,'A'))
    # Authenticate the all-size outer state, not only these two concrete cycles.
    cert=atlas.t.c.certificate('A',0,12,True);atoms=cert['words'][0]
    pair=[['single',9],['head',9]]
    need(any(atoms[i:i+2]==pair for i in range(len(atoms)-1)),'constant 9--16 seam for N=50t+12')
    C,rec=glue(A,B);hashes=set()
    for word in product('01',repeat=6):
        bits=''.join(word[:3])+'0'*22+''.join(word[3:])+'0'*22
        actual=C.walk(bits);hashes.add(digest(actual))
        # Pull back both retained edge sets and restore their old cut edges.
        es=C.edges(bits);a,b=rec['source_scales'];va={a*v for v in A.vertices}
        left={edge(u//a,v//a) for u,v in es if u in va and v in va}|{A.port}
        right={edge(u//b,v//b) for u,v in es if u not in va and v not in va}|{B.port}
        need(left==A.edges(bits[:25]) and right==B.edges(bits[25:]),'exact source and bit recovery')
    need(len(hashes)==64,'complete six-bit product injective')
    D,rec2=glue(C,C)
    tests=[]
    for salt in range(4):
        bits=''.join(str(hashlib.sha256(f'{salt}:{i}'.encode()).digest()[0]%2) for i in range(len(D.moves)))
        tests.append(digest(D.walk(bits)))
    # Failure boundary: an ordinary square-sum edge need not have square endpoints.
    wrong=Closed(A.vertices,A.base,A.moves,(1,3));rejected=False
    try:wrong.check_port()
    except RuntimeError:rejected=True
    need(rejected,'square endpoints are a genuine type guard')
    need(6*6>16 and 6*6*16==8*8*9,'coprimality cannot be dropped from separation')
    print(json.dumps({'status':'passed','base_progression':{'N':'50t+12','t_min':56,'port':[9,16],
          'symbolic_atom_pair':pair},'first_glue':rec,'second_glue':rec2,
          'complete_free_product_bits':6,'distinct_complete_cycles':len(hashes),
          'second_glue_full_code_sha256':tests,'hostile_non_square_endpoint_rejected':rejected,
          'noncoprime_collision':{'legs':[6,8],'source_labels':[16,9],'common_scaled_label':576},
          'scope':'Exact Hamiltonian product for the selected subcubic carrier on a sparse union of square-scaled vertex sets. Not an initial interval or the full induced square-sum graph.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
