"""Exact paired-port observation/lattice certificate; no stochastic tests."""
from pathlib import Path
from itertools import combinations
from fractions import Fraction
from collections import Counter
from functools import reduce
from math import gcd, comb
import json,gzip,hashlib,sys
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
SOURCE=ROOT/'05-knowledge/results/session290_core_descent_20261113/descent/profiles.json.gz'
def need(b,why):
    if not b: raise RuntimeError(why)
def det(a):
    a=[list(r) for r in a];n=len(a);sg=1;prev=1
    for k in range(n-1):
        if a[k][k]==0:
            p=next((j for j in range(k+1,n) if a[j][k]),None)
            if p is None:return 0
            a[k],a[p]=a[p],a[k];sg=-sg
        piv=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                v=a[i][j]*piv-a[i][k]*a[k][j]
                need(v%prev==0,'Bareiss exact division')
                a[i][j]=v//prev
            a[i][k]=0
        prev=piv
    return sg*a[-1][-1]
def run():
    data=json.loads(gzip.decompress(SOURCE.read_bytes()));r=data['even_records']
    need(len(r)==24,'24 inherited charts')
    bs=sorted({b for q in r for b,c in q['U']});need(bs==[1,2,3,4,6,8,9,12],'fixed basis')
    M=[[dict(q['U']).get(b,0) for q in r] for b in bs]
    D=[[row[2*j+1]-row[2*j] for j in range(12)] for row in M]
    pairs=[]
    for j in range(12):
        a,b=r[2*j:2*j+2]
        need(a['rows']==b['rows'] and a['selected']==b['selected'],'same marked graph')
        need(a['ports']!=b['ports'],'different ports')
        need(sum(bs[i]*D[i][j] for i in range(8))==0,'terminal equality')
        pairs.append({'families':[2*j,2*j+1],'ports':[a['ports'],b['ports']], 'selected':a['selected'],'terminal_count':1+2*sum(c*b for b,c in a['U'])})
    need(len({(tuple(x['rows']),x['selected']) for x in r})==12,'12 marked bases')
    # These are presentations; no assertion that the eight row lists are eight isomorphism classes.
    literal_row_count=len({tuple(x['rows']) for x in r})
    sm=[abs(int(x)) for x in smith_normal_form(sp.Matrix(M),domain=ZZ).diagonal()]
    sd=[abs(int(x)) for x in smith_normal_form(sp.Matrix(D),domain=ZZ).diagonal()]
    need(sm==[1,1,1,2,2,2,18,216],'full profile Smith factors')
    need(sd==[1,1,2,2,6,72,216,0],'difference Smith factors')
    # Independent maximal-minor gcd proves the saturation index in ker(evaluation).
    g=0;nonzero=0
    for cs in combinations(range(12),7):
        z=det([[D[i][j] for j in cs] for i in range(1,8)])
        g=gcd(g,abs(z));nonzero+=bool(z)
    need(g==2**9*3**6,'independent saturation index')
    profiles=[]
    for e in range(4096):
        profiles.append(tuple(sum(D[i][j] for j in range(12) if e>>j&1) for i in range(8)))
    need(len(set(profiles))==4096,'all profile sums distinct')
    residue_images={str(p):len({tuple(x%p for x in v) for v in profiles}) for p in (2,3,5,7,13)}
    need(residue_images=={'2':4,'3':81,'5':1480,'7':3056,'13':4096},'finite coefficient-code images')
    L=sum(abs(x) for row in D for x in row);need(L==36288,'universal coefficient bound')
    N=1
    while L*8**N>=9**N:N+=1
    need(N==90,'exact dominance threshold')
    hist={}
    for m in range(1,N):
        vals=[sum(c*b**m for c,b in zip(v,bs)) for v in profiles]
        h=Counter(Counter(vals).values())
        if m>=3:need(len(set(vals))==4096,'all early scalar sums distinct')
        if m<=3:hist[str(m)]={str(k):v for k,v in sorted(h.items())}
    need(hist['1']=={'4096':1} and hist['2']=={'1':3932,'2':82},'small lengths')
    # Every polynomial vanishes at x=2,y=3. Its integer value therefore lies in (x-2,y-3).
    # Rank7 plus exact saturation index allows only factors2 and3 in the reverse inclusion.
    # Sharp valuations use Delta_7 mod4 and Delta_5,Delta_7 mod9.
    need(D[0][7]%4==2 and all(D[i][7]%4==0 for i in [2,4,5,6,7]),'mod4 witness structure')
    # At m>=2, x=2^m is0 mod4 and y=3^m is odd, so direct polynomial evaluation gives2.
    mon=lambda x,y:[1,x,y,x*x,x*y,x**3,y*y,x*x*y]
    for y in (1,3):need(sum(a*b for a,b in zip([D[i][7] for i in range(8)],mon(0,y)))%4==2,'v2=1')
    # m even gives x=1 mod3; m odd x in{2,5,8} mod9, y=0 mod9.
    need(sum(D[i][5]*mon(1,0)[i] for i in range(8))%3!=0,'even m no factor3')
    for x in (2,5,8):
        residues=[sum(D[i][j]*mon(x,0)[i] for i in range(8))%9 for j in (5,7)]
        need(any(z in (3,6) for z in residues),'odd m exact single3')
    for m in range(2,257):
        v=[sum(D[i][j]*bs[i]**m for i in range(8)) for j in range(12)]
        need(reduce(gcd,v)==gcd(2**m-2,3**m-3),'gcd formula controls')
    code_hist=Counter(2*e.bit_count() for e in range(4096))
    need(code_hist=={2*j:comb(12,j) for j in range(13)},'paired code enumerator')
    J=sp.Matrix([[0,-1],[1,-1]]);G=sp.Matrix([[2,-1],[-1,2]])
    need(J*J+J+sp.eye(2)==sp.zeros(2),'Eisenstein relation')
    need(J.T*G*J==G,'A2 norm invariance')
    roots=[(a,b) for a in range(-2,3) for b in range(-2,3) if 2*a*a-2*a*b+2*b*b==2]
    need(len(roots)==6,'A2 six roots')
    # Ternary Golay/Leech not imported: this computes only canonical paired carrier.
    out={'scope':'Inherited 24 even-chart specialization. Formal aggregate observations, not physical count addition or Chan-Pak coverage.',
      'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'bases':bs,'coefficient_matrix':M,'paired_difference_matrix':D,
      'pairs':pairs,'marked_base_count':12,'literal_adjacency_row_lists':literal_row_count,
      'profile_smith':sm,'difference_smith':sd,'difference_lattice_index':g,'nonzero_maximal_minors':nonzero,
      'profile_rank':8,'paired_difference_rank':7,'binary_profile_count':4096,'coefficient_residue_image_counts':residue_images,'scalar_fiber_histograms':hist,
      'dominance_coefficient_bound':L,'dominance_from':N,'finite_scalar_checked_lengths':[3,N-1],
      'gcd_controls':[2,256],'paired_code_weight_enumerator':dict(sorted(code_hist.items())),
      'paired_code_parameters':[24,12,2],'eisenstein_matrix':list(map(list,J.tolist())),
      'eisenstein_gram':list(map(list,G.tolist())),'a2_12_gram_determinant':3**12,'a2_12_root_count':72,
      'sympy_version':sp.__version__}
    name='receipt'+('.optimized' if not __debug__ else '')+'.json'
    (HERE/name).write_text(json.dumps(out,sort_keys=True,indent=2,default=int)+'\n')
    print(json.dumps({'status':'PASS','profiles':4096,'rank':8,'difference_rank':7,'index':g,'scalar_injective_from':3,'dominance_from':N,'receipt':name},sort_keys=True))
if __name__=='__main__':run()
