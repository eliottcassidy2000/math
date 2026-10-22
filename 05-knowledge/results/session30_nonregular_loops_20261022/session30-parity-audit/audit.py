"""Independent integer basis conversion and direct contraction audit.
No producer-code imports. Reads frozen source and candidate JSON tensors only.
"""
from pathlib import Path
import json,itertools,collections,math,time,hashlib
import sympy as sp
import numpy as np
BASE=Path(__file__).resolve().parent.parent
old=json.loads((BASE/'session30-decorated-engine/symbolic-kernels.json').read_text())
new=json.loads((BASE/'session30-parity-basis/parity-kernels.json').read_text())
# Candidate source variables written in independent G,H,O,P variables.
replacement={0:[(1,1),(2,1)],1:[(1,1),(2,-1)],2:[(1,-1),(3,1),(2,-1)],3:[(0,-1)]}
x=sp.symbols('x');G=-x*x/(1+x*x);H=x*x/(1-x*x);O=x/(1-x*x);P=x/(1+x*x)
rp=x/(1-x);rm=-x/(1+x);qpp=-rp*rp*(1-x*x)/(1+x*x);qmm=-rm*rm*(1-x*x)/(1+x*x);qpm=-rp*rm*(1-x*x)/(1+x*x)
identities=[rp-H-O,rm-H+O,qpp+H-P+O,qmm+H+P-O,qpm+G]
if any(sp.cancel(v)!=0 for v in identities):raise RuntimeError('scalar identity failure')
if len(old)!=31 or len(new)!=31:raise RuntimeError('core count')
# Reversal acts linearly on the ORIGINAL four formal basis symbols.
# This check does not rely on the candidate parity coefficients.
reversal={0:[(1,1)],1:[(0,1)],2:[(0,-1),(1,-1),(2,-1)],3:[(3,1)]}
for core in old:
 rev=collections.defaultdict(int)
 for term in core['terms']:
  for choices in itertools.product(*(reversal[i] for i in term['indices'])):
   rev[tuple(i for i,c in choices)]+=term['numerator']*math.prod(c for i,c in choices)
 rev={k:v for k,v in rev.items() if v}
 assert rev=={tuple(t['indices']):t['numerator'] for t in core['terms']}
summary=[]
for a,b in zip(old,new):
 for key in ['V','edges','excess','leaves','weight','automorphisms','denominator']:
  if a[key]!=b[key]:raise RuntimeError('metadata mismatch')
 transformed=collections.defaultdict(int)
 for term in a['terms']:
  for choices in itertools.product(*(replacement[i] for i in term['indices'])):
   idx=tuple(i for i,c in choices);value=term['numerator']*math.prod(c for i,c in choices)
   transformed[idx]+=value
 transformed={k:v for k,v in transformed.items() if v}
 expected={tuple(t['indices']):t['numerator'] for t in b['terms']}
 if transformed!=expected:raise RuntimeError('exact coefficient mismatch')
 if any(sum(i>=2 for i in k)%2 for k in transformed):raise RuntimeError('odd parity survives')
 if any(any(edge[0]==edge[1] and k[j]>=2 for j,edge in enumerate(a['edges'])) for k in transformed):raise RuntimeError('odd loop propagator survives')
 summary.append({'weight':a['weight'],'edges':a['edges'],'old_terms':len(a['terms']),'new_terms':len(transformed)})

def matrices(M,z):
 Q=M*z;I=np.eye(len(M),dtype=Q.dtype);Q2=Q@Q
 rp=Q@np.linalg.inv(I-Q);rm=-Q@np.linalg.inv(I+Q)
 K=(I-Q2)@np.linalg.inv(I+Q2)
 oldb=[rp,rm,-rp@rp@K,-rp@rm@K]
 newb=[-Q2@np.linalg.inv(I+Q2),Q2@np.linalg.inv(I-Q2),Q@np.linalg.inv(I-Q2),Q@np.linalg.inv(I+Q2)]
 # Verify all matrix substitutions separately from tensor evaluation.
 for i in range(4):
  value=sum(c*newb[j] for j,c in replacement[i])
  if not np.allclose(value,oldb[i],rtol=1e-12,atol=1e-12):raise RuntimeError('matrix identity failure')
 return oldb,newb

def contraction(core,B):
 spec=','.join('abcd'[u]+'abcd'[v] for u,v in core['edges'])+'->'
 trial=[B[0]]*len(core['edges']);path=np.einsum_path(spec,*trial,optimize='greedy')[0]
 return sum(t['numerator']*np.einsum(spec,*[B[k] for k in t['indices']],optimize=path) for t in core['terms'])/core['denominator']

def mat(n,bits):
 S=np.zeros((n,n),int)
 for b,(i,j) in zip(bits,itertools.combinations(range(n),2)):
  S[i,j]=1 if b=='1' else -1;S[j,i]=-S[i,j]
 return S/n
fixtures=[('weighted3',np.array([[0,.1,.2],[-.1,0,-.15],[-.2,.15,0]])),('transitive4',np.sign(np.arange(4)[:,None]-np.arange(4)[None,:])/4),('almostregular6',mat(6,'000110011000000')),('regular7',np.array([[0 if i==j else 1 if (j-i)%7 in(1,2,3) else -1 for j in range(7)] for i in range(7)])/7)]
checks=[];maxerr=0
for name,M in fixtures:
 for z in [1.,.3+.1j]:
  A,B=matrices(M,z)
  errors=[abs(contraction(a,A)-contraction(b,B)) for a,b in zip(old,new)]
  worst=max(errors);maxerr=max(maxerr,worst)
  if worst>2e-11:raise RuntimeError('contraction mismatch')
  checks.append({'fixture':name,'z':str(z),'core_count':31,'max_absolute_difference':float(worst)})
# Identical generic sparse-contraction evaluator; not a whole-engine benchmark.
M=fixtures[-1][1];A,B=matrices(M,1.);timings={'old':[],'parity':[]}
for repeat in range(3):
 for name,cores,basis in [('old',old,A),('parity',new,B)]:
  t=time.perf_counter();v=sum(contraction(c,basis) for c in cores);timings[name].append(time.perf_counter()-t)
receipt={'status':'PASS','source_imports':'No producer-code imports; JSON inputs only','exact_scalar_identities':len(identities),'cores':len(old),'old_nonzero_terms':sum(r['old_terms'] for r in summary),'parity_nonzero_terms':sum(r['new_terms'] for r in summary),'original_basis_reversal_invariants':31,'odd_monomials':0,'odd_loop_terms':0,'contraction_checks':checks,'maximum_absolute_difference':maxerr,'generic_sparse_timings_seconds':timings,'timing_scope':'Same direct sparse einsum evaluator on all31 cores of one regular7 fixture; not a full predictor benchmark','summary':summary}
print(json.dumps(receipt,indent=2))
