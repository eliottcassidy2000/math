import sys,itertools,json,math
sys.path.append('/workspace/scratch/07d5532f2fce/session228-fiber-switch')
from transfer import matrix, S,T,local
from search import stats
import sympy as s
adj=[28,5,16,6,10];H,Q,alpha=stats(adj)
assert (H,Q,alpha)==(13,104,[26,18,14,18,15])
# Independent literal permutation enumeration for H,Q,alpha.
full=set(range(5)); h=q=0; al=[0]*5
for mask in range(32):
 P=[v for v in range(5) if mask>>v&1];R=full-set(P)
 for p in itertools.permutations(P):
  if any(not(adj[u]>>v&1) for u,v in zip(p,p[1:])):continue
  for r in itertools.permutations(R):
   if any(not(adj[u]>>v&1) for u,v in zip(r,r[1:])):continue
   q+=1
   if r:al[r[0]]+=1
   else:h+=1
assert (h,q,al)==(H,Q,alpha)
S0,S1=map(s.Matrix,S);A=S1*S0;M=s.Matrix(matrix(adj,1));G=M-13*S0*S1*S0;F=s.Matrix(T[1])-3*A
zero=s.zeros(6);v=s.Matrix([1,0,0,0,0,0]);ell=s.Matrix([[1,1,1,1,0,0]])
assert A**4-6*A**3+9*A**2-4*A==zero
for j,k in itertools.product(range(4),repeat=2):
 assert int((ell*A**k*S1*G*A**j*v)[0])==4**(j+1)-10+4*(j==0)
# All words, not just a sample: reachable row spaces on right of G,F.
# Verify words up to length7 as finite diagnostics only.
for length in range(1,8):
 for bits in itertools.product([0,1],repeat=length):
  W=s.eye(6)
  for b in bits:W=[S0,S1][b]*W
  for X,Y in itertools.product([F,G],repeat=2):assert X*W*Y==zero
# Arbitrary-word proof: K=image(F)+image(G), its image after one
# singleton lies in common kernel, and span remains invariant.
K=F.row_join(G).columnspace();K=s.Matrix.hstack(*K)
L=(S0*K).row_join(S1*K).columnspace();L=s.Matrix.hstack(*L)
assert F*L==s.zeros(6,L.cols) and G*L==s.zeros(6,L.cols)
assert L.row_join(S0*L).row_join(S1*L).rank()==L.rank()
out={'H':H,'Q':Q,'alpha':alpha,'rank_G':G.rank(),'arbitrary_nonempty_binary_spacer_proved':True,'invariant_space_dimension':L.cols,'formula_grid_verified':True,'example_r7_strong':[13*(4**7+6*7-1)//9+4**(j+1)-10 for j in range(1,5)]}
print(json.dumps(out,indent=2))
