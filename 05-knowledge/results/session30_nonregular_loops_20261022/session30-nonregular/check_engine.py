import numpy as np
from decorated_engine import cores,symbolic_kernel,evaluate
rng=np.random.default_rng(1137)
for n in [4,5,8]:
 A=np.triu(rng.choice([-1,1],(n,n)),1);M=(A-A.T)/n;I=np.eye(n);one=np.ones(n);q=M@M
 G=-q@np.linalg.inv(I+q);H=q@np.linalg.inv(I-q);O=M@np.linalg.inv(I-q);d=np.diag(H-G);g=G@one;h=M@np.linalg.solve(I+q,one);de=M@one
 expected={(-1,2):-de@np.linalg.solve(I+q,de),(0,1):-2*d@g,(-1,3):np.sum(-2*g*h*h-2/3*g**3),(1,0):np.sum(2*np.diag(G)*np.diag(H)-1.5*np.diag(G)**2)+np.sum(G*(H*H-O*O-2/3*G*G))-d@G@d}
 for k,t in expected.items():
  val=sum(evaluate(ed,symbolic_kernel(ed,a),M) for ed,a in cores(*k))
  print(n,k,'actual',val,'expected',t,'error',val-t,flush=True)
  if abs(val-t)>1e-10:raise RuntimeError('mismatch')
