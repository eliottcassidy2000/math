import json
import numpy as np
from check_rank_fan import HERE,A,h,need
LIMIT=2**63-1

def mul(X,Y):
 need(np.min(X)>=0 and np.min(Y)>=0)
 need(X.shape[-1]*int(np.max(X))*int(np.max(Y))<=LIMIT)
 return X@Y

def op(w):
 z=np.eye(len(h),dtype=np.int64)
 for a in w:z=mul(A[a],z)
 return z

def main():
 D=json.loads((HERE/'return_automaton.json').read_text())
 W=[tuple(w)+(0,3,5,2,4) for w in D['words']]
 need(len(W)==len(set(W))==168)
 M=np.stack([op(w) for w in W])
 v=M[:,:,0].T
 f=h.reshape((1,-1))
 r=np.stack([mul(f,X)[0] for X in M])
 S=(int(mul(f,v).min()),int(mul(r,v).min()),min(int(mul(mul(r,X),v).min()) for X in M))
 need(S==(55,27143,13460883))
 print('Exact return-block minima over 1,2,3 blocks:',S)
if __name__=='__main__':main()
