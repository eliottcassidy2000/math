from check_motion import tournament,TT,path_dp,lift
from two_live import mask_shelf,sweep
def verify():
 D=tournament(5,955);assert all(D[u][v]==D[4-v][4-u] for u in range(5) for v in range(5))
 h,q,rows=mask_shelf(D); assert (h,q)==(13,104)
 def v(mask):return {k for k in range(5) if mask>>k&1}
 def rv(mask):return {4-k for k in v(mask)}
 def hi(mask):E,L,*_=rows[mask];return h*E+q*L
 def norm(mask):E,L,J,T,S=rows[mask];return 2*h*J+L**2+T**2
 assert hi(21)==hi(7)==6786 and (norm(21),norm(7))==(4104,4102)
 vals={h*h-norm(m) for m in range(32)}
 assert any(x+2 in vals for x in vals)
 assert all(not (x+2 in vals and x+4 in vals) for x in vals)
 for c,a,b in [(0,1,1),(0,1,2),(1,1,1),(2,1,1)]:
  heights=[]
  for m in [21,7]:
   args=[TT(a),D]+([TT(c)]if c else [])+[D,TT(b)]
   masks=[set(),v(m)]+([set()]if c else [])+[rv(m),set()]
   E,L,J,T,S=rows[m]
   pred=q*q*2**(a+b+c)+(h*E+q*L*2**c)*(2**a-2**b)+h*h-2*h*J-2**c*L*L-(T*T if c==0 else 0)
   H=path_dp(lift(args,masks));assert H==pred;heights.append(H)
  assert heights[1]-heights[0]==(2 if c==0 else 104-117*2**c)
 print("Mirrored shared-factor unit and separator failure verified on 8 literal strong head tournaments.")
if __name__=="__main__":verify()
