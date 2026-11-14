import sympy as s
J=s.Matrix([[1,-1,1,0],[1,-1,0,0],[1,0,1,0],[1,0,0,0]])
O=s.Matrix([[1,0,0,0],[0,0,0,0],[0,0,2,0],[0,0,0,1]])
Z=s.Matrix([[1,-1,1,0],[0,0,0,1],[0,0,2,-1],[0,0,0,1]])
words='e O OO OZ OOO OOZ OZO OZZ OOOZ OZOZ OZZZ'.split()
G={}
for w in words:
 P=J
 for c in w.replace('e',''): P=(O if c=='O' else Z)*P
 G[w]=P[:3,:3]
if __name__=='__main__':
 for w,P in G.items(): print(w,list(P),P.det(),P.eigenvals() if P.rank()<3 else P.charpoly().as_expr())
T=s.Matrix([[1,-1,0],[0,1,0],[0,-1,1]])
assert all(all(x>=0 for x in T*M*T.inv()) for M in G.values())
pairs=[(0,2),(0,1),(1,2)]
wedge=lambda M: s.Matrix([[M.extract(i,j).det() for j in pairs] for i in pairs])
for k in range(1,4):
 assert list(wedge(G['O'*k])[0,:])==[0,2**k,-2**k]
assert list(wedge(G['OZO'])[0,:])==[6,6,-8]
examples=[
 ('O',1,2,1,1),
 ('e e OOZ e e OZO',55,138,7,9),
 ('e OZZZ OOZ e OZO',65,162,51,52),
 ('OOZ OZOZ e e OZO',71,178,41,45),
 ('OZ e OZZ OOOZ O',51,104,65,81)]
for word,h,l,u,v in examples:
 P=s.eye(3)
 for w in word.split():
  P=G[w]*P
 assert P.extract([0,2],[0,2])==s.Matrix([[h*u,h*v],[l*u,l*v]])
 assert list(P.extract([1],[0,2]))==[0,0]
 assert s.gcd(h,l)==1
 # Independently apply the generating state updates as in runterm.c.
 for h0,l0 in ((1,2),(15,18),(173,256)):
  H,E,L=h0,0,l0
  for w in word.split():
   old=H
   H,E,L=H-E+L,H-E,H+L
   for c in w.replace('e',''):
    if c=='O': E,L=0,2*L
    else: H,E,L=H+L-E,old,2*L-old
  assert (H,E,L)==((u*h0+v*l0)*h,0,(u*h0+v*l0)*l)
print('PASS: eleven terminal matrices, cone, compounds and five literal separators')
