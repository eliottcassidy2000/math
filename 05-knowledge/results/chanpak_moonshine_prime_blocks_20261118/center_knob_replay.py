from unit_gain_replay import K,rev,L,R,state,forward,H,bits
examples=[]
for j in range(1,11):
 r=4*j+2;u=(7*16**j-2)//5
 C=(299*16**j+36-10*j)//5
 for V in (0,1,2,3,(C-6)//2,10**25+3):
  counts=[]
  for pos in (r-3,r+1):
   wl=2**r+2**pos+2**(r+2)*V
   wr=u+2**pos+2**(r+2)*V
   hl=9+23*wl-5*K(wl)-(wl&1)
   hr=9+25*wr-5*K(wr)-(wr&1)
   assert hl-hr==C-6-2*V
   if j==1 and V<=1:
    assert H(forward(L,3,bits(wl)))==hl
    assert H(forward(R,3,bits(wr)))==hr
   counts.append((hl,hr))
  delta=345*2**(r-3)+75*V+10-5*(V&1)
  assert counts[1][0]-counts[0][0]==counts[1][1]-counts[0][1]==delta
  if j==1 and V in (0,1,(C-6)//2):examples.append((V,counts,delta))
print({'center_knob':'PASS','first_examples':examples,'largest_r':42})
