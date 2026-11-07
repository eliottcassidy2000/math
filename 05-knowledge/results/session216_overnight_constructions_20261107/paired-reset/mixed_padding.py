import json
hit=json.load(open('RESET_PAIRED_RAYS.json'))[0]
def initial(r):
 h,s,e,L,R,d,C,b,T=r['profile'];return(h+L-s,d,e,2*L-d,R+b,d,2*b)
def move(p,c):
 h,s,e,L,R,d,b=p
 return {'O':(h,s,0,L,2*R,d,2*b),'Z':(h+R-e,s,d,L+b,2*R-d,d,2*b),'P':(h,0,e,2*L,R,d,2*b),'Q':(h+L-s,d,e,2*L-d,R+b,d,2*b)}[c]
def key(p):
 h,s,e,L,R,d,b=p
 if s!=d:return None
 r=d%3;j=(d-r)//3
 return(h+L-d+2*j,L-d,b,R+b-e,d-e)
A={initial(hit['A']):''};B={initial(hit['B']):''};out=[]
for depth in range(11):
 AA={key(p):(p,w) for p,w in A.items() if key(p) is not None};BB={key(p):(p,w) for p,w in B.items() if key(p) is not None}
 common=AA.keys()&BB.keys()
 for k in common:out.append({'depth':depth,'key':k,'A':AA[k],'B':BB[k]})
 print('depth',depth,'states',len(A),len(B),'hits',len(common),flush=True)
 if common:break
 if depth==10:break
 A={move(p,c):w+c for p,w in A.items() for c in 'OZPQ'};B={move(p,c):w+c for p,w in B.items() for c in 'OZPQ'}
json.dump(out,open('MIXED_PADDING_MATCHES.json','w'),indent=2)
