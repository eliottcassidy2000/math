import json
bank={};cases=0;hits={}
for i,line in enumerate(open('RESET_SOURCES.txt')):
 v=list(map(int,line.split()));h,s,e,L,R,d,C,b,T=v[3:];r=d%3;j=(d-r)//3
 p=(h+3*L-s-2*d+2*j,2*(L-d),2*b,R+3*b-e,d-e)
 def walk(p,word):
  global cases
  cases+=1;slots=bank.setdefault(p,[None,None]);slots[r-1]=(i,word)
  if all(x is not None for x in slots):hits[p]=slots[:]
  if len(word)==4:return
  X,a,b,w,t=p
  walk((X+2*a,2*a,2*b,w+2*b,t),word+'Q')
  walk((X+w,a+b,2*b,2*w-2*t,0),word+'Z')
 walk(p,'')
out=[{'key':p,'recipes':z} for p,z in sorted(hits.items())];json.dump({'cases':cases,'keys':len(bank),'hits':out},open('RESET_BINARY_EXTENDED.json','w'),indent=2);print('cases',cases,'keys',len(bank),'hits',len(out));print(out[:3])
