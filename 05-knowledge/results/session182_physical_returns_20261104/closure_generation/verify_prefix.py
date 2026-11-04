import itertools,json,pathlib,collections
D=pathlib.Path(__file__).parent
def birth(r,c):return [v|(((c>>i)&1)<<len(r))for i,v in enumerate(r)]+[((1<<len(r))-1)^c]
def response(r):
 n=len(r);h=p=0
 for perm in itertools.permutations(range(n)):
  ok=[bool(r[u]>>v&1)for u,v in zip(perm,perm[1:])];h+=all(ok)
  for k in range(n+1):p+=all(ok[:max(0,k-1)])and all(ok[k:])
 return h,p
menus={4:{},5:{}}
for a in range(8):
 r=birth([6,4,0],a);menus[4].setdefault(response(r),[]).append([a])
 for b in range(16):menus[5].setdefault(response(birth(r,b)),[]).append([a,b])
if set(menus[4])!={(1,16),(3,24),(5,32)}:raise RuntimeError('Order4 response menu mismatch')
if set(menus[5])!={(1,32),(3,48),(5,64),(9,84),(11,100),(13,104),(15,108),(15,120)}:raise RuntimeError('Order5 response menu mismatch')
z={str(n):[dict(H=h,P2=p,cuts=cuts)for (h,p),cuts in sorted(menu.items())]for n,menu in menus.items()};(D/'prefix_verification.json').write_text(json.dumps(dict(status='PASS',menus=z),indent=2)+'\n');print('PASS: 8 one-birth and128 ordered two-birth prefix extensions, by direct permutations')
