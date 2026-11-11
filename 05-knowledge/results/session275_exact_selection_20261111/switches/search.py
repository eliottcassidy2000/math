from carrier import *
lookup={}
for n in range(2,7):
 count=0
 for T in tournaments(n):
  p=profile(T); hp.cache_clear()
  key=(n,)+tuple(p[k] for k in ('H','h','Q','V','e'))
  if key in lookup:
   B,pb=lookup[key]
   if (p['U'],p['d'])!=(pb['U'],pb['d']):
    print(json.dumps(dict(n=n,A=B,B=T,pA=pb,pB=p)),flush=True);quit()
  else:lookup[key]=(T,p)
  count+=1
 print('finished',n,count,flush=True)
