from carrier import *
D={}
for T in tournaments(6):
 p=profile(T);hp.cache_clear()
 if tuple(p[k] for k in ('h','Q','V','e'))==(5,64,32,0):D.setdefault((p['H'],p['U'],p['d']),T)
print(json.dumps([(a,b) for a,b in sorted(D.items())]))
