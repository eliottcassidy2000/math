from itertools import combinations
def need(x):
 if not x:raise AssertionError
def add(x,y,p):return tuple((a+b)%p for a,b in zip(x,y))
def norm(x,p):return sum(a*a for a in x)%p
def main():
 for p in (3,7,11,19):
  pts=[(a,b) for a in range(p) for b in range(p)]
  choice={norm(z,p):z for z in pts}
  need(len(choice)==p)
  cyc=[choice[-2*i%p] for i in range(p)]
  need(all((norm(y,p)-norm(x,p))%p==p-2 for x,y in zip(cyc,cyc[1:]+cyc[:1])))
  by={}
  for d in pts[1:]:
   j={z for z in pts if (norm(add(z,d,p),p)-norm(z,p))%p==p-2}
   need(len(j)==p and all(add(z,d,p) not in j for z in j))
   by[d]=j
  for d,e in combinations(pts[1:],2):
   if (d[0]*e[1]-d[1]*e[0])%p:
    need(sum(add(z,d,p) in by[e] for z in by[d])==1)
 for p in (5,13,17):
  t=next(t for t in range(p) if (t*t+1)%p==0)
  x=(p-1,0);d=(1,t)
  for i in range(p):
   y=add(x,d,p);need((norm(y,p)-norm(x,p))%p==p-2);x=y
 print("Inert affine-line intersections and full norm cycles, and split one-letter cycles verified.")
if __name__=='__main__':main()
