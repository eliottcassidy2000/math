"""Literal count-preserving module trade. No birth-only closure is claimed."""
from core import inventories, require, adjacency_map

def construction(t,j,m):
 require(t>=1 and j>=1 and m>=1,'positive module sizes')
 names=['a']+[('u',i)for i in range(t)]+['b','x0','x1']+[('z',i)for i in range(j)]+['y0','y1']+[('r',i)for i in range(m)]
 n=len(names);T=[0]*n;cut=t+2;R0=t+j+6
 for x in range(n):
  for y in range(x+1,n):
   if x<cut and y<cut:forward=True
   elif x<cut and y<R0:
    forward=(names[x]=='a'and names[y]in('x0','x1'))or(names[x]=='b'and names[y]in('y0','y1'))
   elif x<cut:forward=False
   else:forward=True
   T[x if forward else y]|=1<<(y if forward else x)
 return T,dict(a=0,U=list(range(1,t+1)),b=t+1,X=[t+2,t+3],Z=list(range(t+4,t+4+j)),Y=[t+4+j,t+5+j],R=list(range(R0,n)),root=n-1)

def hp(t,j,m):return 30*2**(t+j)*3**m+1

def ending(t,j):return 30*2**j-3

def cover(t,j):return (90*3**t+24*2**t-27)*2**j+6*3**t+1

def zcount(t,j):return (90*3**t+114*2**t-57)*2**j+6*3**t+5

def append(T,selected,left=False):
 n=len(T);A=list(T)+[0]
 for v in range(n):
  new_wins=(v not in selected)if left else(v in selected)
  A[n if new_wins else v]|=1<<(v if new_wins else n)
 return A

def delete(T,v):
 keep=[i for i in range(len(T))if i!=v]
 return [sum(((T[x]>>y)&1)<<k for k,y in enumerate(keep))for x in keep],keep

def move_right(t,j,m):
 """Delete the first Z, append a vertex beating only b, retain old root."""
 require(j>=2,'keep Z nonempty');T,p=construction(t,j,m);v=p['Z'][0]
 A,keep=delete(T,v);A=append(A,{keep.index(p['b'])})
 order=[p['a']]+p['U']+[v,p['b']]+p['X']+p['Z'][1:]+p['Y']+p['R']
 f=[order.index(x)for x in keep+[v]]
 B,q=construction(t+1,j-1,m)
 require(adjacency_map(A,B,f),'right replacement graph identity')
 ports=[p['root'],p['b']]+p['X']
 require([f[keep.index(x)]for x in ports]==[q['root'],q['b']]+q['X'],'right persistent ports')
 return A,f

def move_left(t,j,m):
 """Delete the last U, append a vertex losing only to x0,x1."""
 require(t>=2,'keep U nonempty');T,p=construction(t,j,m);v=p['U'][-1]
 A,keep=delete(T,v);A=append(A,{keep.index(x)for x in p['X']},True)
 order=[p['a']]+p['U'][:-1]+[p['b']]+p['X']+[v]+p['Z']+p['Y']+p['R']
 f=[order.index(x)for x in keep+[v]]
 B,q=construction(t-1,j+1,m)
 require(adjacency_map(A,B,f),'left replacement graph identity')
 ports=[p['root'],p['b']]+p['X']
 require([f[keep.index(x)]for x in ports]==[q['root'],q['b']]+q['X'],'left persistent ports')
 return A,f

def cross_edges(T,v):
 O=[u for u in range(len(T))if T[v]>>u&1]
 I=[u for u in range(len(T))if T[u]>>v&1]
 return [(a,b)for a in O for b in I if T[a]>>b&1]

def cover_check(T,v,C):
 return all(a in C or b in C for a,b in cross_edges(T,v))

def strong_cycle(t,j,m,p):return [p['a']]+p['X']+p['Z']+p['U']+[p['b']]+p['Y']+p['R']
