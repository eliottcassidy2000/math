"""Exact C1 inverse by a strictly increasing physical source readout."""
from array import array
from core import constants

ONE=(1,1,1,1)
def readout(s):
 D,E,U,d=s;return D-E+U

def inventory(cap,keep_states=False):
 if cap<1:return {},set(),0
 seen={ONE};queue=[ONE];parents=array('i',[-1]);labels=array('i',[-1]);support={}
 for i,(D,E,U,d) in enumerate(queue):
  h=D-E+U
  if h>1:support.setdefault(h,i)
  successors=[(0,(D,0,2*U,d))]
  for L in range(1,cap//U+1):
   g,Q,B=constants(L)
   successors.append((L,(g*D-(L%2)*E+L*U,D-int(L==1)*E,Q*D+B*U,D)))
  for label,s in successors:
   if readout(s)<=cap and s not in seen:
    seen.add(s);queue.append(s);parents.append(i);labels.append(label)
 def recipe(i):
  out=[]
  while parents[i]>=0:
   L=labels[i];out.append(('O',None) if L==0 else ('M',L));i=parents[i]
  return tuple(reversed(out))
 result={h:dict(source=queue[i],recipe=recipe(i))for h,i in support.items()}
 return result,seen if keep_states else None,len(queue)
