from check_automaton import *
A=[[False,True],[False,False]]
phi=[3,1,2,0]
def star(w):return [phi[x] for x in w[::-1]]
w=[1,2,0,2];control=0
for m in range(1,13):
 h=automaton_count(2,A,w)
 assert h==15*2**m-1
 if not control:w=star(w+[2]);control=1
 else:w=star(w)+[2];control=0
 assert automaton_count(2,A,w)==2*h+1
print({'controlled_converse_plus_digits':12,'passed':True})
init,acc,tr=compile_automaton(2,A)
def counts(prefix,b,n):
 v=init
 for c in prefix:v=step(v,tr[c])
 out=[]
 for _ in range(n):
  out.append(sum(x for x,a in zip(v,acc) if a));v=step(v,tr[b])
 return out
h=counts([1,2,0,2],2,10)
assert [h[i+1]-2*h[i]-1 for i in range(9)]==[0]*9
bad=counts([3,1,2,0],0,10);c=bad[1]-2*bad[0]
res=[bad[i+1]-2*bad[i]-c for i in range(9)]
assert res[0]==0 and res[1]==72
print({'good_cycle_nine_residuals':[0]*9,'false_scalar_return_first_carry':c,'false_cycle_residuals':res})
