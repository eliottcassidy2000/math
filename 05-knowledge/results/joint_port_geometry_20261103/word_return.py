"""Four-edge word repair family and its next-step stopping boundary."""
import json
from return_probe import HERE,old,junction,need,family

def repair(w,u,j):
    m=len(w);T=old.source(w,u);a,b,ell=m,m+1,m-1
    U=list(junction.ear(junction.ear(T,1<<a),1<<(m+2)));x,y=m+2,m+3
    for v,z in ((y,b),(a,ell),(a,x),(x,ell)):
        U[v]^=1<<z;U[z]^=1<<v
    order=(y,a)+tuple(i for i in range(m) if i!=j)+(x,b)
    target=junction.sub(tuple(U),order);uu=(1,0)+u[:j]+u[j+1:];ww=(1,)*(len(uu)-1)+(0,)
    need(target==old.source(ww,uu),'literal word transform u -> 10 delete_j(u)')
    return ww,uu

def main():
    successes=[];tested=0;literal=0
    for p in range(1,10):
        for q in range(10):
            w,u=family(p,q);h=old.gap_h(w,u)
            for j in range(len(u)-1):
                ww=(1,)*len(u)+(0,);uu=(1,0)+u[:j]+u[j+1:];out=old.gap_h(ww,uu);tested+=1
                if p<=3 and q<=2:
                    need(repair(w,u,j)==(ww,uu),'uniform word formula literal test');literal+=1
                if out==3*h+2:successes.append({'p':p,'q':q,'deleted_letter':j,'H':h,'target':out,'u':uu})
    need(successes==[{'p':2,'q':1,'deleted_letter':3,'H':77,'target':233,'u':(1,0,1,1,0,0)}],'complete stated one-hole finite bank')
    w=(1,1,1,1,1,0);u=(1,0,1,1,0,0)
    next_values=[old.gap_h(*repair(w,u,j)) for j in range(5)]
    need(next_values==[571,601,673,673,377] and 701 not in next_values,'same four-edge return family is not closed')
    (HERE/'word_return.json').write_text(json.dumps({'status':'FINITE_EXACT_WITH_UNIFORM_WORD_IDENTITY',
        'one_hole_range':{'p':[1,9],'q':[0,9]},'tested':tested,'literal_checks':literal,'successes':successes,
        'next_source_u':u,'next_required_H':701,'next_target_values':next_values},indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS: uniform four-edge word transform; 77 success; next701 outside this deletion bank',flush=True)
if __name__=='__main__':main()
