"""A count25 core replacement and its entry into the inherited infinite renewal."""
from math import gcd
from atomic_ears import profile,pure_ends
from peeling import induced,strong
from triangle_probe import need,save

def birth(T,r):
    n=len(T);out=list(T)+[((1<<n)-1)^(1<<r)];out[r]|=1<<n;return tuple(out)

def renewal(T,r):
    z=profile(T,r);b,t=z[:2];need(z==(b,t,t,b),'balanced entry')
    T=birth(T,r);r=len(T)-1;children=[]
    for _ in range(3):T=birth(T,r);children.append(len(T)-1)
    selected=children[-1];h_before=sum(pure_ends(T)[-1]);U=list(T)
    for old in children[:-1]:U[selected]^=1<<old;U[old]^=1<<selected
    T=tuple(U);h_after=sum(pure_ends(T)[-1])
    need(h_before==h_after==9*b,'literal neutral endpoint preparation')
    T=birth(T,selected);r=len(T)-1;need(sum(pure_ends(T)[-1])==34*b-t,'selected-leaf birth')
    T=birth(T,r);output=profile(T,r)
    need(output==(61*b-2*t,9*b,9*b,61*b-2*t),'literal renewed interface')
    return T,r,dict(input=z,prepared_count=h_after,output=output)

def run():
    old=(2,8,3,69,15,31,55);raw=(2,36,33,7,15,25)
    # Physical old names after deletion are 0,1,2,3,4,6. The physical root6
    # and incoming port3 are retained. Raw root5 is sent to physical6.
    new=induced(raw,[0,2,3,1,4,5]);names=[0,1,2,3,4,6]
    before=profile(old,6);after=profile(new,5)
    need(before==(25,5,5,25) and after==(25,3,3,27),'exact neutral core replacement')
    need(strong(old) and strong(new),'both cores strong')
    need(old[3]>>6&1 and new[3]>>5&1,'same physical port arc')
    flips=[]
    for i,a in enumerate(names):
        for j,b in enumerate(names[i+1:],i+1):
            if ((old[a]>>b)&1)!=((new[i]>>j)&1):flips.append([a,b])
    T=birth(new,5);r=len(T)-1;T=birth(T,r)
    need(profile(T,r)==(79,25,25,79),'two-birth entry')
    U,s,receipt=renewal(T,r);need(strong(U),'strong renewal')
    b,t=79,25;states=[]
    for j in range(25):
        need(gcd(b,t)==1,'primitive exact recurrence')
        if j:need(t%3==0 and b%3!=0,'invariant renewal guard')
        states.append([b,t]);b,t=61*b-2*t,9*b
    save('neutral_entry',dict(scope='Literal nonisomorphic H25 core replacement, named root/port retention, two-birth H79 entry and one full14-vertex renewal; 25 arithmetic iterates are controls, not the all-order proof',old_rows=old,old_root=6,old_port=3,new_rows=new,new_physical_names=names,deleted_physical_name=5,flipped_pairs=flips,old_profile=before,new_profile=after,H79_rows=T,H79_root=r,renewal=receipt,H4769_rows=U,H4769_root=s,recurrence_states=states))
    print('NEUTRAL',before,'->',after,'flips',flips)
    print('RENEWAL',receipt,'first counts',[v[0] for v in states[:5]])

if __name__=='__main__':run()
