"""Complete one-apex clone family, structural halving and ternary carry."""
from functools import lru_cache
from support import need,save,junction,phases,observer
from observer import records

@lru_cache(None)
def B(a):
    if a==0:return 1
    k,bit=divmod(a,2)
    return B(k)-bit*(k+(k&1))

def word(a):return tuple(1-((a>>i)&1) for i in range(a.bit_length()))

def cut_H(bits):
    h,t,l=1,0,0
    for b in bits:
        if b:h,t,l=h,2*t+l,1
        else:h,t,l=h+2*t+2*l,2*t+l,0
    return h

def tournament(a,m):
    bits=(1,)*m+word(a);n=len(bits)+1;w=n-1
    T=[sum(1<<j for j in range(i+1,n-1)) for i in range(n-1)]+[0]
    for i,b in enumerate(bits):
        if b:T[w]|=1<<i
        else:T[i]|=1<<w
    return tuple(T)

def carry(a):
    # State (Delta,E,F,n,parity,1); read the usual binary expansion left to right.
    d=e=f=-2;n=p=0
    for bit in bin(a)[2:]:
        if bit=='0':d,e,f,n,p=d,d-3*n-p,e,2*n,0
        else:d,e,f,n,p=e-2+4*p,f+3*n+3*p,f-2+2*p,2*n+1,1
    return d,e,f,n,p,1

def main():
    N=65536;offsets={};checks=0
    for a in range(1,N):
        direct=cut_H(word(a))-a+(a&1)
        need(B(a)==direct and B(2*a)==B(a),'word formula and structural halving')
        d,e,f,n,p,c=carry(a)
        need((d,e,f)==(B(3*a)-3*B(a),B(3*a+1)-3*B(a),B(3*a+2)-3*B(a)),'closed binary carry recurrence')
        odd=a//(a&-a)
        need((d<=-4 if odd%4==1 else d==0 if odd==3 else d>=4),'proved sign split')
        need(d!=2 and ((d==4)==(odd==7)),'uniform ternary boundary')
        offsets[d]=offsets.get(d,0)+1;checks+=1
    literal=[]
    for a in range(1,16):
        for m in (1,2,3):
            T=tournament(a,m);z=junction.profile(T,0)
            need(junction.strong(T) and z[0]==a*2**m+B(a),'literal family H')
            if m>=2:need(z[3]-z[0]==1,'literal first clone readiness')
            literal.append({'a':a,'m':m,'profile':z});junction.endpoints.cache_clear()
    births=[]
    for m in (1,2,3,5):
        T=tournament(7,m);n=len(T);w=n-1
        # Old Q0,Q1,Q2 become target Q0,Q2,Q4. New Q1,Q3 interleave them.
        first_mask=((1<<m)-1)|(1<<m)|(1<<w)
        U=junction.ear(T,first_mask)
        second_mask=((1<<m)-1)|(1<<m)|(1<<(m+1))|(1<<w)|(1<<n)
        V=junction.ear(U,second_mask)
        labels=list(range(m))+[m,n,m+1,n+1,m+2,w]
        target=tournament(21,m)
        need(junction.sub(V,tuple(labels))==target,'all old names survive two literal joint-cut births')
        need(junction.H(V)==3*junction.H(T)+4,'literal uniform plus-four branch')
        births.append({'m':m,'source_rows':T,'incoming_cuts':[first_mask,second_mask],'final_rows':V,'target_order':labels,
                       'H_source':junction.H(T),'H_target':junction.H(V),'root_profile':junction.profile(V,0),'native_moves':0,'births':2,'old_deletions':0,
                       'root_phases_before':phases(T,0),'root_phases_after':phases(V,0),'spare_name':m,'spare_phases_after':phases(V,m),
                       'two_port_observer_before':records(observer(T,(0,m))),'two_port_observer_after':records(observer(V,(0,m)))})
        junction.endpoints.cache_clear()
    # The root/spare repair's coefficient has an exact, nonvanishing return debt.
    need(B(1243)==-1061 and B(3729)==-1999 and B(3729)-3*B(1243)==1184,'actual repaired-family carry')
    save('apex_carry.json',{'status':'PROVED_ALL_SIZE_WORD_CARRIER_WITH_EXACT_FINITE_CONTROLS','family':'H(A,m)=A*2^m+B(A), A>=1,m>=1; 1-ready at first clone for m>=2',
      'rooted_base_class_series':'x^3/(1-2x)','B_recurrence':{'B0':1,'even':'B(2k)=B(k)','odd':'B(2k+1)=B(k)-k-(k mod2)'},
      'carry_state':['Delta','E','F','n','n_mod2','1'],
      'append0':['Delta','Delta-3n-parity','E','2n','0','1'],
      'append1':['E-2+4parity','F+3n+3parity','F-2+2parity','2n+1','1','1'],
      'all_size_conclusion':'Delta(A) is never2; Delta(A)=4 iff oddpart(A)=7; Delta(A)=0 iff oddpart(A)=3',
      'checked_A_range':[1,N-1],'checked':checks,'literal_small_controls':literal,'uniform_plus4_births':births,
      'repaired_family':{'A':1243,'B':-1061,'target_A':3729,'target_B':-1999,'desired_plus2_B':-3181,'excess_over_plus2':1182},
      'scope':'uniform identities at the same clone parameter, within all one-apex clone families; no exclusion of a single H value or arbitrary tournaments'})
    print('PASS apex carrier: all-size no uniform +2; +4 exactly oddpart(A)=7; explicit two-birth branch',flush=True)
if __name__=='__main__':main()
