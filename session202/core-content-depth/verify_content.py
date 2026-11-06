import itertools,json
from verify import hamilton,glue,need

def core(st,j):
    h,x,y,d=st;hc=(3,5,5,3)[j]
    return (3*h+(9-3*j)*x+(3-j)*y+(hc-3)*d,12*x+3*j*d,j*d,3*d)
def op(st,c):
    h,x,y,d=st
    if c=='O':return(h,2*x+d,d,d)
    if c=='Z':return(h+x+y,2*x,0,d)
    return core(st,int(c))
# Independently count all seven core response constants by permutations/cuts.
for j in range(4):
    rows=glue((2,4,1),1,(1<<j)-1);n=len(rows)
    H=e=start=R=L=0
    def path(p):return all(rows[a]>>b&1 for a,b in zip(p,p[1:]))
    for perm in itertools.permutations(range(n)):
        pos=perm.index(0)
        if path(perm):
            H+=1;e+=perm[-1]==0;start+=perm[0]==0
        if path(perm[:pos]) and path(perm[pos:]):R+=1
        if path(perm[:pos+1]) and path(perm[pos+1:]):L+=1
    need((H,L-3,R-3,3-start,3-e)==((3,5,5,3)[j],9-3*j,3*j,3-j,j),'literal full core constants')

# At most three actual blocks: 1+8+64+512 literal graphs.
literal=0
for m in range(4):
    for mask in range(1<<(3*m)):
        st=(1,0,0,1)
        for i in range(m):st=core(st,((mask>>(3*i))&7).bit_count())
        need(hamilton(glue((2,4,1),m,mask))==st[0],'literal triangle transfer')
        literal+=1
values=[]
for m in range(3):
    vals=set()
    for js in itertools.product(range(4),repeat=m):
        st=(1,0,0,1)
        for j in js:st=core(st,j)
        vals.add(st[0])
    values.append(sorted(vals))
prime_values=sorted(v for v in set().union(*map(set,values)) if v>1 and all(v%d for d in range(2,int(v**.5)+1)))
need(prime_values==[3,5,31,41,61],'complete prime palette')
checks=0
for length in range(8):
    for word in itertools.product('OZ0123',repeat=length):
        st=(1,0,0,1);m=0
        for c in word:
            st=op(st,c);m+=c not in 'OZ'
        h,x,y,d=st
        need(d%3**m==0,'d grade')
        if m>=1:need(x%3**(m-1)==0 and y%3**(m-1)==0,'xy grade')
        if m>=2:need(h%3**(m-2)==0,'h grade')
        need(h>=d,'deletion monotonicity')
        checks+=1
# Literal one-core escape from the five-prime palette: prefix O, then j=1.
# Vertices root0, singleton1, triangle2,3,4; root beats1 and2.
rows=(6,28,8,16,4)
# Ensure missing root incidences and every tournament pair are correct.
rows=(6,28,8,17,5)
need(hamilton(rows)==13,'binary spacer H13 escape')
st=core(op((1,0,0,1),'O'),1);need(st[0]==13,'escape transfer')
print(json.dumps({'status':'PASS','literal_graphs':literal,'independent_full_core_profiles':4,'complete_transfer_words_length_at_most':7,'transfer_word_count':checks,'value_sets_m_0_1_2':values,'prime_palette':prime_values,'singleton_spacer_escape':{'rows':rows,'H':13}},indent=2,sort_keys=True))
