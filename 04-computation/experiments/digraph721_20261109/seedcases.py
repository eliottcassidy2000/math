"""Check of the seed-lemma case analysis (PROOF.md Lemma 6, cases A1/A2/A3/B1/B3) for 13 <= L(j) <= 15 against the
Gamma data dumped by evenclos (dump mode).  usage: python3 seedcases.py gamma_dump.out"""
import sys
Gam={}
import sys
for line in open(sys.argv[1]):
    p=line.split()
    if not p or p[0]!='G': continue
    c=int(p[1]); Gam[c]=set(int(x) for x in p[2:])
def M(n): return 2**n-n if n%2 else 2**n-n-1
def inE(k,n):
    if k<1 or k%2==0 or k>M(n): return False
    if k in (3,7,13): return False
    if n%2==1 and n>=5 and k==2**n-n-2: return False
    return True
def L(k): return k.bit_length()
cnt={}
bad=0
for J in range(13,16):
    for j in range(2**(J-1)+1,2**J,2):
        if j>=3*2**(J-2):
            x=2**J-j
            if inE(x,J-2): case='A1'; ok=True
            elif x in (3,7,13): case='A2'; ok=(J+1) in Gam[x]
            else:
                i=2**(J-2)-x; case='A3'; ok=(1<=i<=J) and ((J-1) in Gam[-i])
        else:
            y=2**J-j; z=y+2
            if inE(z,J-1): case='B1'; ok=True
            else:
                i=j-2**(J-1); case='B3'; ok=(1<=i<=J+3) and (J in Gam[-i])
        cnt[case]=cnt.get(case,0)+1
        if not ok: bad+=1; print('FAIL',J,j,case)
print(cnt,'bad',bad)
