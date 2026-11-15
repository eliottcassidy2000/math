from itertools import combinations
from replay_slot_spectral import counts,covers,strong

def main():
    prof=set()
    checked=0
    for n in range(3,6):
        ed=list(combinations(range(n),2))
        for bits in range(1<<len(ed)):
            r=[0]*n
            for j,(a,b) in enumerate(ed):
                if bits>>j&1:r[a]|=1<<b
                else:r[b]|=1<<a
            if not strong(r):continue
            checked+=1
            h=counts(r)
            if h[-1]==9:prof.add((n,covers(h,2)))
    if prof!={(5,84)}:raise AssertionError(prof)
    print('strong labelled tournaments through order 5:',checked)
    print('fixed H=9, strong profiles:',sorted(prof),'; two triangles give 144')
    print('fixed H=9, exact reduced odd part bank:',[9,21])
if __name__=='__main__':main()
