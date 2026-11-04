"""Exhaustive small literal controls for the feedback-word transfer engine."""
from itertools import product
import subprocess
from support import HERE,save,need,junction
from feedback_search import tournament

def main():
    exe=HERE/'two_feedback.exe'
    subprocess.run(['g++','-std=c++20','-O3',str(HERE/'two_feedback.cpp'),'-o',str(exe)],check=True)
    cases=[(e,w) for n in range(6) for w in product(range(4),repeat=n) for e in range(2)]
    request=''.join(str(e)+' '+(''.join(map(str,w)) or '-')+'\n' for e,w in cases)
    p=subprocess.run([str(exe),'-1'],input=request,text=True,capture_output=True,check=True)
    checks=0;parity=0;strong=0
    for (e,w),line in zip(cases,p.stdout.splitlines(),strict=True):
        T=tournament(w,e,m=0);h=junction.H(T)
        need(int(line)==h,'literal subset DP vs 36-state word transfer')
        need(h%2==1,'Boolean tournament path parity');parity+=1
        strong+=junction.strong(T);checks+=1
        junction.endpoints.cache_clear()
    save('audit_feedback.json',{'status':'FINITE_EXACT_INDEPENDENT_CONTROL','all_four_letter_words_through_length':5,'feedback_edge_orientations':2,'literal_cases':checks,'parity_checks':parity,'strong_cases':strong,
        'scope':'independent subset DP on actual Boolean relations; search pruning separately justified by canonical insertion of a new vertex'})
    print('PASS feedback automaton',checks,'literal cases',flush=True)
if __name__=='__main__':main()
