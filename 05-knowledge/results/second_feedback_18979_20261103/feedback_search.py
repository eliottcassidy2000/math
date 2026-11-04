"""Two-feedback Boolean word search for exact dyadic coefficients."""
import subprocess
from support import HERE,save,need,junction,clone,matching

def tournament(word,edge,m=1):
    letters=(1,)*m+tuple(word);L=len(letters);T=[sum(1<<j for j in range(i+1,L)) for i in range(L)]+[0,0]
    T[L if edge else L+1]|=1<<(L+1 if edge else L)
    for i,letter in enumerate(letters):
        for j in range(2):
            if letter>>j&1:T[L+j]|=1<<i
            else:T[i]|=1<<(L+j)
    return tuple(T)

def main():
    exe=HERE/'two_feedback.exe';subprocess.run(['g++','-std=c++20','-O3',str(HERE/'two_feedback.cpp'),'-o',str(exe)],check=True)
    p=subprocess.run([str(exe),'13','100'],capture_output=True,text=True,check=True)
    hits=[]
    for line in p.stdout.splitlines():
        edge,w=line.split();edge=int(edge);word=tuple(map(int,w));T=tournament(word,edge);U=tournament(word,edge,2)
        need((junction.H(T),junction.H(U))==(4277,11735),'independent literal pair of coefficients')
        strong=junction.strong(T);M,C=matching(T,0);need(len(M)<=1,'unique in-port clone boundary')
        profiles=[junction.profile(T,r) for r in range(len(T))]
        hits.append({'word':word,'edge':edge,'rows':T,'strong':strong,'matching':M,'cover':C,'profiles':profiles})
        junction.endpoints.cache_clear()
    save('feedback_search.json',{'status':'FINITE_EXACT_BOUNDED_WORD_SEARCH','word_bound':13,'hit_cap':100,'counter':p.stderr.strip(),'hits':hits,
      'scope':'two feedback vertices, uniform clone letter1, all four-letter tails up to13; stops at100 coefficient hits; all hits checked as actual tournaments'})
    print('PASS two-feedback coefficient search',p.stderr.strip(),'strong',sum(d['strong'] for d in hits),flush=True)
if __name__=='__main__':main()
