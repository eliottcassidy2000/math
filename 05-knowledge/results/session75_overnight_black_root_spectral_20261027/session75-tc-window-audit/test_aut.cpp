#define main scan_main
#include "tc_window_scan_portable.cpp"
#undef main
#include <algorithm>
static long long brute(){vector<int> p(n);for(int i=0;i<n;i++)p[i]=i;long long ans=0;do{bool ok=true;for(int i=0;i<n&&ok;i++)for(int j=0;j<n;j++)if(((outm[i]>>j)&1)!=((outm[p[i]]>>p[j])&1)){ok=false;break;}ans+=ok;}while(std::next_permutation(p.begin(),p.end()));return ans;}
static void need(bool ok,const char *s){if(!ok){fprintf(stderr,"FAIL %s\n",s);exit(2);}}
int main(){long long total=0;for(n=1;n<=6;n++){int pairs=n*(n-1)/2;for(int bits=0;bits<(1<<pairs);bits++){memset(outm,0,sizeof outm);int b=0;for(int i=0;i<n;i++)for(int j=i+1;j<n;j++,b++){if(bits>>b&1)outm[i]|=1u<<j;else outm[j]|=1u<<i;}need(autOrder()==brute(),"small brute aut mismatch");total++;}}
 for(int m:{3,5,7}){n=m;memset(outm,0,sizeof outm);for(int i=0;i<n;i++)for(int d=1;d<n;d++)if((m==7?(d==1||d==2||d==4):(d<=m/2)))outm[i]|=1u<<((i+d)%n);need(autOrder()==(m==7?21:m),"cyclic block aut");}
 for(int outer:{3,5}){int inner=outer==3?5:3;n=outer*inner;memset(outm,0,sizeof outm);for(int u=0;u<n;u++)for(int v=0;v<n;v++){int a=u/inner,b=v/inner,x=u%inner,y=v%inner;bool arc=a==b?((y-x+inner)%inner>0&&(y-x+inner)%inner<=inner/2):((b-a+outer)%outer<=outer/2);if(u!=v&&arc)outm[u]|=1u<<v;}long long expected=outer;for(int i=0;i<outer;i++)expected*=inner;need(autOrder()==expected,"wreath aut");}
 printf("PASS: exact orbit-stabilizer equals brute permutation enumeration on %lld labelled tournaments through n=6; C3/R5/P7 and C3[R5]/R5[C3] controls pass.\n",total);
}
