#include "dependencies/counting.hpp"
int main(){I graphs=0,roots=0,pairs=0,eq=0,sharp=0,flipped=0,basic=0;for(int n=1;n<=6;n++)for(int code=0;code<(1<<(n*(n-1)/2));code++){
 G a=graph(n,code);Counts c(a);graphs++;
 for(int x=0;x<n;x++)for(int y=x+1;y<n;y++){
  I covers=0;for(int S=1;S<c.F;S++){int A=c.F^S;covers+=c.st[S*n+x]*c.en[A*n+y]+c.st[S*n+y]*c.en[A*n+x];}
  need(covers>0,"two-terminal cover lemma");pairs++;
 }
 for(int v=0;v<n;v++){
  auto[h,s,e,L,R,d,b,B,D]=c.profile(v);int A=c.F^(1<<v),Iset=A^a[v],degree=__builtin_popcount((unsigned)Iset);I starts=0;
  for(int x=0;x<n;x++)if(Iset>>x&1)starts+=c.st[A*n+x];
  need(L-h>=starts+degree*(degree-1)/2,"degree lower bound");
  bool exact=degree==0||(degree==1&&starts==0);need((L==h)==exact,"exact equality classification");if(degree<=1)need(L-h==starts,"low degree exact formula");
  if(L==h)eq++;
  if(n>=2&&s>0&&D==s){
   sharp++;need(__builtin_popcount((unsigned)a[v])==1,"sharp unique outneighbor");int w=__builtin_ctz((unsigned)a[v]),T=A^(1<<w),incoming=A^(a[w]&A)^(1<<w);
   need(__builtin_popcount((unsigned)incoming)<=1,"sharp normal form degree");
   if(incoming){int x=__builtin_ctz((unsigned)incoming);need(c.st[T*n+x]==0&&c.en[T*n+x]==0,"sharp interior vertex");}
   need(T&&c.H[T]==s&&strong(a),"sharp residual and strongness");
  }
  roots++;
 }
 G u=cyclic(a);Counts cu(u);auto pu=cu.profile(0);need(pu[1]==c.h()&&pu[8]==pu[1]&&strong(u),"basic equality construction");basic++;
 for(int x=0;x<n;x++)if(c.st[c.F*n+x]==0&&c.en[c.F*n+x]==0){G z=u;z[1]^=1<<(x+2);z[x+2]|=2;Counts cz(z);auto pz=cz.profile(0);need(pz[1]==c.h()&&pz[8]==pz[1]&&strong(z),"flipped equality construction");
 I h=c.h(),C=c.C(c.F),L=c.L(x,c.F),R=c.R(x,c.F),E3=0;for(int S=0;S<=c.F;S++)E3+=c.en[S*n+x]*c.C(c.F^S);
 need(pz[0]==h+C+L-R&&pz[2]==h+L&&pz[5]==h+L&&pz[6]==2*C-2*R+2*E3&&pz[7]+h==3*C-3*R+2*E3,"flipped full readout");
 G reach=a;for(int i=0;i<n;i++)reach[i]|=1<<i;for(int k=0;k<n;k++)for(int i=0;i<n;i++)if(reach[i]>>k&1)reach[i]|=reach[k];int P=0,M=0,Q=0;for(int i=0;i<n;i++){if((reach[x]>>i&1)&&(reach[i]>>x&1))M|=1<<i;else if(reach[i]>>x&1)P|=1<<i;else Q|=1<<i;}
 need(P&&Q,"interior marked component");I CP3=0,EM3=0;for(int S=P;;S=(S-1)&P){CP3+=c.H[S]*c.C(P^S);if(!S)break;}for(int S=M;;S=(S-1)&M){EM3+=c.en[S*n+x]*c.C(M^S);if(!S)break;}
 need(C==c.C(P)*c.C(M)*c.C(Q)&&R==c.H[P]*c.R(x,M)*c.C(Q)&&E3==CP3*EM3*c.C(Q),"marked SCC factorization");
 flipped++;}
 }
 cout<<"{\n  \"status\": \"PASS\",\n  \"labeled_tournaments\": "<<graphs<<",\n  \"rooted_inputs\": "<<roots<<",\n  \"terminal_pairs\": "<<pairs<<",\n  \"cover_equalities\": "<<eq<<",\n  \"sharp_equalities\": "<<sharp<<",\n  \"basic_sources\": "<<basic<<",\n  \"flipped_sources\": "<<flipped<<"\n}\n";
}
