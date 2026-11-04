#define main inherited_funnel_main
#include "funnel_bank.cpp"
#undef main
#include <bit>
int main(){int n,J,id;while(cin>>n>>J>>id){vector<unsigned>T(n),pred(n);for(auto&x:T)cin>>x;int full=(1<<n)-1;
 for(int i=0;i<n;i++)pred[i]=full^T[i]^(1<<i);
 vector<I>e((full+1)*n),s((full+1)*n),h(full+1),E(full+1),S(full+1);h[0]=1;
 for(int m=1;m<=full;m++)for(unsigned bits=m;bits;bits&=bits-1){int v=countr_zero(bits),r=m^(1<<v);I a=!r,b=!r;
  for(unsigned z=r&pred[v];z;z&=z-1)a+=e[r*n+countr_zero(z)];
  for(unsigned z=r&T[v];z;z&=z-1)b+=s[r*n+countr_zero(z)];
  e[m*n+v]=a;s[m*n+v]=b;h[m]+=a;if(J>>v&1){E[m]+=a;S[m]+=b;}}
 I P2=0,F=0,G=0,C=0;for(int m=0;m<=full;m++){P2+=h[m]*h[full^m];F+=E[m]*h[full^m];G+=S[m]*h[full^m];C+=E[m]*S[full^m];}
 if(S[full]!=2)throw runtime_error("core start budget");
 for(int t=0;t<=17-n;t++){int words=t?1<<(t-1):1;
  for(int w=0;w<words;w++){int bits=t?w|(1<<(t-1)):0;auto c1=tail(t,bits,1);auto c2=tail(t,bits,2);
   I parent=c2[0]*G+c1[2]*h[full]-c1[1]*E[full]-C;I Fparent=c2[0]*P2-F;I child=parent+Fparent-2;
   if(child==18977||child==18979)cout<<"{\"core_id\":"<<id<<",\"n\":"<<n+t+1<<",\"target\":"<<child<<",\"h\":"<<parent<<",\"F\":"<<Fparent<<",\"tail_length\":"<<t<<",\"tail_bits\":"<<bits<<",\"joint_C\":"<<C<<"}\n";
  }
 }
 } }
