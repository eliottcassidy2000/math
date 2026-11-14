#include <array>
#include <vector>
#include <set>
#include <iostream>
#include <string>
using namespace std;
using V=array<long long,36>;
vector<int> go[4][36];
V push(V x,int d){V y{};for(int i=0;i<36;i++)if(x[i])for(int j:go[d][i])y[j]+=x[i];return y;}
V pull(V x,int d){V y{};for(int i=0;i<36;i++)for(int j:go[d][i])y[i]+=x[j];return y;}
long long dot(V a,V b){long long ret=0;for(int i=0;i<36;i++)ret+=a[i]*b[i];return ret;}
V endpoint(int e){V y{};for(int r=0;r<2;r++)
  for(int a=0;a<2;a++)for(int c=0;c<2;c++){
   y[18*r+6*a+2+c]=1;
   if(e!=r)y[18*r+6*a+c]=1;
  }
 return y;
}
#include <algorithm>
#include <map>
#include <cstdint>
#include <cassert>
using U=uint64_t;
int target,m_limit;
bool intercyclic=false;
int hits[128]{},bests[128]{};
string w;
V f;
long long visited=0,pruned=0;
struct Best{int beta=-1,alpha=-1;string w;array<int,4> anchors{};};
Best best;
bool strong(){
 int m=w.size(),n=m+2; if(n<=2)return false;
 vector<U> row(n,0);
 for(int i=0;i<m;i++){
  for(int j=i+1;j<m;j++)row[i]|=U(1)<<j;
  for(int h=0;h<2;h++)if((w[i]-'0')>>(1-h)&1)row[m+h]|=U(1)<<i;else row[i]|=U(1)<<(m+h);
 }
 row[m]|=U(1)<<(m+1);
 auto reach=[&](int start){U have=U(1)<<start,prev=0;while(prev!=have){prev=have;for(int i=0;i<n;i++)if(have>>i&1)have|=row[i];}return have;};
 U all=(U(1)<<n)-1;
 return reach(0)==all&&([&]{for(int i=1;i<n;i++)if(!(reach(i)&1))return false;return true;})();
}
bool allowed(){
 int m=w.size(),n=m+2;
 auto arc=[&](int i,int j){
  if(i<m&&j<m)return i<j;
  if(i>=m&&j>=m)return i<j;
  if(i>=m)return ((w[j]-'0')>>(m+1-i)&1)!=0;
  return ((w[i]-'0')>>(m+1-j)&1)==0;
 };
 vector<U> tris;
 for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)for(int k=j+1;k<n;k++){
   if(arc(i,j)!=arc(j,k)||arc(i,j)!=arc(k,i))continue;
   U tri=(U(1)<<i)|(U(1)<<j)|(U(1)<<k);
   for(U t:tris)if(!(t&tri))return false;
   tris.push_back(tri);
 }
 if(tris.empty())return false;
 U ints=(U(1)<<n)-1;
 for(U t:tris)ints&=t;
 return ints==0;
}
void evaluate(){
 if(!strong()||(intercyclic&&!allowed()))return;
 int m=w.size();hits[m]++;int local=-1;
 for(int i=0;i<m;i++)for(int j=i+1;j<m;j++){
   int jl=w[i]-'0',jr=w[j]-'0';
   if(!((jl&2&&jr&1)||(jl&1&&jr&2)))continue;
   for(int k=j+1;k<m;k++)for(int l=k+1;l<m;l++){
      int kl=3-(w[k]-'0'),kr=3-(w[l]-'0');
      if(!((kl&2&&kr&1)||(kl&1&&kr&2)))continue;
      int beta=k-j-1, alpha=j-i-1+l-k-1;
      local=max(local,beta);
      if(beta>best.beta||(beta==best.beta&&alpha>best.alpha)){
       best={beta,alpha,w,{i,j,k,l}};
      }
   }
 }
 bests[m]=max(bests[m],local);
}
void dfs(V p){
 visited++;
 auto h=dot(p,f);
 if(h>target){pruned++;return;}
 if(h==target)evaluate();
 if((int)w.size()==m_limit)return;
 for(int d=0;d<4;d++){w+=char('0'+d);dfs(push(p,d));w.pop_back();}
}
int main(int argc,char**argv){
 target=argc>1?stoi(argv[1]):385;
 intercyclic=argc>2 && string(argv[2])=="intercyclic";
 assert(target>0); 
 m_limit=0;for(U u=target;u>1;u>>=1)m_limit++;
 m_limit+=3;
 assert(m_limit<40);
 for(int d=0;d<4;d++)for(int r=0;r<2;r++){
  int a=(d>>(1-r))&1,b=(d>>r)&1;
  for(int p=0;p<3;p++)for(int q=0;q<3;q++)for(int s=0;s<2;s++){
   auto&v=go[d][18*r+6*p+2*q+s];
   v.push_back(18*r+6*(a?2:1)+2*q+s);
   if(q||a)v.push_back(18*r+6*p+2*(b?2:1)+s);
   if(s||b)v.push_back(18*r+6*p+2*q+1);
  }
 }
 for(int&m:bests)m=-1;
 f=endpoint(1);V p{};p[0]=p[18]=1;dfs(p);
 cout<<"target "<<target<<" intercyclic "<<intercyclic<<" cap "<<m_limit<<" visited "<<visited<<" pruned "<<pruned<<endl;
 for(int m=0;m<=m_limit;m++)if(hits[m])cout<<"spine length "<<m<<" strong labelled words "<<hits[m]<<" best proper overlap "<<bests[m]<<endl;
 if(best.beta>=0)cout<<"best beta "<<best.beta<<" alpha "<<best.alpha<<" word "<<best.w<<" anchors "<<best.anchors[0]<<" "<<best.anchors[1]<<" "<<best.anchors[2]<<" "<<best.anchors[3]<<endl;
}
