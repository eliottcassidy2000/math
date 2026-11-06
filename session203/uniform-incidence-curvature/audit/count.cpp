#include <iostream>
#include <vector>
#include <functional>
#include <stdexcept>
using namespace std; using Z=long long;
void need(bool b,const char*s){if(!b)throw runtime_error(s);}
vector<int> literal(const vector<int>&core,int k,int m,int polarity,int uniforms,int word){
 int c=core.size(),n=c+m;vector<vector<int>> edge(n,vector<int>(n));
 for(int i=0;i<c;i++)for(int j=0;j<c;j++)edge[i][j]=(core[i]>>j)&1;
 for(int j=0;j<m;j++)for(int i=0;i<c;i++){
  int b=i<k?(((word>>j)&1)!=((polarity>>i)&1)):((uniforms>>(i-k))&1);
  edge[i][c+j]=b;edge[c+j][i]=1-b;
 }
 for(int i=0;i<m;i++)for(int j=i+1;j<m;j++)edge[c+i][c+j]=1;
 vector<int> rows(n);for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(edge[i][j])rows[i]|=1<<j;
 for(int i=0;i<n;i++){need(!edge[i][i],"loop");for(int j=i+1;j<n;j++)need(edge[i][j]+edge[j][i]==1,"not tournament");}
 return rows;
}
Z count_starting(const vector<int>&a){int n=a.size();if(n==0)return 1;vector<Z>d((1<<n)*n);for(int s=1;s<(1<<n);s++)for(int v=0;v<n;v++)if(s>>v&1){int rest=s^(1<<v);Z z=rest?0:1;for(int w=0;w<n;w++)if((rest>>w&1)&&(a[v]>>w&1))z+=d[rest*n+w];d[s*n+v]=z;}Z ans=0;for(int v=0;v<n;v++)ans+=d[((1<<n)-1)*n+v];return ans;}
Z literal_dfs(const vector<int>&a){int n=a.size();Z z=0;function<void(int,int)>visit=[&](int used,int last){if(used==(1<<n)-1){z++;return;}for(int v=0;v<n;v++)if(!(used>>v&1)&&(a[last]>>v&1))visit(used|(1<<v),v);};for(int v=0;v<n;v++)visit(1<<v,v);return n?z:1;}
int main(){try{int cases;cin>>cases;for(int t=0;t<cases;t++){int c,k,m,p,u;cin>>c>>k>>m>>p>>u;vector<int>core(c);for(int&i:core)cin>>i;cout<<t;for(int b=0;b<(1<<m);b++){auto a=literal(core,k,m,p,u,b);Z h=count_starting(a);if(a.size()<=8)need(h==literal_dfs(a),"DFS mismatch");need(h%2==1,"HP parity failure");cout<<" "<<h;}cout<<"\n";}}catch(exception&e){cerr<<e.what()<<endl;return 1;}}
