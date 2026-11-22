#include <bits/stdc++.h>
using namespace std;
struct Stats{
 int h=0;
 int before[10][10]={};
 vector<int> edges;
 explicit Stats(int n):edges(n){}
};
Stats read_paths(const vector<int>&a){
 int n=a.size(),all=(1<<n)-1;
 vector<int> order(n);Stats s(n);
 auto walk=[&](auto&& self,int used,int k)->void{
  if(k==n){
   ++s.h;
   for(int i=0;i<n;i++)for(int j=i+1;j<n;j++)++s.before[order[i]][order[j]];
   for(int i=0;i+1<n;i++)s.edges[order[i]]|=1<<order[i+1];
   return;
  }
  for(int E=(k?a[order[k-1]]:all)&~used;E;E&=E-1){
   int v=__builtin_ctz(E);order[k]=v;self(self,used|(1<<v),k+1);
  }
 };
 walk(walk,0,0);
 return s;
}
int main(int argc,char**argv){
 string root=argc>1?argv[1]:"session212/cycle-levels";
 const int expected[]={1,1,1,2,4,12,56,456,6880};
 for(int n=3;n<=8;n++){
  ifstream in(root+"/classes"+to_string(n)+".rows");
  if(!in)throw runtime_error("class list unavailable");
  vector<int>a(n);int classes=0, occurrence=0;
  while(in>>a[0]){
   for(int i=1;i<n;i++)if(!(in>>a[i]))throw runtime_error("partial row list");
   Stats d=read_paths(a), s=read_paths(d.edges);
   if(!(s.h&1)||s.h!=d.h||s.edges!=d.edges)throw runtime_error("support data");
   auto T=s.edges;
   for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
    bool visible=(T[i]>>j&1)||(T[j]>>i&1);
    if(!visible){
     ++occurrence;
     if(s.before[i][j]&1) T[i]|=1<<j;
     else T[j]|=1<<i;
     int reverse_count=(a[i]>>j&1)?s.before[j][i]:s.before[i][j];
     if(reverse_count%4)throw runtime_error("reverse-order divisibility");
    }
   }
   if(T!=a)throw runtime_error("failed reconstruction");
   ++classes;
  }
  if(classes!=expected[n])throw runtime_error("wrong deck length");
  cout<<n<<" "<<classes<<" "<<occurrence<<"\n";
 }
}
