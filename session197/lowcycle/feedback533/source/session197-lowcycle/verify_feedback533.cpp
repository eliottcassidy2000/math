#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdexcept>
using namespace std;
int main(int argc,char**argv){if(argc!=2)return 2;ifstream f(argv[1]);string line;long long checked=0;while(getline(f,line)){replace(line.begin(),line.end(),',',' ');istringstream in(line);int n,c;in>>n>>c;vector<int>a(n);for(int&x:a)in>>x;if(c>2)continue;checked++;
 vector<int>d((1<<n)*n);for(int v=0;v<n;v++)d[(1<<v)*n+v]=1;
 for(int S=1;S<(1<<n);S++)for(int u=0;u<n;u++)if(d[S*n+u]){int q=a[u]&~S;while(q){int b=q&-q;q-=b;int v=__builtin_ctz((unsigned)b);int&x=d[(S|b)*n+v];x=min(534,x+d[S*n+u]);}}
 int h=0;for(int u=0;u<n;u++)h=min(534,h+d[(((1<<n)-1)*n)+u]);
 if(h==533){cout<<n<<','<<c;for(int x:a)cout<<','<<x;cout<<'\n';}
 }cerr<<checked<<'\n';}
