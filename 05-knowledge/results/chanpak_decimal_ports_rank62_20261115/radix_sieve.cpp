#include <iostream>
#include <array>
#include <vector>
using namespace std;
constexpr size_t N=8;
using M=array<array<int,N>,N>;
constexpr int P=101;
M U[13];
int inv[P];
long seen=0, hasradix=0, hasboth=0;
vector<int> w;
M mul(const M &A, const M &B){
    M C{};
    for(int i=0;i<N;i++)for(int k=0;k<N;k++)for(int j=0;j<N;j++)C[i][j]+=A[i][k]*B[k][j];
    for(auto& row:C)for(auto& x:row)x%=P;
    return C;
}
bool eigen(M a,int x){
    for(int i=0;i<N;i++)a[i][i]=(a[i][i]-x+P)%P;
    for(int j=0;j<N;j++){
        int i=j;while(i<N && a[i][j]==0)++i;
        if(i==N)return true;
        swap(a[i],a[j]);const int d=inv[a[j][j]];
        for(i=j+1;i<N;i++){
            const int t=a[i][j]*d%P;
            for(int k=j;k<N;k++)a[i][k]=(a[i][k]+P*P-t*a[j][k])%P;
        }
    }
    return false;
}
void go(M a,int n){
    if(n){
        ++seen;
        if(eigen(a,1000%P)){
            hasradix++;
            if(eigen(a,1)){
                hasboth++;
                for(int j:w)cout<<j<<" ";cout<<"\n";
            }
        }
    }
    for(int j=0;j<13;j++){
        const int b=n+(j<4?2:3);
        if(b>12)continue;
        w.push_back(j);go(mul(U[j],a),b);w.pop_back();
    }
}
int main(){
    for(int i=1;i<P;i++)for(int j=1;j<P;j++)if(i*j%P==1)inv[i]=j;
    for(M &m:U)for(auto &r:m)for(auto &x:r)if(!(cin>>x))return 1;
    M I{};for(int j=0;j<N;j++)I[j][j]=1;
    go(I,0);
    cerr<<seen<<" "<<hasradix<<" "<<hasboth<<"\n";
}
