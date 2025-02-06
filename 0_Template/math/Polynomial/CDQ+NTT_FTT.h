#include<bits/stdc++.h>
using namespace std;
#define ll long long
#define inf 0x3f3f3f3f
#define mod 998244353
const int N=1<<22;
ll ge=3,gi;
int main()
{
    ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    int t =1;
    void solve();
    // cin>>t;
    while(t--)solve();
}
//CDQ+NTT/FFT
//可用于处理卷积
//若我们求得了左区间，即可求出左区间对于右区间的贡献
//而这个贡献可通过NTT求出
//这就是所谓的分治FFT（用NTT精度更高）
ll f[N],g[N];
int rev[N]={0};
ll ksm(ll a,ll b)
{
    ll ans=1;
    while(b)
    {
        if(b&1)ans=ans*a%mod;
        a=a*a%mod;
        b>>=1;
    }
    return ans;
}
void chang(ll *A,int n)
{
    for(int i=0;i<n;i++)
    {
        rev[i]=rev[i>>1]>>1;
        rev[i]|=(i&1)?(n>>1):0;
    }
    for(int i=0;i<n;i++)
    {
        if(i<rev[i])swap(A[i],A[rev[i]]);
    }
}
void NTT(ll *A,int n,int opt)
{
    chang(A,n);
    for(int mid=1;mid<n;mid<<=1)
    {
        ll g1=ksm((opt==1)?ge:gi,(mod-1)/(mid<<1));
        for(int R=mid<<1,j=0;j<n;j+=R)
        {
            ll gk=1;
            for(int k=0;k<mid;k++,gk=gk*g1%mod)
            {
                ll x=A[j+k],y=A[j+k+mid]*gk%mod;
                A[j+k]=(x+y)%mod;
                A[j+k+mid]=(x-y+mod)%mod;
            }
        }
    }
}
void mul(ll *A,ll *B,ll n)  //用于处理多项式乘法，A、B均为系数式，且A为返回的系数式
{
    NTT(A,n,1);NTT(B,n,1);

    for(int i=0;i<n;i++)A[i]=A[i]*B[i]%mod;
    
    NTT(A,n,-1);

    ll inv=ksm(n,mod-2);

    for(int i=0;i<n;i++)A[i]=A[i]*inv%mod;
}
ll ta[N],tb[N];
void CDQ(int l,int r)
{
    int mid=(l+r)>>1;
    if(l==r)return ;
    CDQ(l,mid);
    ll mx=1;
    while(mx<(mid-l+r-l)+1)mx<<=1;  //取不小于区间长度的 二的幂

    for(int i=0;i<mx;i++)ta[i]=tb[i]=0;  //初始化

    for(int i=l;i<=mid;i++)ta[i-l]=f[i];  

    for(int i=1;i<=r-l;i++)tb[i]=g[i];

    mul(ta,tb,mx);  //ta乘tb

    for(int i=mid+1;i<=r;i++)f[i]=(f[i]+ta[i-l])%mod;//算贡献
    
    CDQ(mid+1,r);
}
void solve()
{
    gi=ksm(ge,mod-2);
    int n;
    cin>>n;
    for(int i=1;i<n;i++)
    {
        cin>>g[i];
        f[i]=0;
    }
    f[0]=1;  //这是f的边界
    g[0]=0;
    CDQ(0,n-1);  //分治
    for(int i=0;i<n;i++)cout<<f[i]<<" ";
    cout<<"\n";
}