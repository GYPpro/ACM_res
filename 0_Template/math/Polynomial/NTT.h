#include<bits/stdc++.h>
using namespace std;
#define ll long long
#define ld long double
#define inf 0x3f3f3f3f
#define mod 998244353
const ll N=1<<22;  //注意这个N
int main()
{
    ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    int t =1;
    void solve();
    // cin>>t;
    while(t--)solve();
}
//NTT最基础的板子
//FFT涉及三角函数和复数，浮点计算导致运算的复杂度大，精度低
//由此诞生了快速数论变化NTT
//换个根就行了
ll g=1,gi;  //g是原根，一般取3，gi是g的乘法逆元
ll tmp1[N],tmp2[N];
ll mx;
int rev[N]={0};

ll ksm(ll a,ll b)  //快速幂
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
void chang(ll *f,int n)  //蝴蝶变换
{
    for(int i=0;i<n;i++)
    {
        rev[i]=rev[i>>1]>>1;
        rev[i]|=(i&1)?(n>>1):0;
    }
    for(int i=0;i<n;i++)
     if(i<rev[i])swap(f[i],f[rev[i]]);
}

void NTT(ll *f,int n,int opt) //快速数论变换
{
    chang(f,n);//先变换
    for(int mid=1;mid<n;mid<<=1)  //枚举半周期
    {
        ll g1=ksm((opt==1)?g:gi,(mod-1)/(mid*2));  //取单位根

        for(int R=mid<<1,j=0;j<n;j+=R)
        {
            
            ll w=1;  //变化的自变量
            
            for(int k=0;k<mid;k++,w=(w*g1)%mod)
            {
                ll a1=f[j+k],a2=f[j+k+mid]*w%mod;
                
                f[j+k]=(a1+a2)%mod;
                
                f[j+k+mid]=(a1-a2+mod)%mod;
            }
        }
    }
}
void solve()
{
    int n,m;
    cin>>n>>m;

    g=3;gi=ksm(g,mod-2);  //原根 和 它的倒数 

    for(int i=0;i<=n;i++)cin>>tmp1[i];
    for(int i=0;i<=m;i++)cin>>tmp2[i];
    
    mx=1;

    while(mx<n+m+1)mx<<=1;  //取不小于最高次数的二的幂
    
    for(int i=n+1;i<mx;i++)tmp1[i]=0;  //高位补0
    for(int i=m+1;i<mx;i++)tmp2[i]=0;
    
    NTT(tmp1,mx,1);  //分别求点值式
    NTT(tmp2,mx,1);
    
    for(int i=0;i<mx;i++)tmp1[i]=tmp1[i]*tmp2[i]%mod;  //值相乘
    
    NTT(tmp1,mx,-1);  //点值式转系数式

    ll inv=ksm(mx,(mod-2));  //取mx的逆元
    for(int i=0;i<n+m+1;i++)
    {
        cout<<tmp1[i]*inv%mod<<" ";  //这里不要用除法，用乘法逆元
    }
    cout<<"\n";
}