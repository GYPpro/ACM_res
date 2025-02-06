#include<bits/stdc++.h>
using namespace std;
#define ll long long
#define mod 998244353
const int N=(1<<22)|10;
int gen=3,gi;
int ksm(int a,int b)
{
    int ans=1;
    while(b)
    {
        if(b&1)ans=(ll)ans*a%mod;
        a=(ll)a*a%mod;
        b>>=1;
    }
    return ans;
}
signed main()
{
    void solve();
    ios::sync_with_stdio(0);cin.tie(0);cout.tie(0);
    gi=ksm(gen,mod-2);
    solve();
    return 0;
}
//多项式求逆元，只能用NTT
int a[N],rev[N]={0},b[N],tmp[N];
void chang(int *A,int n)  //蝴蝶变换
{
    for(int i=0;i<n;i++)
    {
        rev[i]=(rev[i>>1]>>1)|((i&1)?(n>>1):0);
    }
    for(int i=0;i<n;i++)
    {
        if(i<rev[i])swap(A[i],A[rev[i]]);
    }
}
void NTT(int *A,int n,int opt)
{
    chang(A,n);  //先变换
    for(int mid=1;mid<n;mid<<=1)
    {
        int g1=ksm((opt==1)?gen:gi,(mod-1)/(mid<<1));
        for(int j=0;j<n;j+=(mid<<1))
        {
            int gk=1;
            for(int k=0;k<mid;k++,gk=(ll)gk*g1%mod)
            {
                int x=A[j+k],y=(ll)gk*A[j+k+mid]%mod;
                A[j+k]=((ll)x+(ll)y)%mod;
                A[j+k+mid]=((ll)x-(ll)y+mod)%mod;
            }
        }
    }
    if(opt==1)return ;  //如果是系数式求点值式，到这里即可

    int inv=ksm(n,mod-2);         //否则这里直接处理出答案
    for(int i=0;i<n;i++)A[i]=(ll)inv*A[i]%mod;  
}
void INV(int n,int *A,int *B)  //传入系数个数n、原系数式A、用于返回答案的系数式B
{
    //本质是一种分治的思想，n的答案可由n/2求出，而n==1的答案易得
    stack<int> stk;
    int mx=1;
    while(n!=1){stk.push(n);n=(n+1)>>1;}  //这里在用栈模拟递归
    stk.push(1);
    while(!stk.empty())
    {
        n=stk.top();
        stk.pop();  //取栈顶，先

        if(n==1){B[0]=ksm(A[0],mod-2);continue;}  //0次即为常数的逆元

        while(mx<(n<<1))mx<<=1;  //处理出不小于n*2的 二的幂

        for(int i=0;i<n;i++)tmp[i]=A[i];  //低位复制
        for(int i=n;i<mx;i++)tmp[i]=0;  //高位补0

        NTT(tmp,mx,1);NTT(B,mx,1);  //求点值式

        for(int i=0;i<mx;i++)B[i]=(2ll-(ll)tmp[i]*B[i]%mod+mod)%mod*B[i]%mod;  //套公式计算答案值

        NTT(B,mx,-1);  //点值式转系数式

        for(int i=n;i<mx;i++)B[i]=0;  //高位补0
    }
    
}
void solve()
{
    int n;
    cin>>n;
    for(int i=0;i<n;i++)cin>>a[i];
    INV(n,a,b);
    for(int i=0;i<n;i++)cout<<b[i]<<" ";
    cout<<"\n";
}