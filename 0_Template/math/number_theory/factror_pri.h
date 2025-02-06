int n, m, p, b[10000005], prime[1000005], t, min_prime[10000005];
void euler_Prime(int n)
{ // 用欧拉筛求出1~n中每个数的最小质因数的编号是多少，保存在min_prime中
    for (int i = 2; i <= n; i++)
    {
        if (b[i] == 0)
        {
            prime[++t] = i;
            min_prime[i] = t;
        }
        for (int j = 1; j <= t && i * prime[j] <= n; j++)
        {
            b[prime[j] * i] = 1;
            min_prime[prime[j] * i] = j;
            if (i % prime[j] == 0)
                break;
        }
    }
}
long long c(int n, int m, int p)
{ // 计算C(n,m)%p的值
    euler_Prime(n);
    int a[t + 5]; // t代表1~n中质数的个数 ，a[i]代表编号为i的质数在答案中出现的次数
    for (int i = 1; i <= t; i++)
        a[i] = 0; // 注意清0，一开始是随机数
    for (int i = n; i >= n - m + 1; i--)
    { // 处理分子
        int x = i;
        while (x != 1)
        {
            a[min_prime[x]]++; // 注意min_prime中保存的是这个数的最小质因数的编号（1~t）
            x /= prime[min_prime[x]];
        }
    }
    for (int i = 1; i <= m; i++)
    { // 处理分母
        int x = i;
        while (x != 1)
        {
            a[min_prime[x]]--;
            x /= prime[min_prime[x]];
        }
    }
    long long ans = 1;
    for (int i = 1; i <= t; i++)
    { // 枚举质数的编号，看它出现了几次
        while (a[i] > 0)
        {
            ans = ans * prime[i] % p;
            a[i]--;
        }
    }
    return ans;
}
int main()
{
    cin >> n >> m;
    m = min(m, n - m); // 小优化
    cout << c(n, m, MOD);
}
