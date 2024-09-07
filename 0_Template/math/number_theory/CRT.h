
int CRT(vector<int> &r, vector<int> &a)
{ // % r === a
    int n = a.size();
    __int128 k = 1, ans = 0;
    for (int i = 0; i < n; i++) k *= r[i];
    for (int i = 0; i < n; i++)
    {
        __int128 m = k / r[i];
        int b, y;
        exgcd(m, r[i], b, y); // b * m mod r[i] = 1
        ans = (ans + a[i] * m * b % k) % k;
    }
    return (ans % k + k) % k;
}



int mul(int a, int b, int m) {
    return (__int128)a * b % m;
}

int exgcd (int a,int b,int &x,int &y) {
    if (b == 0) { x = 1, y = 0; return a; }
    int g = exgcd(b, a % b, x, y), tp = x;
    x = y, y = tp - a / b * y;
    return g;
};

int EXCRT(vector<int> &a,vector<int> &r) { // % r == a 
    int x, y, k;
    int n = r.size();
    int M = a[0], ans = r[0];
    for (int i = 1; i < n; ++ i) {
        int ca = M, cb = a[i], cc = (r[i] - ans % cb + cb) % cb;
        int gcd = exgcd(ca, cb, x, y), bg = cb / gcd;
        if (cc % gcd != 0) return -1;
        x = mul(x, cc / gcd, bg);
        ans += x * M;
        M *= bg;
        ans = (ans % M + M) % M;
    }
    return (ans % M + M) % M;
}