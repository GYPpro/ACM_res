
template <class T>
T mypow(T n, int k, T r = 1)
{
    for (; k; k >>= 1, n *= n)
    {
        if (k & 1)
            r *= n;
    }
    return r;
}
// #undef int long long
template <int MOD>
struct mInt
{
    int x;
    mInt(int x = 0) : x(norm(x % MOD)) {}
    int val() const { return x; }
    int norm(int x) const { return (x + MOD) % MOD; }
    mInt operator-() const
    {
        int val = norm(MOD - x);
        return mInt(val);
    }
    mInt &operator++() { return x = norm(x + 1), *this; }
    mInt operator++(signed)
    {
        mInt res = *this;
        ++*this;
        return res;
    }
    mInt &operator--() { return x = norm(x - 1), *this; }
    mInt operator--(signed)
    {
        mInt res = *this;
        --*this;
        return res;
    }
    mInt inv() const
    {
        assert(x != 0);
        return mypow(*this, MOD - 2);
    }
    mInt &operator*=(const mInt &i) { return x = x * i.x % MOD, *this; }
    mInt &operator+=(const mInt &i) { return x = norm(x + i.x), *this; }
    mInt &operator-=(const mInt &i) { return x = norm(x - i.x), *this; }
    mInt &operator/=(const mInt &i) { return *this *= i.inv(); }
    mInt &operator%=(const int &i) { return x %= i, *this; }
    friend mInt operator*(const mInt &i, const mInt &j)
    {
        mInt res = i;
        res *= j;
        return res;
    }
    friend mInt operator+(const mInt &i, const mInt &j)
    {
        mInt res = i;
        res += j;
        return res;
    }
    friend mInt operator-(const mInt &i, const mInt &j)
    {
        mInt res = i;
        res -= j;
        return res;
    }
    friend mInt operator/(const mInt &i, const mInt &j)
    {
        mInt res = i;
        res /= j;
        return res;
    }
    friend mInt operator%(const mInt &i, const int &j)
    {
        mInt res = i;
        res %= j;
        return res;
    }
    friend auto &operator>>(istream &it, mInt &j)
    {
        int v;
        it >> v;
        j = mInt(v);
        return it;
    }
    friend auto &operator<<(ostream &o, const mInt &j) { return o << j.x; }
    bool operator<(const mInt &i) const { return x < i.x; }
    bool operator>(const mInt &i) const { return x > i.x; }
    bool operator==(const mInt &i) const { return x == i.x; }
    bool operator!=(const mInt &i) const { return x != i.x; }
};
#define int long long
// using Z = Zmod<998244353>;

const int N = 1 << 21;
static const int mod1 = 1E9 + 7, base1 = 127;
static const int mod2 = 1E9 + 9, base2 = 131;
using U = mInt<mod1>;
using V = mInt<mod2>;
vector<U> val1;
vector<V> val2;
void init(int n = N)
{
    val1.resize(n + 1), val2.resize(n + 2);
    val1[0] = 1, val2[0] = 1;
    for (int i = 1; i <= n; i++)
    {
        val1[i] = val1[i - 1] * base1;
        val2[i] = val2[i - 1] * base2;
    }
}
struct String
{
    vector<U> hash1;
    vector<V> hash2;
    string s;
    String(string s_) : s(s_), hash1{1}, hash2{1}
    {
        for (auto it : s)
        {
            hash1.push_back(hash1.back() * base1 + it);
            hash2.push_back(hash2.back() * base2 + it);
        }
    }
    pair<U, V> get()
    { // 输出整串的哈希值
        return {hash1.back(), hash2.back()};
    }
    pair<U, V> substring(int l, int r)
    { // 输出子串的哈希值
        if (l > r)
            swap(l, r);
        U ans1 = hash1[r + 1] - hash1[l] * val1[r - l + 1];
        V ans2 = hash2[r + 1] - hash2[l] * val2[r - l + 1];
        return {ans1, ans2};
    }
    pair<U, V> modify(int idx, char x)
    { // 修改 idx 位为 x
        int n = s.size() - 1;
        U ans1 = hash1.back() + val1[n - idx] * (x - s[idx]);
        V ans2 = hash2.back() + val2[n - idx] * (x - s[idx]);
        return {ans1, ans2};
    }
};

