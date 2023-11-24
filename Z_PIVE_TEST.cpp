// #include <basic\subsearch.h>
// #include <math\other\binpow.h>
// #include <string\trie_Tree.h>
// #include <ds\segTree.h>
// #include <string\KMP.h>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <string.h>
#include <math.h>
#include <set>
#include <algorithm>
#include <iostream>
#include <queue>

using namespace std;

#define ll long long
#define ld long double
const ll int maxn = 1E5+10;
const ll int mod1 = 998244353;
const ll int mod2 = 1E9+7;

class KMP
{
public:
    string s;
    string inis;
    vector<int> pi;
    KMP(string _s)
    {
        s = _s;
        inis = s;
    }
    void prefix_function()
    {
        pi.clear();
        int n = (int)s.length();
        pi.resize(n);
        for (int i = 1; i < n; i++)
        {
            int j = pi[i - 1];
            while (j > 0 && s[i] != s[j])
                j = pi[j - 1];
            if (s[i] == s[j])
                j++;
            pi[i] = j;
        }
        return;
    }
    vector<int> find_occr(string p)
    {
        s = inis;
        s =p+"#"+s;
        prefix_function();
        vector<int> v;
        for(int i = p.size() +1;i < s.size();i ++) if(pi[i] == p.size())v.push_back(i-2*p.size());
        return v;
    }
};

int main()
{
    string s;
    cin >> s;
    string t;
    cin >> t;
    KMP kmp(s);
    vector<int> vec = kmp.find_occr(t);
    for(int i = 0;i < vec.size();i ++) cout << vec[i] +1<< "\n";
    for(int i = 0;i < t.size();i ++) cout << kmp.pi[i] << " ";
    cout << "\n";
    // system("pause");
};