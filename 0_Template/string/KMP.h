#include <template_overAll.h>

class KMP
{
private:
    string s;
    string inis;
public:
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
        s = p + "#" + s;
        prefix_function();
        vector<int> v;
        for (int i = p.size() + 1; i < s.size(); i++)
            if (pi[i] == p.size())
                v.push_back(i - 2 * p.size());
        return v;
    }
};
 