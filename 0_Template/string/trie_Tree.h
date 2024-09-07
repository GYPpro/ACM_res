#include <template_overAll.h>

class Trie//AC
{
public:
    vector<map<char, int>> t;
    int root = 0;
    Trie()
    {
        t.resize(1);
    }
    void addedge(string _s)
    {
        int pvidx = root;
        _s.push_back('-');
        for (int i = 0; i < _s.size(); i++)
        {
            if (t[pvidx].find(_s[i]) != t[pvidx].end())
            {
                pvidx = t[pvidx][_s[i]];
            }
            else
            {
                t[pvidx][_s[i]] = t.size();
                t.push_back(map<char, int>());
                pvidx = t[pvidx][_s[i]];
            }
        }
    }
    bool ifcmp(string &s)
    {
        int pvidx = root;
        for(int i = 0;i < s.size();i ++)
        {
            if(t[pvidx].find(s[i]) != t[pvidx].end()) pvidx = t[pvidx][s[i]];
            else return 0;
        }
        return t[pvidx].find('-') != t[pvidx].end();
    }
};