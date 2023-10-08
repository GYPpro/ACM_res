#include <basic\subsearch.h>
#include <iostream>

int main()
{
    for (;;)
    {
        vector<ll int> a = testReadLine();
        for (;;)
        {
            string s;
            cin >> s;
            ll int que = str2int(s);
            if (que == INT64_MAX)
            {
                getchar(); break;
            }
            cout << subsearch_3<ll int>(que, a) << endl;
        }
    }
}