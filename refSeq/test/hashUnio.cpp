#include <fstream>  //ifstream读文件，ofstream写文件，fstream读写文件
#include <string>   //文本对象，储存读取的内容
#include <iostream> //屏幕输出cout
#include <cstdlib>  //调用system("pause");
#include <map>
#include <windows.h> //用于函数SetConsoleOutputCP(65001);更改cmd编码为utf8
using namespace std;

int main()
{
    SetConsoleOutputCP(65001);
    ifstream in("./refSeq/test/test");
    string line;

    if (in) // 有该文件
    {
        while (getline(in, line)) // line中不包括每行的换行符
        {
            cout << line << endl;
        }
    }
    else // 没有该文件
    {
        cout << "no such file" << endl;
    }

        return 0;
}