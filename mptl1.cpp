#include <string>
#include "Python.h"

using namespace std;

int main()
{
   Py_Initialize();
   int x[5] = {0, 1, 2, 3, 4};
   int y[5] = {5, 1, 7, 5, 1};
   string command = "pylab.plot([";
   for(int i = 0; i < 4; i++) {
       command += x[i];
       command += ", ";
   }
   command += x[4];
   command += "], [";
   for(int i = 0; i < 4; i++) {
       command += y[i];
       command += ", ";
   }
   command += y[4];
   command += "])";
   PyRun_SimpleString("import pylab");
   PyRun_SimpleString(command.c_str());
   PyRun_SimpleString("pylab.show()");
   Py_Exit(0);
   return 0;
}
