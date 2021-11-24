The `matrix<T>` data type emulates matrices from mathematics, including retaining its dimensions as `dim1` and `dim2`.
The contents of the matrix are stored in an `std::vector<T>` named `mat` and if no data type is specified upon construction, it will default to type `double`.

```c++
#include "calclib.h"


int main(){
  matrix y(3,3);
  y = {2, 9, 9, 7, 9, 2, 4, 5, 8};

  calculus::printMatrixf(y);
  
  return 0;
```
This code will print out:
```c++
2.000000 9.000000 9.000000
7.000000 7.000000 7.000000
4.000000 5.000000 8.000000
```
