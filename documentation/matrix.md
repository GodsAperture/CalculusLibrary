The `matrix<T>` data type emulates matrices from mathematics, including retaining its dimensions as `dim1` and `dim2`.
The contents of the matrix are stored in an `std::vector<T>` named `mat` and if no data type is specified upon construction, it will default to type `double`.

```c++
#include "calclib.h"


int main(){
  matrix y(3,2);
  y = {1, 2, 3, 4, 5, 6};

  calculus::printMatrixf(y);
  
  return 0;
```
