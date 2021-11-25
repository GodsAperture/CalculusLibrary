The `matrix<T>` data type emulates matrices from mathematics, including retaining its dimensions as `dim1` and `dim2`.
The contents of the matrix are stored in an `std::vector<T>` named `mat` and if no data type is specified upon construction, it will default to type `double`.

The `matrix` struct has been templated to hold any defined data type the user wishes to have. The operators that need to be overloaded so that `matrix` use its operator overloads are `/`, `*`, `-`, and `+`. The rest are templated to make overloading easier for the users preferred data type.

All matrices can be constructed with or without dimension specifications.
```c++
matrix A;
matrix b(1,3);
```
Matrices that are generated without dimension specifications default to a 1x1 matrix.

Matrices can then be assigned values through a generated list.
```c++
b.mat = {1, 2, 3};
A = b;
```
`A` will automatically be resized to be the same dimensions as `b`.

It is entirely possible to type cast a copy of a matrix to a new type.
```c++
matrix A; //A is a matrix of type double.

(matrix<int>) A;
//This is now a matrix whose entries have all been cast to type int.
//However, A is still defined as a matrix of type double everywhere else.
```
