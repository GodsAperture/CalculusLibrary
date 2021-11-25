There are many overloaded operators specifically made to make working with the `matrix` struct much easier.

Typecasting using `(matrix<type>)` is possible for any datatype that exists, so long as the data type within the matrix can be cast to the new type.
```c++
matrix A(2, 2);
A.mat = {3.14, 2.718, 1.414, 3.0};

(matrix<int>) A //The contents of the local copy of A are now {3, 2, 1, 3}.
```

The `=` operator will resize the variable on the left to fit whatever the matrix size is on the right.
```c++
matrix A;  //A is currently a 1x1 matrix, since no size was given open construction.
matrix x(3,3);

x.mat = {2, 9, 9, 7, 9, 2, 4, 5, 8};

A = 2 * x; //A will now be a 3x3 matrix with every element multiplied by 2 from x.
```
A will resize to be the same dimensions as `b` and contain every element of `b`.

The `~` will transpose a local matrix.
```c++
matrix A(1,3);
A.mat = {1,2,3};
~A //This matrix will now be a 3x1 matrix with elements moved to the proper position to match the transpose.
```

The `<=` will concatenate matrices along the rows. 
```c++
matrix A(2,3);
matrix b(5,3);

A <= b //This will return a new matrix of dimensions 7x3 where the top 2x3 will be A and the bottom 5x3 will be b.
```

The `<<` operator will concatenate matrices along their columns.
```c++
matrix A(2,3);
matrix b(2,7);

A << b //This will return a new matrix of dimensions 2x10 where the first 3 columns will contain A and the last 7 columns will contain b.
```

The `[]` operator will return the indexed row from the matrix and return a matrix containing those elements.
```c++
matrix A(2, 2);
A.mat = {1, 2, 3, 4};

A[1] //This will return a 1x2 matrix containing {3, 4}.
```

The `()` operator will return the indexed column from the matrix and return a matrix containing those elements.
```c++
matrix A(2, 2);
A.mat = {1, 2, 3, 4};

A(1) //This will return a 2x1 matrix containing {2, 4}.
```

The `*` operator will return the product of two matrices. If one of the matrices is a 1x1,
then it will multiply every element in the other matrix by the single element in the 1x1 matrix.
```c++
matrix A(3, 3);
matrix x(3, 1);
matrix b;

A.mat = {2, 9, 9, 7, 9, 2, 4, 5, 8};
x.mat = {1, 2, 3};
b.mat = {7};

A * x //This will return the product of the two matrices, as dictated by the definition of matrix multiplication.
A * b //Since b is a 1x1 matrix, every element in A will be multiplied by the lone element in b.
```

The `/` operator will perform solving similar to MatLab's implementation. If the matrix on the right is a 1x1,
then every element in the left matrix will be divided by the lone value in the matrix on the right.
```c++
matrix A(3, 3);
matrix b(3, 1);
matrix c;
A.mat = {2, 9, 9, 7, 9, 2, 4, 5, 8};
b.mat = {1, 2, 3};
c.mat = {3};

b / A //This is equivalent to b * A^-1 or solving the system x * A = b for x. The returned value will be a 3x1 matrix containing {0.717, 0.323, -0.306}.
A / c //This will return a new matrix whose contents are the elements of A divided by the lone element in c.
```

The `+` operator will return a matrix that has the index wise sum of every element between two matrices.
```c++
matrix A(2, 2);
matrix x(2, 2);
A.mat = {1, 2, 3, 4};
x.mat = {5, 6, 7, 8};

A + x //This will return a matrix whose contents are {6, 8, 10, 12};
```

The `-` operator has two functions. It will either make a local matrix negative or it will subtract two matrices.
```c++
matrix A;
matrix x;

-A;//This will be the negative of every element in A.
A - x;//This will subtract every index wise element in b from A.
```

The `|` operator behaves similarly to MatLab's `\` operator, solving a system of equations.
```c++
matrix A(3, 3);
matrix b(3, 1);
A.mat = {2, 9, 9, 7, 9, 2, 4, 5, 8};
b.mat = {1, 2, 3};

A | b //This is equivalent to A^-1 * b or solving the system A * x = b for x. The returned value will be a 3x1 matrix containing {-0.510, 2.117, -0.306}.
```

The `^` operator behaves like MatLab's `^` operator. It will create an `MxN` matrix from a `Mx1`
matrix as the base and a `1xN` matrix as the exponent or making a different `MxN` matrix from
a `1xN` matrix as the base with a `Mx1` matrix as the exponent.
```c++
matrix A(1, 3);
matrix x(3, 1);
A.mat = {-1, 0, 1};
x.mat = {0, 1, 2};

A ^ x //This will return a 3x3 matrix whose contents are {(-1)^0, 0^0, 1^0, (-1)^1, 0^1, 1^1, (-1)^2, 0^2, 1^2}.
      //After simplifying, that will become {1, 1, 1, -1, 0, 1, 1, 0, 1}.
x ^ A //This will return a 3x3 matrix whose contents are {0^-1, 0^0, 0^1, 1^-1, 1^0, 1^-1, 2^-1, 2^0, 2^1}.
      //After simplifying, that will become {NaN, 1, 0, 1, 1, 1, 0.5, 1, 2}.
```

The `*=` operator is like MatLab's `.*` operator in that it will just perform index wise multiplication.

The `^=` operator is like MatLab's `.^` operator in that it will just perform index wise exponentiation.
```c++
matrix A(2, 2);
matrix x(2, 2);
A.mat = {1, 2, 3, 4};
x.mat = {1, 0, 0, 1};

A ^= x //This will return {1^1, 2^0, 3^0, 4^1}.
```
