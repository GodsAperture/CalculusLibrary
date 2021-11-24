/* ////DISCLAIMER////
 * THIS LIBRARY WAS BUILT FOR CALCULUS RELATED TASKS.
 * BECAUSE OF THE NUMERICAL METHODS EMPLOYED, IT IS SUGGESTED TO PAY ATTENTION TO YOUR DATA.
 * THE METHODS EMPLOYED ARE APPROXIMATIONS AND MAY CAUSE SEVERE DETERIORATION TO DATA IF USED IMPROPERLY.
 */

//This namespace will only contain commonly used constants, and will be updated as I find more. Constants will be Pascal case.

#include <complex>
#include <iostream>
#include <vector>
#include <cassert>

namespace constant{
//If you don't know what Pi is, God help you.
    const long double Pi = 3.1415926535897932385;

//If you don't know what E is, that's okay.
    const long double E = 2.7182818284590452354;
}

template <typename T = double>
struct matrix{
    int dim1, dim2;
    std::vector<T> mat;

    matrix<T> row(const std::vector<T> &arr){
        matrix<T> temp(1, arr.size());

        for(int i = 0; i < arr.size(); i++){
            temp.mat[i] = arr[i];
        }

        return temp;

    }

    matrix<T> column(const std::vector<T> &arr){
        matrix<T> temp(arr.size(), 1);

        for(int i = 0; i < arr.size(); i++){
            temp.mat[i] = arr[i];
        }

        return temp;

    }

    template<typename T1>
    explicit operator matrix<T1>() {
        matrix<T1> temporaryMatrix(dim1, dim2);

        for(int i = 0; i < dim1 * dim2; i++){
            temporaryMatrix.mat[i] = (T1) mat[i];
        }

        return temporaryMatrix;

    };

    explicit matrix(int dimension1 = 1, int dimension2 = 1): dim1(dimension1), dim2(dimension2), mat(dimension1 * dimension2){}

    void resize(const int newDim1, const int newDim2) {
        this->dim1 = newDim1;
        this->dim2 = newDim2;
        this->mat.resize(this->dim1 * this->dim2);
    }

    matrix<int> size(){
        matrix<int> temporaryMatrix(1, 2);

        temporaryMatrix.mat[0] = dim1;
        temporaryMatrix.mat[1] = dim2;

        return temporaryMatrix;

    }

    matrix<T> index(int dimension1, int dimension2){
        matrix<T> temp;

        temp.mat[0] = mat[dimension1 * dim2 + dimension2];

        return temp;

    }



    //Redefine a matrix, this will automatically resize the matrix on the left as necessary.
    matrix<T>& operator=(const matrix<T> &Matrix){
        resize(Matrix.dim1, Matrix.dim2);

        for(int i = 0; i < Matrix.dim1 * Matrix.dim2; i++){
            mat[i] = Matrix.mat[i];
        }

        return *this;

    }

    //Allows you to make every entry in a matrix negative.
    matrix<T> operator-(){
        matrix<T> temp(dim1, dim2);

        for(int i = 0; i < dim1 * dim2; i++){
            temp.mat[i] = -mat[i];
        }

        return temp;

    }

    //Transpose the matrix to the right. Used as ~Matrix.
    matrix<T> operator~(){
        matrix<T> temp(dim2, dim1);

        for(int i = 0; i < dim2; i++){
            for(int k = 0; k < dim1; k++){

                temp.mat[k + i * dim1] = mat[i + k * dim2];

            }
        }

        return temp;

    }

    //Add a new row or set of rows to the bottom of the matrix on the left.
    matrix<T> operator<=(const matrix<T> &Matrix){
        if(dim2 != Matrix.dim2){
            printf("The matrices being concatenated do not have the same dim2.\n");
            assert(0);
        }
        matrix<T> temp(dim1 + Matrix.dim1, dim2);
        resize(dim1 + Matrix.dim1, dim2);

        for(int i = 0; i < dim1 * dim2; i++){
            temp.mat[i] = mat[i];
        }

        for(int i = dim1 * dim2; i < (dim1 + Matrix.dim1) * dim2; i++){
            temp.mat[i] = Matrix.mat[i];
        }

        return temp;

    }

    //Add a new column or set of columns to the right of the matrix on the left.
    matrix<T> operator<<(const matrix<T> &Matrix){
        if(dim2 != Matrix.dim2){
            printf("The matrices being concatenated do not have the same dim1.\n");
            assert(0);
        }
        matrix<T> temporaryMatrix(dim1, dim2 + Matrix.dim2);

        //This for loop controls which row is being edited.
        for(int i = 0; i < dim1; i++){

            //This for loop adds the elements of *this to temporaryMatrix.
            for(int k = 0; k < dim2; k++){
                temporaryMatrix.mat[i * (dim2 + Matrix.dim2) + k] = mat[i * dim2 + k];
            }

            //This for loop adds the elements of &Matrix to temporaryMatrix.
            for(int k = 0; k < Matrix.dim2; k++){
                temporaryMatrix.mat[i * (dim2 + Matrix.dim2) + k + dim2] = Matrix.mat[i * Matrix.dim2 + k];
            }

        }

        return temporaryMatrix;

    }

    //Index and return a particular row of a matrix.
    matrix<T> operator[](int row){
        matrix<T> temp(1, dim2);

        for (int i = 0; i < dim2; i++) {
            temp.mat[i] = mat[i + dim2 * row];
        }

        return temp;

    }

    //Index and return a particular column of a matrix.
    matrix<T> operator()(int column){
        matrix<T> other(dim1, 1);

        for(int i = 0; i < dim1; i++){
            other.mat[i] = mat[i * dim2 + column];
        }

        return other;

    }



    //Index wise multiplication of two matrices.
    matrix<T> operator*=(const matrix<T> &Matrix){
        if((dim1 != Matrix.dim1) || (dim2 != Matrix.dim2)){
            printf("The matrices using index wise multiplication are not of the same dimensions.\n");
            assert(0);
        }
        matrix<T> temp(dim1, dim2);

        for(int i = 0; i < dim1 * dim2; i++){
            temp.mat[i] = mat[i] * Matrix.mat[i];
        }

        return temp;

    }

    //Multiply the matrix on the left either by a matrix where dim2 == Matrix.dim1 or multiply each entry by a 1x1 matrix.
    matrix<T> operator*(const matrix<T> &Matrix){
        if((dim2 != Matrix.dim1) && ((Matrix.dim1 != 1) && (Matrix.dim2 != 1))){
            printf("The matrices being multiplied either do not share the proper dimensions or the matrix on the right is not 1x1.\n");
            assert(0);
        }
        matrix<T> temp;


        if(dim2 == Matrix.dim1) {
            temp.resize(dim1, Matrix.dim2);
            T sum = (T) 0;

            for (int i = 0; i < dim1; i++) {
                for (int j = 0; j < Matrix.dim2; j++) {
                    for (int k = 0; k < dim2; k++) {

                        sum = sum + mat[k + i * dim2] * Matrix.mat[j + k * Matrix.dim2];

                    }

                    temp.mat[i * Matrix.dim2 + j] = sum;
                    sum = (T) 0;

                }

            }
        }
        else{
            temp.resize(dim1, dim2);

            for(int i = 0; i < dim1 * dim2; i++){

                temp.mat[i] = mat[i] * Matrix.mat[0];

            }
        }

        return temp;
    }

    //Multiply each entry in the matrix on the left by a constant (int, float, double).
    matrix<T> operator*(const double val){
        matrix<T> temp(dim1, dim2);

        for(int i = 0; i < dim1 * dim2; i++){
            temp.mat[i] = mat[i] * val;
        }

        return temp;

    }

    //Index wise division of the matrix on the left by the matrix on the right.
    matrix<T> operator/=(const matrix<T> &Matrix){
        if((dim1 != Matrix.dim1) || (dim2 != Matrix.dim2)){
            printf("The matrices using index wise division are not of the same dimensions.\n");
            assert(0);
        }
        matrix<T> temp(dim1, dim2);

        for(int i = 0; i < dim1 * dim2; i++){
            temp.mat[i] = mat[i] / Matrix.mat[i];
        }

        return temp;

    }

    //Divide each entry in the matrix on the left by a 1x1 matrix.
    //For the equation x * A = b, x is the result of b * A^-1 or equally b / A.
    matrix<T> operator/(matrix<T> Matrix){
        //This is to ensure that the matrix to the right is either a 1x1 or isn't underdetermined.
        if(((Matrix.dim1 != 1) || (Matrix.dim2 != 1)) && ((Matrix.dim1 < Matrix.dim2) || (Matrix.dim1 != dim1))){
            printf("Either the two matrices are of improper sizes or the matrix on the right is not 1x1.\n");
            assert(0);
        }

        //Dividing by a 1x1 matrix will just divide everything in the matrix to the left by the value in the 1x1.
        if((Matrix.dim1 == 1) && (Matrix.dim2 == 1)) {
            matrix<T> temp = *this;

            for (int i = 0; i < dim1 * dim2; i++) {
                temp.mat[i] = temp.mat[i] / Matrix.mat[0];
            }

            return temp;
        }
        else{

            T upperSum;
            T lowerSum;

            if(Matrix.dim1 == Matrix.dim2){
                matrix<T> temp = *this;
                matrix<T> lowerTriangular(Matrix.dim1, Matrix.dim2);//Temporary lower triangular matrix
                matrix<T> upperTriangular(Matrix.dim1, Matrix.dim2);//Temporary upper triangular matrix
                int swappedRowNumber = -1;

                //This if statement exists only to swap rows and do a quick assessment on invertibility.
                if(Matrix.mat[0] == (T) 0){
                    T temp2;

                    for(int i = 1; i < Matrix.dim1; i++){
                        if(Matrix.mat[i * Matrix.dim2] != (T) 0){
                            swappedRowNumber = i;
                            break;
                        }
                    }

                    if(swappedRowNumber == -1){
                        printf("The given matrix is not invertible.\n");
                        assert(0);
                    }

                    for(int i = 0; i < Matrix.dim2; i++){
                        temp2 = Matrix.mat[i];
                        Matrix.mat[i] = Matrix.mat[i + Matrix.dim2 * swappedRowNumber];
                        Matrix.mat[i + Matrix.dim2 * swappedRowNumber] = temp2;
                    }

                    for(int i = 0; i < temp.dim2; i++){
                        temp2 = temp.mat[i];
                        temp.mat[i] = temp.mat[i + Matrix.dim2 * swappedRowNumber];
                        temp.mat[i + temp.dim2 * swappedRowNumber] = temp2;
                    }
                }

                //This for loop turns the lowerTriangularMatrix into an identity matrix.
                for(int i = 0; i < Matrix.dim2; i++){
                    lowerTriangular.mat[i + i * Matrix.dim2] = (T) 1;
                }



                //The following for loop generates the lowerTriangular matrix and the upperTriangular matrix.
                for(int iterator = 0; iterator < Matrix.dim2; iterator++){

                    //This for loop generates the upperTriangular matrix.
                    for(int whichColumn = iterator; whichColumn < Matrix.dim2; whichColumn++){
                        upperSum = 0;

                        for(int element = 0; element < iterator; element++){
                            upperSum = upperSum + upperTriangular.mat[whichColumn + element * Matrix.dim2] * lowerTriangular.mat[iterator * Matrix.dim2 + element];
                        }

                        upperTriangular.mat[whichColumn + iterator * Matrix.dim2] = Matrix.mat[iterator * Matrix.dim2 + whichColumn] - upperSum;
                        if(upperTriangular.mat[iterator + iterator * Matrix.dim2] == 0){
                            printf("\nThe given matrix is not invertible.\n");
                            assert(0);
                        }

                    }

                    //This for loop generates the lowerTriangular matrix.
                    for(int whichRow = iterator + 1; whichRow < Matrix.dim2; whichRow++){
                        lowerSum = 0;

                        for(int element = 0; element < iterator; element++){
                            lowerSum = lowerSum + upperTriangular.mat[whichRow + element * Matrix.dim2 - 1] * lowerTriangular.mat[whichRow * Matrix.dim2 + element];
                        }

                        lowerTriangular.mat[whichRow * Matrix.dim2 + iterator] = (Matrix.mat[iterator + whichRow * Matrix.dim2] - lowerSum) / upperTriangular.mat[iterator * Matrix.dim2 + iterator];

                    }

                }

                //The following for loop begins the solving process by substitution with the lower triangular matrix.
                for(int column = 0; column < Matrix.dim2; column++){
                    for(int column2 = 0; column2 < dim2; column2++){
                        for(int element = column + 1; element < Matrix.dim2; element++) {

                            temp.mat[element * dim2 + column2] = temp.mat[element * dim2 + column2] - lowerTriangular.mat[Matrix.dim2 * element + column] * temp.mat[column * dim2 + column2];

                        }
                    }
                }

                //The following for loop begins the solving process by substitution with the upper triangular matrix.
                for(int column = Matrix.dim2 - 1; column >= 0; column--){
                    for(int column2 = 0; column2 < dim2; column2++){
                        for(int element = 0 ; element < column; element++){

                            temp.mat[element * dim2 + column2] = (temp.mat[element * dim2 + column2] - upperTriangular.mat[Matrix.dim2 * element + column] / upperTriangular.mat[column * Matrix.dim2 + column] * temp.mat[column * dim2 + column2]);

                        }
                    }
                }

                //The following for loop finishes the work by dividing each row by each entry along the diagonal of upperTriangular.
                for(int i = 0; i < Matrix.dim2; i++){
                    for(int k = 0; k < Matrix.dim2; k++){
                        temp.mat[i * Matrix.dim2 + k] = temp.mat[i * Matrix.dim2 + k] / upperTriangular.mat[i + i * Matrix.dim2];
                    }
                }



                return temp;


            }else{
                matrix<T> temp = *this;

                //This code basically multiplies the matrix by its transpose; A = transpose(A) * A.
                temp = temp * ~Matrix;
                Matrix = Matrix * ~Matrix;

                matrix<T> lowerTriangular(Matrix.dim1, Matrix.dim2);//Temporary lower triangular matrix
                matrix<T> upperTriangular(Matrix.dim1, Matrix.dim2);//Temporary upper triangular matrix

                //This for loop turns the lowerTriangularMatrix into an identity matrix.
                for(int i = 0; i < Matrix.dim2; i++){
                    lowerTriangular.mat[i + i * Matrix.dim2] = (T) 1;
                }



                //The following for loop generates the lowerTriangular matrix and the upperTriangular matrix.
                for(int iterator = 0; iterator < Matrix.dim2; iterator++){

                    //This for loop generates the upperTriangular matrix.
                    for(int whichColumn = iterator; whichColumn < Matrix.dim2; whichColumn++){
                        upperSum = 0;

                        for(int element = 0; element < iterator; element++){
                            upperSum = upperSum + upperTriangular.mat[whichColumn + element * Matrix.dim2] * lowerTriangular.mat[iterator * Matrix.dim2 + element];
                        }

                        upperTriangular.mat[whichColumn + iterator * Matrix.dim2] = Matrix.mat[iterator * Matrix.dim2 + whichColumn] - upperSum;
                        if(upperTriangular.mat[iterator + iterator * Matrix.dim2] == 0){
                            printf("\nThe given matrix is not invertible.\n");
                            assert(0);
                        }

                    }

                    //This for loop generates the lowerTriangular matrix.
                    for(int whichRow = iterator + 1; whichRow < Matrix.dim2; whichRow++){
                        lowerSum = 0;

                        for(int element = 0; element < iterator; element++){
                            lowerSum = lowerSum + upperTriangular.mat[whichRow + element * Matrix.dim2 - 1] * lowerTriangular.mat[whichRow * Matrix.dim2 + element];
                        }

                        lowerTriangular.mat[whichRow * Matrix.dim2 + iterator] = (Matrix.mat[iterator + whichRow * Matrix.dim2] - lowerSum) / upperTriangular.mat[iterator * Matrix.dim2 + iterator];

                    }

                }

                //The following for loop begins the solving process by substitution with the lower triangular matrix.
                for(int column = 0; column < temp.dim2; column++){
                    for(int column2 = 0; column2 < temp.dim2; column2++){
                        for(int element = column + 1; element < temp.dim2; element++){

                            temp.mat[element * temp.dim2 + column2] = temp.mat[element * temp.dim2 + column2] - lowerTriangular.mat[temp.dim2 * element + column] * temp.mat[column * temp.dim2 + column2];

                        }
                    }
                }

                //The following for loop begins the solving process by substitution with the upper triangular matrix.
                for(int column = temp.dim2 - 1; column >= 0; column--){
                    for(int column2 = 0; column2 < temp.dim2; column2++){
                        for(int element = 0 ; element < column; element++){

                            temp.mat[element * temp.dim2 + column2] = (temp.mat[element * temp.dim2 + column2] - upperTriangular.mat[Matrix.dim2 * element + column] / upperTriangular.mat[column * Matrix.dim2 + column] * temp.mat[column * temp.dim2 + column2]);

                        }
                    }
                }

                //The following for loop finishes the work by dividing each row by each entry along the diagonal of upperTriangular.
                for(int i = 0; i < temp.dim2; i++){
                    for(int k = 0; k < temp.dim2; k++){
                        temp.mat[i * temp.dim2 + k] = temp.mat[i * temp.dim2 + k] / upperTriangular.mat[i + i * Matrix.dim2];
                    }
                }



                return temp;

            }


        }

    }

    //Divide each entry in the matrix on the left by a constant (int, float, double).
    matrix<T> operator/(const double val){
        matrix<T> temp(dim1, dim2);

        for(int i = 0; i < dim1 * dim2; i++){

            temp.mat[i] = mat[i] / val;

        }

        return temp;

    }

    //Add two matrices together.
    matrix<T> operator+(const matrix<T> &Matrix){
        if((dim1 != Matrix.dim1) || (dim2 != Matrix.dim2)){
            printf("Matrices being added are of different dimensions.\n");
            assert(0);
        }
        matrix<T> temp(dim1, dim2);

        for(int i = 0; i < dim1 * dim2; i++){
            temp.mat[i] = mat[i] + Matrix.mat[i];
        }

        return temp;

    }

    //Find the difference between two matrices of the same size.
    matrix<T> operator-(const matrix<T> &Matrix){
        if((dim1 != Matrix.dim1) || (dim2 != Matrix.dim2)){
            printf("Matrices being subtracted are of different dimensions.\n");
            assert(0);
        }
        matrix<T> temp(dim1, dim2);

        for(int i = 0; i < dim1 * dim2; i++){
            temp.mat[i] = mat[i] - Matrix.mat[i];
        }

        return temp;

    }

    //For the equation A * x = b, x is the result of A^-1 * b.
    matrix<T> operator|(matrix<T> Matrix){


        T upperSum;
        T lowerSum;

        if(dim1 == dim2){
            int swappedRowNumber = -1;
            matrix<T> lowerTriangular(dim2, dim2);//Temporary lower triangular matrix
            matrix<T> upperTriangular(dim2, dim2);//Temporary upper triangular matrix


            //This for loop turns the lowerTriangularMatrix into an identity matrix.
            for(int i = 0; i < dim2; i++){
                lowerTriangular.mat[i + i * dim2] = (T) 1;
            }

            //This if statement exists only to swap rows and do a quick assessment on invertibility.
            if(mat[0] == (T) 0){
                T temp;

                for(int i = 1; i < dim1; i++){
                    if(mat[i * dim2] != 0){
                        swappedRowNumber = i;
                        break;
                    }
                }

                if(swappedRowNumber == -1){
                    assert(0);
                }

                for(int i = 0; i < dim2; i++){
                    temp = mat[i];
                    mat[i] = mat[i + dim2 * swappedRowNumber];
                    mat[i + dim2 * swappedRowNumber] = temp;
                }

                for(int i = 0; i < dim2; i++){
                    temp = Matrix.mat[i];
                    Matrix.mat[i] = Matrix.mat[i + dim2 * swappedRowNumber];
                    Matrix.mat[i + dim2 * swappedRowNumber] = temp;
                }
            }

            //The following for loop generates the lowerTriangular matrix and the upperTriangular matrix.
            for(int iterator = 0; iterator < dim2; iterator++){

                //This for loop generates the upperTriangular matrix.
                for(int whichColumn = iterator; whichColumn < dim2; whichColumn++){
                    upperSum = 0;

                    for(int element = 0; element < iterator; element++){
                        upperSum = upperSum + upperTriangular.mat[whichColumn + element * dim2] * lowerTriangular.mat[iterator * dim2 + element];
                    }

                    upperTriangular.mat[whichColumn + iterator * dim2] = mat[iterator * dim2 + whichColumn] - upperSum;
                    if(upperTriangular.mat[iterator + iterator * dim2] == 0){
                        printf("\nThe given matrix is not invertible.\n");
                        assert(0);
                    }

                }

                //This for loop generates the lowerTriangular matrix.
                for(int whichRow = iterator + 1; whichRow < dim2; whichRow++){
                    lowerSum = 0;

                    for(int element = 0; element < iterator; element++){
                        lowerSum = lowerSum + upperTriangular.mat[whichRow - 1 + element * dim2] * lowerTriangular.mat[whichRow * dim2 + element];
                    }

                    lowerTriangular.mat[whichRow * dim2 + iterator] = (mat[iterator + whichRow * dim2] - lowerSum) / upperTriangular.mat[iterator * dim2 + iterator];

                }

            }

            //The following for loop begins the solving process by substitution with the lower triangular matrix.
            for(int column = 0; column < dim2; column++){
                for(int column2 = 0; column2 < Matrix.dim2; column2++){
                    for(int element = column + 1; element < dim2; element++){

                        Matrix.mat[element * Matrix.dim2 + column2] = Matrix.mat[element * Matrix.dim2 + column2] - lowerTriangular.mat[dim2 * element + column] * Matrix.mat[column * Matrix.dim2 + column2];

                    }
                }
            }

            //The following for loop begins the solving process by substitution with the upper triangular matrix.
            for(int column = dim2 - 1; column >= 0; column--){
                for(int column2 = 0; column2 < dim2; column2++){
                    for(int element = 0; element < column; element++){

                        Matrix.mat[element * Matrix.dim2 + column2] = (Matrix.mat[element * Matrix.dim2 + column2] - upperTriangular.mat[dim2 * element + column] / upperTriangular.mat[column * dim2 + column] * Matrix.mat[column * Matrix.dim2 + column2]);

                    }
                }
            }

            //The following for loop finishes the work by dividing each row by each entry along the diagonal of upperTriangular.
            for(int i = 0; i < dim2; i++){
                for(int k = 0; k < dim2; k++){
                    Matrix.mat[i * dim2 + k] = Matrix.mat[i * dim2 + k] / upperTriangular.mat[i + i * dim2];
                }
            }



            //This for loop is to swap the pivoted row back, if pivoting occurred.
            if(swappedRowNumber != -1){
                T temp;

                for(int i = 0; i < dim1; i++){
                    temp = mat[i];
                    mat[i] = mat[i + dim2 * swappedRowNumber];
                    mat[i + dim2 * swappedRowNumber] = temp;
                }

            }



            return Matrix;

        }else{
            int swappedRowNumber = -1;
            matrix<T> temp = *this;

            //This code basically multiplies the matrix by its transpose; A = transpose(A) * A.
            Matrix = ~temp * Matrix;
            temp = ~temp * temp;

            matrix<T> lowerTriangular(dim2, dim2);//Temporary lower triangular matrix
            matrix<T> upperTriangular(dim2, dim2);//Temporary upper triangular matrix

            //This if statement exists only to swap rows and do a quick assessment on invertibility.
            if(temp.mat[0] == (T) 0){
                T temp2;

                for(int i = 1; i < temp.dim1; i++){
                    if(temp.mat[i * temp.dim2] != 0){
                        swappedRowNumber = i;
                        break;
                    }
                }

                if(swappedRowNumber == -1){
                    assert(0);
                }

                for(int i = 0; i < temp.dim2; i++){
                    temp2 = temp.mat[i];
                    temp.mat[i] = temp.mat[i + temp.dim2 * swappedRowNumber];
                    temp.mat[i + temp.dim2 * swappedRowNumber] = temp2;
                }

                for(int i = 0; i < temp.dim2; i++){
                    temp2 = Matrix.mat[i];
                    Matrix.mat[i] = Matrix.mat[i + temp.dim2 * swappedRowNumber];
                    Matrix.mat[i + temp.dim2 * swappedRowNumber] = temp2;
                }
            }



            //This for loop turns the lowerTriangularMatrix into an identity matrix.
            for(int i = 0; i < temp.dim2; i++){
                lowerTriangular.mat[i + i * temp.dim2] = (T) 1;
            }

            //The following for loop generates the lowerTriangular matrix and the upperTriangular matrix.
            for(int iterator = 0; iterator < temp.dim2; iterator++){

                //This for loop generates the upperTriangular matrix.
                for(int whichColumn = iterator; whichColumn < temp.dim2; whichColumn++){
                    upperSum = 0;

                    for(int element = 0; element < iterator; element++){
                        upperSum = upperSum + upperTriangular.mat[whichColumn + element * temp.dim2] * lowerTriangular.mat[iterator * temp.dim2 + element];
                    }

                    upperTriangular.mat[whichColumn + iterator * temp.dim2] = temp.mat[iterator * temp.dim2 + whichColumn] - upperSum;
                    if(upperTriangular.mat[iterator + iterator * temp.dim2] == 0){
                        printf("\nThe given matrix is not invertible.\n");
                        assert(0);
                    }

                }

                //This for loop generates the lowerTriangular matrix.
                for(int whichRow = iterator + 1; whichRow < temp.dim2; whichRow++){
                    lowerSum = 0;

                    for(int element = 0; element < iterator; element++){
                        lowerSum = lowerSum + upperTriangular.mat[whichRow - 1 + element * temp.dim2] * lowerTriangular.mat[whichRow * temp.dim2 + element];
                    }

                    lowerTriangular.mat[whichRow * temp.dim2 + iterator] = (temp.mat[iterator + whichRow * temp.dim2] - lowerSum) / upperTriangular.mat[iterator * temp.dim2 + iterator];

                }

            }

            //The following for loop begins the solving process by substitution with the lower triangular matrix.
            for(int column = 0; column < temp.dim2; column++){
                for(int column2 = 0; column2 < Matrix.dim2; column2++){
                    for(int element = column + 1; element < temp.dim2; element++){

                        Matrix.mat[element * Matrix.dim2 + column2] = Matrix.mat[element * Matrix.dim2 + column2] - lowerTriangular.mat[temp.dim2 * element + column] * Matrix.mat[column * Matrix.dim2 + column2];

                    }
                }
            }

            //The following for loop begins the solving process by substitution with the upper triangular matrix.
            for(int column = temp.dim2 - 1; column >= 0; column--){
                for(int column2 = 0; column2 < temp.dim2; column2++){
                    for(int element = 0; element < column; element++){

                        Matrix.mat[element * Matrix.dim2 + column2] = (Matrix.mat[element * Matrix.dim2 + column2] - upperTriangular.mat[temp.dim2 * element + column] / upperTriangular.mat[column * temp.dim2 + column] * Matrix.mat[column * Matrix.dim2 + column2]);

                    }
                }
            }

            //The following for loop finishes the work by dividing each row by each entry along the diagonal of upperTriangular.
            for(int i = 0; i < temp.dim2; i++){
                for(int k = 0; k < temp.dim2; k++){
                    Matrix.mat[i * temp.dim2 + k] = Matrix.mat[i * temp.dim2 + k] / upperTriangular.mat[i + i * temp.dim2];
                }
            }



            return Matrix;

        }

    }

    //Exponentiate a row/column vector with a column/row vector to obtain a rectangular matrix of size (dim1, Matrix.dim2) or (Matrix.dim1, dim2).
    matrix<T> operator^(matrix<T> Matrix){
        if(((dim1 != 1) || (Matrix.dim2 != 1)) && ((Matrix.dim1 != 1) || (dim2 != 1))){
            printf("The given matrices are not the same dimensions, they are not compatible, or the matrix on the right is not 1x1.\n");
            assert(0);
        }

        if((dim1 == 1) && (Matrix.dim2 == 1)){
            matrix<T> other(Matrix.dim1, dim2);

            for(int i = 0; i < Matrix.dim1; i++){
                for(int k = 0; k < dim2; k++){

                    other.mat[i * dim2 + k] = std::pow(mat[k], Matrix.mat[i]);

                }
            }

            return other;

        }
        else if ((dim2 == 1) && (Matrix.dim1 == 1)){
            matrix<T> other(dim1, Matrix.dim2);

            for(int i = 0; i < dim1; i++){
                for(int k = 0; k < Matrix.dim2; k++){

                    other.mat[i * Matrix.dim2 + k] = std::pow(mat[i], Matrix.mat[k]);

                }
            }

            return other;

        }
        else{
            matrix<T> other(dim1, dim2);
            
            for(int i = 0; i < dim1 * dim2; i++){
                
                other.mat[i] = std::pow(other.mat[i], Matrix.mat[0]);
                
            }
            
            return other;
            
        }



    }

    //Index wise exponentiation of the matrix on the left by the matrix on the right.
    matrix<T> operator^=(const matrix<T> &Matrix){
        if((dim1 != Matrix.dim1) || (dim2 != Matrix.dim2)){
            printf("The matrices being used for exponentiation do not share the proper dimensions.\n");
            assert(0);
        }
        matrix<T> temp;

        for(int i = 0; i < dim1 * dim2; i++){

            temp.mat[i] = std::pow(mat[i], Matrix.mat[i]);

        }

        return temp;

    }




};

namespace calculus{

    ////vector stuff

    template<typename T>
    void array(T (*fun)(T), std::vector<T> &arr, const T delta = 0.1, const T start = 0){

        for(int i = 0; i < arr.size(); i++){
            arr[i] = fun(start + i * delta);
        }

    }

    template<typename T1, typename T2>
    void function(T1 (*fun)(T1), std::vector<T2> &arr){

        for(int i = 0; i < arr.size(); i++){
            arr[i] = fun(arr[i]);
        }

    }

    template<typename T1, typename T2>
    matrix<double> row(T1 start, T2 end, double delta = 1.0){
        double iterations = ((end - start) / delta) + 1;
        matrix temporaryMatrix(1, (int) iterations);

        for(int i = 0; i < (int) iterations; i++){
            temporaryMatrix.mat[i] = start + i * delta;
        }

        return temporaryMatrix;

    }

    template<typename T1, typename T2>
    matrix<double> column(T1 start, T2 end, double delta = 1.0){
        double iterations = ((end - start) / delta) + 1;
        matrix temporaryMatrix(((int) iterations), 1);

        for(int i = 0; i < (int) iterations; i++){
            temporaryMatrix.mat[i] = start + i * delta;
        }

        return temporaryMatrix;

    }

    template<typename T>
    std::vector<T> difference(std::vector<T> &arr, int depth = 1){

        for(int i = 0; i <= depth; i++) {
            for (int k = 0; k < arr.size() - 1; k++) {

                arr[k] = arr[k + 1] - arr[k];

            }
            arr.resize(arr.size() - 1);
        }
    }

    template<typename T>
    T dotProduct(std::vector<T> &arr1, std::vector<T> &arr2){
        assert(arr1.size() == arr2.size());
        T sum = 0;

        for(int i = 0; i < arr1.size(); i++){
            sum += arr1[i] * arr2[i];
        }

        return sum;

    }

    template<typename T>
    T area(std::vector<T> &arr1, std::vector<T> &arr2){
        assert(arr1.size() == arr2.size());
        T sum = 0;
        T val;

        for(int i = 0; i < arr1.size(); i++){
            for(int k = i + 1; k < arr2.size(); k++){
                val = arr1[i] * arr2[k] - arr1[k] * arr2[i];
                sum += val * val;
            }
        }

        return sqrtl(sum);

    }

    template<typename T>
    T magnitude(std::vector<T> &arr){
        T sum = 0;

        for(int i = 0; i < arr.size(); i++){
            sum += arr[i] * arr[i];
        }

        return sqrtl(sum);

    }




    ////Differentiation functions

        ////function based derivatives

    template<typename T1, typename T2>
    T1 derivative2_1(T1 (*fun)(T1), const T2 point, const T1 delta = 0.1) {

        return (fun(point + 0.5 * delta) - fun(point - 0.5 * delta)) / delta;

    }

    template<typename T1, typename T2>
    T1 derivative2_2(T1 (*fun)(T1), const T2 point, const T1 delta = 0.1) {

        return (fun(point + delta) - 2 * fun(point) + fun(point - delta)) / (delta * delta);

    }

    template<typename T1, typename T2>
    T1 derivative2_3(T1 (*fun)(T1), const T2 point, const T1 delta = 0.1) {

        return (fun(point + 1.5 * delta) - 3 * fun(point + 0.5 * delta) + 3 * fun(point - 0.5 * delta) - fun(point - 1.5 * delta)) / (delta * delta * delta);

    }

    template<typename T1, typename T2>
    T1 derivative2_4(T1 (*fun)(T1), const T2 point, const T1 delta = 0.1) {

        return (fun(point + 2 * delta) - 4 * fun(point + delta) + 6 * fun(point) - 4 * fun(point - delta) + fun(point - 2 * delta)) / (delta * delta * delta * delta);

    }

////In its current form, derivative2_N is very ineffective for derivatives of `num` greater than two. I really need a better method.
    template<typename T1, typename T2>
    T1 derivative2_N(T1 (*fun)(T1), const T2 point, const int num, const T1 delta = 0.1){
        T1 temp[num + 1];

        for(int i = 0; i < num + 1; i++){
            temp[i] = fun(point - 0.5 * num * delta - 0.5 * delta + i * delta);
        }

        for(int i = num; i > 0; i--){

            for(int k = 0; k < i + 1; k++){
                printf("%d\n", k + 1);
                temp[k] = (temp[k + 1] - temp[k]) / delta;
            }

        }

        return temp[0];

    }

    template<typename T1, typename T2>
    T1 derivative3_1(T1 (*fun)(T1), const T2 point, const T1 delta = 0.1) {

        //I modified the proof to Simpson's rule so that it may be extended to derivatives.
        //These coefficients will produce the approximate derivative of a 3rd degree interpolation.

        const T1 c1 = T1 (1) / T1 (24);
        const T1 c2 = T1 (9) / T1 (8);

        return (-c1 * fun(point + 1.5 * delta) + c2 * fun(point + 0.5 * delta) - c2 * fun(point - 0.5 * delta) + c1 * fun(point - 1.5 * delta)) / delta;
    }

    template<typename T1, typename T2>
    T1 derivative3_2(T1 (*fun)(T1), const T2 point, const T1 delta = 0.1) {

        //These coefficients are the same as nesting derivative3_1 two times.

        const T1 c1 = T1 (1) / T1(576);
        const T1 c2 = T1 (3) / T1 (32);
        const T1 c3 = T1 (87) / T1 (64);
        const T1 c4 = T1 (365) / T1 (144);

        return (c1 * fun(point + 3 * delta) - c2 * fun(point + 2 * delta) + c3 * fun(point + delta) - c4 * fun(point) + c3 * fun(point - delta) - c2 * fun(point - 2 * delta) + c1 * fun(point - 3 * delta)) / (delta * delta);
    }

    template<typename T1, typename T2>
    T1 derivative3_3(T1 (*fun)(T1), const T2 point, const T1 delta = 0.1) {

        //These coefficients are the same as nesting der1 three times.

        const T1 c1 = T1 (1) / T1 (13824);
        const T1 c2 = T1 (3) / T1 (512);
        const T1 c3 = T1 (21) / T1 (128);
        const T1 c4 = T1 (2005) / T1 (1152);
        const T1 c5 = T1 (1137) / T1 (256);

        return (-c1 * fun(point + 4.5 * delta) + c2 * fun(point + 3.5 * delta) - c3 * fun(point + 2.5 * delta) + c4 * fun(point + 1.5 * delta) - c5 * fun(point + 0.5 * delta) + c5 * fun(point - 0.5 * delta) - c4 * fun(point - 1.5 * delta) + c3 * fun(point - 2.5 * delta) - c2 * fun(point - 3.5 * delta) + c1 * fun(point - 4.5 * delta)) / (delta * delta * delta);
    }

    template<typename T1, typename T2>
    T1 derivative3_4(T1 (*fun)(T1), const T2 point, const T1 delta = 0.1) {

        //These coefficients are the same as nesting der1 four times.

        const T1 c1 = T1 (1) / T1 (331776);
        const T1 c2 = T1 (1) / T1 (3072);
        const T1 c3 = T1 (83) / T1 (6144);
        const T1 c4 = T1 (21871) / T1 (82944);
        const T1 c5 = T1 (9535) / T1 (4096);
        const T1 c6 = T1 (3659) / T1 (512);
        const T1 c7 = T1 (280301) / T1 (27648);

        return (c1 * fun(point + 6 * delta) - c2 * fun(point + 5 * delta) + c3 * fun(point + 4 * delta) - c4 * fun(point + 3 * delta) + c5 * fun(point + 2 * delta) - c6 * fun(point + delta) + c7 * fun(point) - c6 * fun(point - delta) + c5 * fun(point - 2 * delta) - c4 * fun(point - 3 * delta) + c3 * fun(point - 4 * delta) - c2 * fun(point - 5 * delta) + c1 * fun(point - 6 * delta)) / (delta * delta * delta * delta);
    }

////In its current form, derivative3_N is ineffective for derivatives of 'num' greater than six. I need a new method.
    template<typename T1, typename T2>
    T1 derivative3_N(T1 (*fun)(T1), const T2 point, const int num, const T1 delta = 0.1){
        T1 temp[3 * num + 1];
        const T1 c1 = T1 (1) / T1 (24);
        const T1 c2 = T1 (9) / T1 (8);

        for(int i = 0; i < 3 * num + 1; i++){
            temp[i] = fun(point - (1.5 * num) * delta + i * delta);
        }

        printf("\n");

        for(int i = num; i > 0; i--){

            for(int k = 0; k < 3 * i - 2; k++){
                temp[k] = (-c1 * temp[k + 3] + c2 * temp[k + 2] - c2 * temp[k + 1] + c1 * temp[k]) / delta;
            }


        }

        return temp[0];

    }



    ////Integration functions

        ////function based integrals

    template<typename T1, typename T2, typename T3>
    T1 integral2_1(T1 (*fun)(T1), const T2 low, const T3 up, const T1 delta = 0.1){
        T1 sum = 0;
        const T1 val = up - low;
        const int ratio = val / delta;
        const T1 diff = val - ratio * delta;

        //These constants are derived from Simpson's rule, when slightly modified so that points on a parabola are spaced by `delta`
        const T1 c1 = T1 (1) / T1 (6);
        const T1 c2 = T1 (2) / T1 (3);

        T1 old1 = c1 * fun(low);

        for(int i = 0; i < ratio; i++){

            sum += old1 + c2 * fun(low + (i + 0.5) * delta);

            old1 = c1 * fun(low + (i + 1) * delta);

            sum += old1;

        }

        return delta * sum + diff * (c1 * fun(up - diff) + c2 * fun(up - 0.5 * diff) + c1 * fun(up));

    }

    template<typename T1, typename T2, typename T3>
    T1 integral2_2(T1 (*fun)(T1), const T2 low, const T3 up, const T1 delta = 0.1){
        T1 sum = 0;
        const T1 val = up - low;
        const int ratio = val / delta;
        const T1 diff = val - ratio * delta;

        const T1 c1 = T1 (1) / T1 (6);
        const T1 c2 = T1 (2) / T1 (3);

        T1 old1 = val * c1 * fun(low);

        for(int i = 0; i < ratio; i ++) {

            sum += old1 + c2 * (val - (i + 0.5) * delta) * fun(low + (i + 0.5) * delta);

            old1 = c1 * (val - (i + 1) * delta) * fun(low + (i + 1) * delta);

            sum += old1;

        }

        return delta * sum + diff * (c1 * diff * fun(up - diff) + T1 (1) / T1 (3) * diff * fun(up - 0.5 * diff));

    }

    template<typename T1, typename T2, typename T3>
    T1 integral2_3(T1 (*fun)(T1), const T2 low, const T3 up, const T1 delta = 0.1){
        T1 sum = 0;
        const T1 val = up - low;
        const int ratio = val / delta;
        const T1 diff = val - ratio * delta;

        T1 old1 = val * 0.08333333333333333 * fun(low);

        for(int i = 0; i < ratio; i ++) {

            sum += old1 + 0.3333333333333333 * (val - (i + 0.5) * delta) * (val - (i + 0.5) * delta) * fun(low + (i + 0.5) * delta);

            old1 = 0.08333333333333333 * (val - (i + 1) * delta) * (val - (i + 1) * delta) * fun(low + (i + 1) * delta);

            sum += old1;

        }

        return delta * sum + diff * (0.08333333333333333 * diff * diff * fun(up - diff) + 0.1666666666666667 * diff * diff * fun(up - 0.5 * diff));

    }



    template<typename T1, typename T2, typename T3>
    T1 integral3_1(T1 (*fun)(T1), const T2 low, const T3 up, const T1 delta = 0.1){
        T1 sum = 0;
        const int length = (up - low) / delta;
        const T1 diff = up - low - length * delta;

        T1 old1 = 0.125 * fun(low);

        for(int i = 0; i < length; i++){

            sum += old1 + 0.375 * fun(low + (i + 0.3333333333333333) * delta) + 0.375 * fun(low + (i + 0.6666666666666667) * delta);

            old1 = 0.125 * fun(low + (i + 1) * delta);

            sum += old1;

        }

        return delta * sum + diff * (0.125 * fun(up - diff) + 0.375 * fun(up - 0.6666666666666667 * diff) + 0.375 * fun(up - 0.3333333333333333 * diff) + 0.125 * fun(up));

    }

    template<typename T1, typename T2, typename T3>
    T1 integral3_2(T1 (*fun)(T1), const T2 low, const T3 up, const T1 delta = 0.1){
        T1 sum = 0;
        const T1 val = up - low;
        const int ratio = val / delta;
        const T1 diff = val - ratio * delta;

        T1 old1 = val * 0.125 * fun(low);

        for(int i = 0; i < ratio; i ++) {

            sum += old1 + 0.375 * (val - (i + 0.3333333333333333) * delta) * fun(low + (i + 0.3333333333333333) * delta) + 0.375 * (val - (i + 0.6666666666666667) * delta) * fun(low + (i + 0.6666666666666667) * delta);

            old1 = 0.125 * (val - (i + 1) * delta) * fun(low + (i + 1) * delta);

            sum += old1;

        }

        return delta * sum + diff * (0.125 * diff * fun(up - diff) + 0.25 * diff * fun(up - 0.6666666666666667 * diff) + 0.125 * diff * fun(up - 0.3333333333333333 * diff));

    }

    template<typename T1, typename T2, typename T3>
    T1 integral3_3(T1 (*fun)(T1), const T2 low, const T3 up, const T1 delta = 0.1){
        T1 sum = 0;
        const T1 val = up - low;
        const int ratio = val / delta;
        const T1 diff = val - ratio * delta;

        T1 coe1;
        T1 coe2;
        T1 coe3;

        T1 old1 = val * 0.0625 * fun(low);

        for(int i = 0; i < ratio; i ++) {

            coe1 = (val - (i + 0.3333333333333333) * delta);
            coe2 = (val - (i + 0.3333333333333333) * delta);
            coe3 = (val - (i + 1) * delta);

            sum += old1 +
                   0.1875 * coe1 * coe1 * fun(low + (i + 0.3333333333333333) * delta) +
                   0.1875 * coe2 * coe2 * fun(low + (i + 0.6666666666666667) * delta);

            old1 = 0.0625 * coe3 * coe3 * fun(low + (i + 1) * delta);

            sum += old1;

        }

        return delta * sum + diff * (0.0625 * diff * diff * fun(up - diff) + 0.125 * diff * diff * fun(up - 0.6666666666666667 * diff) + 0.0625 * diff * diff * fun(up - 0.3333333333333333 * diff));

    }



    ////ODE solvers (non-DSE)
    //These basic implementations are not as effective as I thought, I need a new method for these.

    template<typename T1, typename T2>
    void ode1_1(void (*fun)(std::vector<T1>, T1), std::vector<T1> initialConditions, T2 end, double delta = 0.1){
        // initialConditions[3] = {y(t0), y'(t0), t0}
        T1 change;
        long int iterations = std::abs(end - initialConditions[2]) / delta;

        if(initialConditions[2] < end){
            change = std::abs(delta);
        }
        else{
            change = -std::abs(delta);
        }

        for(int i = 0; i <= iterations ; i++){
            initialConditions[0] += change * initialConditions[1];
            initialConditions[2] += change;
            fun(initialConditions, change);
        }

    }

    template<typename T1, typename T2>
    void ode1_2(void (*fun)(std::vector<T1>, T1), std::vector<T1> initialConditions, T2 end, double delta = 0.1){
        // initialConditions[4] = {y(t0), y'(t0), y''(t0), t0}
        T1 change;
        long int iterations = std::abs(end - initialConditions[3]) / delta;

        if(initialConditions[3] < end){
            change = std::abs(delta);
        }
        else{
            change = -std::abs(delta);
        }

        for(int i = 0; i <= iterations ; i++){
            initialConditions[0] += change * initialConditions[1];
            initialConditions[1] += change * initialConditions[2];
            initialConditions[3] += change;
            fun(initialConditions, change);
        }

    }

    template<typename T1, typename T2>
    void ode2_2(void (*fun)(std::vector<T1>, T1), std::vector<T1> initialConditions, T2 end, double delta = 0.1){
        // initialConditions[4] = {y(t0), y'(t0), y''(t0), t0}
        T1 change;
        long int iterations = std::abs(end - initialConditions[3]) / delta;

        if(initialConditions[3] < end){
            change = std::abs(delta);
        }
        else{
            change = -std::abs(delta);
        }

        T1 del[2] = {change, change * change / T1(2)};

        for(int i = 0; i <= iterations; i++){
            initialConditions[0] += change * initialConditions[1] + del[1] * initialConditions[2];
            initialConditions[1] += change * initialConditions[2];
            initialConditions[3] += change;
            fun(initialConditions, change);
        }

    }

    template<typename T1, typename T2>
    void ode1_3(void (*fun)(std::vector<T1>, T1), std::vector<T1> initialConditions, T2 end, double delta = 0.1){
        // initialConditions[5] = {y(t0), y'(t0), y''(t0), y'''(t0), t0}
        T1 change;
        long int iterations = std::abs(end - initialConditions[4]) / delta;

        if(initialConditions[4] < end){
            change = std::abs(delta);
        }
        else{
            change = -std::abs(delta);
        }

        for(int i = 0; i <= iterations; i++){
            initialConditions[0] += change * initialConditions[1];
            initialConditions[1] += change * initialConditions[2];
            initialConditions[2] += change * initialConditions[3];
            initialConditions[4] += change;
            fun(initialConditions, change);
        }

    }

    template<typename T1, typename T2>
    void ode2_3(void (*fun)(std::vector<T1>, T1), std::vector<T1> initialConditions, T2 end, double delta = 0.1){
        // initialConditions[5] = {y(t0), y'(t0), y''(t0), y'''(t0), t0}
        T1 change;
        long int iterations = std::abs(end - initialConditions[4]) / delta;

        if(initialConditions[4] < end){
            change = std::abs(delta);
        }
        else{
            change = -std::abs(delta);
        }

        T1 del[2] = {change, change * change / T1(2)};

        for(int i = 0; i <= iterations; i++){
            initialConditions[0] += del[0] * initialConditions[1] + del[1] * initialConditions[2];
            initialConditions[1] += del[0] * initialConditions[2] + del[1] * initialConditions[3];
            initialConditions[2] += del[0] * initialConditions[3];
            initialConditions[4] += change;
            fun(initialConditions, change);
        }

    }

    template<typename T1, typename T2>
    void ode3_3(void (*fun)(std::vector<T1>, T1), std::vector<T1> initialConditions, T2 end, double delta = 0.1){
        // initialConditions[5] = {y(t0), y'(t0), y''(t0), y'''(t0), t0}
        T1 change;
        long int iterations = std::abs(end - initialConditions[4]) / delta;

        if(initialConditions[4] < end){
            change = std::abs(delta);
        }
        else{
            change = -std::abs(delta);
        }

        T1 del[3] = {change, change * change / T1(2), change * change * change / T1(6)};

        for(int i = 0; i <= iterations; i++){
            initialConditions[0] += del[0] * initialConditions[1] + del[1] * initialConditions[2] + del[2] * initialConditions[3];
            initialConditions[1] += del[0] * initialConditions[2] + del[1] * initialConditions[3];
            initialConditions[2] += del[0] * initialConditions[3];
            initialConditions[4] += del[0];
            fun(initialConditions, del[0]);
        }

    }





    ////matrix stuff

    matrix<double> identityMatrix(int size){
        matrix temp(size, size);

        for(int i = 0; i < size; i++){
            for(int j = 0; j < size; j++){
                temp.mat[i + j * size] = (i == j);
            }
        }

        return temp;

    }

    template<typename T>
    void printMatrixf(matrix<T> Matrix){

        for(int i = 0; i < Matrix.dim1; i++){
            for(int k = 0; k < Matrix.dim2; k++) {

                std::cout << std::fixed << Matrix.mat[k + i * Matrix.dim2] << " ";

            }
            printf("\n");
        }

    }

    template<typename T>
    void printMatrixs(matrix<T> Matrix){

        for (int i = 0; i < Matrix.dim1; i++) {
            for (int k = 0; k < Matrix.dim2; k++) {

                std::cout << std::scientific << Matrix.mat[k + i * Matrix.dim2] << " ";

            }
            printf("\n");
        }

    }

    template<typename T>
    matrix<T> function(T (*fun)(T), matrix<T> Matrix){

        for(int i = 0; i < Matrix.dim1 * Matrix.dim2; i++){
            Matrix.mat[i] = fun(Matrix.mat[i]);
        }

        return Matrix;
    }

    template<typename T>
    void LUPDecomposition(matrix<T> Matrix, matrix<T> &lowerTriangular, matrix<T> &upperTriangular, matrix<T> &permutation){
        assert(Matrix.dim1 == Matrix.dim2);
        if( (lowerTriangular.dim1 != Matrix.dim1) || (lowerTriangular.dim2 != Matrix.dim2) ){
            lowerTriangular.resize(Matrix.dim1, Matrix.dim2);
        }
        if( (upperTriangular.dim1 != Matrix.dim1) || (upperTriangular.dim2 != Matrix.dim2) ){
            upperTriangular.resize(Matrix.dim1, Matrix.dim2);
        }
        if( (permutation.dim1 != Matrix.dim1) || (permutation.dim2 != Matrix.dim2) ){
            permutation.resize(Matrix.dim1, Matrix.dim2);
        }

        int swappedRowNumber = -1;

        for(int i = 0; i < Matrix.dim2; i++){
            lowerTriangular.mat[i + i * Matrix.dim2] = (T) 1;
        }

        //The following if statement and for loop check for invertibility and pivots only as necessary.
        if(Matrix.mat[0] == (T) 0){
            T temp2;

            for(int i = 1; i < Matrix.dim1; i++){
                if(Matrix.mat[i * Matrix.dim2] != (T) 0){
                    swappedRowNumber = i;
                    break;
                }
            }

            if(swappedRowNumber == -1){
                printf("The given matrix is not invertible.\n");
                assert(0);
            }

            for(int i = 0; i < Matrix.dim2; i++){
                temp2 = Matrix.mat[i];
                Matrix.mat[i] = Matrix.mat[i + Matrix.dim2 * swappedRowNumber];
                Matrix.mat[i + Matrix.dim2 * swappedRowNumber] = temp2;
            }
        }

        T upperSum;
        T lowerSum;

        //The following for loop generates the lowerTriangular matrix and the upperTriangular matrix.
        for(int iterator = 0; iterator < Matrix.dim2; iterator++){

            //This for loop generates the upperTriangular matrix.
            for(int whichColumn = iterator; whichColumn < Matrix.dim2; whichColumn++){
                upperSum = 0;

                for(int element = 0; element < iterator; element++){
                    upperSum = upperSum + upperTriangular.mat[whichColumn + element * Matrix.dim2] * lowerTriangular.mat[iterator * Matrix.dim2 + element];
                }

                upperTriangular.mat[whichColumn + iterator * Matrix.dim2] = Matrix.mat[iterator * Matrix.dim2 + whichColumn] - upperSum;
                if(upperTriangular.mat[iterator + iterator * Matrix.dim2] == 0){
                    printf("\nThe given matrix is not invertible.\n");
                    assert(0);
                }

            }

            //This for loop generates the lowerTriangular matrix.
            for(int whichRow = iterator + 1; whichRow < Matrix.dim2; whichRow++){
                lowerSum = 0;

                for(int element = 0; element < iterator; element++){
                    lowerSum = lowerSum + upperTriangular.mat[whichRow + element * Matrix.dim2 - 1] * lowerTriangular.mat[whichRow * Matrix.dim2 + element];
                }

                lowerTriangular.mat[whichRow * Matrix.dim2 + iterator] = (Matrix.mat[iterator + whichRow * Matrix.dim2] - lowerSum) / upperTriangular.mat[iterator * Matrix.dim2 + iterator];

            }

        }

        for(int i = 0; i < Matrix.dim2; i++){
            permutation.mat[i * Matrix.dim2 + i] = (T) 1;
        }

        if(swappedRowNumber != -1) {
            permutation.mat[swappedRowNumber * Matrix.dim2 + swappedRowNumber] = (T) 0;
            permutation.mat[swappedRowNumber * Matrix.dim2] = (T) 1;
            permutation.mat[0] = (T) 0;
            permutation.mat[swappedRowNumber] = (T) 1;
        }
    }

    template<typename T>
    matrix<T> inverse(matrix<T> Matrix){
        if(Matrix.dim1 != Matrix.dim2){
            printf("The given matrix is not square.");
            assert(0);
        }

        return calculus::identityMatrix(Matrix.dim1) / Matrix;

    }

    template<typename T>
    T determinant(matrix<T> Matrix){
        T temporaryProduct = (T) 1;
        matrix<T> l(Matrix.dim1, Matrix.dim1), u(Matrix.dim1, Matrix.dim1), p(Matrix.dim1, Matrix.dim1);

        calculus::LUPDecomposition(Matrix, l, u, p);

        for(int i = 0; i < Matrix.dim1; i++){

            temporaryProduct = temporaryProduct * u.mat[i * Matrix.dim1 + i];

        }

        return temporaryProduct;

    }





    ////other stuff

    long int factorial(long int x){
        long int product = 1;

        for(int i = 1; i < x + 1; i++){

            product *= i;

        }

        return product;

    }

}



template<typename T>
matrix<T> sin(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::sin(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> cos(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::cos(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> tan(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::tan(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> csc(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = 1 / std::sin(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> sec(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = 1 / std::cos(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> cot(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = 1 / std::tan(Matrix.mat[i]);

    }

    return Matrix;

}



template<typename T>
matrix<T> sinh(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::sinh(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> cosh(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::cosh(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> tanh(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::tanh(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> csch(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = 1 / std::sinh(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> sech(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = 1 / std::cosh(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> coth(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = 1 / std::tanh(Matrix.mat[i]);

    }

    return Matrix;

}



template<typename T>
matrix<T> asin(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::asin(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> acos(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::acos(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> atan(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::atan(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> acsc(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::asin(1 / Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> asec(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::acos(1 / Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> acot(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::atan(1 / Matrix.mat[i]);

    }

    return Matrix;

}



template<typename T>
matrix<T> asinh(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::asinh(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> acosh(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::acosh(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> atanh(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::atanh(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> acsch(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::asinh(1 / Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> asech(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::acosh(1 / Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> acoth(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::atanh(1 / Matrix.mat[i]);

    }

    return Matrix;

}



template<typename T>
matrix<T> sqrt(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::sqrt(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> log(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::log(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> log10(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::log10(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> cbrt(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::cbrt(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> exp(matrix<T> Matrix) {

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::exp(Matrix.mat[i]);

    }

    return Matrix;

}

template<typename T>
matrix<T> pow(matrix<T> Matrix, double val){

    for (int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {

        Matrix.mat[i] = std::pow(Matrix.mat[i], val);

    }

    return Matrix;
}



//Multiply the matrix on the right by a constant on the left (int, float, double).
template<typename T1>
matrix<T1> operator*(double val, matrix<T1> Matrix){

    for(int i = 0; i < Matrix.dim1 * Matrix.dim2; i++){
        Matrix.mat[i] = Matrix.mat[i] * val;
    }

    return Matrix;

}

//Apply the function pointer 'fun' to each entry in the matrix on the right.
template<typename T>
matrix<T> operator&&(T (*fun)(T), matrix<T> Matrix){

    for(int i = 0; i < Matrix.dim1 * Matrix.dim2; i++) {
        Matrix.mat[i] = fun(Matrix.mat[i]);
    }

    return Matrix;

}

//Apply the function pointer 'fun' to each row of the matrix on the right.
template<typename T1, typename T2>
matrix<T2> operator&(matrix<T2> (*fun)(matrix<T1>), matrix<T1> Matrix){
    matrix<T2> temporaryMatrix(0, fun(Matrix[0]).dim2);

    for(int i = 0; i < Matrix.dim1; i++){

        temporaryMatrix = temporaryMatrix <= fun(Matrix[i]);

    }

    return temporaryMatrix;

}
