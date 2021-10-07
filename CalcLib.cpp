//Special definitions will go here later on that make coding more easy with inspiration from MatLab and Wolfram.

#include <complex>
#include <iostream>
#include <vector>
#include <functional>
#include <cassert>

/* ////DISCLAIMER////
 * THIS LIBRARY WAS BUILT FOR CALCULUS RELATED TASKS.
 * BECAUSE OF THE NUMERICAL METHODS EMPLOYED, IT IS SUGGESTED TO PAY ATTENTION TO YOUR DATA.
 * THE METHODS EMPLOYED ARE APPROXIMATIONS AND MAY CAUSE SEVERE DETERIORATION TO DATA IF USED IMPROPERLY.
 *
 *
 *
 * void array(T (*fun)(T), const int size, T arr[], const T delta = 0.1, const T start = 0)
 *      This function generates an array using the function `fun` evaluated at evenly spaced intervals from `start` until `size` * `delta`.
 *
 *
 * void function(T1 (*fun)(T1), const int size, T2 arr[])
 *      This function applies the function head (*fun) to the first `size` elements of the array `arr`.
 *
 *
 * void linspace(T1 arr[], const int size, T2 start, T3 end)
 *      This function fills the array `arr` with evenly spaced points starting from `start` to `end`.
 *
 *
 * void difference(T arr[], const int size) **OR** void difference(T arr[], const int size, T arr2[])
 *      This function takes in an array `arr` and finds the difference between the ith element and the (i - 1)th element.
 *      The second function differs in that the results will be placed into the array `arr2` instead of `arr1`.
 *
 *
 * T dotProduct(T arr1[], const int size, T arr2[])
 *      This function takes in the two arrays `arr1` and `arr2` and finds the dot product up to `size` elements in the array.
 *
 *
 * T crossProduct(T arr1[], const int size, T arr2[])
 *      This function does ***NOT*** behave the same as the normal cross product in that it does *not* return a vector,
 *      this instead returns the AREA of the PARALLELOGRAM between the two vectors `arr1` and `arr2`.
 *
 *
 * T magnitude(T arr[], const int size)
 *      This function takes in the array `arr` and returns the magnitude of the vector up to `size` elements.
 *
 *
 * void printMatrix(T arr[], const int dim1, const int dim2 = 1)
 *      This function takes in the array `arr` and prints it out as a matrix of dimensions `dim1`x`dim2`.
 *
 *
 * void transpose(T arr[], const int dim1, const int dim2) **OR** transpose(T arr1[], const int dim1, const int dim2, T arr2[])
 *      This function takes in the array `arr` and rearranges the elements so that the matrix of dimensions `dim1`x`dim2` is
 *      now a transposed matrix of dimensions `dim2`x`dim1`.
 *      The second function differs in that instead of overwriting `arr` it is transferred to `arr2`.
 *
 *
 * void multiply(T arr1[], const int dim1, const int dim2, const int dim3, T arr2[]) **OR** multiply(T arr1[], const int dim1, const int dim2, const int dim3, T arr2[], T arr3[])
 *      This function multiplies the two array matrices `arr1` and `arr2` and stores the newly multiplied matrix to `arr1`.
 *      The second function differs in that it stores the resulting matrix to `arr3`.
 *
 *
 * T index(T arr1[], const int dim1, const int dim2, const int dim3, T arr2[], const int index1, const int index2)
 *      This function takes in two array matrices `arr1` and `arr2` and return the index (index1, index2) of the resulting matrix.
 */

//This namespace will only contain commonly used constants, and will be updated as I find more. Constants will be Pascal case.
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

    matrix<T> row(std::vector<T> arr){
        matrix<T> temp(1, arr.size());

        for(int i = 0; i < arr.size(); i++){
            temp.mat[i] = arr[i];
        }

        return temp;

    }

    matrix<T> col(std::vector<T> arr){
        matrix<T> temp(arr.size(), 1);

        for(int i = 0; i < arr.size(); i++){
            temp.mat[i] = arr[i];
        }

        return temp;

    }

    explicit matrix(int x = 1, int y = 1): dim1(x), dim2(y), mat(x * y) {}

    matrix<T> transpose(matrix<T> arr){
        T temp[arr.dim1 * arr.dim2];

        for(int i = 0; i < arr.dim2; i++){
            for(int k = 0; k < arr.dim1; k++){

                temp[k + i * arr.dim1] = arr.mat[i + k * arr.dim2];

            }
        }



        for(int i = 0; i < arr.dim1 * arr.dim2; i++){
            arr.mat[i] = temp[i];
        }

        arr.resize(arr.dim2, arr.dim1);

        return arr;

    }

    void resize(int newDim1, int newDim2) {
        this->dim1 = newDim1;
        this->dim2 = newDim2;
        this->mat.resize(this->dim1 * this->dim2);
    }

    matrix<T> operator*(const matrix<T> other){
        assert((dim2 == other.dim1) || ((other.dim1 == 1) && (other.dim2 == 1)));
        matrix<T> temp(dim1, other.dim2);
        T sum = 0;

        for (int i = 0; i < dim1; i++) {
            for (int j = 0; j < other.dim2; j++) {
                for (int k = 0; k < dim2; k++) {

                    sum += mat[k + i * dim2] * other.mat[j + k * other.dim2];

                }

                temp.mat[i * other.dim2 + j] = sum;
                sum = 0;

            }

        }

        return temp;
    }
    template<typename T2>
    matrix<T> operator*(const T2 val){

        for(int i = 0; i < dim1 * dim2; i++){
                mat[i] *= val;
            }

        return *this;
    }

    matrix<T> operator+(const matrix<T> other){
        assert((dim1 == other.dim1) && (dim2 == other.dim2));

        matrix<T> temp(dim1, dim2);

        for(int i = 0; i < dim1 * dim2; i++){
            temp.mat[i] = mat[i] + other.mat[i];
        }

        return temp;

    }

    matrix<T> operator=(const matrix<T> other){
        resize(other.dim1, other.dim2);

        for(int i = 0; i < other.dim1 * other.dim2; i++){
            mat[i] = other.mat[i];
        }

        return other;

    }

    matrix<T> operator-(const matrix<T> other){
        assert((dim1 == other.dim1) && (dim2 == other.dim2));

        matrix<T> temp(dim1, dim2);

        for(int i = 0; i < dim1 * dim2; i++){
            temp.mat[i] = mat[i] - other.mat[i];
        }

        return temp;

    }

    matrix<T> operator-(){
        matrix<T> temp(dim1, dim2);

        for(int i = 0; i < dim1 * dim2; i++){
            temp.mat[i] = -mat[i];
        }

        return temp;

    }

    matrix<T> operator*=(matrix<T> other){

        *this = *this * other;

        return *this;

    }

    matrix<T> operator+=(matrix<T> other){

        *this = *this + other;

        return *this;

    }

    matrix<T> operator-=(matrix<T> other){

        *this = *this - other;

        return *this;

    }

    matrix<T> operator~(){

        return transpose(*this);

    }

    matrix<T> operator|(matrix<T> other){
        assert((dim1 == other.dim1) && (dim2 == other.dim2));

        for(int i = 0; i < dim1 * dim2; i++) {
            mat[i] /= other.mat[i];
        }

        return *this;

    }

    matrix<T> operator^(matrix<T> other){
        assert(((dim1 == other.dim1) && (other.dim2 == 1)) || ((other.dim1 == 1) && (dim2 == other.dim2)));

        if((other.dim1 == 1) && (dim2 == other.dim2)){

            for(int i = 0; i < dim1; i++){
                for(int k = 0; k < dim2; k++){

                    mat[i * dim2 + k] = std::pow(mat[i * dim2 + k], other.mat[k]);

                }
            }

            return *this;

        }
        else{

            for(int i = 0; i < dim1; i++){
                for(int k = 0; k < dim2; k++){

                    mat[i * dim2 + k] = std::pow(mat[i * dim2 + k], other.mat[i]);

                }
            }

                return *this;
        }


    }

    matrix<T> operator[](int val){
        matrix<T> temp(1, dim2);

        for (int i = 0; i < dim2; i++) {
            temp.mat[i] = mat[i + dim2 * val];
        }

        return temp;

    }

    T operator()(int val){

        return mat[val];

    }

};

namespace calculus{

    ////array stuff

    ////basic assets

    template<typename T>
    void array(T (*fun)(T), const int size, T arr[], const T delta = 0.1, const T start = 0){

        for(int i = 0; i < size; i++){
            arr[i] = fun(start + i * delta);
        }

    }

    template<typename T1, typename T2>
    void function(T1 (*fun)(T1), const int size, T2 arr[]){

        for(int i = 0; i < size; i++){
            arr[i] = fun(arr[i]);
        }

    }

    template<typename T1, typename T2, typename T3>
    void linspace(T1 arr[], const int size, T2 start, T3 end){
        T1 spacing = (end - start) / (size - 1);

        for(int i = 0; i < size; i++){
            arr[i] = start + i * spacing;
        }

    }

    template<typename T>
    void difference(T arr[], const int size){

        for(int i = 0; i < size - 1; i++){
            arr[i] = arr[i + 1] - arr[i];
        }

    }

    template<typename T>
    void difference(T arr[], const int size, T arr2[]){

        for(int i = 0; i < size - 1; i++){
            arr2[i] = arr[i + 1] - arr[i];
        }

    }

    template<typename T>
    T dotProduct(T arr1[], const int size, T arr2[]){
        T sum = 0;

        for(int i = 0; i < size; i++){
            sum += arr1[i] * arr2[i];
        }

        return sum;

    }

    template<typename T>
    T area(T arr1[], const int size, T arr2[]){
        T sum = 0;
        T val;

        for(int i = 0; i < size; i++){
            for(int k = i + 1; k < size; k++){
                val = arr1[i] * arr2[k] - arr1[k] * arr2[i];
                sum += val * val;
            }
        }

        return sqrtl(sum);

    }

    template<typename T>
    T magnitude(T arr[], const int size){
        T sum = 0;

        for(int i = 0; i < size; i++){
            sum += arr[i] * arr[i];
        }

        return sqrtl(sum);

    }


//    template<typename T>
//    void rootsToCoefficients(T arr[], const int size){
//        T temporaryArray[size + 1];
//        temporaryArray[size] = 1;
//
//        for(int i = 0; i < size; i++){
//            for(int k = size; k > size - i; k++){
//
//                temporaryArray[k] =
//
//            }
//        }
//
//    }

    ////matrix related functions

    template<typename T>
    void printMatrixf(T arr[], const int dim1, const int dim2 = 1){

        for(int i = 0; i < dim1; i++){
            for(int k = 0; k < dim2; k++){
                if(arr[k + i * dim2] >= 0) {
                    std::cout << std::fixed << arr[k + i * dim2] << "  ";
                }
                else{
                    std::cout << std::fixed << arr[k + i * dim2] << " ";
                }
            }
            printf("\n");
        }

    }

    template<typename T>
    void printMatrixs(T arr[], const int dim1, const int dim2 = 1){

            for (int i = 0; i < dim1; i++) {
                for (int k = 0; k < dim2; k++) {
                    if (arr[k + i * dim2] >= 0) {
                        std::cout << std::scientific << arr[k + i * dim2] << "  ";
                    } else {
                        std::cout << std::scientific << arr[k + i * dim2] << " ";
                    }
                }
                printf("\n");
            }



    }

    template<typename T>
    void transpose(T arr[], const int dim1, const int dim2){
        T temp[dim1 * dim2];

        for(int i = 0; i < dim2; i++){
            for(int k = 0; k < dim1 + 1; k++){

                temp[k + i * dim1] = arr[i + k * dim2];

            }
        }

        for(int i = 0; i < dim1 * dim2; i++){
            arr[i] = temp[i];
        }

    }

    template<typename T>
    void transpose(T arr[], const int dim1, const int dim2, T arr2[]){

        for(int i = 0; i < dim2; i++){
            for(int k = 0; k < dim1 + 1; k++){

                arr2[k + i * dim1] = arr[i + k * dim2];

            }
        }

    }

    template<typename T>
    void multiply(T arr1[], const int dim1, const int dim2, const int dim3, T arr2[]){
        T temp[dim1 * dim3];
        T sum = 0;

        for(int i = 0; i < dim1; i++) {
            for (int j = 0; j < dim3; j++) {
                for (int k = 0; k < dim2; k++) {

                    sum += arr1[k + i * dim2] * arr2[j + k * dim3];

                }

                temp[i * dim3 + j] = sum;
                sum = 0;

            }

        }

        for(int i = 0; i < dim1 * dim3; i++){
            arr1[i] = temp[i];
        }

    }

    template<typename T>
    void multiply(T arr1[], const int dim1, const int dim2, const int dim3, T arr2[], T arr3[]) {
        T sum = 0;

        for (int i = 0; i < dim1; i++) {
            for (int j = 0; j < dim3; j++) {
                for (int k = 0; k < dim2; k++) {

                    sum += arr1[k + i * dim2] * arr2[j + k * dim3];

                }

                arr3[i * dim3 + j] = sum;
                sum = 0;

            }
        }
    }

////Perfect if you only need a few indices from the resulting matrix of `arr1` x `arr2`. Saves time and always works.
    template <typename T>
    T index(T arr1[], const int dim1, const int dim2, const int dim3, T arr2[], const int index1, const int index2){
        T sum = 0;

        for(int i = 0; i < dim1; i++){

            printf("%lf, %lf\n", arr1[index1 * dim2 + i], arr2[index2 * dim3 + i * dim3]);
            sum += arr1[index1 * dim2 + i] * arr2[index2 * dim3 + i * dim3];

        }

        return sum;

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

            printf("\n");

        }

        return temp[0];

    }



        ////array based derivatives
////IT IS ADVISED TO USE ARRAY BASED DERIVATIVES SPARINGLY TO PREVENT DATA CORRUPTION
////I'll eventually add more of these, but that doesn't mean you should actually use them.
    template<typename T1, typename T2>
    void derivative3_1(const T1 arr[], int size, T2 arr2[], const T1 delta = 0.1){
        const T2 coe = 1 / delta;

        for(int i = 1; i < size - 4; i++){

            arr2[i] = coe * (-T1 (1) / T1 (3) * arr[i - 1] - T1 (1) / T1 (2) * arr[i] + arr[i + 1] - T1 (1) / T1 (6) * arr[i + 2]);

        }

        arr2[0] = coe * (-1.833333333333333 * arr[0] + 3 * arr[1] - 1.5 * arr[2] + T1 (1) / T1 (3) * arr[3]);

        arr2[size - 4] = coe * (0.1666666666666667 * arr[size - 6] - arr[size - 5] + 0.5 * arr[size - 4] + T1 (1) / T1 (3) * arr[size - 3]);
        arr2[size - 3] = coe * (0.1666666666666667 * arr[size - 5] - arr[size - 4] + 0.5 * arr[size - 3] + T1 (1) / T1 (3) * arr[size - 2]);
        arr2[size - 2] = coe * (0.1666666666666667 * arr[size - 4] - arr[size - 3] + 0.5 * arr[size - 2] + T1 (1) / T1 (3) * arr[size - 1]);
        arr2[size - 1] = coe * (-T1 (1) / T1 (3) * arr[size - 4] + T1 (3) / T1 (2) * arr[size - 3] - 3 * arr[size - 2] + 1.8333333333333333 * arr[size - 1]);

    }

    template<typename T>
    void derivative3_1(T arr[], int size, const T delta = 0.1){
        const T coe = 1 / delta;
        T temp = arr[0];
        T atemp;

        T temp1 = arr[1];
        T temp2 = arr[2];
        T temp3 = arr[3];
        T temp4 = arr[size - 6];
        T temp5 = arr[size - 5];
        T temp6 = arr[size - 1];


        for(int i = 1; i < size - 4; i++){
            atemp = arr[i];

            arr[i] = coe * (-0.3333333333333333 * temp - 0.5 * arr[i] + arr[i + 1] - 0.1666666666666667 * arr[i + 2]);

            temp = atemp;
        }

        arr[0] = coe * (-1.833333333333333 * arr[0] + 3 * temp1 - 1.5 * temp2 + 0.3333333333333333 * temp3);

        temp1 = temp4;
        temp2 = temp5;
        temp3 = arr[size - 4];
        temp4 = arr[size - 3];
        temp5 = arr[size - 2];

        arr[size - 4] = coe * (0.1666666666666667 * temp1 - temp2 + 0.5 * temp3 + 0.3333333333333333 * temp4);
        arr[size - 3] = coe * (0.1666666666666667 * temp2 - temp3 + 0.5 * temp4 + 0.3333333333333333 * temp5);
        arr[size - 2] = coe * (0.1666666666666667 * temp3 - temp4 + 0.5 * temp5 + 0.3333333333333333 * temp6);
        arr[size - 1] = coe * (-0.3333333333333333 * temp3 + 1.5 * temp4 - 3 * temp5 + 1.8333333333333333 * arr[size - 1]);

    }



        ////function to array based differentiation

    template<typename T1, typename T2, typename T3, typename T4>
    void derivative2_1(T1 (*fun)(T1), const int size, T2 arr[], const T3 delta = 0.1, const T4 start = 0) {
        const T2 coe = 1 / delta;
        const T2 temp1 = fun(start - 0.5 * delta);
        const T2 temp2 = fun(start + 0.5 * delta);

        for(int i = 0; i < size; i++){
            arr[i] = coe * (temp1 - temp2);

            temp1 = temp2;
            temp2 = fun(start + (i + 1.5) * delta);
        }

    }

    template<typename T>
    void derivative2_2(T (*fun)(T), const int size, T arr[], const T delta = 0.1, const T start = 0) {
        const T coe = 1 / delta;
        const T coe2 = coe * coe;
        T temp1 = fun(start - delta);
        T temp2 = fun(start);
        T temp3 = fun(start + delta);

        for(int i = 0; i < size; i++){
            arr[i] = coe2 * (temp3 - 2 * temp2 + temp1);

            temp1 = temp2;
            temp2 = temp3;
            temp3 = fun(start + (i + 2) * delta);
        }

    }



    template<typename T>
    void derivative3_1(T (*fun)(T), const int size, T arr[], const T delta = 0.1, const T start = 0){
        const T coe = 1 / delta;
        const T c1 = T(1) / T(24);
        const T c2 = T(9) / T(8);

        T temp1 = fun(start - 1.5 * delta);
        T temp2 = fun(start - 0.5 * delta);
        T temp3 = fun(start + 0.5 * delta);
        T temp4 = fun(start + 1.5 * delta);

        for(int i = 0; i < size; i++){
            arr[i] = coe * (c1 * temp1 - c2 * temp2 + c2 * temp3 - c1 * temp4);

            temp1 = temp2;
            temp2 = temp3;
            temp3 = temp4;
            temp4 = fun(start + (i + 2.5) * delta);
        }

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
    void ode1_1(void (*fun)(T1 [3], T1), T1 initialConditions[3], T2 end, double delta = 0.1){
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
    void ode1_2(void (*fun)(T1 [4], T1), T1 initialConditions[4], T2 end, double delta = 0.1){
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
    void ode2_2(void (*fun)(T1 [4], T1), T1 initialConditions[4], T2 end, double delta = 0.1){
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
    void ode1_3(void (*fun)(T1 [5], T1), T1 initialConditions[5], T2 end, double delta = 0.1){
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
    void ode2_3(void (*fun)(T1 [5], T1), T1 initialConditions[5], T2 end, double delta = 0.1){
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
    void ode3_3(void (*fun)(T1 [5], T1), T1 initialConditions[5], T2 end, double delta = 0.1){
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

    template<typename T1, typename T2, typename T3>
    void linspace(std::vector<T1> &arr, T2 start, T3 end){
        T1 spacing = (end - start) / (arr.size() - 1);

        for(int i = 0; i < arr.size(); i++){
            arr[i] = start + i * spacing;
        }

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
    T dotProduct(std::vector<T> arr1, std::vector<T> arr2){
        assert(arr1.size() == arr2.size());
        T sum = 0;

        for(int i = 0; i < arr1.size(); i++){
            sum += arr1[i] * arr2[i];
        }

        return sum;

    }

    template<typename T>
    T area(std::vector<T> arr1, std::vector<T> arr2){
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
    T magnitude(std::vector<T> arr){
        T sum = 0;

        for(int i = 0; i < arr.size(); i++){
            sum += arr[i] * arr[i];
        }

        return sqrtl(sum);

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
    void printMatrixf(matrix<T> arr){

        for(int i = 0; i < arr.dim1; i++){
            for(int k = 0; k < arr.dim2; k++){
                if(arr.mat[k + i * arr.dim2] >= 0) {
                    std::cout << std::fixed << arr.mat[k + i * arr.dim2] << "  ";
                }
                else{
                    std::cout << std::fixed << arr.mat[k + i * arr.dim2] << " ";
                }
            }
            printf("\n");
        }

    }

    template<typename T>
    void printMatrixs(matrix<T> arr){

        for (int i = 0; i < arr.dim1; i++) {
            for (int k = 0; k < arr.dim2; k++) {
                if (arr.mat[k + i * arr.dim2] >= 0) {
                    std::cout << std::scientific << arr.mat[k + i * arr.dim2] << "  ";
                } else {
                    std::cout << std::scientific << arr.mat[k + i * arr.dim2] << " ";
                }
            }
            printf("\n");
        }

    }

    template <typename T>
    T index(matrix<T> arr1, matrix<T> arr2, const int index1, const int index2){
        T sum = 0;

        for(int i = 0; i < arr1.dim1; i++){

            sum += arr1[index1 * arr1.dim2 + i] * arr2[index2 * arr2.dim2 + i * arr2.dim2];

        }

        return sum;

    }

    template<typename T1, typename T2>
    matrix<T2> function(T1 (*fun)(T1), matrix<T2> arr){

        for(int i = 0; i < arr.dim1 * arr.dim2; i++){
            arr.mat[i] = fun(arr.mat[i]);
        }

        return arr;
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

template<typename T1, typename T2>
    matrix<T1> operator*(T2 val, matrix<T1> arr){

    for(int i = 0; i < arr.dim1 * arr.dim2; i++){
        arr.mat[i] *= val;
    }

    return arr;

}

    template<typename T1, typename T2>
    matrix<T1> operator&(T2 (*fun)(T2), matrix<T1> other){

        for(int i = 0; i < other.dim1 * other.dim2; i++) {
            other.mat[i] = fun(other.mat[i]);
        }

        return other;

}


    //This namespace is for lazier versions of functions that the computer guesses `delta` for a good, working answer.
    //These will adjust `delta` automatically for a decent answer or just break entirely due to machine precision limitations. Use these wisely.
namespace lazy{

    template<typename T1, typename T2>
    T1 derivative2_1(T1 (*fun)(T1), const T2 point, const int precision = 2) {
        T1 delta = 0.1;

        T1 val1;
        T1 val2 = (fun(point + 0.5 * delta) - fun(point - 0.5 * delta)) / delta;


        do{

            val2 = val1;

            delta *= 0.1;

            val1 = (fun(point + 0.5 * delta) - fun(point - 0.5 * delta)) / delta;

        }while(std::abs(val1 - val2) > std::pow(10, -precision));

        return val1;

    }

    //THESE ARE UNFINISHED, THEY ARE NOT LAZY ENOUGH
    template<typename T1, typename T2>
    T1 derivative2_2(T1 (*fun)(T1), const T2 point, const int precision = 2) {
        T1 delta = 1.0;
        T1 val = fun(point);

        T1 val1 = (fun(point + delta) - 2 * val + fun(point - delta)) / (delta * delta);
        T1 val2;

        printf("%lf\n", val1);


        do{

            val2 = val1;

            delta *= T1 (std::sqrt(0.1));

            val1 = (fun(point + delta) - 2 * val + fun(point - delta)) / (delta * delta);

        }while(std::abs(val1 - val2) > std::pow(10, -precision));

        return val1;

    }

    template<typename T1, typename T2>
    T1 derivative2_3(T1 (*fun)(T1), const T2 point, const int precision = 2) {
        T1 delta = 0.1;

        T1 val1 = (fun(point + 1.5 * delta) - 3 * fun(point + 0.5 * delta) + 3 * fun(point - 0.5 * delta) - fun(point - 1.5 * delta)) / (delta * delta * delta);
        T1 val2;

        do{

            val2 = val1;

            delta *= T1 (std::pow(0.10, 1.0 / 3.0));

            val1 = (fun(point + 1.5 * delta) - 3 * fun(point + 0.5 * delta) + 3 * fun(point - 0.5 * delta) - fun(point - 1.5 * delta)) / (delta * delta * delta);

        }while(std::abs(val1 - val2) > std::pow(10, -precision));

        return val1;

    }

//    template<typename T1, typename T2>
//    T1 derivative2_4(T1 (*fun)(T1), const T2 point, const T1 delta = 0.1) {
//
//        return (fun(point + 2 * delta) - 4 * fun(point + delta) + 6 * fun(point) - 4 * fun(point - delta) + fun(point - 2 * delta)) / (delta * delta * delta * delta);
//
//    }
//
//
//
//    template<typename T1, typename T2>
//    T1 derivative3_1(T1 (*fun)(T1), const T2 point, const T1 delta = 0.1) {
//
//        //I modified the proof to Simpson's rule so that it may be extended to derivatives.
//        //These coefficients will produce the approximate derivative of a 3rd degree interpolation.
//        //Fun fact, summing up the coefficients c1, c2, c3, etc. on the `return` line for any derivativeM_N will always add to 0.
//
//        const T1 c1 = T1 (1) / T1 (24);
//        const T1 c2 = T1 (9) / T1 (8);
//
//        return (-c1 * fun(point + 1.5 * delta) + c2 * fun(point + 0.5 * delta) - c2 * fun(point - 0.5 * delta) + c1 * fun(point - 1.5 * delta)) / delta;
//    }
//
//    template<typename T1, typename T2>
//    T1 derivative3_2(T1 (*fun)(T1), const T2 point, const T1 delta = 0.1) {
//
//        //These coefficients are the same as nesting derivative3_1 two times.
//
//        const T1 c1 = T1 (1) / T1(576);
//        const T1 c2 = T1 (3) / T1 (32);
//        const T1 c3 = T1 (87) / T1 (64);
//        const T1 c4 = T1 (365) / T1 (144);
//
//        return (c1 * fun(point + 3 * delta) - c2 * fun(point + 2 * delta) + c3 * fun(point + delta) - c4 * fun(point) + c3 * fun(point - delta) - c2 * fun(point - 2 * delta) + c1 * fun(point - 3 * delta)) / (delta * delta);
//    }
//
//    template<typename T1, typename T2>
//    T1 derivative3_3(T1 (*fun)(T1), const T2 point, const T1 delta = 0.1) {
//
//        //These coefficients are the same as nesting der1 three times.
//
//        const T1 c1 = T1 (1) / T1 (13824);
//        const T1 c2 = T1 (3) / T1 (512);
//        const T1 c3 = T1 (21) / T1 (128);
//        const T1 c4 = T1 (2005) / T1 (1152);
//        const T1 c5 = T1 (1137) / T1 (256);
//
//        return (-c1 * fun(point + 4.5 * delta) + c2 * fun(point + 3.5 * delta) - c3 * fun(point + 2.5 * delta) + c4 * fun(point + 1.5 * delta) - c5 * fun(point + 0.5 * delta) + c5 * fun(point - 0.5 * delta) - c4 * fun(point - 1.5 * delta) + c3 * fun(point - 2.5 * delta) - c2 * fun(point - 3.5 * delta) + c1 * fun(point - 4.5 * delta)) / (delta * delta * delta);
//    }
//
//    template<typename T1, typename T2>
//    T1 derivative3_4(T1 (*fun)(T1), const T2 point, const T1 delta = 0.1) {
//
//        //These coefficients are the same as nesting der1 four times.
//
//        const T1 c1 = T1 (1) / T1 (331776);
//        const T1 c2 = T1 (1) / T1 (3072);
//        const T1 c3 = T1 (83) / T1 (6144);
//        const T1 c4 = T1 (21871) / T1 (82944);
//        const T1 c5 = T1 (9535) / T1 (4096);
//        const T1 c6 = T1 (3659) / T1 (512);
//        const T1 c7 = T1 (280301) / T1 (27648);
//
//        return (c1 * fun(point + 6 * delta) - c2 * fun(point + 5 * delta) + c3 * fun(point + 4 * delta) - c4 * fun(point + 3 * delta) + c5 * fun(point + 2 * delta) - c6 * fun(point + delta) + c7 * fun(point) - c6 * fun(point - delta) + c5 * fun(point - 2 * delta) - c4 * fun(point - 3 * delta) + c3 * fun(point - 4 * delta) - c2 * fun(point - 5 * delta) + c1 * fun(point - 6 * delta)) / (delta * delta * delta * delta);
//    }

}
