#include <complex>
#include <iostream>
#include <cstring>
//#include <vector>

/* THIS LIBRARY WAS BUILT FOR CALCULUS RELATED TASKS.
 * BECAUSE OF THE NUMERICAL METHODS EMPLOYED, IT IS SUGGESTED TO PAY ATTENTION TO YOU DATA.
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
 *
 */

//This namespace will only contain commonly used constants, and will be updated as I find more. Constants will be Pascal case.
namespace constant{

    const long double Pi = 3.1415926535897932385;

    const long double E = 2.7182818284590452354;
}

namespace calculus{
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
    T crossProduct(T arr1[], const int size, T arr2[]){
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


    long int factorial(long int x){
        long int product = 1;

        for(int i = 1; i < x + 1; i++){
            product *= i;
        }

        return product;

    }

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



    ////Differentiation functions

        ////function based derivatives

    template<typename T>
    T derivative2_1(T (*fun)(T), const T point, const T delta = 0.1) {

        return (fun(point + 0.5 * delta) - fun(point - 0.5 * delta))/ delta;

    }

    template<typename T>
    T derivative2_2(T (*fun)(T), const T point, const T delta = 0.1) {

        return (fun(point + delta) - 2 * fun(point) + fun(point - delta)) / (delta * delta);

    }

    template<typename T>
    T derivative2_3(T (*fun)(T), const T point, const T delta = 0.1) {

        return (fun(point + 1.5 * delta) - 3 * fun(point + 0.5 * delta) + 3 * fun(point - 0.5 * delta) - fun(point - 1.5 * delta)) / (delta * delta * delta);

    }

    template<typename T>
    T derivative2_4(T (*fun)(T), const T point, const T delta = 0.1) {

        return (fun(point + 2 * delta) - 4 * fun(point + delta) + 6 * fun(point) - 4 * fun(point - delta) + fun(point - 2 * delta)) / (delta * delta * delta * delta);

    }

    template<typename T>
    T derivative3_1(T (*fun)(T), const T point, const T delta = 0.1) {

        //I modified the proof to Simpson's rule so that it may be extended to derivatives.
        //These coefficients will produce the approximate derivative of a 3rd degree interpolation.

        return (-0.04166666667 * fun(point + 1.5 * delta) + 1.125 * fun(point + 0.5 * delta) - 1.125 * fun(point - 0.5 * delta) + 0.04166666667 * fun(point - 1.5 * delta)) / delta;
    }

    template<typename T>
    T derivative3_2(T (*fun)(T), const T point, const T delta = 0.1) {

        //These coefficients are the same as nesting derivative3_1 two times.

        return (0.001736111111 * fun(point + 3 * delta) - 0.09375 * fun(point + 2 * delta) + 1.359375 * fun(point + delta) - 2.534722222 * fun(point) - 1.359375 * fun(point - delta) - 0.09375 * fun(point - 2 * delta) + 0.001736111111 * fun(point - 3 * delta)) / (delta * delta);
    }

    template<typename T>
    T derivative3_3(T (*fun)(T), const T point, const T delta = 0.1) {

        //These coefficients are the same as nesting der1 three times.

        return (-0.00007233796296 * fun(point + 4.5 * delta) + 0.005859375 * fun(point + 3.5 * delta) - 0.1640625 * fun(point + 2.5 * delta) + 1.740451389 * fun(point + 1.5 * delta) - 4.44140625 * fun(point + 0.5 * delta) + 4.44140625 * fun(point - 0.5 * delta) - 1.740451389 * fun(point - 1.5 * delta) + 0.1640625 * fun(point - 2.5 * delta) - 0.005859375 * fun(point - 3.5 * delta) + 0.00007233796296 * fun(point - 4.5 * delta)) / (delta * delta * delta);
    }



        ////array based derivatives
////IT IS ADVISED TO USE ARRAY BASED DERIVATIVES SPARINGLY TO PREVENT DATA CORRUPTION
    template<typename T1, typename T2, typename T3>
    void derivative3_1(const T1 arr[], int size, T2 arr2[], const T3 delta = 0.1){

        const T2 coe = 1 / delta;

        for(int i = 1; i < size - 4; i++){

            arr2[i] = coe * (-0.3333333333333333 * arr[i - 1] - 0.5 * arr[i] + arr[i + 1] - 0.1666666666666667 * arr[i + 2]);

        }

        arr2[0] = coe * (-1.833333333333333 * arr[0] + 3 * arr[1] - 1.5 * arr[2] + 0.3333333333333333 * arr[3]);

        arr2[size - 4] = coe * (0.1666666666666667 * arr[size - 6] - arr[size - 5] + 0.5 * arr[size - 4] + 0.3333333333333333 * arr[size - 3]);
        arr2[size - 3] = coe * (0.1666666666666667 * arr[size - 5] - arr[size - 4] + 0.5 * arr[size - 3] + 0.3333333333333333 * arr[size - 2]);
        arr2[size - 2] = coe * (0.1666666666666667 * arr[size - 4] - arr[size - 3] + 0.5 * arr[size - 2] + 0.3333333333333333 * arr[size - 1]);
        arr2[size - 1] = coe * (-0.3333333333333333 * arr[size - 4] + 1.5 * arr[size - 3] - 3 * arr[size - 2] + 1.8333333333333333 * arr[size - 1]);

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

            ////spanned differentiation

    template<typename T>
    void derivative2_1(T (*fun)(T), T arr[], const int size, const T delta = 0.1, const T start = 0, const T span = 0.1){
        const T coe = 1 / delta;
        const T space = span + delta;
        const T val1 = 0.5 * delta;

        for(int i = 0; i < size; i++){
            arr[i] = coe * (fun(start + val1 + i * space) - fun(start - val1 + i * space));
        }

    }

    template<typename T>
    void derivative2_2(T (*fun)(T), T arr[], const int size, const T delta = 0.1, const T start = 0, const T span = 0.1){
        const T coe = 1 / delta;
        const T space = span + delta;
        const T coe2 = coe * coe;

        for(int i = 0; i < size; i++){
            arr[i] = coe2 * (fun(start + delta + i * space) - 2 * fun(start + i * space) + fun(start - delta + i * space));
        }

    }

    template<typename T>
    void derivative2_3(T (*fun)(T), T arr[], const int size, const T delta = 0.1, const T start = 0, const T span = 0.1){
        const T coe = 1 / delta;
        const T space = span + delta;
        const T coe3 = coe * coe * coe;
        const T val1 = 0.5 * delta;
        const T val2 = 1.5 * delta;


        for(int i = 0; i < size; i++){
            arr[i] = coe3 * (fun(start + val2 + i * space) - 3 * fun(start + val1 + i * space) + 3 * fun(start - val1 + i * space) - fun(start - val2 + i * space));
        }

    }

            ////spanless differentiaton

    template<typename T>
    void derivative2_1(T (*fun)(T), const int size, T arr[], const T delta = 0.1, const T start = 0){
        const T coe = 1 / delta;
        const T temp1 = fun(start - 0.5 * delta);
        const T temp2 = fun(start + 0.5 * delta);

        for(int i = 0; i < size; i++){
            arr[i] = coe * (temp1 - temp2);

            temp1 = temp2;
            temp2 = fun(start + (i + 1.5) * delta);
        }

    }

    template<typename T>
    void derivative2_2(T (*fun)(T), const int size, T arr[], const T delta = 0.1, const T start = 0){
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
        T temp1 = fun(start - 1.5 * delta);
        T temp2 = fun(start - 0.5 * delta);
        T temp3 = fun(start + 0.5 * delta);
        T temp4 = fun(start + 1.5 * delta);

        for(int i = 0; i < size; i++){
            arr[i] = coe * (0.0416666666666667 * temp1 - 1.125 * temp2 + 1.125 * temp3 - 0.0416666666666667 * temp4);

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

        T1 old1 = 0.1666666666666667 * fun(low);

        for(int i = 0; i < ratio; i++){

            sum += old1 + 0.6666666666666667 * fun(low + (i + 0.5) * delta);

            old1 = 0.1666666666666667 * fun(low + (i + 1) * delta);

            sum += old1;

        }

        return delta * sum + diff * (0.1666666666666667 * fun(up - diff) + 0.6666666666666667 * fun(up - 0.5 * diff) + 0.1666666666666667 * fun(up));

    }

    template<typename T1, typename T2, typename T3>
    T1 integral2_2(T1 (*fun)(T1), const T2 low, const T3 up, const T1 delta = 0.1){
        T1 sum = 0;
        const T1 val = up - low;
        const int ratio = val / delta;
        const T1 diff = val - ratio * delta;

        T1 old1 = val * 0.1666666666666667 * fun(low);

        for(int i = 0; i < ratio; i ++) {

            sum += old1 + 0.6666666666666667 * (val - (i + 0.5) * delta) * fun(low + (i + 0.5) * delta);

            old1 = 0.1666666666666667 * (val - (i + 1) * delta) * fun(low + (i + 1) * delta);

            sum += old1;

        }

        return delta * sum + diff * (0.1666666666666667 * diff * fun(up - diff) + 0.3333333333333333 * diff * fun(up - 0.5 * diff));

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



    ////ODE solvers (non-DSAE)
//Coming...whenever I make it.

};
