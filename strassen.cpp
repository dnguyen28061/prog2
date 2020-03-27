#include <stdio.h>
#include<string.h>
#include <iostream> 
#include <cassert> 
#include <fstream> 
#include <string>
#include<chrono> 
#include<math.h> 
#include <cmath>

int* createNewArray(int);
struct Matrix{
    public: 
    int* firstNum; 
    int rowLength; 
    int matrixLength;  
    // point to this dummy element when referring to an entry in padding 
    int* zeroForPadding; 
    // the start of padding in the x and y direction of the matrix  
    int startOfHorizontalPadding; 
    int startOfVerticalPadding; 
    // cached pointers to the ABCD submatrices of this larger matrix. 
    Matrix* blocks[4] = {nullptr, nullptr, nullptr, nullptr};

    Matrix(int* numbers = nullptr, int blockRows = 0, int matrixRows = 0){ 
        this->firstNum = numbers;
        this->rowLength = blockRows; 
        this->startOfVerticalPadding = blockRows; 
        this->startOfHorizontalPadding = blockRows; 
        this-> matrixLength = matrixRows; 
        this-> zeroForPadding = new int; 
        *zeroForPadding = 0; 
    };
    ~Matrix(){
        for(int i = 0; i < 4; ++i){ 
            if (blocks[i] != nullptr){
                delete blocks[i]; 
            }
        }
        delete zeroForPadding; 
        // delete this; 
    };
    
    // For a matrix of 
    private: 
    int* elt(int n){ 
        int offset = (n / rowLength * matrixLength) + (n % rowLength);
        // make sure element is in bounds of array. 
        assert(firstNum + offset < firstNum + (matrixLength * matrixLength));
        return (firstNum + offset); 
    };
    public: 
    int* rowCol(int i, int j){
        // assert((0 <= i) && (0 <= j) && (i < rowLength) && (j < rowLength));
        // if i, j are in the padding or all of the matrix is in the padding, then return 0. 
        if (firstNum == zeroForPadding)
        {
            return zeroForPadding; 
        }
        else if (i >= startOfVerticalPadding || j >= startOfHorizontalPadding){ 
            return zeroForPadding;
        }
        else{ 
            return elt((rowLength * i) + j);
        }
    };
    private: 
    int getPadding(int firstRow, int StartOfPadding){ 
        // if padding is in the elements of the matrix. 
        if (firstRow + (rowLength / 2) >= StartOfPadding){ 
            return StartOfPadding - firstRow; 
        } 
        // if no elements are in the padding. 
        else{ 
            return rowLength / 2; 
        }
    }; 
    public: 
    Matrix* block(int n){
        assert((0 <= n) &&  (n < 4)); 
        assert(rowLength % 2 == 0);
        if (blocks[n] != nullptr){ 
            return blocks[n]; 
        }
        int* firstNumber; 
        int startOfHorizPadding; 
        int startOfVertPadding; 

        if (n == 0){ 
            firstNumber = this->rowCol(0, 0); 
            startOfHorizPadding = getPadding(0, startOfHorizontalPadding); 
            startOfVertPadding = getPadding(0, startOfVerticalPadding); 
        }
        else if (n == 1){ 
            firstNumber = this->rowCol(0, rowLength / 2); 
            startOfHorizPadding = getPadding(rowLength/2, startOfHorizontalPadding); 
            startOfVertPadding = getPadding(0, startOfVerticalPadding); 
        }
        else if (n == 2){ 
            firstNumber = this->rowCol(rowLength / 2, 0); 
            startOfHorizPadding = getPadding(0, startOfHorizontalPadding); 
            startOfVertPadding = getPadding(rowLength/2, startOfVerticalPadding); 
        }
        // take the first num of second block and add rowLength / 2 
        else if (n == 3){ 
            firstNumber = this->rowCol(rowLength / 2, rowLength / 2);
            startOfHorizPadding = getPadding(rowLength/2, startOfHorizontalPadding); 
            startOfVertPadding = getPadding(rowLength/2, startOfVerticalPadding);  
        }
        blocks[n] = new Matrix(firstNumber, rowLength / 2, matrixLength);
        blocks[n]->startOfHorizontalPadding = startOfHorizPadding; 
        blocks[n]->startOfVerticalPadding = startOfVertPadding; 

        if(firstNumber == zeroForPadding){ 
            blocks[n]->firstNum = blocks[n]->zeroForPadding; 
        }

        return blocks[n];
    }

    void printElts(){ 
        for (int i = 0; i < rowLength; ++i){ 
            for (int j = 0; j < rowLength; ++j){
                std::cout << *rowCol(i, j) << " ";
                if (j == rowLength - 1){
                    std::cout << "\n";
                } 
            } 
        }
    }
    void pad(){ 
        assert(rowLength % 2 == 1); 
        rowLength += 1; 
    }
    void unpad(){ 
        assert(rowLength % 2 == 0); 
        rowLength -= 1; 
    }

    Matrix* copy(){
        int* copyArray = createNewArray(rowLength);
        memcpy(copyArray, firstNum, rowLength * rowLength);
        return new Matrix(copyArray, rowLength, matrixLength);
    }
}; 

//creates dynamically sized array of size n, initializes it to 0, and returns a pointer to that object 

int* createNewArray(int dimension){ 
    int* newArray = new int[dimension * dimension]; 
    int* newArrayPointer = newArray; 

    for (int i = 0; i < dimension * dimension; ++i){ 
        *newArrayPointer = 0; 
        newArrayPointer++; 
    }
    return newArray; 
};
// adds two matrices and stores them in resMatrix. Allocates new matrix if it does not exist. 
Matrix* addMatrix(Matrix* m1, Matrix* m2, bool isAddition = true, Matrix* resultMatrix = NULL){ 
    if (resultMatrix == NULL){ 
        int* resultArray = createNewArray(m1->matrixLength); 
        resultMatrix = new Matrix(resultArray, m1->rowLength, m1->rowLength); 
    }
    for (int i = 0; i < m1->rowLength; ++i){ 
        for(int j = 0; j < m1->rowLength; ++j){ 
            int sum;
            if (isAddition){ 
                sum = *(m1->rowCol(i, j)) + *(m2->rowCol(i, j)); 
            }
            else{ 
                sum = *(m1->rowCol(i, j)) - *(m2->rowCol(i, j)); 
            }
            *resultMatrix->rowCol(i, j) = sum; 
        }
    }
    return resultMatrix;
};

// multiplies two matrices and stores the result in a new matrix 
int* multConv(Matrix*, Matrix*); 
Matrix* multiplyStrassen(Matrix* m1, Matrix* m2, int n0){
    assert(m1->rowLength = m2->rowLength); 
    // base case 
    int* resultArray; 
    if (m1->rowLength <= n0){ 
        // multiply two matrices and store product in m1 
        resultArray = multConv(m1, m2); 
        return new Matrix(resultArray, m1->rowLength, m1->rowLength); 
    } 
    else{ 
        resultArray = createNewArray(m1->matrixLength); 
        Matrix* resultMatrix = new Matrix(resultArray, m1->rowLength, m1->matrixLength); 
        bool padded = false; 
        if (m1->rowLength % 2 == 1){ 
            m1->pad(); 
            m2->pad(); 
            resultMatrix->pad(); 
        }
        Matrix* FH = addMatrix(m2->block(1), m2->block(3), false); 
        Matrix* AB = addMatrix(m1->block(0), m1->block(1)); 
        Matrix* CD = addMatrix(m1->block(2), m1->block(3)); 
        Matrix* GE = addMatrix(m2->block(2), m2->block(0), false); 
        Matrix* AD = addMatrix(m1->block(0), m1->block(3)); 
        Matrix* EH = addMatrix(m2->block(0), m2->block(3)); 
        Matrix* BD = addMatrix(m1->block(1), m1->block(3), false); 
        Matrix* GH = addMatrix(m2->block(2), m2->block(3)); 
        Matrix* AC = addMatrix(m1->block(0), m1->block(2), false); 
        Matrix* EF = addMatrix(m2->block(0), m2->block(1)); 
        Matrix* intermediates[7]; 
        intermediates[0] = multiplyStrassen(m1->block(0), FH, n0); 
        intermediates[1] = multiplyStrassen(AB, m2->block(3), n0);
        intermediates[2] = multiplyStrassen(CD, m2->block(0), n0);
        intermediates[3] = multiplyStrassen(m1->block(3), GE, n0);
        intermediates[4] = multiplyStrassen(AD, EH, n0);
        intermediates[5] = multiplyStrassen(BD, GH, n0);
        intermediates[6] = multiplyStrassen(AC, EF, n0);

        addMatrix(intermediates[4], intermediates[3], true, resultMatrix->block(0)); 
        addMatrix(resultMatrix->block(0), intermediates[1], false, resultMatrix->block(0)); 
        addMatrix(resultMatrix->block(0), intermediates[5], true, resultMatrix->block(0));
        
        addMatrix(intermediates[0], intermediates[1], true, resultMatrix->block(1)); 
        addMatrix(intermediates[2], intermediates[3], true, resultMatrix->block(2));

        addMatrix(intermediates[4], intermediates[0], true, resultMatrix->block(3)); 
        addMatrix(resultMatrix->block(3), intermediates[2], false, resultMatrix->block(3)); 
        addMatrix(resultMatrix->block(3), intermediates[6], false, resultMatrix->block(3)); 
        if (padded){ 
            m1->unpad(); 
            m2->unpad(); 
            resultMatrix->unpad(); 
        }
        Matrix* matricesToDelete[17] = {FH, AB, CD, GE, AD, EH, BD, GH, AC, EF}; 
        for(int j = 10; j < 17; ++j){ 
            matricesToDelete[j] = intermediates[j-10]; 
        }
        for (int i = 0; i < 17; ++i){ 
            delete[] matricesToDelete[i]->firstNum; 
            delete matricesToDelete[i]; 
        }
        return resultMatrix; 
    }
};


int* returnOneDArray(int** arr, int dimension){ 
    int* flattenedArray = createNewArray(dimension); 
    int currElt = 0; 
    for (int i = 0; i < dimension; ++i){ 
        for (int j = 0; j < dimension; ++j){ 
            flattenedArray[currElt] = arr[i][j]; 
            currElt += 1; 
        }
    }
    return flattenedArray; 

};
// multiplies matrices using conventional method and returns product matrix
int* multConv(Matrix* m1, Matrix* m2){
    int N = m1->rowLength;
    int** resMatrix = 0;
    resMatrix = new int*[N];
    for (int i = 0; i < N; i++){
        resMatrix[i] = new int[N];
    }
    for (int i = 0; i < N; ++i){ 
        for(int j = 0; j < N; ++j){ 
            resMatrix[i][j] = 0; 
        }
    }
    int i, j, k; 
    for (i = 0; i < N; i++){ 
        for (k = 0; k < N; k++){ 
            for (j = 0; j < N; j++){
                resMatrix[i][j] += *(m1->rowCol(i,k)) * *(m2->rowCol(k,j));
            }
        }
    }
    int* oneDimMatrix = returnOneDArray(resMatrix, N); 
    for (int i = 0; i < N; ++i){ 
        delete[] resMatrix[i]; 
    }
    delete[] resMatrix; 
    return oneDimMatrix;
}

// Creates a random graph of edges with probability p
Matrix* createRandGraph(double p, int dim){
    int* coords = new int[dim * dim];
    int* coordsPtr = coords;
    Matrix* randGraph = new Matrix(coordsPtr, dim, dim);
    int edgeProb;
    for (int i = 0; i < dim * dim; ++i){ 
        coords[i] = -1; 
    }
    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            if (*randGraph->rowCol(j,i) != -1){
                *randGraph->rowCol(i,j) = *randGraph->rowCol(j,i);
            }
            else{
                edgeProb = (rand() % 100) + 1;
                if (edgeProb <= p * 100.){
                    *randGraph->rowCol(i,j) = 1;
                }
                else{
                    *randGraph->rowCol(i,j) = 0;
                }
            }
        }
    }
    return randGraph;
}

int numOfTriangles(double p, int dim){
    Matrix* A = createRandGraph(p, dim);
    // A->printElts();
    Matrix* Acopy = A->copy();
    Matrix* A2 = multiplyStrassen(A, Acopy, 15);
    // A2->printElts();
    Matrix* A3 = multiplyStrassen(A2, A, 15);
    // A3->printElts();
    int numTriangles = 0;
    for (int j = 0; j < dim; j++){
        numTriangles += *A3->rowCol(j,j);
    }
    numTriangles /= 6;
    std::cout << "Number of triangles for p = " << p << ": " << numTriangles << "\n";
    return numTriangles;
}
 
// Returns factorial of n 
int fact(int n) 
{ 
    int res = 1; 
    for (int i = 2; i <= n; i++){
        res = res * i;
    }
    return res; 
}

void readFileIntoArray(std::ifstream* file, int dimension, Matrix* matrix){ 
    for (int i = 0; i < dimension; ++i){ 
        for (int j = 0; j < dimension; ++j){ 
            std::string numb; 
            std::getline(*file, numb); 
            *matrix->rowCol(i, j) = atoi(numb.c_str()); 
        }
    }
}

int main(int argc, char** argv){ 
    if (argc != 4) {
        throw std::invalid_argument("Usage: ./strassen 0 dimension inputfile");
    }
    int dimension = atoi(argv[2]); 
    int* matrix_1 = createNewArray(dimension); 
    int* matrix_2 = createNewArray(dimension); 
    Matrix* matrixStruct = new Matrix(matrix_1, dimension, dimension);
    Matrix* matrixStruct2 = new Matrix(matrix_2, dimension, dimension); 
    std::ifstream file(argv[3], std::ios::in);
    readFileIntoArray(&file, dimension, matrixStruct); 
    readFileIntoArray(&file, dimension, matrixStruct2);
    file.close(); 
    auto start = std::chrono::high_resolution_clock::now(); 
    Matrix* resMatrix = multiplyStrassen(matrixStruct, matrixStruct2, 15); 
    auto end = std::chrono::high_resolution_clock::now(); 
    int* convRes = multConv(matrixStruct, matrixStruct2);
    auto endConventional = std::chrono::high_resolution_clock::now();
    std::cout << "Time for Strassen: " << (std::chrono::duration_cast<std::chrono::microseconds>(end - start)).count() << "\n"; 
    std::cout << "Time for Conventional: " << (std::chrono::duration_cast<std::chrono::microseconds>(endConventional - end)).count() << "\n"; 
    // resMatrix->printElts(); 
    // int i = dimension;
    // int prod = 1;
    // while (i >= dimension - 3){
    //     prod *= i;
    //     i--;
    // }
    // int comb = prod / 6;
    // // Compute number of triangles for each prob p
    // // for (int i = 1; i < 6; i++){
    // const double p = 1 / 100.;
    // int numTri = numOfTriangles(p, dimension);
    // int exp = comb * pow(p,3);
    // std::cout << "Expected # of Triangles for p = " << p << ": " << exp << "\n";
    // }
    delete[] matrix_1;
    delete[] matrix_2; 
    delete matrixStruct; 
    delete matrixStruct2;
    delete[] resMatrix->firstNum; 
    delete resMatrix; 
    delete[] convRes;
 

}
