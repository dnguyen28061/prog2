#include <stdio.h>
#include <iostream> 
#include <cassert> 
#include <fstream> 
#include <string>
#include<chrono> 

struct Matrix{
    public: 
    int* firstNum; 
    int rowLength; 
    int matrixLength; 
    int* numsInPadding; 
    bool padded; 
    Matrix(int* numbers = nullptr, int blockRows = 0, int matrixRows = 0){ 
        this->firstNum = numbers;
        this->rowLength = blockRows; 
        this-> matrixLength = matrixRows; 
    }
    
    // For a matrix of 
    int* elt(int n){ 
        // assert((0 <= n) && (n < rowLength * rowLength)); 
        int offset = (n / rowLength * matrixLength) + (n % rowLength);
        return (firstNum + offset); 
    }
    Matrix block(int n){
        assert((0 <= n) &&  (n < 4)); 
        int* firstNumber; 
        if (n == 0){ 
            firstNumber = this->elt(0); 
        }
        else if (n == 1){ 
            firstNumber = this->elt(rowLength / 2); 
        }
        else if (n == 2){ 
            firstNumber = this->elt(rowLength * rowLength / 2); 
        }
        // take the first num of second block and add rowLength / 2 
        else if (n == 3){ 
            firstNumber = this->elt(rowLength * rowLength / 2) + rowLength / 2; 
        }
        return Matrix(firstNumber, rowLength / 2, matrixLength); 
    }

    // returns element at ith row and jth column
    int* rowCol(int i, int j){
        // assert((0 <= i) && (0 <= j) && (i < rowLength) && (j < rowLength));
        return elt((rowLength * i) + j);
    }

    void pad(){ 
        assert(rowLength % 2 == 1); 
        rowLength += 1; 
        padded = true; 
        numsInPadding = new int[rowLength + rowLength]; 
        int currElt = 0; 
        for (int i = 0; i < rowLength; ++i){ 
            std::cout << i << " " << rowLength << "\n"; 
            std::cout << *(rowCol(i, rowLength)) << "\n"; 
            numsInPadding[currElt] = *(rowCol(i, rowLength)); 
            currElt++; 
            *(rowCol(i, rowLength)) = 0; 
            numsInPadding[currElt] = *(rowCol(rowLength, i)); 
            currElt++; 
            *(rowCol(rowLength, i)) = 0; 
        }
    }

    void unpad(){ 
        // std::cout << "rowlength: " << rowLength << "\n";
        padded = false; 
        assert(rowLength % 2 == 0); 
        int currElt = 0; 
        for (int i = 0; i < rowLength; ++i){ 
            *(rowCol(i, rowLength)) = numsInPadding[currElt]; 
            currElt++; 
            *(rowCol(rowLength, i)) = numsInPadding[currElt]; 
            currElt++; 
        }
        delete[] numsInPadding; 
        rowLength -= 1; 
    }
    void printElts(){ 
        for (int i = 0; i < rowLength * rowLength; ++i){ 
            std::cout << *elt(i) << " ";
            if (i % rowLength == rowLength - 1){
                std::cout << "\n";
            } 
        }
    }
}; 

//creates dynamically sized array of size n, initializes it to 0, and returns a pointer to that object 

int* createNewArray(int dimension){ 
    if (dimension % 2 == 1){ 
        dimension++; 
    }
    int* newArray = new int[dimension * dimension]; 
    int* newArrayPointer = newArray; 

    for (int i = 0; i < dimension * dimension; ++i){ 
        *newArrayPointer = 0; 
        newArrayPointer++; 
    }
    return newArray; 
};
// adds two matrices and stores them in resMatrix. Allocates new matrix if it does not exist. 
Matrix addMatrix(Matrix m1, Matrix m2, bool isAddition = true, Matrix resultMatrix = Matrix(nullptr, 0, 0)){ 
    if (resultMatrix.firstNum == nullptr){ 
        int* resultArray = createNewArray(m1.matrixLength); 
        resultMatrix = Matrix(resultArray, m1.rowLength, m1.rowLength); 
    }
    for (int i = 0; i < m1.rowLength * m1.rowLength; ++i){ 
        int sum; 
        if (isAddition){ 
            sum = *(m1.elt(i)) + *(m2.elt(i)); 
        }
        else{ 
            sum = *(m1.elt(i)) - *(m2.elt(i)); 
        }
        *resultMatrix.elt(i) = sum; 
    }
    return resultMatrix;
};

// multiplies two matrices and stores the result in a new matrix 
int* multConv(Matrix m1, Matrix m2); 
Matrix multiplyStrassen(Matrix m1, Matrix m2, int n0){
    assert(m1.rowLength = m2.rowLength); 
    // base case 
    Matrix* matrix; 
    int* resultArray; 
    if (m1.rowLength <= n0){ 
        // multiply two matrices and store product in m1 
        resultArray = multConv(m1, m2); 
        return Matrix(resultArray, m1.rowLength, m1.rowLength); 
    } 
    else{ 
        if(m1.rowLength % 2 == 1){ 
            m1.pad(); 
            m2.pad(); 
            m1.printElts(); 
            m2.printElts(); 
            std::cout << "padded"; 
        }
        resultArray = createNewArray(m1.matrixLength); 
        Matrix FH = addMatrix(m2.block(1), m2.block(3), false); 
        Matrix AB = addMatrix(m1.block(0), m1.block(1)); 
        Matrix CD = addMatrix(m1.block(2), m1.block(3)); 
        Matrix GE = addMatrix(m2.block(2), m2.block(0), false); 
        Matrix AD = addMatrix(m1.block(0), m1.block(3)); 
        Matrix EH = addMatrix(m2.block(0), m2.block(3)); 
        Matrix BD = addMatrix(m1.block(1), m1.block(3), false); 
        Matrix GH = addMatrix(m2.block(2), m2.block(3)); 
        Matrix AC = addMatrix(m1.block(0), m1.block(2), false); 
        Matrix EF = addMatrix(m2.block(0), m2.block(1)); 
        Matrix intermediates[7]; 
        intermediates[0] = multiplyStrassen(m1.block(0), FH, n0); 
        intermediates[1] = multiplyStrassen(AB, m2.block(3), n0);
        intermediates[2] = multiplyStrassen(CD, m2.block(0), n0);
        intermediates[3] = multiplyStrassen(m1.block(3), GE, n0);
        intermediates[4] = multiplyStrassen(AD, EH, n0);
        intermediates[5] = multiplyStrassen(BD, GH, n0);
        intermediates[6] = multiplyStrassen(AC, EF, n0);

        Matrix resultMatrix = Matrix(resultArray, m1.rowLength, m1.matrixLength); 

        addMatrix(intermediates[4], intermediates[3], true, resultMatrix.block(0)); 
        addMatrix(resultMatrix.block(0), intermediates[1], false, resultMatrix.block(0)); 
        addMatrix(resultMatrix.block(0), intermediates[5], true, resultMatrix.block(0));
        
        addMatrix(intermediates[0], intermediates[1], true, resultMatrix.block(1)); 
        addMatrix(intermediates[2], intermediates[3], true, resultMatrix.block(2));

        addMatrix(intermediates[4], intermediates[0], true, resultMatrix.block(3)); 
        addMatrix(resultMatrix.block(3), intermediates[2], false, resultMatrix.block(3)); 
        addMatrix(resultMatrix.block(3), intermediates[6], false, resultMatrix.block(3)); 

        // 
        Matrix matricesToDelete[17] = {FH, AB, CD, GE, AD, EH, BD, GH, AC, EF}; 
        for(int j = 10; j < 17; ++j){ 
            matricesToDelete[j] = intermediates[j-10]; 
        }
        for (int i = 0; i < 17; ++i){ 
            delete[] matricesToDelete[i].firstNum; 
        }
        if (m1.padded == true && m2.padded == true){ 
            m1.unpad();
            m2.unpad(); 
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

}
// multiplies matrices using conventional method and returns product matrix
int* multConv(Matrix m1, Matrix m2){
    int N = m1.rowLength;
    int** resMatrix = 0;
    resMatrix = new int*[N];
    int i, j, k; 
    for (i = 0; i < N; i++){ 
        resMatrix[i] = new int[N];
        for (j = 0; j < N; j++){ 
            resMatrix[i][j] = 0;
            for (k = 0; k < N; k++){
                resMatrix[i][j] += *(m1.rowCol(i,k)) * *(m2.rowCol(k,j));
            }
        }
    }
    int* oneDimMatrix = returnOneDArray(resMatrix, N); 
    // for (int i = 0; i < N * N; ++i){ 
    //     std::cout << oneDimMatrix[i] << "\n"; 
    // }
    for (int i = 0; i < N; ++i){ 
        delete[] resMatrix[i]; 
    }
    delete[] resMatrix; 
    return oneDimMatrix;
};

// Creates a random graph of edges with probability p
// Matrix createRandGraph(int p){
//     int* coords = new int[1024];
//     int* coordsPtr = coords;
//     for (int i = 0; i < 1024; i++){
//         *coordsPtr = 
//     }
// }

void readFileIntoArray(std::ifstream* file, int dimension, Matrix matrix, bool isOdd){ 
    for (int i = 0; i < dimension; ++i){ 
        for (int j = 0; j < dimension; ++j){ 
            if (isOdd && (i == dimension - 1 || j == dimension - 1)){
                *matrix.rowCol(i, j) = 0; 
            }
            else{
                std::string numb; 
                std::getline(*file, numb); 
                *matrix.rowCol(i, j) = atoi(numb.c_str()); 
            }
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
    bool isOddDimension = false; 
    if(dimension % 2 == 1){ 
        dimension += 1; 
        isOddDimension = true; 
    }
    Matrix matrixStruct = Matrix(matrix_1, dimension, dimension);
    Matrix matrixStruct2 = Matrix(matrix_2, dimension, dimension); 

    std::ifstream file(argv[3], std::ios::in);
    readFileIntoArray(&file, dimension, matrixStruct, isOddDimension);
    readFileIntoArray(&file, dimension, matrixStruct2, isOddDimension);

    auto start = std::chrono::high_resolution_clock::now(); 
    Matrix resMatrix = multiplyStrassen(matrixStruct, matrixStruct2, 1); 
    auto end = std::chrono::high_resolution_clock::now(); 
    int* convRes = multConv(matrixStruct, matrixStruct2);
    auto endConventional = std::chrono::high_resolution_clock::now(); 
    std::cout << "Time for Strassen: " << (std::chrono::duration_cast<std::chrono::microseconds>(end - start)).count() << "\n"; 
    std::cout << "Time for Conventional: " << (std::chrono::duration_cast<std::chrono::microseconds>(endConventional - end)).count() << "\n"; 
    resMatrix.printElts();
    delete[] matrix_1; 
    delete[] matrix_2;
    delete[] resMatrix.firstNum; 
    delete[] convRes;
    file.close();  
}
