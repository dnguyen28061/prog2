#include <stdio.h>
#include <iostream> 
#include <cassert> 

struct Matrix{
    public: 
    int* firstNum; 
    int rowLength; 
    int matrixLength; 
    Matrix(int* numbers = nullptr, int blockRows = 0, int matrixRows = 0){ 
        this->firstNum = numbers;
        this->rowLength = blockRows; 
        this-> matrixLength = matrixRows; 
    }
    
    // For a matrix of 
    int* elt(int n){ 
        assert((0 <= n) && (n < rowLength * rowLength)); 
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

    void printElts(){ 
        for (int i = 0; i < rowLength * rowLength; ++i){ 
            std::cout << *elt(i) << "\n"; 
        }
    }
}; 

//creates dynamically sized array of size n, initializes it to 0, and returns a pointer to that object 
int* createNewArray(int n){ 
    int* newArray = new int[n]; 
    int* newArrayPointer = newArray; 
    for (int i = 0; i < n; ++i){ 
        *newArrayPointer = 0; 
        newArrayPointer++; 
    }
    return newArray; 
};
// adds two matrices and stores them in resMatrix. Allocates new matrix if it does not exist. 
Matrix addMatrix(Matrix m1, Matrix m2, bool isAddition = true, Matrix resultMatrix = Matrix(nullptr, 0, 0)){ 
    if (resultMatrix.firstNum == nullptr){ 
        int* resultArray = createNewArray(m1.matrixLength * m1.matrixLength); 
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
Matrix multiplyStrassen(Matrix m1, Matrix m2){
    // base case 
    Matrix* matrix; 
    int* resultArray = createNewArray(m1.matrixLength * m1.matrixLength); 
    if (m1.rowLength == 1 && m2.rowLength == 1){ 
        // multiply two matrices and store product in m1 
        *(resultArray) = *(m1.firstNum) * *(m2.firstNum); 
        return Matrix(resultArray, m1.rowLength, m1.matrixLength); 
    } 
    else{ 
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
        intermediates[0] = multiplyStrassen(m1.block(0), FH); 
        intermediates[1] = multiplyStrassen(AB, m2.block(3));
        intermediates[2] = multiplyStrassen(CD, m2.block(0));
        intermediates[3] = multiplyStrassen(m1.block(3), GE);
        intermediates[4] = multiplyStrassen(AD, EH);
        intermediates[5] = multiplyStrassen(BD, GH);
        intermediates[6] = multiplyStrassen(AC, EF);

        Matrix resultMatrix = Matrix(resultArray, m1.rowLength, m1.matrixLength); 

        addMatrix(intermediates[4], intermediates[3], true, resultMatrix.block(0)); 
        addMatrix(resultMatrix.block(0), intermediates[1], false, resultMatrix.block(0)); 
        addMatrix(resultMatrix.block(0), intermediates[5], true, resultMatrix.block(0));
        
        addMatrix(intermediates[0], intermediates[1], true, resultMatrix.block(1)); 
        addMatrix(intermediates[2], intermediates[3], true, resultMatrix.block(2));

        addMatrix(intermediates[4], intermediates[0], true, resultMatrix.block(3)); 
        addMatrix(resultMatrix.block(3), intermediates[2], false, resultMatrix.block(3)); 
        addMatrix(resultMatrix.block(3), intermediates[6], false, resultMatrix.block(3)); 

        Matrix matricesToDelete[17] = {FH, AB, CD, GE, AD, EH, BD, GH, AC, EF}; 
        for(int j = 10; j < 17; ++j){ 
            matricesToDelete[j] = intermediates[j-10]; 
        }
        for (int i = 0; i < 17; ++i){ 
            delete[] matricesToDelete[i].firstNum; 
        }
        return resultMatrix; 
    }
}; 

// multiplies matrices using conventional method and returns product matrix
int** multConv(Matrix m1, Matrix m2){
    int N = m1.matrixLength;
    int** resMatrix = 0;
    resMatrix = new int*[N];
    int i, j, k; 
    for (i = 0; i < N; i++){ 
        resMatrix[i] = new int[N];
        for (j = 0; j < N; j++){ 
            resMatrix[i][j] = 0;
            for (k = 0; k < N; k++){
                resMatrix[i][j] += m1.rowCol(i,k) * m2.rowCol(k,j);
            }
        }
    }
    return resMatrix;
};

int main(int argc, char* argv){ 
    if (argc != 4) {
        throw std::invalid_argument("Usage: ./strassen 0 dimension inputfile");
    }
    int dimension = atoi(argv[2]); 
    int* matrix1 = new int[dimension * dimension]; 
    int* matrix2 = new int[dimension * dimension]; 
    

    int* matrix = new int[16];
    int* matrixPtr = matrix; 
    for(int i = 0; i < 16; ++i){ 
        *matrixPtr = rand() % 10; 
        matrixPtr++; 
    } 
    int* matrix2 = new int[16];
    int* matrixPtr2 = matrix2; 
    for(int i = 0; i < 16; ++i){ 
        *matrixPtr2 = rand() % 10; 
        matrixPtr2++; 
    } 

    Matrix matrixStruct = Matrix(matrix, 4, 4);
    Matrix matrixStruct2 = Matrix(matrix2, 4, 4); 
    matrixStruct.printElts(); 
    std::cout << "\n"; 
    matrixStruct2.printElts(); 
    Matrix resMatrix = multiplyStrassen(matrixStruct, matrixStruct2); 
    std::cout << "\n";
    resMatrix.printElts(); 
    delete[] matrix; 
    delete[] matrix2;
    delete[] resMatrix.firstNum; 
}

