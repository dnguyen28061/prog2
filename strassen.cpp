#include <stdio.h>
#include <iostream> 

struct Matrix{
    public: 
    int* firstNum; 
    int rowLength; 
    int matrixLength; 
    Matrix(int* numbers, int blockRows, int matrixRows){ 
        this->rowLength = blockRows; 
        this->firstNum = numbers;
        this-> matrixLength = matrixRows; 
        }
    
    // For a matrix of 
    Matrix block(int n){
        assert((0 <= n) &&  (n < 4)); 
        int* firstNumber; 
        if (n == 0){ 
            firstNumber = firstNum; 
        }
        else if (n == 1){ 
            firstNumber = firstNum + rowLength / 2; 
        }
        else if (n == 2){ 
            firstNumber = firstNum + rowLength * rowLength / 2; 
        }
        // take the first num of second block and add rowLength / 2 
        else if (n == 3){ 
            firstNumber = firstNum + (rowLength * rowLength + rowLength) / 2;  
        }
        return Matrix(firstNumber, rowLength / 2, matrixLength); 
    }
    int elt(int n){ 
        // std::cout << n << " , rowLength = " << this->rowLength << "\n"; 
        assert((0 <= n) && (n < rowLength * rowLength)); 
        int offset = (n / rowLength * matrixLength) + (n % rowLength);
        return *(firstNum + offset); 
    }
}; 

// multiplies two matrices and stores the result in a new matrix 
Matrix multiplyStrassen(Matrix m1, Matrix m2){
    // base case 
    if (m1.rowLength == 1 && m2.rowLength == 1){ 
        int numsArray[1] = {(*(m1.firstNum) * *(m2.firstNum))}; 
        return Matrix(numsArray, m1.rowLength, m1.matrixLength); 
    } 
    else{ 
        
    }

    return Matrix(dummy, 1, 1);
}; 

// adds two matrices and stores them in a 
Matrix addMatrix(Matrix matrices[], bool isAddition){ 

    int* dummy = new int[1]; 
    return Matrix(dummy, 1, 1);
};

Matrix multConv(Matrix m1, Matrix m2){
    int* dummy = new int[1];
    return Matrix(dummy, 1, 1);
}

int main(){ 
    int* matrix = new int[64];
    int* matrixPtr = matrix; 
    for(int i = 0; i < 64; ++i){ 
        *matrixPtr = i; 
        matrixPtr++; 
    } 
    Matrix matrixStruct = Matrix(matrix, 8, 8); 
    Matrix matrices[2] = {matrixStruct, matrixStruct}; 
    addMatrix(matrices, true); 
    // std::cout << *(matrixStruct.block(1).firstNum) << "\n"; 
    for (int i = 0; i < 4; ++i){ 
        std::cout << matrixStruct.block(0).block(1).elt(i) << "\n"; 
    }
    delete[] matrix; 
}

