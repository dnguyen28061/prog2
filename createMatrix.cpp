#include<stdio.h> 
#include<iostream> 
#include<random> 
#include<fstream> 

int main(int argc, char** argv){ 
    if (argc != 5){ 
        std::cout << "Usage: ./createMatrix lowerBound upperBound dimension inputfile" << "\n"; 
    }
    char* file = argv[4]; 
    int dimension = atoi(argv[3]); 
    int lowerBound = atoi(argv[1]); 
    int upperBound = atoi(argv[2]); 

    std::fstream fileStream; 
    fileStream.open(file, std::ios::out); 
    std::default_random_engine generator; 
    std::uniform_int_distribution<int> distribution(lowerBound, upperBound); 
    for(int i = 0; i < dimension * dimension + dimension * dimension; ++i){ 
        fileStream << distribution(generator) << "\n"; 
    }
    fileStream.close(); 
}