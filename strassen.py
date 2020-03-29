import random

def addMatrix(m1, m2):
    return 0

def multiplyStrassen(m1, m2, n0):
    return 0

def conventionalMult(m1, m2):
    rows, cols = (len(m1), len(m1[0]))
    # initializes product matrix
    prod = [[0 for i in range(cols)] for j in range(rows)]
    # finds product entries
    for i in range(rows):
        for k in range(cols):
            for j in range(cols):
                prod[i][j] += m1[i][k] * m2[k][j]
    # prints product matrix
    for i in range(len(prod)):  
        for j in range(len(prod)):  
            print(prod[i][j],end=" ") 
        print("\n",end="")
    return prod
    # weird zip method from piazza
    # [[sum(ab for a,b in zip(X_row,Y_col))for Y_col in zip(Y)] for X_row in X]

def createRandomGraph(p, dim):
    randGraph = [[-1 for i in range(dim)] for j in range(dim)]
    for i in range(dim):
        for j in range(dim):
            if randGraph[i][j] != -1:
                continue
            else:
                edgeProb = random.randint(1, 101)
                if edgeProb <= p * 100:
                    randGraph[i][j] = 1
                    randGraph[j][i] = 1
                else:
                    randGraph[i][j] = 0
                    randGraph[j][i] = 0
    return randGraph

def numOfTriangles(p, dim):
    A = createRandomGraph(p, dim)
    A2 = multiplyStrassen(A, A, 15)
    A3 = multiplyStrassen(A2, A, 15)
    numTriangles = 0
    for i in range(dim):
        numTriangles += A3[i][i]
    return numTriangles

def factorial(n):
    res = 1
    for i in range(2, n + 1):
        res *= i
    return res