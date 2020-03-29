# import statements we'll probably need

def addMatrix(m1, m2):
    return 0

def pad(matrix): 
    for i in range(len(matrix)): 
        matrix[i].append(0)
    matrix.append([0 for _ in range(len(m1) + 1)])
    
def unpad(matrix): 
    matrix.pop()
    for i in range(len(matrix)): 
        matrix[i].pop() 
        
def multiplyStrassen(m1, m2, n0):
    assert(len(m1) == len(m2)) 
    assert(len(m1[0]) == len(m2[0]))
    if len(m1) < n0: 
        return conventionalMult(m1, m2)
    else: 

    return 0

def conventionalMult(m1, m2):
    return 0

def createRandomGraph(p, dim):
    return 0

def numOfTriangles(p, dim):
    return 0

def factorial(n):
    return 0