from numpy import argmax

def transpose(a):
    return [list(t) for t in zip(*a)]    

def column(a, j):
    return [a[i][j] for i in range(len(a))]

def multiply(a, b):
    if len(a) != len(b):
        raise ValueError('Impossible to multiply A and B.')
    return sum([a[i] * b[i] for i in range(len(a))])
    
def gauss(a, b):
    for k in range(len(a) - 1):
        # a[k][k] is the pivot
        if a[k][k] == 0:
            j = argmax([abs(a[i][k]) for i in range(k + 1, len(a))]) + k + 1
            a[k], a[j] = a[j], a[k] # swap a lines
            b[k], b[j] = b[j], b[k] # swap b lines
        
        for i in range(k + 1, len(a)):
            m = a[i][k] / a[k][k]
            
            for j in range(k + 1, len(a)):
                a[i][j] -= m * a[k][j]
            b[i] -= m * b[k]
            a[i][k] = 0.0
    return _solve(a, b)

def _solve(a, b):
    x = [0 for i in range(len(a))]
    x[-1] = b[-1] / a[-1][-1]
    
    for i in reversed(range(len(a))):
        acc = sum([a[i][j] * x[j] for j in range(i + 1, len(a))])
        x[i] = (b[i] - acc) / a[i][i]    
    return x