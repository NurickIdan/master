import numpy

class Coordinate:
    def __init__(self, line):
        line = line.split()
        self.i = int(line[0])
        self.j = int(line[1])
        self.Mij = float(line[2])

    def __str__(self):
        return "%s %s %s\n " % (self.i, self.j, self.Mij)

    def __repr__(self):
        return self.__str__()

def matrix2contracted(matrix, badrows, badcols):
    ##This gets rid of rows and columns whose value is zero. But watch out, it doesn't copy the original matrix (for speed) but can muck it up as a result.

    rows, cols = matrix.shape
    contracted = matrix[:rows - len(badrows), :cols - len(badcols)]
    goodrows = [i for i in range(rows) if i not in badrows]
    goodcols = [i for i in range(cols) if i not in badcols]
    for i in range(len(goodrows)):
        for j in range(len(goodcols)):
            contracted[i, j] = matrix[goodrows[i], goodcols[j]]

    return contracted

def vector2expanded(vector,badrows):
##This restores rows and columns removed by matrix2contracted into the eigenvector.

        from numpy import zeros
        newvector=zeros(len(vector)+len(badrows),dtype=float)
        goodrows=[i for i in range(len(vector)+len(badrows)) if i not in badrows]
        for i in range(len(goodrows)):
                newvector[goodrows[i]]=vector[i]

        return newvector


def matrix2expanded(matrix,badrows):
##This restores rows and columns removed by matrix2contracted into the eigenvector.

        from numpy import zeros
        dim = len(matrix)+len(badrows)
        newmatrix=zeros((dim,dim), dtype=float)
        goodrows=[i for i in range(len(matrix)+len(badrows)) if i not in badrows]
        for i in range(len(goodrows)):
            for j in range(len(goodrows)):
                newmatrix[goodrows[i],goodrows[j]]=matrix[i,j]
        return newmatrix

#################################################
def matrix2zeroindex(matrix):
    rows, cols = matrix.shape
    zerorows = []
    zerocols = []
    for i in range(rows):
        if cols == len([1 for j in range(cols) if matrix[i, j] == 0]):
            zerorows.append(i)
    for j in range(cols):
        if rows == len([1 for i in range(rows) if matrix[i, j] == 0]):
            zerocols.append(j)

    return zerorows, zerocols

"""
the function gets the obs\exp normed Hi-C matrix and
returns the PCA results and the correlation matrix
eigs is returned with bad rows that contained only 0
 """
def PCA(normmatrix):
    from numpy import corrcoef
    from scipy import cov
    from scipy import linalg
    from numpy import matrix
    from copy import deepcopy

    badrows, badcols = matrix2zeroindex(normmatrix)
    newnormmatrix = matrix2contracted(normmatrix, badrows, badcols)
    corrmatrix = corrcoef(newnormmatrix)
    means = numpy.matrix([numpy.ones(len(corrmatrix[row]))*(float(sum(corrmatrix[row, :])) / corrmatrix.shape[0])\
                          for row in range(len(corrmatrix))])
    copymatrix = deepcopy(corrmatrix) - means
    print "calculating eigvecs"
    covmatrix = cov(copymatrix)
    evals, eigs = linalg.eig(covmatrix)
    evreal = [eval.real for eval in evals]
    eigs = [vector2expanded(eigs[:,i],badrows) for i in range(10)]
    copymatrix = matrix2expanded(copymatrix,badrows)
    return evreal, eigs, copymatrix

def ABFromPCA(corrmatrix, eigens):
    from numpy import dot
    result = [numpy.zeros(len(eigen)) - ((dot(corrmatrix,eigen))>0) for eigen in eigens]
    return result

def compareAB(real, test):
    minlen = min(len(real),len(test))
    re = numpy.zeros(minlen) + (real[:minlen] > 0)
    te = numpy.zeros(minlen) - (test[:minlen] > 0)
    plt.bar(range(len(re)), list(re))
    plt.bar(range(len(te)), list(te), edgecolor="blue")
    plt.show()
