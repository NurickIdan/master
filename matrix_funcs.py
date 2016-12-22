class Coordinate:
    def __init__(self,line):
        line = line.split()
        self.i = int(line[0])
        self.j = int(line[1])
        self.Mij = float(line[2])
        
    def __str__(self):
        return "%s %s %s\n"%(self.i,self.j,self.Mij)

    def __repr__(self):
        return self.__str__()


def matrixFromFile(filename):
    import numpy
    data = file(filename,"rb").readlines()
    result = numpy.zeros((len(data)-2,len(data)-2))
    for i in range(2,len(data)):
        vals = map(float,data[i].split()[1:])
        for j in range(len(vals)):
            result[i-2,j]=vals[j]
    return result

def matrixFromSparseFile(filename, res = 10**5,exp = None):
    data = map(Coordinate, file(filename,"rb").readlines())
    dim = max(max([x.i for x in data]), max([x.j for x in data]))
    dim = dim/res + 1
    mat = numpy.zeros((dim,dim))
    for d in data:
        if d.Mij > 0 or d.Mij < 0:
            mat[d.i/res, d.j/res] = d.Mij
        if exp:
            mat[d.i/res, d.j/res] /= exp[abs(d.i-d.j)/res]
    return mat

"""res is the resolution of the HIC (we mostly use 5 KB resolution)"""
def normalizeIntraMatrix(cells, norm,exp, res = 5000):
    for k in range(len(cells)):
        cells[k].Mij = cells[k].Mij / \
                       (norm[cells[k].i/res] * norm[cells[k].j/res])
        distance = abs(cells[k].i - cells[k].j)
        cells[k].Mij = cells[k].Mij / exp[distance/res]
    return cells

def normalizeInterMatrix(cells, norm1, norm2, res = 5000):
    for k in range(len(cells)):
        cells[k].Mij = cells[k].Mij / \
                       (norm1[cells[k].i/res] * norm2[cells[k].j/res])
    return cells

def normIntraFile(filename,res = 1000, expected = True):
    
    f = file(filename + "RAWobserved","rb").read().splitlines()
    print "reading file, contains %d coordinate"%(len(f))
    result = file(filename + "normalized","wb")
    if not expected:
        result = file(filename + "normalized_wo_expected","wb")
    print "create result file"
    norm = map(float,file(filename + "KRnorm","rb").read().splitlines())
    print "reading norm file", len(norm)
    exp = map(float,file(filename + "KRexpected","rb").read().splitlines())
    print "reading exp file", len(exp)
    for k in range(len(f)):
        sp = f[k].split()
        if len(sp)<3:
            break
        Mij = float(sp[-1])
        i = int(sp[0])
        j = int(sp[1])
        Mij = Mij / \
                       (norm[i/res] * norm[j/res])
        distance = abs(i-j)
        # checking for NaN values
        try:
            if expected:
                Mij = Mij / exp[distance/res]
        except:
                Mij = Mij / exp[distance / res -1]
        result.write("%s %s %s\n"%(i,j,round(Mij,2)))        
        if k % 10000000 == 0:
            print k
    result.close()
    del f
    del norm
    del exp
    del result
    return 0

def normInterFile(filename, chr1, chr2, outputfile = r"E:\Users\idannuri\interchromosomal\chr",res = 100000):
    
    f = file(filename + "%s_%s_100kb.RAWobserved"%(chr1,chr2),"rb").read().splitlines()
    print "reading file, contains %d coordinate"%(len(f))
    result = file(outputfile + "%s_%s_100kb.normalized"%(chr1,chr2),"wb")
    print "create result file"    
    norm1 = map(float,file(filename + "%s_100kb.KRnorm"%(chr1),"rb").read().splitlines())
    norm2 = map(float,file(filename + "%s_100kb.KRnorm"%(chr2),"rb").read().splitlines())
    for k in range(len(f)):
        sp = f[k].split()
        Mij = float(sp[-1])
        i = int(sp[0])
        j = int(sp[1])
        Mij = Mij / \
                       (norm1[i/res] * norm2[j/res])
        result.write("%s %s %s\n"%(i,j,Mij)) 
        if k % 10000000 == 0:
            print k
    result.close()
    del f
    del norm1
    del norm2
    del result
    return 0
    
ranges = [1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000]

def create_diffs(chri):
    diffs = []
    for i in range(len(chri)):
            c = chr1[i].split()
            try:
                    val = int(math.log(float(c[2]))*10)
                    d = abs(int(c[0])-int(c[1]))
                    diffs.append((d,val))
            except:
                    pass
            if i%1000000 == 0:
                    print i,
    return diffs

def create_results():
    for i in range(len(diff)):
	if diff[i][0] > ranges[-1]:
		continue
	for j in range(len(ranges)):
		if diff[i][0] <= ranges[j]:
			for k in range(j,len(ranges)):
				results[k].append(diff[i][1])
			break
	if i % 1000000 == 0:
		print i


