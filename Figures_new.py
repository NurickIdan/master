import math
import matplotlib.pyplot as plt

try:
    a = hun
except:
    print "loading Figures data"
    hun = load("./ABcompartmentsFromHiC/Final/dic_100KB_sym_KRnorm.pkl")
    pairs = load("/home/gaga/data-scratch/idannuri/ChipSeq/pipeline_input/pairs")
    expressions = load("./Data/gene_expression/rna_seq/dic_FPKM_new.pkl")
    enhs = load("./Data/Enhancer/dic_enhancers_new.pkl")
    dic_chia = load("./Data/Chiapet/dic_chia.pkl")
    dic_chia_full = load("./Data/Chiapet/dic_chia_full.pkl")
    print "Finished loading Figures data"

def generateFigure1(full = True, verbose = False):
    """
    Figure 1 compares genes expression between AB and BA for all pairs in keys
    """
    import numpy as np
    import math
    import scipy.stats
    maxi = 0
    #keys = ['HMEC', 'T47D', 'MCF7', 'IMR90', 'NHEK', 'GM12878','HUVEC','K562']
    keys = ['GM12878','T47D', 'MCF7','HUVEC','HMEC', 'IMR90', 'NHEK', 'K562']
    data = np.zeros((len(keys), len(keys)))
    for i in range(data.shape[0]):
        for j in range(i):
            if verbose:
                print;print;print keys[i],keys[j]
            ratios = expressionCompartmentAnalysis(keys[i],keys[j],ret=1, verbose = verbose)
            real_ratio = [[x[0]/x[1] for x in r] for r in ratios]
            data[i,j] = scipy.stats.ks_2samp(real_ratio[1],real_ratio[2])[1]
            data[i,j] = -1*math.log(data[i,j],10)
            if data[i,j] > maxi:
                maxi = data[i,j]
            if full:
                data[j,i] = data[i,j]
    for i in range(data.shape[0]):
        data[i,i] = 0
    return data

def Figure0_compartmentVsExpression(cls = None, ylabel = "log2 FPKM", pval = 0.05, thresh = 10):
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as GS
    import math
    if cls == None:
        cls = ['HUVEC', 'NHEK', 'K562', 'HMEC', 'IMR90', 'MCF7', 'GM12878']
    tmp = {}
    for cl in cls:
        tmp[cl] = {}
        t = expressionCompartmentAnalysis(cl, cl, AA=True)
        for x in t:
            tmp[cl][x]  = 1
    nexp =[[math.log(expressions[cl][x.ID]+1,2) for x in TSSs if x not in tmp[cl] and x.ID in expressions[cl]\
            and expressions[cl][x.ID] != None] for cl in cls]

    exp = [[math.log(expressions[cl][x.ID]+1,2) for x in tmp[cl] if x.ID in expressions[cl]]\
            for cl in cls]

    pval = [round(math.log(scipy.stats.ks_2samp(exp[i],nexp[i])[1],10),1) for i in range(len(exp))]

    fig,ax = plt.subplots(1,1)
    bp = ax.boxplot(sum_list([[[],exp[i],nexp[i],[]] for i in range(len(exp))]), patch_artist = True)
    for whisker in bp['whiskers']:
        whisker.set(color='#7570b3', linewidth=2)
    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set( linewidth=4)
    for cap in bp['caps']:
        cap.set(color='#7570b3', linewidth=2)
    for flier in bp['fliers']:
        flier.set(marker='', color='#e7298a', alpha=0.5)
    
    ax.set_xticklabels(sum_list([["",'%s\n A\nlog p-val = %s'%(cl,pval[cls.index(cl)]), '\n B',""]  for cl in cls]), fontsize = 24)
    ax.set_ylabel(ylabel, fontsize = 34)
    ax.set_yticklabels(range(0,14,2),fontsize = 24)
    ax.set_xlabel("Cell lines", fontsize = 34)
    ax.set_ylim(-0.2,12)
    #ax.set_title("Figure 0\n",  fontsize = 30)
    for i in range(len(bp['boxes'])):
        # change outline color
        bp['boxes'][i].set( color='#7570b3', linewidth=2)
        # change fill color
        bp['boxes'][i].set( facecolor = 'rbymcgw'[i/2])
    #return bp
    plt.show()

def Figure0b_expressionCompartmentComparisson(cl1,cl2):
    ratios = expressionCompartmentAnalysis(cl1,cl2,ret  = 1)
    ratios = [[math.log(x[0]/x[1],2) for x in r] for r in ratios]
    ratios = [ratios[-1],ratios[1],ratios[2]]
    plt.boxplot(ratios,1,'',patch_artist = True)
    plt.show()

def plotHeatmap(data, keys, cmap =plt.cm.Oranges,tri = True, maxval = None):
    from matplotlib.collections import LineCollection
    import matplotlib.pyplot as plt
    import matplotlib

    maxValueInMap = data.max()
    if maxval:
        maxValueInMap = maxval
    
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            data[i,j] = min(data[i,j],maxValueInMap)
    x = np.arange(maxValueInMap+1)
    ys = [x+i for i in x]
    line_segments = LineCollection([zip(x, y) for y in ys],
                    linewidths=(0.5, 3, 1.5, 2),
                    linestyles='solid', cmap = cmap)

    line_segments.set_array(x)
    
    fig, ax = plt.subplots()
    if tri:
        mask = 1-numpy.tri(data.shape[0])
        data = numpy.ma.array(data,mask=mask)
    heatmap = ax.pcolor(data, cmap=cmap)
    a = ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
    b = ax.set_yticks(np.arange(data.shape[0])+0.5, minor=False)
    c = ax.set_xticklabels(keys, minor=False)
    d = ax.set_yticklabels(keys, minor=False)
    axcb = fig.colorbar(line_segments)
    plt.show()
    return 0

class bin:
    def __init__(self, ch, start, end):
        self.start = start
        self.end = end
        self.chromosome = ch
        self.TF = {}

    def __repr__(self):
        return "bin in chr %s from %s to %s"%(self.chromosome,self.start,self.end)

def generateTF2Pval(cl,path):
    tf2pval = {}
    tfs = glob.glob(path+  "/*")
    for tf in tfs:
        try:
            f = file(tf + "/ABenrichment.tsv").readlines()
            line = [l for l in f if cl in l][0]
        except:
            continue
        tf2pval[tf.split("/")[-1]] = float(line.split()[-1])
    return tf2pval

def Figure1C(cl, tf2pval, family = "TF"):
    import statsmodels
    import statsmodels.sandbox
    import statsmodels.sandbox.stats.multicomp
    multipletests = statsmodels.sandbox.stats.multicomp.multipletests
    try:
        fdr = multipletests(tf2pval.values(), method = "fdr_bh")
    except:
        fdr = multipletests([10**(-1*x[1]) for x in tf2pval.values()], method = "fdr_bh")
    print "FDR-BH < 0.05 for %s out of %s"%(sum(fdr[0]),len(tf2pval))
    sortedpval = sorted(tf2pval.values())
    k = tf2pval.keys()
    sortedkeys = []
    for i in range(len(sortedpval)):
        for j in range(len(k)):
            if tf2pval[k[j]] == sortedpval[i]:
                sortedkeys.append(k[j])
                break
        a = k.pop(j)
    fig, ax = plt.subplots(1, 1, figsize = (20,15))
    #a = ax.plot([math.log(-1*math.log(tf2pval[x] + 10**(-300),10),10) for x in sortedkeys])
    #a = ax.plot([math.log(2,10) for i in range(len(sortedkeys))])
    a = ax.plot([-1*tf2pval[x] for x in sortedkeys], c = 'b', linewidth = 5)
    a = ax.plot([2 for i in range(len(sortedkeys))], c = 'r')
    ax2 = ax.twinx()
    width = .4
    enrichments = []
    for x in sortedkeys:
        f = file("/home/gaga/idannuri/Data/%s/all/%s/ABenrichment.tsv"%(family,x)).readlines()
        rel = [line for line in f if cl in line][0].splitlines()[0]
        params = rel.split("\t")
        A = params[1]
        B = params[2]
        en = int(A)*1.0/int(B)
        a = ax2.bar(sortedkeys.index(x) + width/2, math.log(en,2),width=width, color = "black")
    a = ax.set_xticks(range(len(sortedkeys)))
    a = ax.set_ylabel("-log(pvalue)", fontsize = 30, color = "b")
    ax2.set_ylabel("AB Density Factor (log)",fontsize = 30, color = "black")
    """ADJUST following 4 lines according to specific values each time you use function"""
    a = ax.set_yticks(range(0,350,50))
    a = ax2.set_yticks(range(-1,5))
    a = ax.set_yticklabels([x for x in range(0,350,50)],fontsize = 25)
    a = ax2.set_yticklabels([x for x in range(-1,5)],fontsize = 25)
    
    ax.grid(True)
    a = ax.set_xticklabels(sortedkeys, fontsize = 12, rotation = "vertical")
    a = ax2.set_xticklabels(sortedkeys, fontsize = 12, rotation = "vertical")
    plt.xticks(rotation = "vertical")
    fig.tight_layout()
    plt.savefig("/home/gaga/idannuri/workdays/%sABAuto_%s"%(family,cl))
    
def Figure2b(tf2pval, T = 0, CL1 = "GM12878", CL2 = "K562", ret_dic = 0,family = "TF"):
    #import statsmodels
    #import statsmodels.sandbox
    #import statsmodels.sandbox.stats.multicomp
    #from statsmodels.sandbox.stats.multicomp import *
    #fdr = multipletests(tf2pval.values(), method = "fdr_bh")
    #print "FDR-BH < 0.05 for %s out of %s"%(sum(fdr[0]),len(tf2pval))
    sortedpval = sorted(tf2pval.values())
    k = tf2pval.keys()
    sortedkeys = []
    for i in range(len(sortedpval)):
        for j in range(len(k)):
            if tf2pval[k[j]] == sortedpval[i]:
                sortedkeys.append(k[j])
                break
        a = k.pop(j)
    fig, ax = plt.subplots(1, 1)
    ax.set_title("%s %s"%(CL1,CL2))
    #a = ax.plot([math.log(-1*math.log(tf2pval[x] + 10**(-300),10),10) for x in sortedkeys])
    #a = ax.plot([math.log(2,10) for i in range(len(sortedkeys))])
    try:
        a = ax.plot([-1*math.log(tf2pval[x] + 10**(-300),10) for x in sortedkeys], c = 'b', linewidth = 5)
    except:
        a = ax.plot([tf2pval[x] for x in sortedkeys], c = 'b', linewidth = 5)
    #a = ax.plot([2 for i in range(len(sortedkeys))], c = 'r')
    ax2 = ax.twinx()
    width = .4
    enrichments = []
    for x in sortedkeys:
        f = file("/home/gaga/idannuri/Data/%s/all/%s/result.tsv"%(family,x)).readlines()
        try:
            index = f.index("%s %s\n"%(CL1, CL2))
        except:
            index = f.index("%s %s\n"%(CL2, CL1))
        cl1 = f[index+2].split()
        cl2 = f[index+3].split()
        cl3 = f[index+4].split()
        if T==1:
            en = (int(cl1[2])*1.0/int(cl1[3]))/(int(cl3[2])*1.0/int(cl3[3]))
        elif T==2:
            en = (int(cl2[3])*1.0/int(cl2[2]))/(int(cl3[3])*1.0/int(cl3[2]))
        else:
            en = (int(cl1[2]) + int(cl2[3]))*1.0/(int(cl1[3]) + int(cl2[2]))
        enrichments.append(math.log(en,2))
        a = ax2.bar(sortedkeys.index(x) + width/2,enrichments[-1],width=width, color = "dimgray")
    if ret_dic:
        #if you want to use Figure4a(exp2pval)
        d = {}
        for i in range(len(sortedkeys)):
            key = sortedkeys[i]
            d["%s_%s_%s"%(CL1[:2],CL2[:2],key)] = [tf2pval[key], enrichments[i]]
        return d
    
    a = ax.set_xticks(range(len(sortedkeys)))
    a = ax.set_ylabel("-log(pvalue)", fontsize = 30, color = "b")
    if T:
        ax2.set_ylabel("Transition Enrichment Ratio T%s"%T,fontsize = 30, color = "dimgray")
    else:
        ax2.set_ylabel("Occypancy Enrichment Ratio",fontsize = 30, color = "dimgray")
    """ADJUST following 2 lines according to specific values each time you use function"""
    #a = ax.set_yticklabels([x for x in range(0,100,10)],fontsize = 25)
    #a = ax2.set_yticklabels([x for x in range(0,7)],fontsize = 25)
    
    ax.grid(True)
    a = ax.set_xticklabels(sortedkeys, fontsize = 16, rotation = "vertical")
    a = ax2.set_xticklabels(sortedkeys, fontsize = 16, rotation = "vertical")
    plt.xticks(rotation = "vertical")
    #plt.show()
    plt.tight_layout()
    plt.savefig("/home/gaga/idannuri/workdays/Figure2/%s_%s"%(CL1,CL2))
    

def generateFigure2b_T1T2(T = 1, CL1 = "GM12878", CL2 = "K562", rv = False,family = "TF"):
    tfs = [x.split("/")[-1] for x in glob.glob("/home/gaga/idannuri/Data/%s/all/*"%(family))]
    tdic = {}
    for x in tfs:
        try:
            f = file("/home/gaga/idannuri/Data/%s/all/%s/result.tsv"%(family,x)).readlines()
	    index = f.index("%s %s\n"%(CL1,CL2))
	except:
		try:
		    index = f.index("%s %s\n"%(CL2,CL1))
		except:
			continue
	    
	cl1 = f[index+2].split()
	cl2 = f[index+3].split()
	cl3 = f[index+4].split()
	AB1 = int(cl1[2]);BA1 = int(cl1[3])
	AB2 = int(cl2[2]);BA2 = int(cl2[3])
	AB3 = int(cl3[2]);BA3 = int(cl3[3])
	if T==1:
            pval = chisq([AB1,BA1],[AB3,BA3])
        if T == 2:
            pval = chisq([AB2,BA2],[AB3,BA3])
        else:
            pval = chisq([AB1,BA1],[AB2,BA2])
	tdic[x] = pval
	
    if rv:
        return tdic
    Figure2b(tdic,T,CL1,CL2)

def updateDicBins(TFdic, TFname, keys = None):
    if not keys:
        keys = TFdic.keys()
    for key in dic_bins:
        if TFname not in dic_bins[key].TF:
            dic_bins[key].TF[TFname] = []
    for key in keys:
        for chip in TFdic[key]:
            if chip.chromosome == 'chrY':
                continue
            if chip.chromosome == 'chrX':
                chip.chromosome = 23
            start = chip.start/5000*5000
            tmp = bin(chip.chromosome, start, start+5000)
            dic_bins[str(tmp)].TF[TFname].append(key)
    print "Do you want to save the new version?"
    return 0

def TFABratio(TFs, cl1,cl2):
    groups = [0,0,0,0]
    for tf in TFs:
        A1 = False
        A2 = False
        try:
            start = tf.TSS
            end = tf.TSS
        except:
            start = tf.start
            end = tf.end
        for a in [d for d in hun[cl1][0] if d.chromosome == tf.chromosome]:
            if start >= a.start and end <= a.end:
                A1 = True
                break
        for a in [d for d in hun[cl2][0] if d.chromosome == tf.chromosome]:
            if start >= a.start and end <= a.end:
                A2 = True
                break
        groups[3 - (A1*2+A2)] += 1
    return groups
    
def generateChipSeqTable(path, indexes = None, cl = None, title = None):
    """
    function compares chip seq of TF according to A\B compartments and return beatuiful table
    :param path: path of different chip seq files, filename need to be the cell line
    :return: saves result file in the given path
    """		
    files = [x for x in glob.glob(path + "/*") if "." not in x and "cell" not in x]
    result = ""
    for k1 in range(len(files)):
        for k2 in range(k1):
            key1 = files[k1].split("/")[-1]
            key2 = files[k2].split("/")[-1]
            if cl:
                """ in case we want to return the pvalue of a certain pair
                for example for figure 2b
                """
                if not ((cl[0] == key1 and cl[1] == key2) or (cl[0] == key2 and cl[1] == key1)):
                    continue
            print key1, key2
            result += "%s %s\n"%(key1,key2)
            result += "\tAA\tAB\tBA\tBB\ttotal\t\tAA%\tAB%\tBA%\tBB%\t\tR\tEnrichment\t\tp-value\n"
            os.system("bedtools intersect -wa -wb -a %s -b %s -v > %s/cell1_only"\
                      % (files[k1], files[k2],path))
            os.system("bedtools intersect -wa -wb -a %s -b %s -v > %s/cell2_only" \
                      % (files[k2], files[k1],path))
            os.system("bedtools intersect -wa -wb -a %s -b %s > %s/cell_common" \
                      % (files[k1], files[k2],path))

            temp = [path+"/cell1_only",path+"/cell2_only",path+"/cell_common"]
            temp = [map(TF, file(f).readlines()) for f in temp]
            ABdiv = [TFABratio(t, key1, key2) for t in temp]
            totals = [sum(abdiv) for abdiv in ABdiv]
            try:
                ABratio = [map(lambda x: round(x * 1.0 / totals[i], 2), ABdiv[i]) \
                       for i in range(len(ABdiv))]
                R = [ABratio[0][1]/ABratio[0][2],ABratio[1][2]/ABratio[1][1]]
            except:
                return ABdiv, totals
            pval = chisq(ABdiv[0][1:3], ABdiv[1][1:3])
            if cl:
                print "idan"
                return pval
            result += "Cell1_only_BSs\t" + '\t'.join(map(str, ABdiv[0]))+"\t%s\t\t"%(totals[0]) + \
                      '\t'.join(map(str, ABratio[0])) + "\t\t%s\t%s\t\t%s\n"%(R[0],R[0]*R[1],pval)
            result += "Cell2_only_BSs\t" + '\t'.join(map(str, ABdiv[1])) + "\t%s\t\t" % (totals[1]) + \
                      '\t'.join(map(str, ABratio[1])) + "\t\t%s\n"%(R[1])
            result += "Common_BSs\t" + '\t'.join(map(str, ABdiv[2])) + "\t%s\t\t" % (totals[2]) + \
                      '\t'.join(map(str, ABratio[2])) + "\n"

    rfile = file(path+"/result.tsv",'wb')
    rfile.write(result)
    rfile.close()
    os.system("rm %s/cell*" % (path))
    return 0

def generateFigure1c_ABautoTFBSenrichment(path):
    """
    function checks enrichment for A or B in TFBS
    :param path: path of different chip seq files, filename need to be the cell line
    :return: saves result file in the given path
    """		
    files = [x for x in glob.glob(path + "/*") if "." not in x and "cell" not in x]
    result = "Cell line\tA\tB\ttotal\t\tA%\tB%\t\tp-value\n"
    for k1 in range(len(files)):
        key1 = files[k1].split("/")[-1]
        print key1
        tfbs = map(TF,file(files[k1]).readlines())
        print len(tfbs)
        AB = check_regexp_AB(hun[key1][0],hun[key1][1],tfbs)
        print AB
        lenA = sum([x.end-x.start for x in hun[key1][0]])
        lenB = sum([x.end-x.start for x in hun[key1][1]])
        pval = round(math.log(chisq([lenA,lenB], AB)+10**(-300),10),2)
        result += "%s\t%s\t%s\t%s\t\t%s\t%s\t\t%s\n"\
                  %(key1, AB[0],AB[1], sum(AB), sum(AB)*1.0/AB[0],sum(AB)*1.0/AB[1],\
                    pval)
    rfile = file(path+"/ABenrichment.tsv",'wb')
    rfile.write(result)
    rfile.close()
    return 0

def generateFigure3b_domains(cl1, cl2, thresh, sm = 3):
    import matplotlib.pyplot as plt
    aa = expressionCompartmentAnalysis(cl1, cl2, AA=True)
    cl1Highcl2Low = [x for x in aa if  expressions[cl1][x.ID] > thresh and expressions[cl2][x.ID] <= 1]
    cl2Highcl1Low = [x for x in aa if  expressions[cl1][x.ID] <= 1 and expressions[cl2][x.ID] > thresh]
    relative1 = findRelativeLocationTest(cl1Highcl2Low, domains[cl1])
    relative1 = [x*1.0/sum(relative1) for x in relative1]
    relative2 = findRelativeLocationTest(cl1Highcl2Low, domains[cl2])
    relative2 = [x * 1.0 / sum(relative2) for x in relative2]
    relative3 = findRelativeLocationTest(cl2Highcl1Low, domains[cl1])
    relative3 = [x * 1.0 / sum(relative3) for x in relative3]
    relative4 = findRelativeLocationTest(cl2Highcl1Low, domains[cl2])
    relative4 = [x * 1.0 / sum(relative4) for x in relative4]
    fig,ax = plt.subplots(2,1)
    ax[0].plot([sum(relative1[i - sm:i + sm]) for i in range(sm, len(relative1) - sm)])
    ax[0].plot([sum(relative2[i - sm:i + sm]) for i in range(sm, len(relative2) - sm)])
    ax[0].legend([cl1,cl2])
    ax[0].set_xticklabels(np.arange(-20,130,20))
    ax[1].plot([sum(relative3[i - sm:i + sm]) for i in range(sm, len(relative3) - sm)])
    ax[1].plot([sum(relative4[i - sm:i + sm]) for i in range(sm, len(relative4) - sm)])
    ax[1].legend([cl1, cl2])
    ax[1].set_xticklabels(np.arange(-20, 130, 20))
    plt.savefig("/home/gaga/idannuri/workdays/Figure 3/Figure3D/%s_%s"%(cl1,cl2))
    return relative1,relative2,relative3,relative4

def generateFigure3a_HiC(cl,thresh = 1):
    aa = expressionCompartmentAnalysis(cl, cl, AA=True)
    expressed = [x.ID for x in aa if\
                 expressions[cl][x.ID] > thresh]
    nexpressed = [x.ID for x in aa if\
                 expressions[cl][x.ID] <= 1]
    enexp = sum_list([map(Enhancer,enhs[cl][x][1:]) for x in expressed\
                       if x in enhs[cl]])
    ennexp = sum_list([map(Enhancer,enhs[cl][x][1:]) for x in nexpressed\
                       if x in enhs[cl]])
    return expressed, nexpressed,enexp,ennexp

def Figure3a_HiC_1(cls, ylabel = "FPKM", pval = 0.05):
    import matplotlib.pyplot as plt
    import math
    tmp = {}
    for cl in cls:
        tmp[cl] = expressionCompartmentAnalysis(cl, cl, AA=True)
    aa = [sum_list([[expressions[cl][x.ID] for x in tmp[cl] if x.ID in enhs[cl]\
           and len([en for en in enhs[cl][x.ID] if en.pval < pval])==i]\
                    for cl in cls]) for i in range(1,5)]

    aa.append(sum_list([[expressions[cl][x.ID] for x in tmp[cl] if x.ID in enhs[cl]\
           and len([en for en in enhs[cl][x.ID] if en.pval < pval])>4]\
                    for cl in cls]))

    aa = [map(lambda x:math.log(x+1,2), li) for li in aa]
    print map(len,aa)
    plt.boxplot(aa, 0 , '')
    plt.ylabel(ylabel)
    plt.show()

def Figure3a_HiC_2(cls, ylabel = "number of interactions", pval = 0.05, thresh = 10):
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as GS
    import math
    tmp = {}
    for cl in cls:
        tmp[cl] = expressionCompartmentAnalysis(cl, cl, AA=True)
    nexp = [[len([en for en in enhs[cl][x.ID] if en.pval < pval]) for x in tmp[cl] if x.ID in enhs[cl]\
           and expressions[cl][x.ID] <= 1] for cl in cls]

    exp = [[len([en for en in enhs[cl][x.ID] if en.pval < pval]) for x in tmp[cl] if x.ID in enhs[cl]\
           and expressions[cl][x.ID] > thresh] for cl in cls]

    print map(len, nexp)
    print map(len, exp)
    fig = plt.figure(figsize=(18, 12))
    fig.subplots_adjust(wspace=1,hspace=0.5)
    gs = GS.GridSpec(4,4)
    for i in range(len(cls)):
        print i
        if i%2:
            ax = fig.add_subplot(gs[i/2,2:])
        else:
            ax = fig.add_subplot(gs[i/2,:2])
        bp = ax.boxplot([nexp[i],exp[i]], patch_artist = True)
        for box in bp['boxes']:
            # change outline color
            box.set( color='#7570b3', linewidth=2)
            # change fill color
            box.set( facecolor = '#1b9e77' )
        for whisker in bp['whiskers']:
            whisker.set(color='#7570b3', linewidth=2)
        ## change color and linewidth of the medians
        for median in bp['medians']:
            median.set( linewidth=4)
        for cap in bp['caps']:
            cap.set(color='#7570b3', linewidth=2)
        #for flier in bp['fliers']:
        #    flier.set(marker='o', color='#e7298a', alpha=0.5)
        ax.set_xticklabels(['FPKM < 1', 'FPKM > %s'%thresh])
        ax.set_ylabel(ylabel)
        ax.set_xlabel("expression %s"%cls[i])
        ax.set_ylim(0,5)
    plt.show()

def Figure3a_HiC_3(cls, ylabel = "number of interactions", pval = 0.05, thresh = 10):
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as GS
    import math
    tmp = {}
    for cl in cls:
        tmp[cl] = expressionCompartmentAnalysis(cl, cl, AA=True)
    nexp = [[len([en for en in enhs[cl][x.ID] if en.pval < pval]) for x in tmp[cl] if x.ID in enhs[cl]\
           and expressions[cl][x.ID] <= 1] for cl in cls]

    exp = [[len([en for en in enhs[cl][x.ID] if en.pval < pval]) for x in tmp[cl] if x.ID in enhs[cl]\
           and expressions[cl][x.ID] > thresh] for cl in cls]

    print map(len, nexp)
    print map(len, exp)
    fig,ax = plt.subplots(1,1)
    bp = ax.boxplot(nexp + [[],[]] + exp, patch_artist = True)
    for whisker in bp['whiskers']:
        whisker.set(color='#7570b3', linewidth=2)
    ## change color and linewidth of the medians
    for median in bp['medians']:
        median.set( linewidth=4)
    for cap in bp['caps']:
        cap.set(color='#7570b3', linewidth=2)
    for flier in bp['fliers']:
        flier.set(marker='o', color='#e7298a', alpha=0.5)
    ax.set_xticklabels(['%s\n FPKM < 1'%cl for cl in cls] +["",""] +  \
                       ['%s\n FPKM > %s'%(cl,thresh) for cl in cls])
    ax.set_ylabel(ylabel)
    ax.set_xlabel("expression")
    ax.set_ylim(0,5)
    for i in range(len(bp['boxes'])):
        # change outline color
        bp['boxes'][i].set( color='#7570b3', linewidth=2)
        # change fill color
        bp['boxes'][i].set( facecolor = 'rbymcgw'[i%5])
    plt.show()

def Figure3a_pdf(cl,limits = [1,10,100], bar = 1,top=0.8, division = 10):
    aa = [x.ID for x in expressionCompartmentAnalysis(cl, cl, AA=True)]
    axes = plt.axes()
    axes.grid()
    fig,ax = plt.subplots(1,1)
    ax.set_xlabel("number of contacts")
    ax.set_ylabel("ratio")
    nexp = map(len,[enhs[cl][x] for x in aa if x in  enhs[cl]\
                    and expressions[cl][x] <=1])
    if bar:
        bins = numpy.linspace(0,10,10)
        ratios = [[nexp]]
        ax.set_yticklabels([x*1.0/division for x in range(10)], fontsize = 26)
        ax.set_xticklabels([x*1.0 for x in range(10)], fontsize = 26)
        for l in limits:
            exp = map(len,[enhs[cl][x] for x in aa if x in enhs[cl] \
                            and expressions[cl][x] >l])
            ratios.append(exp)
        ax.hist(ratios,normed = 1)
        ax.set_title("%s"%cl, fontsize = 35)	
	ax.set_ylabel("fraction", fontsize = 30)
	pval = round(math.log(scipy.stats.ks_2samp(exp,nexp)[1],10),2)
        ax.set_xlabel("number of interactions\nlog10 p-val=%s"%pval, fontsize = 24)
        ax.legend(['FPKM <  1'] + ["FPKM > %s"%l for l in limits], fontsize = 26)
        ax.set_ylim(bottom=0,top = top)
    else:
        ax.plot(range(10),[nexp.count(i)*1.0/len(nexp) for i in range(10)])
        for l in limits:
            exp = map(len,[enhs[cl][x] for x in aa if x in enhs[cl]\
                        and expressions[cl][x] >l])
            ax.plot(range(10),[exp.count(i)*1.0/len(exp) for i in range(10)])
        ax.legend(['FPKM<1'] + ["FPKM > %s"%l for l in limits])
    plt.show()

def Figure3a_cdf(cl,limits = [1,10,100]):
    axes = plt.axes()
    axes.grid(True)
    plt.xlabel("number of contacts")
    plt.ylabel("ratio")
    nexp = map(len,[enhs[cl][x] for x in enhs[cl]\
                    if expressions[cl][x] <1])
    ratios = [nexp.count(i)*1.0/len(nexp) for i in range(10)]
    plt.plot(range(11),[sum(ratios[:i]) for i in range(11)])
    for l in limits:
        exp = map(len,[enhs[cl][x] for x in enhs[cl]\
                    if expressions[cl][x] >l])
        ratios = [exp.count(i) * 1.0 / len(exp) for i in range(10)]
        plt.plot(range(11), [sum(ratios[:i]) for i in range(11)])
    plt.legend(['not expressed'] + ["FPKM > %s"%l for l in limits])
    plt.show()

def Figure3a_pdf_all(cls = None, limits = [1,10,100], bar = 1,top = 0.8):
    if cls == None:
        cls = ['HUVEC', 'K562', 'HMEC', 'IMR90', 'GM12878']
    print cls
    fig,ax = plt.subplots(5,1)
    fig.subplots_adjust(wspace=0.5,hspace=1)
    res = []
    for i in range(5):
		cl = cls[i]
		aa = [x.ID for x in expressionCompartmentAnalysis(cl, cl, AA=True)]
		ax[i].grid(True)
		
		ax[i].set_xticks(range(10))
		ax[i].set_xticklabels(range(10),fontsize = 14)
		ax[i].set_xlim(0,5)
		ax[i].set_title("%s"%cl, fontsize = 25)	
		ax[i].set_ylabel("fraction", fontsize = 20)
		nexp = map(len,[enhs[cl][x] for x in aa if x in enhs[cl] and expressions[cl][x] <1])
		if bar:
                    bins = numpy.linspace(0,10,10)
                    ratios = [[nexp]]
                    ax[i].set_yticklabels([x*1.0/10 for x in range(10)], fontsize = 16)
                    for l in limits:
                        exp = map(len,[enhs[cl][x] for x in aa if x in enhs[cl] \
                                        and expressions[cl][x] >l])
                        ratios.append(exp)
                    res.append(ax[i].hist(ratios,normed = 1))
                    ax[i].set_ylim(top = top)
                else:
                    ratios = [nexp.count(k)*1.0/len(nexp) for k in range(10)]
                    ax[i].plot(range(10),ratios)
                    ax[i].set_ylim(top = 0.7)
                    for l in limits:
                            exp = map(len,[enhs[cl][x] for x in aa if x in enhs[cl] \
                                        and expressions[cl][x] >l])
                            ratios = [exp.count(k) * 1.0 / len(exp) for k in range(10)]
                            ax[i].plot(range(10), ratios)
                pval = round(math.log(scipy.stats.ks_2samp(exp,nexp)[1],10),2)
                ax[i].set_xlabel("number of interactions\nlog10 p-val=%s"%pval, fontsize = 20)
                ax[i].legend(['FPKM <  1'] + ["FPKM > %s"%l for l in limits], fontsize = 16)
    plt.show()


def Figure3a_inverse_pdf_all(cls = None, limits = [1,100]):

    """
    boxplot like figure 22 in the paper
    """
    
    if cls == None:
        cls = ['HUVEC', 'K562', 'HMEC', 'IMR90', 'GM12878',"NHEK"]

    print cls
    fig,ax = plt.subplots(3,2,figsize = (8*0.9,20*0.9))
    fig.subplots_adjust(wspace=1,hspace=1)
    for i in range(3):
	for j in range(2):
		cl = cls[i*2+j]
		aa = [x.ID for x in expressionCompartmentAnalysis(cl, cl, AA=True)]
		ax[i][j].grid(True)
		
		ax[i][j].set_title("%s"%cl, fontsize = 25)	
		ax[i][j].set_ylabel("log2 FPKM", fontsize = 20)
		ax[i][j].set_yticklabels(range(0,16,2), fontsize = 14)
		ninteract = [math.log(expressions[cl][x]+1,2) for x in aa if x in enhs[cl] and len(enhs[cl][x]) == 0] 
		#oneinteract = [math.log(expressions[cl][x]+1,2) for x in aa if x in enhs[cl] and len(enhs[cl][x]) > 0 and len(enhs[cl][x])<limits[0]]
		manyinteract = [math.log(expressions[cl][x]+1,2) for x in aa if x in enhs[cl] and len(enhs[cl][x]) > limits[0]]
		pval = round(math.log(scipy.stats.mannwhitneyu(ninteract,manyinteract)[1],10),2)
		ax[i][j].set_xlabel("number of interactions\n log10 pval = %s"%pval, fontsize = 20)

		ax[i][j].boxplot([ninteract, manyinteract],1,patch_artist=True,widths=0.9)
                #ax[i][j].set_xticklabels(["0","<=%s"%limits[0],">%s"%limits[0]], fontsize = 20)
		ax[i][j].set_xticklabels(["0",">=%s"%limits[0]], fontsize = 20)
    plt.show()

def Figure3a_cdf_all(cls = None, limits = [1,10,100]):
    if cls == None:
        cls = ['HUVEC', 'K562', 'HMEC', 'IMR90', 'GM12878',"NHEK"]

    print cls
    fig,ax = plt.subplots(3,2)
    fig.subplots_adjust(wspace=0.5,hspace=1)
    for i in range(3):
	for j in range(2):
		cl = cls[i*2+j]
		aa = [x.ID for x in expressionCompartmentAnalysis(cl, cl, AA=True)]
		ax[i][j].grid(True)
		ax[i][j].set_title("%s"%cl)
		ax[i][j].set_xlabel("number of interactions", fontsize = 20)
		ax[i][j].set_ylabel("percentage", fontsize = 20)
		nexp = map(len,[enhs[cl][x] for x in aa if x in enhs[cl] and expressions[cl][x] <1])
		ratios = [nexp.count(k)*1.0/len(nexp) for k in range(10)]
		ax[i][j].plot(range(11),[sum(ratios[:k+1]) for k in range(11)])
		ax[i][j].set_ylim(0,1.02)
		for l in limits:
			exp = map(len,[enhs[cl][x] for x in aa if x in enhs[cl] \
				    and expressions[cl][x] >l])
			ratios = [exp.count(k) * 1.0 / len(exp) for k in range(10)]
			ax[i][j].plot(range(11), [sum(ratios[:k+1]) for k in range(11)])
	ax[i][j].legend(['FPKM < 1'] + ["FPKM > %s"%l for l in limits])
    plt.show()

def Figure3a_pdf_chia(cls = None,limits = [1,10,100], top = 1, division = 10):
    if cls == None:
        cls = ["GM12878","MCF7","K562"]

    print cls
    fig,ax = plt.subplots(3,1)
    fig.subplots_adjust(wspace=0.5,hspace=1)
    for i in range(3):
        cl = cls[i]
        aa = [x.ID for x in expressionCompartmentAnalysis(cl, cl, AA=True)]
        ax[i].grid(True)
        ax[i].set_title("%s"%cl, fontsize = 30)
        
        ax[i].set_ylabel("fraction", fontsize = 30)
        nexp = map(len,[dic_chia[cl][x] for x in aa if x in dic_chia[cl] and expressions[cl][x] <1])
        ratios = [[nexp]]
        #ratios = [nexp.count(k)*1.0/len(nexp) for k in range(10)]
        #ax[i].plot(range(10),ratios)
        ax[i].set_ylim(0,top = top)
        ax[i].set_xlim(0,10)
        for l in limits:
            exp = map(len,[dic_chia[cl][x] for x in aa if x in dic_chia[cl] \
                        and expressions[cl][x] >l])
            ratios.append(exp)
            #ratios = [exp.count(k) * 1.0 / len(exp) for k in range(10)]
            #ax[i].plot(range(10), ratios)
        bins = numpy.linspace(0,10,10)
        ax[i].hist(ratios,bins = bins,  normed = 1)
        pval = round(math.log(scipy.stats.mannwhitneyu(nexp,exp)[1]+10**(-300),10),2)
        ax[i].set_xlabel("number of interactions\n log10 pval < %s"%pval, fontsize = 30)
        ax[i].legend(['FPKM < 1'] + ["FPKM > %s"%l for l in limits],fontsize = 16)
        ax[i].set_yticklabels([x*1.0/division for x in range(10)], fontsize = 16)
        ax[i].set_xticks(range(9))
        ax[i].set_xticklabels(range(9), fontsize = 20)
    #return ratios
    plt.show()

def Figure3a_inverse_pdf_chia(cls = None, limits = [1,100], family = None):
    if cls == None:
        cls = ["MCF7",'K562', 'GM12878']

    print cls
    fig,ax = plt.subplots(3,1)
    fig.subplots_adjust(wspace=0.5,hspace=1)
    for i in range(3):
	    cl = cls[i]
            aa = [x.ID for x in expressionCompartmentAnalysis(cl, cl, AA=True)]
            if family:
                #let say we only want housekeeping genes
                #family = load("/home/gaga/idannuri/HK/HK_TSSIds.pkl")
                aa = [x for x in aa if x in family]
            ax[i].grid(True)
            
            ax[i].set_title("%s"%cl, fontsize = 25)	
            ax[i].set_ylabel("log2 FPKM", fontsize = 20)
            ninteract = [math.log(expressions[cl][x]+1,2) for x in aa if x in dic_chia[cl] and len(dic_chia[cl][x]) == 0] 
            oneinteract = [math.log(expressions[cl][x]+1,2) for x in aa if x in dic_chia[cl] and len(dic_chia[cl][x]) == limits[0]]
            manyinteract = [math.log(expressions[cl][x]+1,2) for x in aa if x in dic_chia[cl] and len(dic_chia[cl][x]) > limits[0] and len(dic_chia[cl][x])<limits[1]]
            pval = round(math.log(scipy.stats.mannwhitneyu(ninteract,oneinteract+manyinteract)[1]+10**(-300),10),2)
            ax[i].set_xlabel("number of interactions\n log10 pval < %s"%pval, fontsize = 20)

            ax[i].boxplot([ninteract, oneinteract, manyinteract],1,patch_artist=True)
            ax[i].set_xticklabels(["0","1",">1"], fontsize = 20)
            ax[i].set_aspect(0.1)
    plt.show()


def Figure3b_HiC(cls, pval = 0.05,keys = [0,1,3,100], log = True,bins = 15,\
                 output = "", colorbar = None, to_R = False):
    import numpy as np
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import AxesGrid
    import math

    if colorbar == None:
        colorbar = cm.Blues
    all_data = [[[[],[]] for j in range(len(keys)-1)]\
                for i in range(len(keys)-1)]
    for i in range(len(cls)):
        for j in range(i):
            cl1 = cls[i]
            cl2 = cls[j]
            aa = expressionCompartmentAnalysis(cl1, cl2, AA=True)
            data = []
            for i1 in range(len(keys)-1):
                data.append([])
                for i2 in range(len(keys)-1):
                    gid = [x.ID for x in aa if x.ID in enhs[cl1] and x. ID in enhs[cl2] and\
                           len([en for en in enhs[cl1][x.ID] if en.pval < pval]) >= keys[i1] and\
                           len([en for en in enhs[cl1][x.ID] if en.pval < pval]) < keys[i1+1]\
                           and len([en for en in enhs[cl2][x.ID] if en.pval < pval]) >= keys[i2]\
                           and len([en for en in enhs[cl2][x.ID] if en.pval < pval]) < keys[i2+1]]

                    gid = [x for x in gid if not\
                           (expressions[cl1][x] < 1 and expressions[cl2][x] < 1)]

                    d = [[max(expressions[cl1][x],1) for x in gid],\
                         [max(expressions[cl2][x],1) for x in gid]]

                    if log:
                        d = [map(lambda x: math.log(x+1,1.1),p) for p in d]
                    data[-1].append([d[1],d[0]])
            for i1 in range(len(keys)-1):
                for i2 in range(len(keys)-1):
                    all_data[i1][i2][0] += data[i1][i2][0]
                    all_data[i1][i2][1] += data[i1][i2][1]

    if to_R:
        return all_data
    fig = plt.figure()

    grid = AxesGrid(fig, 111,
                    nrows_ncols=(len(keys)-1, len(keys)-1),
                    axes_pad=0.05,
                    share_all=True,
                    label_mode="L",
                    cbar_location="right",
                    cbar_mode="single",
                    )
    l = len(keys)-1
    for i1 in range(len(keys)-1):
        for j1 in range(len(keys)-1):
            ax = grid[(l-1-i1)*l + j1]
            ax.grid(False)
            ax.plot([0,2],[0,2],":")
            print map(len, all_data[i1][j1])
            #im = ax.scatter(data[i1][j1][0],data[i1][j1][1])
            x = all_data[i1][j1][0]
            y = all_data[i1][j1][1]
            hist,xedges,yedges = numpy.histogram2d(x,y,bins=bins)
            extent = [0, 1.99, 0, 1.99]
            
            im = ax.imshow(hist.T,extent=extent,origin='lower',\
                           cmap = colorbar)

    grid.cbar_axes[0].colorbar(im)

    for cax in grid.cbar_axes:
        cax.toggle_label(False)

    if output != "":
        fig.savefig(output)
        return 0
    plt.show()


def Figure3b_Chiapet(cls, keys=[0, 1, 3, 10000], log=True, bins=30, \
                 output="", colorbar=None, bp = 0, united = 0, parts = 0,to_R = False):
    """
    PLEASE DOCUMENT
    """
    import numpy as np
    import matplotlib.cm as cm
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import AxesGrid
    import math
    from scipy.stats import gaussian_kde

    if parts:
        cl1 = cls[0]
        cl2 = cls[1]
        sort_cl1 = sorted(map(len,dic_chia[cl1].values()))
        sort_cl2 = sorted(map(len,dic_chia[cl2].values()))
        keys_cl1 = [0, sort_cl1[len(sort_cl1)/3]+0.5,sort_cl1[len(sort_cl1)*2/3],10000]
        keys_cl2 = [0, sort_cl1[len(sort_cl2)/3]+0.5,sort_cl2[len(sort_cl1)*2/3],10000]
        #keys_cl1 = [0, sort_cl1[len(sort_cl1)/2]+0.5,10000]
        #keys_cl2 = [0, sort_cl1[len(sort_cl2)/2]+0.5,10000]
        print keys_cl1, keys_cl2
    else:
        keys_cl1 = keys
        keys_cl2 = keys

    if colorbar == None:
        colorbar = cm.Blues
    all_data = [[[[], []] for j in range(len(keys_cl1) - 1)] \
                for i in range(len(keys_cl2) - 1)]
    for i in range(len(cls)):
        for j in range(i):
            cl1 = cls[i]
            cl2 = cls[j]
            aa = expressionCompartmentAnalysis(cl1, cl2, AA=True)
            data = []
            for i1 in range(len(keys_cl1) - 1):
                data.append([])
                for i2 in range(len(keys_cl2) - 1):
                    gid = [x.ID for x in aa if x.ID in dic_chia[cl1] and x.ID in dic_chia[cl2] and \
                           len([en for en in dic_chia[cl1][x.ID]]) >= keys_cl1[i1] and \
                           len([en for en in dic_chia[cl1][x.ID] ]) < keys_cl1[i1 + 1] \
                           and len([en for en in dic_chia[cl2][x.ID] ]) >= keys_cl2[i2] \
                           and len([en for en in dic_chia[cl2][x.ID] ]) < keys_cl2[i2 + 1]]

                    gid = [x for x in gid if not \
                        (expressions[cl1][x] < 1 and expressions[cl2][x] < 1)]

                    d = [[max(expressions[cl1][x], 1) for x in gid], \
                         [max(expressions[cl2][x], 1) for x in gid]]

                    if log:
                        d = [map(lambda x: math.log(x, 2), p) for p in d]
                    data[-1].append([d[1], d[0]])
            for i1 in range(len(keys_cl1) - 1):
                for i2 in range(len(keys_cl2) - 1):
                    all_data[i1][i2][0] += data[i1][i2][0]
                    all_data[i1][i2][1] += data[i1][i2][1]
    if to_R:
        return all_data
    if united:
        #return all_data
        fig,ax = plt.subplots(1,1)
        cmaps = [cm.Pastel1, cm.Reds, cm.Greys, cm.Pastel1]
        ps = [[0,0],[1,1],[0,1],[1,0]]
        for p in ps:
            i = p[0]
            j = p[1]
            x = all_data[i][j][0]
            y = all_data[i][j][1]
            xy = np.vstack([x,y])
            z = gaussian_kde(xy)(xy)
            ax.scatter(x, y, c=z,s = 1, edgecolor='', cmap = cmaps[i*2+j])
        plt.show()
        return 7
    
    if bp == 1:
        ratios = []
        for i in range(len(all_data)):
            for j in range(len(all_data[i])):
                try:
                    ratios.append([math.log((all_data[i][j][0][k]+1)*1.0/(all_data[i][j][1][k]+1),2)\
                               for k in range(len(all_data[i][j][0]))])
                except:
                    return all_data
	fig,ax = plt.subplots(1,1)
	ax.boxplot(ratios,1,patch_artist = True)
	ax.set_ylim(-8,10)
	ax.set_yticklabels(range(-8,10,2),fontsize = 40)
	ax.set_xticklabels(["0 in CL1\n in CL2",">1 in CL1\n0 in CL2","0 in CL1\n>1 in CL2"\
                           ,">1 in CL1\n>1 in CL2"], fontsize = 34)
	pval = round(math.log(scipy.stats.mannwhitneyu(ratios[1],ratios[2])[1],10),2)
	print pval
	ax.set_xlabel("Promoter interactions\npval=%s"%pval, fontsize = 36)
	ax.set_ylabel("log2 FPKM(CL1)/FPKM(CL2)", fontsize = 40)
	
	plt.show()
	return 0
    
    fig = plt.figure()
    grid = AxesGrid(fig, 111,
                    nrows_ncols=(len(keys) - 1, len(keys) - 1),
                    axes_pad=0.05,
                    share_all=True,
                    label_mode="L",
                    cbar_location="right",
                    cbar_mode="single",
                    )
    l = len(keys) - 1
    for i1 in range(len(keys) - 1):
        for j1 in range(len(keys) - 1):
            ax = grid[(l - 1 - i1) * l + j1]
            ax.grid(False)
            ax.plot([0, 15], [0, 15], ":")
            ax.set_xlabel('log FPKM %s'%cls[0],fontsize = 25)
            ax.set_ylabel('log FPKM %s'%cls[1], fontsize = 25)
            
            ax.set_xticklabels(range(0,12,2),fontsize = 20)
            ax.set_yticklabels(range(0,12,2),fontsize = 20)
            """
            print map(len, all_data[i1][j1])
            # im = ax.scatter(data[i1][j1][0],data[i1][j1][1])
            x = all_data[i1][j1][0]
            y = all_data[i1][j1][1]
            hist, xedges, yedges = numpy.histogram2d(x, y, bins=bins)
            extent = [0, 1.99, 0, 1.99]
            im = ax.imshow(hist.T, extent=extent, origin='lower', \
                       cmap=colorbar)
            """
            x = all_data[i1][j1][0]
            y = all_data[i1][j1][1]
            xy = np.vstack([x,y])
            z = gaussian_kde(xy)(xy)
            ax.set_xlim(0,12)
            ax.set_ylim(0,12)
            sc = ax.scatter(x, y, c=z,s = 50, edgecolor='', cmap = cm.Reds)

    fig.colorbar(sc)
    #grid.cbar_axes[0].colorbar(im)

    #for cax in grid.cbar_axes:
     #   cax.toggle_label(False)

    if output != "":
        fig.savefig(output)
        return 0
    plt.show()


                    
def processChiapetFile_TSS(filename, key):
    f = file(filename).readlines()
    chpet = [[line for line in f if line.count("chr%s:"%i) ==  2]\
             for i in range(1,23) + ["X"]]
    chpet = [map(Chiapet, x) for x in chpet]
    print map(len,chpet)
    dump(chpet, "./Data/Chiapet/%s_intra"%(key))
    dic_chia[key] = {}
    for t in ids:
        dic_chia[key][t] = []
    for i in range(len(TSSs)):
	tss = TSSs[i]
	if tss.chromosome == "X":
		ind = 22
	elif type(tss.chromosome) == type(1):
		ind = tss.chromosome -1
	else:
		continue
	for con in chpet[ind]:
		if tss.TSS >= con.start1 and tss.TSS <= con.end1:
			dic_chia[key][tss.ID].append(con)
			continue
		if tss.TSS >= con.start2 and tss.TSS <= con.end2:
			dic_chia[key][tss.ID].append(con)
	if i%1000 == 0:
		print i,
    print "don't forget to save dic_chia"

def processChiapetFile_Transcripts(filename, key):
    f = file(filename).readlines()
    chpet = [[line for line in f if line.count("chr%s:"%i) ==  2]\
             for i in range(1,23) + ["X"]]
    chpet = [map(Chiapet, x) for x in chpet]
    print map(len,chpet)
    dump(chpet, "/home/gaga/idannuri/Data/Chiapet/%s_intra"%(key))
    dic_chia_full[key] = {}
    for t in transcripts.keys():
        # transcripts have multiple keys for each trans we only use one
        if "." not in t: 
            dic_chia_full[key][t] = []
    for i in range(len(transcripts)):
        # transcripts have multiple keys for each trans we only use one
        if "." in transcripts.keys()[i]:
            continue
	trans = transcripts.values()[i]
	if trans.chromosome == "chrX":#UNTIL we fix it to be "X"
		ind = 22
	elif type(trans.chromosome) == type(1):
		ind = tra97ns.chromosome -1
	else:
		continue
	for con in chpet[ind]:
                mid = trans.start + (trans.end-trans.start)/2
		if mid >= con.start1 and mid <= con.end1:
			dic_chia_full[key][transcripts.keys()[i]].append(con)
			continue
		if mid >= con.start2 and mid <= con.end2:
			dic_chia_full[key][transcripts.keys()[i]].append(con)
	if i%1000 == 0:
		print i,
    print "don't forget to save dic_chia_full"

def generateFigure3a_chiapet(cl,thresh = 1):
    aa = expressionCompartmentAnalysis(cl, cl, AA=True)
    expressed = [x.ID for x in aa if\
                 expressions[cl][x.ID] > thresh]
    nexpressed = [x.ID for x in aa if\
                 expressions[cl][x.ID] <= 1]
    conexp = sum_list([dic_chia[cl][x] for x in expressed])
    connexp = sum_list([dic_chia[cl][x] for x in nexpressed])
    return expressed, nexpressed,conexp, connexp

def generateFigure3b_rawHiCstep1(cl1, cl2, res = 100*1000, thresh = 10):
    import matplotlib.pyplot as plt
    aa = expressionCompartmentAnalysis(cl1, cl2, AA=True)
    cl1Highcl2Low = [x for x in aa if \
                     expressions[cl1][x.ID]/max(expressions[cl2][x.ID],1) > thresh]
    cl2Highcl1Low = [x for x in aa if \
                     expressions[cl2][x.ID] / max(expressions[cl1][x.ID],1) > thresh]
    print len(cl1Highcl2Low),len(cl2Highcl1Low)
    hic1 = [getHicValuesForTSS(t, cl1, res) for t in cl1Highcl2Low]
    print "finishied hic1"
    hic2 = [getHicValuesForTSS(t, cl2, res) for t in cl1Highcl2Low]
    print "finishied hic1"
    hic3 = [getHicValuesForTSS(t, cl1, res) for t in cl2Highcl1Low]
    print "finishied hic1"
    hic4 = [getHicValuesForTSS(t, cl2, res) for t in cl2Highcl1Low]
    print "finishied hic1"
    return hic1, hic2, hic3, hic4

def generateFigure3b_rawHiCstep2(cl1,cl2):
    aa = expressionCompartmentAnalysis(cl1, cl2, AA=True)
    cl1Highcl2Low = [x for x in aa if \
                     expressions[cl1][x.ID]/max(expressions[cl2][x.ID],1) > thresh]
    cl2Highcl1Low = [x for x in aa if \
                     expressions[cl2][x.ID] / max(expressions[cl1][x.ID],1) > thresh]
    res = load("./workdays/Figure 3/%s_%s"%(cl1,cl2))
    print "finished preprocessing"
    result1 = []
    control = []
    for tssi in range(len(cl1Highcl2Low)):
        tss = cl1Highcl2Low[tssi]
        en = [enhs[cl1][tss.ID][1:],enhs[cl2][tss.ID][1:]]
        if en[0] and en[1]:
            en = [map(TF, eni) for eni in en]
            start = [min([x.start for x in eni]) for eni in en]
            end = [max([x.start for x in eni]) for eni in en]
            r1 = [x.Mij for x in res[0][tssi] if\
                  ((x.i > start[0]) and (x.i < end[0]))\
                  or ((x.j > start[0]) and (x.j < end[0]))]
            if not r1:
                r1 =[0]
            r2 = [x.Mij for x in res[1][tssi] if\
                  ((x.i > start[1]) and (x.i < end[1]))\
                  or ((x.j > start[1]) and (x.j < end[1]))]
            if not r2:
                r2 =[0]
            c1 = [x.Mij for x in res[0][tssi]]
            c2 = [x.Mij for x in res[1][tssi]]
            result1.append([sum(r1)/len(r1), sum(r2)/len(r2)])
            control.append([sum(c1)/len(c1), sum(c2)/len(c2)])
    return result1, control
    
def generateFigure3b_HiC(cl1,cl2, thresh = 10):
    aa = expressionCompartmentAnalysis(cl1, cl2, AA=True)
    cl1Highcl2Low = [x for x in aa if \
                     expressions[cl1][x.ID]/max(expressions[cl2][x.ID],1) > thresh]
    cl2Highcl1Low = [x for x in aa if \
                     expressions[cl2][x.ID] / max(expressions[cl1][x.ID],1) > thresh]
    cl11en = sum_list([map(Enhancer,enhs[cl1][tss.ID][1:]) for tss in cl1Highcl2Low\
                       if tss.ID in enhs[cl1]])
    cl21en = sum_list([map(Enhancer,enhs[cl2][tss.ID][1:]) for tss in cl1Highcl2Low\
                       if tss.ID in enhs[cl2]])
    cl12en = sum_list([map(Enhancer,enhs[cl1][tss.ID][1:]) for tss in cl2Highcl1Low\
                       if tss.ID in enhs[cl1]])
    cl22en = sum_list([map(Enhancer,enhs[cl2][tss.ID][1:]) for tss in cl2Highcl1Low\
                       if tss.ID in enhs[cl2]])
    return cl11en, cl21en, cl12en, cl22en
    
def generateFigure3b_Chiapet(cl1,cl2, thresh = 10):
    aa = expressionCompartmentAnalysis(cl1, cl2, AA=True)
    cl1Highcl2Low = [x for x in aa if \
                     expressions[cl1][x.ID]/max(expressions[cl2][x.ID],1) > thresh]
    cl2Highcl1Low = [x for x in aa if \
                     expressions[cl2][x.ID] / max(expressions[cl1][x.ID],1) > thresh]
    cl11en = sum_list([dic_chia[cl1][tss.ID] for tss in cl1Highcl2Low\
                       if tss.ID in dic_chia[cl1]])
    cl21en = sum_list([dic_chia[cl2][tss.ID] for tss in cl1Highcl2Low\
                       if tss.ID in dic_chia[cl2]])
    cl12en = sum_list([dic_chia[cl1][tss.ID] for tss in cl2Highcl1Low\
                       if tss.ID in dic_chia[cl1]])
    cl22en = sum_list([dic_chia[cl2][tss.ID] for tss in cl2Highcl1Low\
                       if tss.ID in dic_chia[cl2]])
    return cl11en, cl21en, cl12en, cl22en


def Figure3c(cls = None, limits=[10, 100]):
    if cls == None:
        cls = ["GM12878", "K562","HUVEC","HMEC","NHEK","IMR90"]
    nrows = max(len(cls)/2 + len(cls)%2,2)
    fig,ax = plt.subplots(nrows,2)
    fig.subplots_adjust(wspace=1, hspace=0.5)
    for i in range(nrows):
        ncols = min(len(cls) - 2*i,2)
        for j in range(ncols):
            cl = cls[i*2+j]
            aa = [x.ID for x in expressionCompartmentAnalysis(cl, cl, AA=True)]
            print len(aa), len(expressions[cl].keys())
            nexp = [x for x in aa if x in expressions[cl] and expressions[cl][x] < 1]
            rel = findRelativeLocationTest([ids[x] for x in nexp], domains[cl])
            norm = [x * 1.0 / sum(rel) for x in rel]
            ax[i][j].plot(norm)
            for l in limits:
                exp = [x for x in aa if x in expressions[cl] and expressions[cl][x] > l]
                rel = findRelativeLocationTest([ids[x] for x in exp], domains[cl])
                norm = [x * 1.0 / sum(rel) for x in rel]
                ax[i][j].plot(norm)
            ax[i][j].legend(['not expressed'] + ["FPKM > %s" % l for l in limits], fontsize = 24)
            ax[i][j].set_xlabel("relative location",fontsize = 24)
            ax[i][j].set_ylabel("ratio",fontsize = 24)
            #ax[i][j].set_yticks(range(8))
            ax[i][j].set_yticks([0.005*k for k in range(8)])
            #ax[i][j].set_xticks(np.arange(8))
            ax[i][j].set_xticklabels(np.arange(-20,130,20),fontsize = 24)
            ax[i][j].set_title("%s"%cl,fontsize = 30)
    plt.show()
    return 0

def generateFigure4a_HiC(expressedIds, cl, ret = 0):
    res = sum_list([map(Enhancer,enhs[cl][tss][1:]) for tss in expressedIds\
                       if tss in enhs[cl]])
    if ret == 1:
        return res
    bg = sum(map(len,enhs[cl].values())) - len(enhs[cl])
    return chisq([len(expressedIds), len(res)],[len(ids), bg])

def generateFigure4a_Chipseq(pair, cl, compare = None, ret = 0):
    peaks_ind = "/home/gaga/data-scratch/idannuri/ChipSeq/Step3/macs2_out/%s_%s_%s_induced_peaks.narrowPeak"%(cl,pair[0],pair[1])
    peaks_rep = "/home/gaga/data-scratch/idannuri/ChipSeq/Step3/macs2_out/%s_%s_%s_repressed_peaks.narrowPeak"%(cl,pair[1],pair[0])
    try:
        peaks_ind = map(Chipseq, file(peaks_ind).readlines())
        peaks_rep = map(Chipseq, file(peaks_rep).readlines())
    except:
        print pair, " no matching files in macs2_out"
        return ""
    
    if compare:
        cl = compare
    ind = check_regexp_AB(hun[cl][0],hun[cl][1], peaks_ind)
    rep = check_regexp_AB(hun[cl][0],hun[cl][1], peaks_rep)
    lenA = sum([x.end-x.start for x in hun[cl][0]])
    lenB = sum([x.end-x.start for x in hun[cl][1]])
    pval = [chisq([x+1 for x in ind],[lenA, lenB]), chisq([x+1 for x in rep],[lenA, lenB])]
    pval = [round(math.log(x+10**(-300),10),2) for x in pval]
    header = '\t'.join(["Cell line", "Control", "Treatment", "induced in A", "induced in B","log p-value",\
                        "repressed in A", "repressed in B", "log p-value"]) + "\n"
    result = '\t'.join(map(str, [cl] + pair + list(ind) + [pval[0]] + list(rep) + [pval[1]])) + "\n"
    print header + result
    if ret:
        return result
    output = file("/home/gaga/data-scratch/idannuri/ChipSeq/Results/%s_%s_%s.txt"%(cl, pair[0],pair[1]), "wb")
    output.write(header + result)
    output.close()
    return 0
    

def mult_generateFigure4a_Chipseq(filename):
    tmp = file(filename).readlines()
    pairs = [x.split("\t") for x in tmp]
    pairs = [[x.split()[0] for x in p] for p in pairs]
    header = '\t'.join(["Cell line", "Control", "Treatment", "induced in A", "induced in B","log p-value",\
                        "repressed in A", "repressed in B", "log p-value"]) + "\n"
    rows = ""
    for pair in pairs:
        p0 = pair[0]
        p1 = pair[1]
        cl = pair[2]
        rows += generateFigure4a_Chipseq([p0,p1], cl, ret = 1)
    rf = filename.split("/")[-1].split(".")[0]
    output = file("/home/gaga/data-scratch/idannuri/ChipSeq/Results/%s_4a_results.txt"%(rf), "wb")
    output.write(header + rows)
    output.close()
    return 0

def generateFigure4a_processedChipseq(indexes, cl, tf, treatment, Debug = False):
    """
    TFBS is a file containing all TFBS and the cell lines that we see them in
    example:
    :param indexes: must be provided as a list of strings! running the function with indexes
    generates all binding site for the cell lines corresponding to the indexes in TFBS file from encode.
    First index is control second is after treatment. 
    :param cl: cell line name of the indexes ("K562" for example)
    :param title: the name of the experiement ("K562_pol2_IFNa" for example, for K562
    with chip seq of pol2 before and after treated with IFNa
    :return: the function returns Figure 4a for the experiement
    """
    title = "%s_%s_%s"%(cl,tf,treatment)
    data_folder = "/home/gaga/data-scratch/idannuri/ChipSeq/Step3/macs2_out"
    result_folder = "/home/gaga/data-scratch/idannuri/ChipSeq/Results"
    path = "/home/gaga/idannuri/Data/TF/induced"
    cfname = path + "/%s_control"%title
    ifname = path + "/%s_induced"%title
    cdic = {}
    if not Debug:
        TFBS = file("/home/gaga/idannuri/Data/TF/TFBS.bed").readlines()
        for line in TFBS:
            a =  set(line.split("\t")[6].split(",")).intersection(set(indexes))
            for c in a:
                if c not in cdic:
                    cdic[c] = []
                cdic[c].append(line)
        #writing induced BS and control BS to files
        control = file(cfname,"wb")
        control.write(''.join(cdic[indexes[0]]))
        control.close()
        induced = file(ifname,"wb")
        induced.write(''.join(cdic[indexes[1]]))
        induced.close()
    # comparing binding sites
    os.system("bedtools intersect -wa -wb -a %s -b %s -v > %s/%s_%s_%s_repressed_peaks.narrowPeak"\
              %(cfname, ifname, data_folder, cl,treatment, tf))
    import time
    os.system("bedtools intersect -wa -wb -a %s -b %s -v > %s/%s_%s_%s_induced_peaks.narrowPeak"\
              %(ifname, cfname ,data_folder, cl,tf, treatment))
    return generateFigure4a_Chipseq([tf,treatment],cl)



def convertFPKMtoTF(locus, name):
    ch = locus.split(":")[0]
    start, end = locus.split(":")[1].split("-")
    return TF("%s\t%s\t%s\t%s\n"%(ch,start,end,name))

def FPKMfile2dic(filename):
    """
    parsing FPKM file generated from cufflinks to a convenient dictionary
    """
    fpkm = file(filename).readlines()
    # we use x.split() for each object to remove white spaces from header
    headers = [x.split()[0] for x in fpkm[0].split("\t") if x.split()]
    fpkm_dic = {}
    gene_id = headers.index("gene_id")
    expression = headers.index("FPKM")
    #locus = headers.index("locus")
    #name = headers.index("gene_short_name")
    for line in fpkm[1:]:
        table = line.split("\t")
        idd = table[gene_id].split(".")[0]
        exp = max(1,float(table[expression]))
        #gname = table[name].split(".")[0]
        #gene = convertFPKMtoTF(table[locus], gname)
        fpkm_dic[idd] = exp
    return fpkm_dic


def generateFigure4a_Chiapet(expressedIds, cl,ret = 0):
    res = sum_list([dic_chia[cl][x] for x in expressedIds])
    if ret == 1:
        return res
    bg = sum(map(len, dic_chia[cl].values()))
    return chisq([len(expressedIds), len(res)],[len(ids), bg])

def Figure4a(exp2pval, name =""):
    """
    function recieves a dictionary with keys as the name of the experiment
    for example - "MCF7 ERa antibody TNFa" and the value is pval, enrichment factor
    """
    sortedpval = sorted(exp2pval.values())
    k = exp2pval.keys()
    sortedkeys = []
    for i in range(len(sortedpval)):
        for j in range(len(k)):
            if exp2pval[k[j]] == sortedpval[i]:
                sortedkeys.append(k[j])
                break
        a = k.pop(j)
    fig, ax = plt.subplots(1, 1, figsize = (20,10))
    #a = ax.plot([math.log(-1*math.log(exp2pval[x] + 10**(-300),10),10) for x in sortedkeys])
    #a = ax.plot([math.log(2,10) for i in range(len(sortedkeys))])
    a = ax.plot([exp2pval[x][0] for x in sortedkeys], c = 'b', linewidth = 5)
    a = ax.plot([2 for i in range(len(sortedkeys))], color = 'r')
    ax2 = ax.twinx()
    width = .4
    #enrichments = []
    for x in sortedkeys:
        en = exp2pval[x][1]
        #en = math.log(exp2pval[x][1],2)
        a = ax2.bar(sortedkeys.index(x)+ width/2, en,width=width, color = 'dimgray')
    a = ax.set_xticks(range(len(sortedkeys)))
    a = ax.set_ylabel("-log(pvalue)", fontsize = 30, color = 'b')
    ax2.set_ylabel("Enrichment Factor (log)",fontsize = 30, color = 'dimgray')
    """
    for induced TFBS
    """
    #a = ax.set_yticklabels([x for x in range(0,350,50)],fontsize = 25)
    """
    for induced genes
    """
    #ax2.set_yticks(range(-1,7))
    #ax.set_yticks(range(0,70,10))
    #a = ax.set_yticklabels([x for x in range(0,70,10)],fontsize = 25)
    #ax2.set_yticks(range(18))
    #a = ax2.set_yticklabels([x/4. for x in range(-4,4)],fontsize = 25)
    #ax.grid(True)
    ax2.grid(True)
    ax.set_xlabel("Experiments", fontsize = 30, color = "k")
    a = ax.set_xticklabels([x for x in sortedkeys], fontsize = 8, rotation = "vertical")
    #a = ax2.set_xticklabels(range(1, len(exp2pval)+1), fontsize = 10)
    #plt.xticks(rotation = "vertical")
    plt.tight_layout()
    plt.savefig("/home/gaga/idannuri/workdays/Figure4/revised_" + name)

def generateFigure4b_HiC(expressedIds1,cl1,expressedIds2,cl2, ret = 0):
    c11 = sum_list([map(Enhancer,enhs[cl1][tss][1:]) for tss in expressedIds1\
                       if tss in enhs[cl1]])
    c21 = sum_list([map(Enhancer, enhs[cl2][tss][1:]) for tss in expressedIds1 \
                    if tss in enhs[cl2]])
    c12 = sum_list([map(Enhancer, enhs[cl1][tss][1:]) for tss in expressedIds2 \
                    if tss in enhs[cl1]])
    c22 = sum_list([map(Enhancer, enhs[cl2][tss][1:]) for tss in expressedIds2 \
                    if tss in enhs[cl2]])
    if ret:
        return c11,c21,c12,c22
    return chisq([len(c11),len(c21)],[len(c12),len(c22)])


def generateFigure4b_Chipseq(induced1,induced2, path, title):
    """
    function compares chip seq of TF according to A\B compartments and return beatuiful table
    :param induced1: list of chip seq induced BS in cell line 1
    :param induced2: list of chip seq induced BS in cell line 2
    :return: saves result file in the given path
    """
    result = ""
    key1 = induced1[0].split("/")[-1].split("_")[0]
    key2 = induced2[0].split("/")[-1].split("_")[0]
    result += "%s %s\n"%(key1,key2)
    for i in range(len(induced1)):
        f1 = induced1[i]
        f2 = induced2[i]
        srr1 = '_'.join(f1.split("/")[-1].split("_")[1:3])
        srr2 = '_'.join(f2.split("/")[-1].split("_")[1:3])
        result += "%s %s\n"%(srr1,srr2)
        result += "\tAA\tAB\tBA\tBB\ttotal\t\tAA%\tAB%\tBA%\tBB%\t\tR\tEnrichment\t\tp-value\n"
        os.system("bedtools intersect -wa -wb -a %s -b %s -v > %s/cell1_only"\
                  % (f1,f2,path))
        os.system("bedtools intersect -wa -wb -a %s -b %s -v > %s/cell2_only" \
                  % (f2,f1,path))
        os.system("bedtools intersect -wa -wb -a %s -b %s > %s/cell_common" \
                  % (f1,f2,path))

        temp = [path+"/cell1_only",path+"/cell2_only",path+"/cell_common"]
        temp = [map(TF, file(f).readlines()) for f in temp]
        print "groups: ", map(len,temp)
        ABdiv = [TFABratio(t, key1, key2) for t in temp]
        print "ABdiv: ",ABdiv
        totals = [sum(abdiv) for abdiv in ABdiv]
        
        try:
            ABratio = [map(lambda x: round(x * 1.0 / totals[i], 2), ABdiv[i]) \
                   for i in range(len(ABdiv)) if totals[i]]
            ABratio += [[0,0,0,0] for i in range(3-len(ABratio))]
            R = [ABratio[0][1]/ABratio[0][2],ABratio[1][2]/ABratio[1][1]]
            R = map(lambda x:round(x,2),R)
        except:
            print "Error in calculating ABratio or R"
            os.system("rm %s/cell*" % (path))
        enrichment = (ABdiv[0][1]+ABdiv[1][2])*1.0/(ABdiv[0][2] + ABdiv[1][1])
        pval = chisq(ABdiv[0][1:3], ABdiv[1][1:3])

        result += "Cell1_only_BSs\t" + '\t'.join(map(str, ABdiv[0]))+"\t%s\t\t"%(totals[0]) + \
                  '\t'.join(map(str, ABratio[0])) + "\t\t%s\t%s\t\t%s\n"%(R[0],enrichment,pval)
        result += "Cell2_only_BSs\t" + '\t'.join(map(str, ABdiv[1])) + "\t%s\t\t" % (totals[1]) + \
                  '\t'.join(map(str, ABratio[1])) + "\t\t%s\n"%(R[1])
        result += "Common_BSs\t" + '\t'.join(map(str, ABdiv[2])) + "\t%s\t\t" % (totals[2]) + \
                  '\t'.join(map(str, ABratio[2])) + "\n"

    rfile = file(path+"/%s_4b_result.tsv"%title,'wb')
    rfile.write(result)
    rfile.close()
    os.system("rm %s/cell*" % (path))
    return 0

def generateFigure5a_FPKM(pair, cl, fiveB = False, bg = 1):
    background = load("/home/gaga/idannuri/Data/gene_expression/rna_seq/dic_5a_background")
    template = "/home/gaga/data-scratch/idannuri/ChipSeq/Step3/FPKM"
    rna = [template + "/%s.FPKM"%p for p in pair]
    control_dic = FPKMfile2dic(rna[0])
    treated_dic = FPKMfile2dic(rna[1])
    induced = []
    repressed = []
    count = 0
    #return control_dic,treated_dic
    print len(list(set(control_dic.keys()).intersection(treated_dic.keys())))
    for gene in list(set(control_dic.keys()).intersection(treated_dic.keys())):
        if gene not in transcripts:
            continue
        count += 1
        ratio = treated_dic[gene]*1.0/ control_dic[gene]
        #if  >= 1.5:
        if ratio >= 1.5 and control_dic[gene] <=2:
            induced.append(gene)
        if ratio <= 0.66:
            repressed.append(gene)
    print "compared genes: ", count
    if fiveB:
        return induced, repressed
    induced = [transcripts[g] for g in induced]
    repressed = [transcripts[g] for g in repressed]
    ind = check_regexp_AB(hun[cl][0],hun[cl][1], induced)
    rep = check_regexp_AB(hun[cl][0],hun[cl][1], repressed)
    lenA = sum([x.end-x.start for x in hun[cl][0]])
    lenB = sum([x.end-x.start for x in hun[cl][1]])
    if bg:
        lenA, lenB = background[cl]
        if bg == 2:
            basal = [transcripts[x] for x in control_dic if \
                           control_dic[x] > 1]
            lenA, lenB = check_regexp_AB(hun[cl][0],hun[cl][1],basal)
    print "background: %s %s"%(lenA,lenB)
    pval = [chisq([x+1 for x in ind],[lenA, lenB]), chisq([x+1 for x in rep],[lenA, lenB])]    
    pval = [-1*round(math.log(x+10**(-300),10),2) for x in pval]
    result = '\t'.join(map(str, [cl] + pair + list(ind) + [pval[0]] + list(rep) + [pval[1]])) + "\n"
    print result
    tot = [ind[0]+rep[0],ind[1]+rep[1]]
    print round(math.log(chisq(tot,[lenA,lenB])+10**(-300),10),2)
    return result

def generateFigure5a_RNAseq(filename,pair,cl, fiveB = False, bg = 1):
    background = load("/home/gaga/idannuri/Data/gene_expression/rna_seq/dic_5a_background_TSSs")
    rna = file(filename).readlines()
    print pair, cl
    # we use x.split() for each object to remove white spaces from header
    headers = [x.split()[0] for x in rna[0].split("\t") if x.split()]
    print headers
    try:
            #WARNING - this is for the case that the file is R template with blank space
            #as first column header, that's why we add 1
            index1 = headers.index(pair[0])+1
            index2 = headers.index(pair[1])+1
            print index1, index2
    except:
            print "invalid sample"
            return -1
    induced = []
    repressed = []
    count = 0
    for line in rna[1:]:
            table = line.split("\t")
            expression_control = max(1,float(table[index1]))
            expression_treated = max(1,float(table[index2]))
            ratio = expression_treated*1.0/ expression_control
            #if ratio >= 2 
            if ratio >= 2 and expression_control <=2:
                try:
                    induced.append(TF(line))
                except: # line is not in format for TF__init__, we will use ids
                    try:
                        idd = line.split("\t")[0]
                        induced.append(ids[idd])
                    except:
                        #print line
                        count += 1
            if ratio <= 0.5:
                try:
                    repressed.append(TF(line))
                except:  # line is not in format for TF__init__, we will use ids
                    try:
                        idd = line.split("\t")[0]
                        repressed.append(ids[idd])
                    except:
                        count += 1
    print "ERROR %s"%count
    if fiveB:
        return induced, repressed
    if cl not in hun:
        print "No Hi-C data for ",cl
        return ""
    ind = check_regexp_AB(hun[cl][0],hun[cl][1], induced)
    rep = check_regexp_AB(hun[cl][0],hun[cl][1], repressed)
    lenA = sum([x.end-x.start for x in hun[cl][0]])
    lenB = sum([x.end-x.start for x in hun[cl][1]])
    if bg:
        lenA, lenB = background[cl]
    print "background: %s %s"%(lenA,lenB)
    pval = [chisq([x+1 for x in ind],[lenA, lenB]), chisq([x+1 for x in rep],[lenA, lenB])]
    pval = [-1*round(math.log(x+10**(-300),10),2) for x in pval]
    Rs = [((ind[0] + 1)*1.0/(ind[1]+1))/(lenA*1.0/lenB), ((rep[0] + 1)*1.0/(rep[1]+1))/(lenA*1.0/lenB)]
    Rs = [round(x,2) for x in Rs]
    result = '\t'.join(map(str, [cl] + pair + list(ind) + [pval[0],Rs[0]] +  list(rep) + [pval[1],Rs[1]])) + "\n"
    print result
    return result

def mult_generateFigure5a_RNAseq(pairs_file, filename="", FPKM = False, rf = None):
    """
    if not FPKM:
        
        #in this case pairs is the ID of the process ran through pipeline
        #and we can use this ID to extract all parameters
        
        tmp = pairs_file
        pairs = "/home/gaga/data-scratch/idannuri/ChipSeq/pipeline_input/%s.txt"%tmp
        filename = "/home/gaga/data-scratch/idannuri/ChipSeq/Step3/FPKM/%s_rpkms_data_qnorm.20q.floor.txt"%tmp
        #return filename
        rf = filename.split("/")[-1].split("_")[0]
    else:
        rf = pairs_file.split("/")[-1].split("_")[0]
    """
    tmp = file(pairs_file).readlines()
    #pairs = [x.split("\t") for x in tmp]
    #pairs = [[x.split()[0] for x in p] for p in pairs]
    pairs = [p.split() for p in tmp]
    header = '\t'.join(["Cell line", "Control", "Treatment", "induced in A", "induced in B","log p-value", "R",\
                        "repressed in A", "repressed in B", "log p-value","R"]) + "\n"
    rows = ""
    for pair in pairs:
        p0 = pair[0]
        p1 = pair[1]
        cl = pair[2]
        if not FPKM:
            print filename, [p0,p1], cl
            rows += generateFigure5a_RNAseq(filename, [p0,p1], cl)
        else:
            rows += generateFigure5a_FPKM([p0,p1], cl)
    
    output = file("/home/gaga/data-scratch/idannuri/ChipSeq/Results/%s_5a_results_new.txt"%(rf), "wb")
    output.write(header + rows)
    output.close()
    return 0

def generateFigure5b_RNAseq(title,pair1, cl1,pair2, cl2,filename1 = "", filename2 = "", FPKM  =False, fiveC = False):
    """
    usage example:
    generateFigure5b_RNAseq("LNCAP_MCF7_TNFa",["SRR3728890","SRR3728892"],"LNCAP",...")
    """
    path = "/home/gaga/data-scratch/idannuri/ChipSeq/Results"
    result = ""
    key1 = cl1
    key2 = cl2
    result += "%s %s\n"%(key1,key2)

    srr1 = "%s_%s"%(pair1[0],pair1[1])
    srr2 = "%s_%s"%(pair2[0],pair2[1])
    result += "%s %s\n"%(srr1,srr2)
    result += "\tAA\tAB\tBA\tBB\ttotal\t\tAA%\tAB%\tBA%\tBB%\t\tR\tEnrichment\t\tp-value\n"
    #result += "\tAA\tAB\tBA\tBB\ttotal\t\tR\tEnrichment\t\tp-value\n"
    if not FPKM:
        if filename1.find("gene_expression") > 0:
            ind1 = [ids[x] for x in load(filename1)]
        else:
            ind1,rep1 = generateFigure5a_RNAseq(filename1,pair1,cl1, fiveB = True)
        if filename2.find("gene_expression") > 0:
            ind2 = [ids[x] for x in load(filename2)]
        else:
            ind2,rep2 = generateFigure5a_RNAseq(filename2,pair2,cl2, fiveB = True)
    else:
        ind1,rep1 = generateFigure5a_FPKM(pair1,cl1, fiveB = True)
        ind2,rep2 = generateFigure5a_FPKM(pair2,cl2, fiveB = True)   
    ind1_only = set(ind1).difference(set(ind2))
    #ind1_only = [transcripts[x] for x in ind1_only]
    #ind1_only = ind1
    ind2_only = set(ind2).difference(set(ind1))
    #ind2_only = [transcripts[x] for x in ind2_only]
    #ind2_only = ind2
    common = set(ind1).intersection(set(ind2))
    #common = [transcripts[x] for x in common]
    #common = ind1+ind2
    if fiveC:
        return ind1, ind2, common
    temp = [ind1_only, ind2_only, common]
    ABdiv = [TFABratio(t, key1, key2) for t in temp]
    totals = [sum(abdiv) for abdiv in ABdiv]
    try:
        ABratio = [map(lambda x: round(x * 1.0 / totals[i], 2), ABdiv[i]) \
               for i in range(len(ABdiv)) if totals[i]]
        ABratio += [[0,0,0,0] for i in range(3-len(ABratio))]
        R = [ABratio[0][1]/ABratio[0][2],ABratio[1][2]/ABratio[1][1]]
        R = map(lambda x:round(x,2),R)
    except:
        print "Error in calculating ABratio or R"

    enrichment = (ABdiv[0][1]+ABdiv[1][2])*1.0/(ABdiv[0][2] + ABdiv[1][1])
    pval = chisq(ABdiv[0][1:3], ABdiv[1][1:3])


    result += "Cell1_only_BSs\t" + '\t'.join(map(str, ABdiv[0]))+"\t%s\t\t"%(totals[0]) + \
              '\t'.join(map(str, ABratio[0])) + "\t\t%s\t%s\t\t%s\n"%(R[0],enrichment,pval)
    result += "Cell2_only_BSs\t" + '\t'.join(map(str, ABdiv[1])) + "\t%s\t\t" % (totals[1]) + \
              '\t'.join(map(str, ABratio[1])) + "\t\t%s\n"%(R[1])
    result += "Common_BSs\t" + '\t'.join(map(str, ABdiv[2])) + "\t%s\t\t" % (totals[2]) + \
              '\t'.join(map(str, ABratio[2])) + "\n"
    """
    return ABdiv, ABratio, totals, R
    result += "Cell2_only_BSs\t" + '\t'.join(["%s (%s%%)"%(ABdiv[1][i],ABratio[1][i]) for i in range(len(ABdiv[0]))]) + \
              "\t%s\t\t" % (totals[1])+ "\t\t%s\n"%(R[1])
    result += "Common_BSs\t" + '\t'.join(["%s (%s%%)"%(ABdiv[1][i],ABratio[1][i]) for i in range(len(ABdiv[0]))]) + \
              "\t%s\t\t" % (totals[2]) + '\t'.join(map(str, ABratio[2])) + "\n"
    """
    rfile = file(path+"/%s_5b_result.tsv"%title,'wb')
    rfile.write(result)
    rfile.close()
    os.system("rm %s/cell*" % (path))
    return 0

def generateFigure5c(data,p1,p2, FPKM  =False, ret = 0):
    ind1,ind2,com = generateFigure5b_RNAseq("", p1.pair,p1.cell_line,p2.pair,p2.cell_line,p1.sample_table, p2.sample_table, FPKM, fiveC= True)
    if data == "chia":
        dd = dic_chia
    elif data == "chia_full":
        dd = dic_chia_full
    elif data == "hic":
        dd = enhs
    elif data == "norm":
        dd = load("/home/gaga/idannuri/Data/Chiapet/dic_chia_norm.pkl")
    else:
        print "Data type can only be chia\chia_full\hic"
        return -1
    if data != "norm" and data!="hic":
        for cl in dd:
            for k in dd[cl]:
                dd[cl][k] = len(dd[cl][k])
    ind1 = [x.ID for x in ind1]
    ind2 = [x.ID for x in ind2]
    inds1 = [x for x in ind1 if x in dd[p1.cell_line] and x in x in dd[p2.cell_line]]
    inds2 = [x for x in ind2 if x in dd[p1.cell_line] and x in x in dd[p2.cell_line]]
    cp11 = [min(dd[p1.cell_line][x],20) for x in inds1]
    cp12 = [min(dd[p1.cell_line][x],20) for x in inds2]
    cp21 = [min(dd[p2.cell_line][x],20) for x in inds1]
    cp22 = [min(dd[p2.cell_line][x],20) for x in inds2]
    print len(inds1),len(inds2)
    print "only CL1 genes - average_cl1:%s, average_cl2:%s"\
          %(sum(cp11)*1.0/len(cp11), sum(cp21)*1.0/len(cp21))
    print "only CL2 genes - average_cl1:%s, average_cl2:%s"\
          %(sum(cp12)*1.0/len(cp12), sum(cp22)*1.0/len(cp22))
    if ret:
        return [cp11,cp21,cp12,cp22]
    cpps = [plt.hist(c,20)[0] for c in [cp11,cp21,cp12,cp22]]
    if ret:
        return cpps


"""
For version 2 of the permutation test where we permute genes with the same expresion
distribution as the original list
"""
def generate_same_expression_random_gene_list(cl, gene_list):

    import random
    exp_list = [expressions[cl][g] for g in gene_list if g in expressions[cl] and expressions[cl][g]]
    mif = [map(int,exp_list).count(i) for i  in range(100)]
    # adding all > 100
    mif.append(len([x for x in exp_list if x >= 100]))
    new = sum_list([random.sample(expressions_by_range[cl][i],mif[i]) for i in range(len(mif))])
    return new



def permutation_test(data, cl, gene_list, it_num = 10*1000, hist = False, v2 = False):
    """
    A very bad name for Figure 7 and Table 4 from the paper
    """
    if data == "chia":
        dd = dic_chia
    elif data == "chia_full":
        dd = dic_chia_full
    elif data == "hic":
        dd = enhs
    elif data == "norm":
        dd = load("/home/gaga/idannuri/Data/Chiapet/dic_chia_norm.pkl")
    else:
        print "Data type can only be chia\chia_full\hic"
        return -1
    
    dd = dd[cl]
    if type(dd.values()[0]) != type(1):
        for k in dd:
            dd[k] = len(dd[k])
            
    print len(gene_list),
    gene_list = [g.ID for g in gene_list if g.ID in dd]
    exp_list = [g for g in gene_list if g in expressions[cl] and expressions[cl][g]]
    print len(gene_list)
    test_sample = [dd[gene] for gene in gene_list]
    test_median = sorted(test_sample)[len(test_sample)/2]
    test_mean = sum(test_sample)*1.0/len(test_sample)
    test_exp = sum([expressions[cl][gene] for gene in exp_list])/len(exp_list)
    print test_median, test_mean,test_exp
    counter_median = 0
    counter_mean = 0
    counter_exps = 0
    vals = dd.values()
    exps = [val for val in expressions[cl].values() if val]
    #if hist:
    means = [test_mean]
    for iteration in range(it_num):
        if not v2:
            perm = numpy.random.permutation(len(dd))[:len(gene_list)]
            sample = [vals[i] for i in perm]
        else:
            perm = generate_same_expression_random_gene_list(cl, gene_list)
            sample = [dd[gene] for gene in perm if gene in dd]
        median = sorted(sample)[len(sample)/2]
        mean = sum(sample)*1.0/len(sample)
        #if hist:
        means.append(mean)
        #exp = sum([exps[i%len(exps)] for i in perm])/len(perm)
        if median >= test_median:
            counter_median += 1
        if mean >= test_mean:
            counter_mean += 1
        #if exp >= test_exp:
        #    counter_exps += 1

        if iteration %1000 == 0:
            print iteration, mean, median

    pval_median = counter_median*1.0/it_num
    pval_mean = counter_mean*1.0/it_num
    #pval_exp = counter_exps*1.0/it_num
    print "mean p-value %s, median p-value %s"\
          %(pval_mean, pval_median)
    if hist:
        return means
    print test_mean,sum(means)*1.0/it_num,pval_mean, test_median, pval_median
    return test_mean,sum(means)*1.0/it_num,pval_mean, test_median, pval_median

def generateFigure6(data, pair,compare = None):
    if compare == None:
        compare = pair.cell_line
    if pair.sample_table.find("gene_expression") > 0:
        ind = [ids[x] for x in load(pair.sample_table)]
        return permutation_test(data,compare, ind)
    else:
        ind,rep = generateFigure5a_RNAseq(pair.sample_table,pair.pair,pair.cell_line, fiveB = True)
    return permutation_test(data,pair.cell_line, ind)#,permutation_test(data,pair.cell_line, rep)
    
    
def mult_generateFigure6(pairs):
    result = "Cell line\tdescription\tdata type\tsamples\tmean\trandom mean\tmean p-val\tmedian\tmedian p-val\n"
    for cl in pairs:
        for k in pairs[cl]:
            if k.startswith("RNA"):
                try:
                    if cl in dic_chia:
                        res = map(lambda x: round(x,2),generateFigure6("chia",pairs[cl][k]))
                        result += "%s\t%s\t%s\t%s_%s\t%s\t%s\t%s\t%s\t%s\n"%\
                                  (cl, pairs[cl][k].desc,"ChIA-PET",pairs[cl][k].pair[0],pairs[cl][k].pair[1],\
                                   res[0],res[1],res[2],res[3],res[4])
                    if cl in enhs:
                        res = map(lambda x: round(x,2),generateFigure6("hic",pairs[cl][k]))
                        result += "%s\t%s\t%s\t%s_%s\t%s\t%s\t%s\t%s\t%s\n"%\
                                  (cl, pairs[cl][k].desc,"Hi-C",pairs[cl][k].pair[0],pairs[cl][k].pair[1],\
                                   res[0],res[1],res[2],res[3],res[4])
                except:
                    print "ERROR in %s %s"%(cl,k)
    
    
    output = file("/home/gaga/data-scratch/idannuri/ChipSeq/Results/Figure6.tsv",'wb')
    output.write(result)
    output.close()
    return 0

"""
Figure 7 in final papaer
"""
def Figure6(data,pair):
    if pair.sample_table.find("gene_expression") > 0:
        ind = [ids[x] for x in load(pair.sample_table)]
        means =  permutation_test(data,compare, ind, hist = True)
    else:
        ind,rep = generateFigure5a_RNAseq(pair.sample_table,pair.pair,pair.cell_line, fiveB = True)
        means = permutation_test(data,pair.cell_line, ind, hist = True)
    bins = numpy.linspace(min(means),max(means),100)
    plt.hist(means,bins, alpha = 0.5)
    #plt.title("")
    plt.xlabel("TSS interactions measured by ChIA-PET", fontsize = 30)
    plt.ylabel("count", fontsize = 30)
    plt.show()
    
    
def FigurePaper1(cl = "", dos = None):
    if not dos:
        dos = domains[cl]
        print "not dos"
    ABstart = [sorted(sum_list([[x.start] for x in dos if x.chromosome == i])) for i in range(1,23)]
    diss = [binarysearch(ABstart[d.chromosome-1], d.start) for d in hun[cl][0] if d.chromosome < 23]
    """
    mif = [0 for i in range(0,9000000,100000)]
    for x in diss:
	mif[x/100000] += 1
    showFunction(mif)
    """
    #print "Printing median and mean"
    return sorted(diss)[len(diss)/2], sum(diss)*1.0/len(diss)

def FullFigurePaper1(cl, samples = 10000):
    print "The real mean is ", FigurePaper1(cl)[1]
    try:
        means = load("/home/gaga/idannuri/workdays/Figure1Paper/%s_means_%s"%(cl,samples))
        medians = load("/home/gaga/idannuri/workdays/Figure1Paper/%s_medians_%s"%(cl,samples))
    except:    
        reses = [FigurePaper1(cl, permute_domain(cl)) for i in range(samples)]
        means = [res[1] for res in reses]
        medians = [res[0] for res in reses]
        dump(means, "/home/gaga/idannuri/workdays/Figure1Paper/%s_means_%s"%(cl,samples))
        dump(medians, "/home/gaga/idannuri/workdays/Figure1Paper/%s_medians_%s"%(cl,samples))
    bins = numpy.linspace(min(means),max(means),max(20,samples/100))
    a = plt.hist(means,bins, alpha = 0.5)
    plt.xlabel("Mean distance between compartment border and TAD border in %s"%cl, fontsize = 30)
    plt.ylabel("count", fontsize = 30)
    plt.show()


def FigurePaper2(filename,pair,cl,plot = 0, rand = 0, relative = True):
    if cl not in domain2TSS:
        #print "not available TADs for %s"%cl
        return []

    if not rand:
        background = load("/home/gaga/idannuri/Data/gene_expression/rna_seq/dic_5a_background_TSSs")
        rna = file(filename).readlines()
        # we use x.split() for each object to remove white spaces from header
        headers = [x.split()[0] for x in rna[0].split("\t") if x.split()]
        try:
                #WARNING - this is for the case that the file is R template with blank space
                #as first column header, that's why we add 1
                index1 = headers.index(pair[0])+1
                index2 = headers.index(pair[1])+1
        except:
                print "invalid sample"
                return -1
        induced = []
        repressed = []
        count = 0
        for line in rna[1:]:
                table = line.split("\t")
                expression_control = max(1,float(table[index1]))
                expression_treated = max(1,float(table[index2]))
                ratio = expression_treated*1.0/ expression_control
#                if ratio > 1.5 and expression_control <=2:
                if ratio > 1.5:
                    try:
                        induced.append(TF(line))
                    except: # line is not in format for TF__init__, we will use ids
                        try:
                            idd = line.split("\t")[0]
                            induced.append(ids[idd])
                        except:
                            count += 1
                if ratio < (1/1.5):
                    try:
                        repressed.append(TF(line))
                    except:  # line is not in format for TF__init__, we will use ids
                        try:
                            idd = line.split("\t")[0]
                            repressed.append(ids[idd])
                        except:
                            count += 1
        #print "Induced and repressed equals to ", len(induced), len(repressed)
    if rand:
        induced = rand[0]; repressed = rand[1]
        allt = induced + repressed
        perm = random_permute(0,len(allt) - 1)
        allt = [allt[i] for i in perm]
        induced = allt[:len(induced)]
        repressed = allt[len(induced):]

    ind_dic = {}
    for x in induced:
        ind_dic[x.ID] = 1
    rep_dic = {}
    for x in repressed:
        rep_dic[x.ID] = 1

    domain_scores = []
    for domain in domain2TSS[cl]:
        if not domain2TSS[cl][domain] or len(domain2TSS[cl][domain]) < 2:
            continue
        score = 0
        count = 0
        for t in domain2TSS[cl][domain]:
            if t.ID in ind_dic:
                score += 1
                count += 1
            if t.ID in rep_dic:
                score -= 1
                count += 1
        if count > 2:
            if relative:
                score = score*1.0/count
            domain_scores.append(score)
    return domain_scores, induced, repressed

def FullFigurePaper2(pairs, filename="", FPKM = False, rand = 0):
    if not rand:
        rand = [0 for i in range(len(pairs))]
    if not FPKM:
        """
        in this case pairs is the ID of the process ran through pipeline
        and we can use this ID to extract all parameters
        """
        tmp = pairs
        pairs = "/home/gaga/data-scratch/idannuri/ChipSeq/pipeline_input/%s.txt"%tmp
        filename = "/home/gaga/data-scratch/idannuri/ChipSeq/Step3/FPKM/%s_rpkms_data_qnorm.20q.floor.txt"%tmp
    tmp = file(pairs).readlines()
    pairs = [x.split("\t") for x in tmp]
    pairs = [[x.split()[0] for x in p] for p in pairs]
    res = []
    #for pair in pairs:
    for i in range(len(pairs)):
        if not rand[0]:
            print pairs[i]
        p0 = pairs[i][0]
        p1 = pairs[i][1]
        cl = pairs[i][2]
        res.append(FigurePaper2(filename, [p0,p1], cl, rand = rand[i]))
    return res

def AnalysisPaper2(title, rand_it = 100, rang = 10):
#ratios = AnalysisPaper2("PRJNA184350_RNA",rand_it = 1000, )

    res = FullFigurePaper2(title, rand=0)
    
    rand = []
    for i in range(len(res)):
        if res[i]:
            rand.append([res[i][1],res[i][2]])
        else:
            rand.append(0)
            
    rand_res = [FullFigurePaper2(title, rand=rand) for i in range(rand_it)]

    #return real_res, rand_res
    ratios = []
    for experiment in range(len(res)):
        print 'exp num ' ,experiment, len(res)
        if not res[experiment]:
            print "Empty exp: ",experiment
            continue
    
        real_res = res[experiment][0]
        print len(real_res)
        real_mif = [0 for i in range(2*rang+1)]
        for s in real_res:
            real_mif[int((s+1)*rang)] +=1

        rand_mif = [0 for i in range(2*rang+1)]
        for r in rand_res:
            #print rand_mif
            for s in r[experiment][0]:
                rand_mif[int((s+1)*rang)] += 1
            
        for cell in range(len(rand_mif)):
            rand_mif[cell] = rand_mif[cell]*1.0/rand_it

        ratio = [(real_mif[i]*1.0+1)/(rand_mif[i]+1) for i in range(len(real_mif))]

        ratios.append(ratio)
        
    return ratios

def ds(res, rang):
    mif = [0 for i in range(2*rang+1)]
    for s in res:
        mif[int((s+1)*rang)] +=1
    showFunction(mif)
    return 0
    

def FigurePaper3(cl1,cl2, rand_it = 100):

    CL1up = []
    CL2up = []
    for gene in expressions[cl1]:
	if gene not in expressions[cl2]:
		continue
	if max(1,expressions[cl1][gene])*1.0/max(1,expressions[cl2][gene])*1.0 > 2:
		CL1up.append(gene)
	if max(1,expressions[cl2][gene])*1.0/max(1,expressions[cl1][gene])*1.0 > 2:
		CL2up.append(gene)

    cl = "GM12878"

    domain_scores = []
    good = []
    for domain in domain2TSS[cl]:
        if not domain2TSS[cl][domain] or len(domain2TSS[cl][domain]) < 2:
            continue
        score = 0
        count = 0
        for t in domain2TSS[cl][domain]:
            if t.ID in CL1up:
                score += 1
                count += 1
            if t.ID in CL2up:
                score -= 1
                count += 1
        if count > 2:
            score = score*1.0/count
            if abs(score) == 1 and count > 4:
                good.append([domain,score])
            domain_scores.append(score)

    rand_domain_scores = []
    for it in range(rand_it):
        if it%(rand_it/5) == 0:
            print it,
        allt = CL1up + CL2up
        perm = random_permute(0,len(allt) - 1)
        allt = [allt[i] for i in perm]
        rand_CL1up = allt[:len(CL1up)]
        rand_CL2up = allt[len(CL1up):]
        
        for domain in domain2TSS[cl]:
            if not domain2TSS[cl][domain] or len(domain2TSS[cl][domain]) < 2:
                continue

        #print len(rand_CL2up), len(rand_CL1up)      
            score = 0
            count = 0
            for t in domain2TSS[cl][domain]:
                if t.ID in rand_CL1up:
                    score += 1
                    count += 1
                if t.ID in rand_CL2up:
                    score -= 1
                    count += 1
            if count > 2:
                
                score = score*1.0/count
                rand_domain_scores.append(score)

    re = plt.hist(domain_scores,10)
    rand = plt.hist(rand_domain_scores,10)
    vals = (re[0]+1)/((rand[0]+1)/rand_it)
    fig,ax = plt.subplots(1,1)
    ax.plot(list(vals))
    ax.set_xticks(np.arange(len(vals)+1), minor=False)
    ax.set_xticklabels(re[1], minor=False)
    ax.set_xlabel("Coordinate Response Score",fontsize = 18)
    ax.set_ylabel("Ratio between real genes\nand permuted genes",fontsize = 12)
    ax.set_title("%s %s"%(cl1, cl2), fontsize = 30)
    plt.savefig("/home/gaga/idannuri/workdays/Coordinate%s_%s_%s"%(cl1,cl2,rand_it))
    return 0
    
def FigurePaper3v2(cl1,cl2, rand_it = 100):
    """
    here instead of randomizing the genes we randomize the domain borders
    """
    TSSList = [[]] + [[x for x in TSSs if x.chromosome == i] for i in range(1,24)]
    CL1up = []
    CL2up = []
    for gene in expressions[cl1]:
	if gene not in expressions[cl2]:
		continue
	if max(1,expressions[cl1][gene])*1.0/max(1,expressions[cl2][gene])*1.0 > 2:
		CL1up.append(gene)
	if max(1,expressions[cl2][gene])*1.0/max(1,expressions[cl1][gene])*1.0 > 2:
		CL2up.append(gene)

    cl = "GM12878"
    real = inter_intra_domains[cl][0]

    """
    generating domain2TSS dictionary for given domains
    """
    d2TSS = {}
    for domain in real:
        d2TSS[str(domain)] = []
        for t in TSSList[domain.chromosome]:
                if t.TSS >= domain.start and t.TSS <= domain.end:
                    d2TSS[str(domain)].append(t)
        
    domain_scores = []
    for domain in d2TSS:
        if not d2TSS[domain] or len(d2TSS[domain]) < 2:
            continue
        score = 0
        count = 0
        for t in d2TSS[domain]:
            if t.ID in CL1up:
                score += 1
                count += 1
            if t.ID in CL2up:
                score -= 1
                count += 1
        if count > 2:
            score = round(score*1.0/count,1)
            domain_scores.append(score)

    rand_domain_scores = []
    for it in range(rand_it):
        if it%(rand_it/5) == 0:
            print it,

        """
        same procedure now with randomized domains
        """
        dos = permute_domain(cl)
        d2TSS = {}
        for domain in dos:
            d2TSS[str(domain)] = []
            for t in TSSList[domain.chromosome]:
                if t.TSS >= domain.start and t.TSS <= domain.end:
                    d2TSS[str(domain)].append(t)

        for domain in d2TSS:
            if not d2TSS[domain] or len(d2TSS[domain]) < 2:
                continue
            score = 0
            count = 0
            for t in d2TSS[domain]:
                if t.ID in CL1up:
                    score += 1
                    count += 1
                if t.ID in CL2up:
                    score -= 1
                    count += 1
            if count > 2:
                score = round(score*1.0/count,1)
                rand_domain_scores.append(score)
        
    re = plt.hist(domain_scores,11)
    rand = plt.hist(rand_domain_scores,11)
    return re, rand
    vals = (re[0]+1)/((rand[0]+1)/rand_it)
    fig,ax = plt.subplots(1,1)
    ax.plot(list(vals))
    ax.set_xticks(np.arange(len(vals)+1), minor=False)
    ax.set_xticklabels(re[1], minor=False)
    ax.set_xlabel("Coordinate Response Score",fontsize = 18)
    ax.set_ylabel("Ratio between real genes\nand permuted genes",fontsize = 12)
    ax.set_title("%s %s"%(cl1, cl2), fontsize = 30)
    plt.savefig("/home/gaga/idannuri/workdays/Coordinate%s_%s_%s_v2"%(cl1,cl2,rand_it))
    return 0

def autoPCAExpCorrelation(cl,bins = [-3,-2,-1,0,1,2,3],rv = False):
    PCValsNorm = load("/home/gaga/idannuri/Data/PCValsNorm.dic")
    pc2exp = {}
    for gene in TSSs:
        if gene.ID in expressions[cl] and gene.chromosome in PCValsNorm[cl]:
            exp = math.log(max(expressions[cl][gene.ID],1),2)
            pc = PCValsNorm[cl][gene.chromosome][gene.TSS/10**5]
            if pc in pc2exp:
                pc2exp[pc].append(exp)
            else:
                pc2exp[pc] = [exp]

    if rv:
        return pc2expsc
    sortpc = sorted(pc2exp.keys())

    expBins = []
    for i in range(10):
        expBin = []
        for val in sortpc[len(sortpc)*i/10:len(sortpc)*(i+1)/10]:
            expBin += pc2exp[val]
        expBins.append(expBin)

    fig,ax = plt.subplots(1,1)
    ax.boxplot(expBins,1,patch_artist=True)

    avgPC = [sum(sortpc[len(sortpc)*i/10:len(sortpc)*(i+1)/10])/(len(sortpc)*0.1)\
             for i in range(10)]
    ax.set_xticklabels(map(lambda x: round(x,1),avgPC), fontsize = 20)

    if rv:
        return expBins, avgPC

    plt.savefig("/home/gaga/idannuri/workdays/AutoPCAExpCorr%s"%(cl))
            

def PCAExpressionCorrelation(cl1, cl2, bins = [-2,-1,-0.5,0,0.5,1,2], rv = 0):
    PCValsNorm = load("/home/gaga/idannuri/Data/PCValsNorm.dic")

    #norm Expressions from Rani
    #expressions = load("/home/gaga/idannuri/Data/gene_expression/rna_seq/dic_norm_expressions")
    
    """
    function gets two cell lines and plots gene expression delta as a function
    of PC delta
    """
    pc2exp = {}

    #This is the version that runs on all TSSs, we also implement
    #a version that runs only on AB\BA genes
    
    for gene in TSSs:
        if not(gene.ID in expressions[cl1] and gene.ID in expressions[cl2]):
            continue
        if gene.chromosome not in range(1,23) + ["X"]:
            continue
        
        exp1 = max(expressions[cl1][gene.ID],1)
        exp2 = max(expressions[cl2][gene.ID],1)

        if exp1 ==1 and exp2 == 1:
            continue
        
        PCVal1 = PCValsNorm[cl1][gene.chromosome][gene.TSS/10**5]
        PCVal2 = PCValsNorm[cl2][gene.chromosome][gene.TSS/10**5]

        DelPC = PCVal1 - PCVal2
     
        expRatio = math.log((exp1/exp2),2)
        #expRatio = numpy.sign(expRatio)*min(abs(expRatio),5)
        if DelPC in pc2exp:
            pc2exp[DelPC].append(expRatio)
        else:
            pc2exp[DelPC] = [expRatio]

    if rv:
        return pc2exp
    
    sortpc = sorted(pc2exp.keys())

    expBins = []
    if bins:
        avgPC = []
        bin_size = []
        for i in range(len(bins)-1):
            pcs = [x for x in sortpc if x>=bins[i] and x<bins[i+1]]
            expBin = sum_list([pc2exp[x] for x in pcs])
            expBins.append(expBin)
            
            avgPC.append(sum(pcs)/len(pcs))
            bin_size.append(len(pcs))
        
    else:    
        for i in range(10):
            expBin = []
            for val in sortpc[len(sortpc)*i/10:len(sortpc)*(i+1)/10]:
                expBin += pc2exp[val]
            expBins.append(expBin)
            
        avgPC = [sum(sortpc[len(sortpc)*i/10:len(sortpc)*(i+1)/10])/(len(sortpc)*0.1)\
             for i in range(10)]


    fig,ax = plt.subplots(1,1)
    ax.boxplot(expBins,1,patch_artist=True)

    avgPC = map(lambda x: round(x,1),avgPC)
    ax.set_xticklabels(["%s\n%s"%(avgPC[i],bin_size[i]) for i in range(len(avgPC))]\
                        , fontsize = 14)

    plt.savefig("/home/gaga/idannuri/workdays/PCAExpressionFigures/PCAExpCorr%s_%s_new"%(cl1, cl2))
    return 0

        
def PCValDist(cl, filename):
    PCValsNorm = load("/home/gaga/idannuri/Data/PCValsNorm.dic")
    locs = file(filename).readlines()
    vals = []
    for loc in locs:
        tf = TF(loc)
        if tf.chromosome == "chrX":
            tf.chromosome = 'X'
        if tf.chromosome not in range(1,23) +["X"]:
            continue
        try:
            vals.append(PCValsNorm[cl][tf.chromosome][int(((tf.end+tf.start)/2.)/10**5)])
        except:
            return cl, tf.chromosome, (tf.end+tf.start)/2.
    return vals
    

def PCdist(cl, family):
    hm = glob.glob("/home/gaga/idannuri/Data/%s/all/*/%s"%(family,cl))
    names = [x.split("/")[-2] for x in hm]
    pcs = [PCValDist(cl, h) for h in hm]
    fig,ax = plt.subplots(1,1, figsize = (20,10))
    a = ax.boxplot(pcs)
    ax.set_xticklabels(names,fontsize = 12, rotation = "vertical")
    plt.savefig("/home/gaga/idannuri/workdays/PCDistLocs_%s_%s"%(family,cl))


def polyVSnonpoly(cl):

    os.system("bedtools intersect -wa -wb -a /home/gaga/idannuri/Data/TF/all/EZH2/%s\
                -b /home/gaga/idannuri/Data/Histone/all/H3k27me3/%s > \
                /home/gaga/idannuri/workdays/Figure7-Polycomb/%s"%(cl,cl,cl))
    os.system("bedtools intersect -wa -wb -a /home/gaga/idannuri/Data/Histone/all/H3k27me3/%s\
                -b /home/gaga/idannuri/Data/TF/all/EZH2/%s -v > \
                /home/gaga/idannuri/workdays/Figure7-Polycomb/%s_noEZH"%(cl,cl,cl))
    poly = PCValDist(cl, "/home/gaga/idannuri/workdays/Figure7-Polycomb/%s"%cl)
    nonpoly = PCValDist(cl, "/home/gaga/idannuri/workdays/Figure7-Polycomb/%s_noEZH"%cl)

    fig,ax = plt.subplots(1,1)
    ax.boxplot([nonpoly,poly])
    ax.set_xticklabels(["All H3K27ME3","H3K27ME3 and EZH2"]\
                        ,fontsize = 12, rotation = "vertical")
    fig.tight_layout()
    plt.savefig("/home/gaga/idannuri/workdays/Figure7-Polycomb/%s_new"%cl)
    

def generate_exp2pval(diff_genes):
    induced = {}
    repressed = {}
    for exp in diff_genes:
        cl = exp.split("_")[0].upper()
        base = check_regexp_AB(hun[cl][0],hun[cl][1], TSSs)
        ind = check_regexp_AB(hun[cl][0],hun[cl][1], diff_genes[exp]["induced"])
        rep = check_regexp_AB(hun[cl][0],hun[cl][1], diff_genes[exp]["repressed"])
        induced[exp] = [-1*math.log(chisq(ind,base),10), \
                        (ind[0]*1.0/ind[1])/(base[0]*1.0/base[1])]
        repressed[exp] = [-1*math.log(chisq(rep,base),10), \
                          (rep[0]*1.0/rep[1])/(base[0]*1.0/base[1])]
    return induced, repressed
        
    
