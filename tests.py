try:
    execfile("/home/rack-shamir2/home/idannuri/Python/permute_domain_test.py")
except:
    print "ERROR in permute_domains_test.py"

def middle(gene):
    return (gene.start + gene.end)/2

def get_matrix_cell_for_pair(g1,g2, resolution=1000,expected = True, res_dir = "1k", verbose = 0):
    try:
        if expected:
            matrix = file("/home/rack-shamir2/home/idannuri/intra"+ res_dir + "/expected/chr"+str(g1.chromosome)+"_"+res_dir + "b.normalized")
        else:

            matrix = file("/home/rack-shamir2/home/idannuri/intra"+ res_dir +"/not_expected/chr"+str(g1.chromosome)+"_"+res_dir +"b.normalized_wo_expected")
    except:
        matrix = file(r"R:\home\idannuri\intra" + res_dir +"\chr"+str(g1.chromosome)+"_"+res_dir +"b.normalized")
    r = resolution
    indice1 = r * int(middle(g1) / r)
    indice2 = r * int(middle(g2) / r)
    search = "\n%s %s "%(str(indice1),str(indice2))
    chunk = matrix.read(1024*1024)
    i = 0
    while chunk:
        if verbose:
            if i%100 == 0:
                print i,
            i+=1
        try:
            index = chunk.index(search)
            value = float(chunk[index:index+1000].splitlines()[1].split()[-1])
            return value
        except:
            chunk = matrix.read(1024*1024)
    matrix.close()
    del chunk
    return 0


""" a lot FASTER

def pairs_dic_find_cofunc_pairs_inDomain(pairs_dic,domains, number = False):
    good_pairs = []
    bad_pairs = []
    doms = []
    score = 0
    for p in pairs_dic.keys():
        midi = middle(p[0])
        domains0 = domains0 = [d for d in domains if InDomain(d,p[0])]
        if p[0].MCLID == p[1].MCLID:
            continue
        midj = middle(p[1])
        common_domains = [d for d in domains0 if InDomain(d,p[1])]
        if common_domains:
            good_pairs.append([p[0],p[1],common_domains[0]])
            doms += common_domains
            if not number:
                score += len(pairs_dic[p])
            else:
                score += pairs_dic[p]
        else:
            bad_pairs.append([p[0],p[1]])
    good_pairs = sort_pairs_by_dis(set_pairs_fast(good_pairs))
    bad_pairs = sort_pairs_by_dis(set_pairs_fast(bad_pairs))
    return len(good_pairs),len(bad_pairs), len(doms), score
def dic_find_cofunc_pairs_inDomain(genes, dic, max_dis = 1000*1000, min_dis = 0, ret = 0):
    good_pairs = []
    bad_pairs = []
    doms = []
    counter = 0
    for gr in genes:
        mids = sorted([middle(g) for g in gr])
        sor_gr = []
        for mid in mids:
            for g in gr:
                if middle(g) == mid:
                    sor_gr.append(g)
                    break
        for i in range(len(sor_gr)):
            midi = middle(sor_gr[i])
            domains0 = dic[str(sor_gr[i])]
            for j in range(i+1, len(sor_gr)):
                if sor_gr[i].chromosome != sor_gr[j].chromosome:
                    continue
                if sor_gr[i].MCLID == sor_gr[j].MCLID:
                    continue
                midj = middle(sor_gr[j])
                dis = abs(midj-midi)
                if (dis > max_dis) or (dis < min_dis):
                    break
                common_domains = [d for d in domains0 if InDomain(d,sor_gr[j])]
                if common_domains:
                    good_pairs.append([sor_gr[i],sor_gr[j],common_domains[0]])
                    doms += common_domains
                else:
                    bad_pairs.append([sor_gr[i],sor_gr[j]])
    #print "finished, now removing dups, this might take a while"
    good_pairs = sort_pairs_by_dis(set_pairs_fast(good_pairs))
    bad_pairs = sort_pairs_by_dis(set_pairs_fast(bad_pairs))
    if ret == 1:
        return list(set(doms))
    return good_pairs,bad_pairs,list(set(doms))
"""
def find_cofunc_pairs_inDomain_fast(genes, max_dis = 1000*1000, min_dis = 0, cell_line = "GM",domains = kbdomains, ret = 0):
    import itertools 
    good_pairs = []
    bad_pairs = []
    counter = 0
    for gr in genes:
        mids = sorted([middle(g) for g in gr])
        sor_gr = []
        for mid in mids:
            for g in gr:
                if middle(g) == mid:
                    sor_gr.append(g)
                    break
        """DEBUG"""
        if len(genes)>10 and counter%(len(genes)/10) == 0:
            print counter,
        counter +=1
        """DEBUG"""
        for i in range(len(sor_gr)):
            midi = middle(sor_gr[i])
            domains0 = [d for d in domains if InDomain(d,sor_gr[i])]
            for j in range(i+1, len(sor_gr)):
                if sor_gr[i].chromosome != sor_gr[j].chromosome:
                            continue
                midj = middle(sor_gr[j])
                dis = abs(midj-midi)
                if (dis > max_dis) or (dis < min_dis):
                    break
                common_domains = [d for d in domains0 if InDomain(d,sor_gr[j])]
                if common_domains:
                    good_pairs.append([sor_gr[i],sor_gr[j],common_domains[0]])
                else:
                    bad_pairs.append([sor_gr[i],sor_gr[j]])
    print "finished, now removing dups, this might take a while"
    good_pairs = sort_pairs_by_dis(set_pairs_fast(good_pairs))
    bad_pairs = sort_pairs_by_dis(set_pairs_fast(bad_pairs))
    return good_pairs,bad_pairs

def tad_cofunc_pairs(genes, cell_line = "GM", doms = domains["GM"]):
    gmap = {}
    for i in range(len(genes)):
        for gene in genes[i]:
            if str(gene) in gmap:
                gmap[str(gene)].append(i)
            else:
                gmap[str(gene)] = [i]
    good_domains = {}
    for d in doms:
        if str(d) not in domain2gene[cell_line]:
            continue
        gs = domain2gene[cell_line][str(d)]
        good_domains[d] = []
        for i in range(len(gs)):
            if str(gs[i]) not in gmap:
                continue
            for j in range(i):
                if str(gs[j]) not in gmap:
                    continue
                if set(gmap[str(gs[i])]).intersection(set(gmap[str(gs[j])])):
                    good_domains[d].append((gs[i],gs[j]))
        if good_domains[d] == []:
            a = good_domains.pop(d)
    return good_domains
    

def genes_domain_dic(genes, domains):
    result = {}
    for domain in domains:
        result[str(domain)] = [[g for g in k if InDomain(domain,g)] for k in new_version]
    return result

def sort_pairs_by_dis(ps):
    d = {}
    for p in ps:
        if not d.has_key(pair_dis(p)):
            d[pair_dis(p)] = [p]
        else:
            d[pair_dis(p)].append(p)
    diss = sorted(d.keys())
    sor_pairs = []
    for dis in diss:
        sor_pairs += d[dis]
    return sor_pairs

def merge_by_dis(l1,l2):
    result = []
    i = 0
    j = 0
    while True:
        while pair_dis(l1[i]) <= pair_dis(l2[j]):
            result.append(l1[i])
            i+=1
            if i == len(l1):
                return result + l2[j:]
        while pair_dis(l1[i]) > pair_dis(l2[j]):
            result.append(l2[j])
            j+=1
            if j == len(l2):
                return result + l1[i:]

def fraction_test_old(pairs, background_pairs, rang = 10, small = True):
    fractions = []
    for p in pairs:
        i = 0
        while pair_dis(p) > pair_dis(background_pairs[i]):
            i +=1
        if small:
            fraction_group = background_pairs[i-rang:i]
        else:
            fraction_group = background_pairs[i:i+rang]
        good = 0; bad = 0
        for f in fraction_group:
            domains0 = [d for d in domainsGM if InDomain(d,f[0])]
            common_domains = [d for d in domains0 if InDomain(d,f[1])]
            if common_domains:
                good += 1
            else:
                bad += 1
        fractions.append([pair_dis(p), good, bad])
    return fractions

def fraction_test(pairs, background_pairs, domains,NN = 5, small = True, ret = 1):
    if ret == 2:
        print "WARNING: make sure you are not printing result to SCREEN"
    fractions = []
    good = []
    bad = []
    for p in pairs:
        low = 0
        high = len(background_pairs)
        i = (high+low)/2
        while high - low > 10:
            if pair_dis(p) > pair_dis(background_pairs[i]):
                low = i
                i = (low+high)/2
            if pair_dis(p) <= pair_dis(background_pairs[i]):
                high = i
                i = (low+high)/2
        i = low
        while pair_dis(p) > pair_dis(background_pairs[i]):
            i+=1
        if small:
            fraction_group = background_pairs[i-NN:i]
        else:
            fraction_group = background_pairs[i:i+NN]
        for f in fraction_group:
            domains0 = [d for d in domains if InDomain(d,f[0])]
            common_domains = [d for d in domains0 if InDomain(d,f[1])]
            if common_domains:
                good.append(f)
            else:
                bad.append(f)
    if ret ==2 :
            return good,bad
    if ret == 1:
        return map(len, [good, bad])
    return map(lambda x:len(set_pairs_fast(x)), [good, bad])


def chisq(kegg, bg):
    import scipy
    import scipy.stats
    B11 = kegg[0]*1.0
    C11 = kegg[1]*1.0
    B12 = bg[0]*1.0
    C12 = bg[1]*1.0
    D11 = B11+C11
    D12 = B12+C12
    B13 = B11+B12
    C13 = C11+C12
    D13 = B11+B12+C11+C12
    B15 = D11*B13/D13
    C15 = D11*C13/D13
    B16 = D12*B13/D13
    C16 = D12*C13/D13
    B18 = (B11-B15)**2/B15
    C18 = (C11-C15)**2/C15
    B19 = (B12-B16)**2/B16
    C19 = (C12-C16)**2/C16
    score = sum([B18,B19,C18,C19])
    return scipy.stats.chi2.sf(score,1)

def sort_by_loc(genes):
    sor = {}
    for g in genes:
        if g.start in sor:
            sor[g.start].append(g)
        else:
            sor[g.start] = [g]
    return sum_list([sor[x] for x in sorted(sor.keys())])

def value_to_pairs(value):
    pairs = []
    for domain in value:
        for path in domain:
            if len(path)>2:
                p = sort_by_loc(path)
                pairs += [[p[i],p[i+1]] for i in range(len(p)-1)]
    return pairs

def count_bad_pairs(value):
    counter = 0
    bypath = [[x[i] for x in value if x[i]] for i in range(len(value[0]))]    
    for path in bypath:
        #DEBUG:
        #print bypath.index(path), len(path), " " 
        for i in range(len(path)):
            for j in range(i):
                if path[i][0].chromosome == path[j].chromosome:
                    counter += max(len(path[i]),len(path[j]))
    return counter

def dic_findRelativeLocationTest(l = [], domains_dic = [], resolution = 2):
    excess = 0
    reloc = []
    for i in range(len(l)):
            for d in domains_dic[str(l[i])]:
                reloc.append((l[i].start - d.start + (l[i].end-l[i].start)/2.)*1.0/(d.end-d.start))
            if i%10000 == 0:
                print i,
    relocround = map(lambda x: round(x+excess, resolution), reloc)
    mifkad = [relocround.count(1.0*i/(10**resolution)) for i in range(10**resolution+int(excess*2*100))]
    import scipy
    #from scipy import stats
    #score = scipy.stats.chisquare(mifkad[3:-3]) # actually comparing to ([sum(mifkad)/100. for i in range(100)])
    return mifkad#, score
    
def findRelativeLocationTest(l = [], domains = [], resolution = 2, isTF = False, TFname = None, excess = 0.2, edge = False):
    if isTF:
        import pickle
        l = pickle.load(file("./TF/%s"%TFname,"rb"))
    reloc = []
    for i in range(len(l)):
        if edge:
            if l[i].strand == "+":
                point = 0
            else:
                point = l[i].end- l[i].start
        else:
            point = (l[i].end - l[i].start)/2
        for d in domains:
            if str(d.chromosome) == str(l[i].chromosome):
                buffer = (d.end - d.start) * excess
                if ((d.start - buffer) < l[i].start) and ((d.end + buffer) > l[i].end):
                    reloc.append((l[i].start - d.start + point)*1.0/(d.end-d.start))
        if i%10000 == 0:
            print i,
    relocround = map(lambda x: round(x+excess, resolution), reloc)
    mifkad = [relocround.count(1.0*i/(10**resolution)) for i in range(10**resolution+int(excess*2*10**resolution))]
    import scipy
    #from scipy import stats
    #score = scipy.stats.chisquare(mifkad[3:-3]) # actually comparing to ([sum(mifkad)/100. for i in range(100)])
    return mifkad#, score

def oneD_threed_test(gene_group, d = 1000*1000, res = 1000, res_dir = "1k", expected = True):
    print "res is ", res
    pairs = find_cofunc_pairs_inDomain(gene_group,d, ret = 1)
    print "finished find_cofunc_pairs_inDomain, number of pairs is ", len(pairs)
    pairs = [[p[0],p[1]] for p in pairs]
    pairs = set_pairs(pairs)
    print "finished performing set(pairs), new len is", len(pairs)
    vals = [get_matrix_cell_for_pair(p[0],p[1],res,expected, res_dir) for p in pairs]
    return vals
    

def rand_set(start, end, size):
    import random
    sett = []
    i = 0
    while len(sett) < size:
        i += 1
        sett.append(random.randint(start, end))
        sett = list(set(sett))
        if i > size*1000:
            print "ERROR randoming group",start,end,size
            return []
    return sett
        

def permute_oneD_threeD_test(tested_domains_pairs_dic, all_genes_domains_pairs_dic, res = 1000, expected = False):
    print "res is ", res
    import random
    pairs = []
    for key in tested_domains_pairs_dic.keys():
        nkey = [d for d in all_genes_domains_pairs_dic.keys() if (str(d.chromosome) == str(key.chromosome)) and (d.start == key.start)][0]
        pair_num = len(tested_domains_pairs_dic[key])
        rang = len(all_genes_domains_pairs_dic[nkey]) - 1
        rand = rand_set(0, rang, pair_num)
        rand_pairs = [all_genes_domains_pairs_dic[nkey][i] for i in rand]
        pairs += rand_pairs
    print len(pairs)
    vals = [get_matrix_cell_for_pair(p[0],p[1],res,expected) for p in pairs]
    print len(vals)
    name = random.randint(0,1000000)
    f = file("/home/rack-shamir2/home/idannuri/1D3D/random_1kb_wo/vals" + str(name), "wb")
    pickle.dump(vals, f)
    return vals
    

def find_gene(gdic, name):
    for gr in gdic.values():
        for gene in gr:
            try:
                if gene.Symbol_from_nomenclature_authority == name:
                    return gr
            except:
                if gene.feature_name == name:
                    return gr
    return -1
    
def coins_test(good1,bad1,good2,bad2, buc_size = 100):
    """all_probs = [0.37583333333333335, 0.35640886965927526, 0.3595505617977528, 0.33945239576851277, 0.32915831663326656, 0.33849751534638994, 0.31776913099870296, 0.27332621082621084, 0.27149724051985047, 0.24017540900657783, 0.18185388845247447, 0.13509335236056605, 0.10662043558348684, 0.06901635567841512, 0.00828888328471981]
    probs = [0.3076923076923077, 0.47619047619047616, 0.43478260869565216, 0.39215686274509803, 0.46511627906976744, 0.45454545454545453, 0.3225806451612903, 0.25, 0.3389830508474576, 0.3508771929824561, 0.19801980198019803, 0.22988505747126436, 0.21505376344086022, 0.10471204188481675, 0.00909090909090909]
    """
    gdists1 = sorted([abs(middle(p[0]) - middle(p[1])) for p in good1])
    bdists1 = [abs(middle(p[0]) - middle(p[1])) for p in bad1]
    gdists2 = [abs(middle(p[0]) - middle(p[1])) for p in good2]
    bdists2 = [abs(middle(p[0]) - middle(p[1])) for p in bad2]
    
    buckets1 = [gdists1[buc_size*i:buc_size*(i+1)] for i in range(len(good1)/buc_size + 1)]
    limits = [max(b) for b in buckets1]
    
    background1 = [len([d for d in bdists1 if d < max(buckets1[i])]) for i in range(len(buckets1))]
    buckets2 = [len([d for d in gdists2 if d < max(buckets1[i])]) for i in range(len(buckets1))]
    background2 = [len([d for d in bdists2 if d < max(buckets1[i])]) for i in range(len(buckets1))]

    buckets1 = map(len,buckets1)
    background1 = [background1[0]] + [background1[i] - background1[i-1] for i in range(1,len(buckets1))]
    buckets2 = [buckets2[0]] + [buckets2[i] - buckets2[i-1] for i in range(1,len(buckets1))]
    background2 = [background2[0]] + [background2[i] - background2[i-1] for i in range(1,len(buckets1))]

    gprobs1 = [buckets1[i]*1.0/(buckets1[i] + background1[i]) for i in range(len(buckets1))]
    gprobs2 = [buckets2[i]*1.0/(buckets2[i] + background2[i]) for i in range(len(buckets1))]
    return gprobs1, gprobs2, limits

                 
def calc_pvalue_isVally(mif, n = 20):
    mid = len(mif)/2
    ends_sum = sum(mif[:n/2] + mif[-n/2:])*1.0/n
    print "ends sum is", ends_sum
    avg_sum = sum(mif[mid-n/2:mid+n/2])*1.0/n
    print "avg sum is", avg_sum
    p = 1./len(mif)
    q = 1-p
    sigma = (2*sum(mif)*p*q/n)**0.5
    print "p,q,sigma", p,q,sigma
    z = abs(ends_sum-avg_sum)/sigma
    print "z",z
    import scipy
    from scipy.stats import norm
    return 2*(1-norm.cdf(z))
    
def ks_isVally(mif, n = 20):
    mid = len(mif)/2
    ends = mif[:n/2] + mif[-n/2:]
    avg = mif[mid-n/2:mid+n/2]
    import scipy
    from scipy.stats import ks_2samp
    return ks_2samp(ends,avg)

def score_isVally(mif,n=40):
    print "Z test is ",calc_pvalue_isVally(mif,n)
    print "KS ",ks_isVally(mif,n)


def pair_dist(p):
    return abs(middle(p[0])-middle(p[1]))


def set_pairs_fast(pairs):

    ids = map(lambda x: x[0]+" "+x[1],[sorted(map(lambda x:str(x.ID), p[:2])) for p in pairs])
    d =  {}
    for p in pairs:
       d[int(p[0].ID)] = p[0]
       d[int(p[1].ID)] = p[1]
    ids = list(set(ids))
    ids = [map(int,x.split()) for x in ids]
    ids = [x for x in ids if x[0] != x[1]]
    res = []
    for x in ids:
    	res.append([d[x[0]],d[x[1]]])
    return res
    
def pairs_bin_search(list1, dis, smaller = True):
    low = 0
    high = len(list1)
    while True:
        if high - low < 10:
            break
        ind = (high+low)/2
        if pair_dis(list1[ind]) >= dis:
            high = ind
        else:
            low = ind
    while pair_dis(list1[low]) < dis:
        low += 1
    if smaller:
        return low-1
    return low

def pair_key(pair):
    ids = sorted(map(lambda x:str(x.ID),pair))
    return ''.join(ids)

def consecutive_pairs_comparison(kegg_res, all_res, rang = 100, saf = 1.3,resfile = None):
    import math
    kegg = {}
    for p in sum_list(kegg_res):
        kegg[pair_key(p)] = 1
    for i in range(rang,len(all_res[0]),rang):
        mindis = pair_dis(all_res[0][i-rang])
        maxdis = pair_dis(all_res[0][i])
        good = all_res[0][pairs_bin_search(all_res[0],mindis,0):pairs_bin_search(all_res[0],maxdis,1)]
        bad = all_res[1][pairs_bin_search(all_res[1],mindis,0):pairs_bin_search(all_res[1],maxdis,1)]
        kgood = [p for p in good if pair_key(p) in kegg]
        kbad = [p for p in bad if pair_key(p) in kegg]
        kg,kb,g,b = len(kgood),len(kbad),len(good),len(bad)
        if kg > 10:
            score =  -1*math.log(chisq([kg,kb],[g-kg,b-kb]),10)
            if kg*1.0/kb < (g-kg)*1.0/(b-kb):
                score = score*(-1)
            average = sum(map(pair_dis, good))/len(good)
            res = "%d %d : %d %d %d %d, mindis - %d, maxdis - %d, score - %f\n"\
                  %(i, average, kg, kb, g-kg,b-kb,mindis,maxdis, score)
            if abs(score) > saf:
                print res
                if kg*1.0/kb > (g-kg)*1.0/(b-kb):
                    print "KEGG wins"
                else:
                    print "ALL wins"
            if resfile:
                resfile.write(res)

def test_cell_line(name ,rangs = [150], NNs = [5], exp = None, kegg_fn ="MCLBased_newKEGGPathways", bg_fn = "MCLBased_Genes", ret = 1, saf = 1.3):
    cell_domains_dic = gene2domain[name]
    cell_domains = domains[name]
    import math
    kegg_genes = pickle.load(file("./FinalData/" + kegg_fn))
    genes = pickle.load(file("./FinalData/" + bg_fn))
    kegg = dic_find_cofunc_pairs_inDomain(kegg_genes, cell_domains_dic)
    if exp:
        ger = pickle.load(file("./GeneExp/kegg_exp_%s"%name))
        genes_kegg_tmp = [genes_kegg[i] for i in range(len(genes_kegg_tmp)) if ger[i] > exp]
        print "testing cell line only for groups with gene_expression above %d, %d groups remains"%(exp, len(genes_kegg_tmp))
        kegg = find_cofunc_pairs_inDomain_fast(genes_kegg_tmp, domains=cell_domains)
    print name, map(len,kegg)
    genes_find = dic_find_cofunc_pairs_inDomain(genes, cell_domains_dic) 
    print name, map(len, genes_find)
    k = {}
    for x in [[p[0].ID, p[1].ID] for p in sum_list(kegg)]:
        k[str(x[0]) + " " + str(x[1])] = True
    bg = genes_find[0]+genes_find[1]
    bg = [p for p in bg if str(p[0].ID)+" " + str(p[1].ID) not in k and p[0].ID!=p[1].ID]
    #pickle.dump(bg, file("./TestResults/find_cofunc_pairs_inDomain_fast(genes, domains = domains%s)_WO_KEGG"%name, "wb"))
    bg = sort_pairs_by_dis(bg)
    gk, bk = kegg
    print len(gk), len(bk), len(bg)
    for rang in rangs:
        """ no overlapping allowed"""
        if ret == 1:
            fn = "./TestResults/final_fractions/%s %s consecutive %s buck size=%d"%(kegg_fn, bg_fn, name, rang)
            resfile = file(fn, "wb")
            cons = consecutive_pairs_comparison(kegg, genes_find, rang,saf, resfile)
            resfile.close()
        if ret == 2:
            jumps =  rang
            for NN in NNs:
                fn = "./TestResults/final_fractions/%s %s %s NN=%d rang=%d jumps=%d"\
                     %(kegg_fn, bg_fn, name,NN, rang, jumps)
                resfile = file(fn,"wb")
                for i in range(0,len(gk)-rang,jumps):
                    mindis = pair_dis(gk[i])
                    maxdis = pair_dis(gk[i+rang])
                    fr = fraction_test(gk[i:i+rang], bg,domains = cell_domains,NN=NN)
                    good, bad = fr
                    bktmp = [x for x in bk if pair_dis(x) > mindis and pair_dis(x) < maxdis]
                    average = sum(map(pair_dis, gk[i:i+rang]))/rang
                    resfile.write("%d [%d,%d,%d,%d]\n"%(average,rang,len(bktmp),good,bad))
                    score =  -1*math.log(chisq([rang, len(bktmp)],[good,bad]),10)
                    if rang*1.0/len(bktmp) < good*1.0/bad:
                            score = score*(-1)
                    if abs(score) > saf:
                        print i, [rang,len(bktmp),good,bad], name, maxdis
                        print score
                        if rang*1.0/len(bktmp) > good*1.0/bad:
                            print "KEGG wins"
                        else:
                            print "ALL wins"
                resfile.close()

def parse_cons_results(kegg_fn = "MCLBased_KEGGPathways"):
    cons = cons = glob.glob(r"r:\home\idannuri\TestResults\final_fractions\*%s*consecutive*"%(kegg_fn))
    dd = {}
    for fn in cons:
	data = file(fn, "rb").read().splitlines()
	cl = fn.split()[2]
	size = int(fn.split()[-1][5:])
	dd[(cl,size)] = data
    for k in dd:
	res = ""
	for x in dd[k]:
		kg,kb,ag = map(float, x.split()[3:6])
		ab = float(x.split()[6][:-1])
		kp = kg/(kg+kb)
		ap = ag/(ag+ab)
		res += "%s,%d,%d,%d,%d,%f,%f,%s\n"\
		       %(x.split()[1],kg,kb,ag,ab,kp,ap,x.split()[-1])
	dd[k] = res
    fs = {}
    for cl in ['HUVEC', 'NHEK', 'K562', 'HMEC', 'IMR90', 'HeLa', 'GM', 'KBM7']:
	fs[cl] = ""
	for size in range(1000,10000,1000):
		fs[cl] += dd[(cl,size)] + "\n"
    return fs

def parse_matching_results(kegg_fn = "MCLBased_KEGGPathways"):
    fns = glob.glob(r"r:\home\idannuri\TestResults\final_fractions\*%s*NN=*"%(kegg_fn))
    dd = {}
    for fn in fns:
	data = file(fn, "rb").read().splitlines()
	cl = fn.split()[1]
	size = int(fn.split()[3][5:])
	nn = int(fn.split()[2][3:])
	dd[(cl,size,nn)] = data
    for k in dd:
	res = ""
	for y in dd[k]:
		x = y.split()[1][1:-1]	
		kg,kb,ag,ab = map(float, x.split(","))
		kp = kg/(kg+kb)
		ap = ag/(ag+ab)
		score =  -1*math.log(chisq([kg,kb],[ag,ab]),10)
		if ap>kp:
                    score = score*(-1)
		res += "%s,%d,%d,%d,%d,%f,%f,%f\n"\
		       %(y.split()[0],kg,kb,ag,ab,kp,ap,score)
	dd[k] = res
    fs = {}
    for cl in ['HUVEC', 'NHEK', 'K562', 'HMEC', 'IMR90', 'HeLa', 'GM', 'KBM7']:
	fs[cl] = ""
	for NN in [1,2,3,5,10]:
		fs[cl] += dd[(cl,150,NN)] + "\n"
    return fs

def liat(cell_line ="GM", dis = 1000*1000*1000, tough = 0, ret = 0, dists = True):
    domainswkegg = pickle.load(file("./domains/mdomainswkegg"))[cell_line]
    wo = pickle.load(file("./FinalData/MCLBased_newKEGGPathways_withoutOlfactory"))
    import random
    good = 0
    bad = 0
    rand_vec = {}
    dvec = {}
    for key in domainswkegg:
        rindex = random.randint(0, len(domainswkegg[key])-1)
        rand_vec[key] = domainswkegg[key][rindex]
    options = {}
    for key in rand_vec:
        g2 = rand_vec[key][0]
        path = [g for g in wo[rand_vec[key][1]] \
                if g.chromosome == g2.chromosome and \
                pair_dis([g2,g]) < dis and \
                g2.MCLID!=g.MCLID and\
                gene2domain[cell_line][str(g)]]
        if len(path) > 1:
            options[key] = path
    for key in [k for k in rand_vec if k not in options]:
        a = rand_vec.pop(key)
        if tough:
            bad +=1
    if ret:
        good = []
    for key in rand_vec:
        g1 = rand_vec[key][0]
        rindex = random.randint(0, len(options[key])-1)
        g2 = options[key][rindex]
        if dists:
            dvec[key] = pair_dis([g1,g2])
        if g2.MCLID in [x.MCLID for x in mdomain2gene[cell_line][key]]:
            if ret:
                good.append(key)
            else:
                good += 1
        else:
            bad+=1
    if dists:
        return good, bad, dvec
    return good, bad
        
def liat2(cell_line=  "GM", dis = 1000*1000*1000, dists =None):
    domainswkegg = pickle.load(file("./domains/mdomainswkegg"))[cell_line]
    wo = pickle.load(file("./FinalData/MCLBased_newKEGGPathways_withoutOlfactory"))
    import random
    rand_vec = {}
    for key in domainswkegg:
        rindex = random.randint(0, len(mdomain2gene[cell_line][key])-1)
        rand_vec[key] = mdomain2gene[cell_line][key][rindex]
    good = 0
    bad = 0
    keys = rand_vec.keys()
    if dists:
        keys = dists.keys()
    for key in keys:
        g1 = rand_vec[key]
        options = genes[int(rand_vec[key].chromosome)-1]
        if not options:
            continue
        for i in range(100):
            rindex = random.randint(0, len(options)-1)
            g2 = options[rindex]
            if dists:
                dis = dists[key]
            if str(g2) not in mgene2domain[cell_line]:
                break
            if pair_dis([g1,g2]) < dis and g1.MCLID != g2.MCLID\
                and mgene2domain[cell_line][str(g2)]:
                break
        if i == 100:
            continue
        if g2.MCLID in [x.MCLID for x in mdomain2gene[cell_line][key]]:
            good += 1
        else:
            bad+=1
    return good, bad
