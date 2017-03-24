def mult_Step1(SRRIds, ID, run = False):
    os.chdir("/home/gaga/data-scratch/idannuri/ChipSeq/Step1")
    filename = "/home/gaga/data-scratch/idannuri/ChipSeq/Step1/SRRstep1_%s"%ID
    step1file = file(filename, "wb")
    for srr in SRRIds:
        line = "%s\t%s\n"%(srr, srr)
        step1file.write(line)
        print line
    step1file.close()
    if run:
        a = os.system("/home/elkon/bash_scripts/multiple_SRA_fastq_dump.sh --gzip %s 0"%(filename))
        return a
    return 0

def single_Step1(SRRId,  run = False):
    os.chdir("/home/gaga/data-scratch/idannuri/ChipSeq/Step1")
    command = "fastq-dump --gzip %s 0"%(SRRId)
    print command
    if run:
        a = os.system(command)
        return a
    return 0

def Step2(title, SRRIds, ID, run = False):
    if type(SRRIds) == type("a"):
        tmp = file(SRRIds).readlines()
        SRRIds = [line.split("\t")[0] for line in tmp]
    os.chdir("/home/gaga/data-scratch/idannuri/ChipSeq/Step2")
    filename = "/home/gaga/data-scratch/idannuri/ChipSeq/Step2/Step2_%s"%ID
    files = sum_list([glob.glob("/home/gaga/data-scratch/idannuri/ChipSeq/Step1/*%s*"%x)\
                      for x in SRRIds])
    step2file = file(filename, "wb")
    for srr in files:
        line = "%s\t%s\n"%(srr, srr.split("/")[-1].split(".")[0])
        step2file.write(line)
        print line
    step2file.close()
    if run:
        a = os.system("/home/elkon/bash_scripts/multiple_bwt2_gAlign.sh %s %s Hs 0 0 1 0"\
                  %(filename, title))
        return a
    return 0

def Step3(pairs, ID, run = False):
    if type(pairs)== type("a"):
        tmp = file(pairs).readlines()
        pairs = [x.split("\t") for x in tmp]
    os.chdir("/home/gaga/data-scratch/idannuri/ChipSeq/Step3")
    filename = "/home/gaga/data-scratch/idannuri/ChipSeq/Step3/SRRstep3_%s"%ID
    step3file = file(filename, "wb")
    for pair in pairs:
        p0 = pair[0].split("/")[-1].split(".")[0]
        p1 = pair[1].split()[0].split("/")[-1].split(".")[0]
        p0_bam = "/home/gaga/data-scratch/idannuri/ChipSeq/Step2/bwt2_gAlign/%s.hg19.filtered.sorted.bam"%p0
        p1_bam = "/home/gaga/data-scratch/idannuri/ChipSeq/Step2/bwt2_gAlign/%s.hg19.filtered.sorted.bam"%p1
        line = "%s\t%s\t%s\n"%(p0_bam, p1_bam, p0+"_"+p1)
        step3file.write(line)
        print line
    step3file.close()
    if run:
        a = os.system("/home/elkon/bash_scripts/multiple_macs_two_samples_bamfiles.sh %s Hs 0 0"%(filename))
        return a
    return 0

def GEOpipeline(title, pairs, ID, run = False):
    SRRIds = list(set(sum_list(pairs)))
    print "running Step1"
    a = mult_Step1(SRRIds, ID, run)
    print "running Step2"
    a = Step2(title, SRRIds, ID, run)
    print "running Step3"
    a = Step3(pairs, ID, run)
    
def generate_Figure4a_Chipseq(pair, cl):
    peaks = "/home/gaga/data-scratch/idannuri/ChipSeq/Step3/macs2_out/%s_%s_peaks.narrowPeak"%(pair[0],pair[1])
    peaks = map(Chipseq, file(peaks).readlines())
    return check_regexp_AB(hun[cl][0],hun[cl][b], peaks)
    
    

    
