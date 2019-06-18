class Run:
    def __init__(self,header, line):
        self.RunAttr = {}
        headers = header.split("\t")
        attr = line.split("\t")
        for i in range(len(headers)):
            self.RunAttr[headers[i]] = attr[i]

    def __str__(self):
        res = ""
        for x in self.RunAttr:
            res +="%s %s\n"%(x, self.RunAttr[x])
        return res
        
    def __repr__(self):
        return str(self)

class Exp:
    """
    file example at
    /home/gaga/data-scratch/idannuri/ChipSeq/SraRunTables/
    """
    def __init__(self,filename):
        self.id = filename.split("/")[-1].split(".")[0]
        self.runs = {}
        runs = file(filename).readlines()
        for run in runs[1:]:
            #exp is indexing the runs by the SRRId
            r = Run(runs[0],run)
            self.runs[r.RunAttr["Run"]] = r
            
    def Step1(self, run = False, force = False, parallel = True):
        os.chdir("/home/gaga/data-scratch/idannuri/ChipSeq/Step1")
        filename = "/home/gaga/data-scratch/idannuri/ChipSeq/Step1/SRRstep1_%s"%self.id
        step1file = file(filename, "wb")
        status = self.Status_Step1()
        for srr in self.runs:
            if not force:
                if srr not in status:
                    continue
            line = "%s\t%s\n"%(srr, srr)
            step1file.write(line)
            print line
        step1file.close()
        if run:
            if parallel:
                a = os.system("/home/elkon/bash_scripts/multiple_SRA_fastq_dump.sh %s 0"%(filename))
            else:
                for srr in SRRIds:
                    print srr
                    os.system("/home/elkon/tools/sratoolkit.2.8.0-ubuntu64/bin/fastq-dump.2.8.0 --gzip %s 0"%srr)
            os.system("rm /specific/a/home/cc/students/cs/idannuri/ncbi/public/sra/*")
            return a
        return 0
       
    def Status_Step1(self, debug = False):
        from subprocess import check_output
        os.chdir("/home/gaga/data-scratch/idannuri/ChipSeq/Step1")
        res = ""
        rv = []
        for srr in self.runs:
            try:
                out = check_output("ls -la ./%s.fastq.gz"%srr, shell = True)
                print out
            except:
                rv.append(srr)
                res += "%s %s File not found\n"%(srr, self.runs[srr].RunAttr["Assay_Type"])
                continue
            size = int(out.split()[4])/(1024*1024)
            exp_size = int(self.runs[srr].RunAttr["MBytes"])
            if size < exp_size:
                rv.append(srr)
                res += "%s %s %s %s Size too small\n"%(srr, size, exp_size, self.runs[srr].RunAttr["Assay_Type"])
        if debug:
            print res
        return rv

    def Step2(self, run = False, force = False, Assay_Type = "RNA-Seq"):
        """
        Step 2 of pipeline
        if Assay_Type = 'ChIP-Seq' - we run for chipseq
        if Assay_Type = 'RNA'-Seq - we run for RNA seq
        """
        os.chdir("/home/gaga/data-scratch/idannuri/ChipSeq/Step2")
        filename = "/home/gaga/data-scratch/idannuri/ChipSeq/Step2/Step2_%s_%s"%(self.id, Assay_Type)
        if not force:
            #we don't want to forcely rerun preprocessing we already did
            SRRids = self.Status_Step2(Assay_Type)
            if SRRids == []:
                print "All SRRIds already processed, run with force = True if reprocessing is needed"
                return 0
        paths = sum_list([glob.glob("/home/gaga/data-scratch/idannuri/ChipSeq/Step1/*%s*"%x)\
                          for x in SRRids])
        step2file = file(filename, "wb")
        for srr in paths:
            line = "%s\t%s\n"%(srr, srr.split("/")[-1].split(".")[0])
            step2file.write(line)
            print line
        step2file.close()
        if run:
            if Assay_Type == 'ChIP-Seq':
                a = os.system("/home/elkon/bash_scripts/multiple_bwt2_gAlign.sh %s %s Hs 0 0 1 0 0"\
                          %(filename, self.id))
            elif Assay_Type == 'RNA-Seq':
                a = os.system("/home/elkon/bash_scripts/multiple_tophat2.sh %s 0 Hs 0 1 0 %s 0 0"\
                          %(filename, self.id))
            else:
                print "Wrong Assay_Type=%s"%data
                return -1
            return a
        return 0
    
    def Status_Step2(self, Assay_Type, debug = False):
        os.chdir("/home/gaga/data-scratch/idannuri/ChipSeq/Step2")
        res = []
        for srr in self.runs:
            if self.runs[srr].RunAttr["Assay_Type"] == Assay_Type:
                if Assay_Type == "ChIP-Seq":
                    if glob.glob("./bwt2_gAlign/%s.hg19.filtered.sorted.bam"%srr)==[]:
                        res.append(srr)
                elif Assay_Type == "RNA-Seq":
                    if glob.glob("./tophat2/%s"%srr)==[]:
                        res.append(srr)
        return res

    def pairs_Step3(self, Assay_Type):
        # making pairs for comparisson - basal vs. treatment
        os.chdir("/home/gaga/data-scratch/idannuri/ChipSeq/Step3")
        filename = "./Step3_%s_%s"%(self.id, Assay_Type) 
        step3file = file(filename,"wb")
        assay = {}
        for x in exp.runs:
            if exp.runs[x].RunAttr["Assay_Type"] == Assay_Type:
                    assay[x] = exp.runs[x]
        for x in assay:
            #finding basal assays
            treat = assay[x].RunAttr["treated_with"]
            if treat and "none" not in treat:
                    continue
            cell_line = assay[x].RunAttr["source_name"]
            #finding treated assays with same cell line and conditions
            tmp = [srr for srr in assay if \
                   assay[srr].RunAttr["source_name"] == cell_line]
            if Assay_Type == "ChIP-Seq":
                antibody  = assay[x].RunAttr["chip_antibody"]
                tmp = [srr for srr in tmp if \
                       assay[srr].RunAttr["chip_antibody"] == antibody]
            tmp = [srr for srr in tmp if assay[srr].RunAttr["treated_with"]]
            #adding pairs to the file
            for treatment in tmp:
                if x != treatment:
                    step3file.write("%s %s %s\n"%(x,treatment,cell_line))
        step3file.close()
        
    def Step3_ChIP(self, Assay_Type = "ChIP-Seq", run = False):
        
        #comparnig treatment with basal
        os.chdir("/home/gaga/data-scratch/idannuri/ChipSeq/Step3")
        pairs = file("./Step3_%s_%s"%(self.id, Assay_Type)).readlines()
        filename = "./Step3_%s_%s_cli"%(self.id, Assay_Type)
        step3file = file(filename, "wb")
        for pair in pairs:
            p0,p1,cl = pair.split()
            p0_bam = "/home/gaga/data-scratch/idannuri/ChipSeq/Step2/bwt2_gAlign/%s.hg19.filtered.sorted.bam"%p0
            p1_bam = "/home/gaga/data-scratch/idannuri/ChipSeq/Step2/bwt2_gAlign/%s.hg19.filtered.sorted.bam"%p1
            line_induced = "%s\t%s\t%s\n"%(p0_bam, p1_bam, cl+"_"+p0+"_"+p1+"_induced")
            line_repressed = "%s\t%s\t%s\n"%(p1_bam, p0_bam, cl+"_"+p1+"_"+p0+"_repressed")
            step3file.write(line_induced)
            step3file.write(line_repressed)
            print line_induced,line_repressed
        step3file.close()
        if run:
            a = os.system("/home/elkon/bash_scripts/multiple_macs_two_samples_bamfiles.sh %s Hs 0 0"%(filename))
            os.system("rm ./macs2_out/*.bdg")
            return a
        return 0

    def Step3_RNA(self, force = False, run = False):
        #RNA-Seq to FPKM
        os.chdir("/home/gaga/data-scratch/idannuri/ChipSeq/Step3")
        if not force and glob.glob("./FPKM/%s_rpkms*"%title)!=[]:
            print "Counts already processed for %s"%title
            return 0
        os.chdir("/home/gaga/data-scratch/idannuri/ChipSeq/Step2")
        line = "/home/elkon/bash_scripts/cnts_data_to_RPKM_qnorm_10qfloor_using_R.sh %s ./%s_HTseq_merged_levels.txt Hs"%(title,title)
        print title, line
        if run:
            a = os.system(line)
            f = "\t" + file("./%s_rpkms_data_qnorm.20q.floor.txt"%title).read()
            new = file("./%s_rpkms_data_qnorm.20q.floor.txt"%title,"wb")
            new.write(f)
            new.close()
            a = os.system("cp ./%s_rpkms_data_qnorm.20q.floor.txt ../Step3/FPKM/"%title)
        return 0
  
"""
DEPRECATED
"""


def mult_Step1(SRRIds, ID, run = False, force = False, parallel = True):
    os.chdir("/home/gaga/data-scratch/idannuri/ChipSeq/Step1")
    filename = "/home/gaga/data-scratch/idannuri/ChipSeq/Step1/SRRstep1_%s"%ID
    step1file = file(filename, "wb")
    for srr in SRRIds:
        if not force:
            if glob.glob("/home/gaga/data-scratch/idannuri/ChipSeq/Step1/*%s*"%srr) != []:
                continue
        line = "%s\t%s\n"%(srr, srr)
        step1file.write(line)
        print line
    step1file.close()
    if run:
        if parallel:
            a = os.system("/home/elkon/bash_scripts/multiple_SRA_fastq_dump.sh %s 0"%(filename))
        else:
            for srr in SRRIds:
                print srr
                os.system("/home/elkon/tools/sratoolkit.2.8.0-ubuntu64/bin/fastq-dump.2.8.0 --gzip %s 0"%srr)
        os.system("rm /specific/a/home/cc/students/cs/idannuri/ncbi/public/sra/*")
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

def Step2(title, SRRIds, ID, run = False, force = False, data = 'Chip'):
    """
    Step 2 of pipeline
    if data = 'Chip' - we run for chipseq
    if data = 'RNA' - we run for RNA seq
    """
    if type(SRRIds) == type("a"):
        tmp = file(SRRIds).readlines()
        SRRIds = [line.split("\t")[0] for line in tmp]
    os.chdir("/home/gaga/data-scratch/idannuri/ChipSeq/Step2")
    filename = "/home/gaga/data-scratch/idannuri/ChipSeq/Step2/Step2_%s"%ID
    if not force:
        #we don't want to forcely rerun preprocessing we already did
        tmp = []
        if data == "Chip":
            for srr in SRRIds:
                if glob.glob("/home/gaga/data-scratch/idannuri/ChipSeq/Step2/bwt2_gAlign/%s.hg19.filtered.sorted.bam"%srr)==[]:
                    tmp.append(srr)
        elif data == "RNA":
            for srr in SRRIds:
                if glob.glob("/home/gaga/data-scratch/idannuri/ChipSeq/Step2/tophat2/%s"%srr)==[]:
                    tmp.append(srr)
        SRRIds = tmp
        if SRRIds == []:
            print "All SRRIds already processed, run with force = True if reprocessing is needed"
            return 0
    files = sum_list([glob.glob("/home/gaga/data-scratch/idannuri/ChipSeq/Step1/*%s*"%x)\
                      for x in SRRIds])
    step2file = file(filename, "wb")
    for srr in files:
        line = "%s\t%s\n"%(srr, srr.split("/")[-1].split(".")[0])
        step2file.write(line)
        print line
    step2file.close()
    if run:
        if data == 'Chip':
            a = os.system("/home/elkon/bash_scripts/multiple_bwt2_gAlign.sh %s %s Hs 0 0 1 0 0"\
                      %(filename, title))
        elif data == 'RNA':
            a = os.system("/home/elkon/bash_scripts/multiple_tophat2.sh %s 0 Hs 0 1 0 %s 0 0"\
                      %(filename, title))
        else:
            print "Wrong flag data=%s"%data
            return -1
        return a
    return 0

def Step3(pairs, ID, run = False):
    if type(pairs)== type("a"):
        tmp = file(pairs).readlines()
        pairs = [x.split("\t") for x in tmp]
        pairs = [[x.split()[0] for x in p] for p in pairs]
    os.chdir("/home/gaga/data-scratch/idannuri/ChipSeq/Step3")
    filename = "/home/gaga/data-scratch/idannuri/ChipSeq/Step3/SRRstep3_%s"%ID
    step3file = file(filename, "wb")
    for pair in pairs:
        p0 = pair[0].split("/")[-1].split(".")[0]
        p1 = pair[1].split("/")[-1].split(".")[0]
        cl = pair[2]
        if p0 == p1:
            continue
        p0_bam = "/home/gaga/data-scratch/idannuri/ChipSeq/Step2/bwt2_gAlign/%s.hg19.filtered.sorted.bam"%p0
        p1_bam = "/home/gaga/data-scratch/idannuri/ChipSeq/Step2/bwt2_gAlign/%s.hg19.filtered.sorted.bam"%p1
        line_induced = "%s\t%s\t%s\n"%(p0_bam, p1_bam, cl+"_"+p0+"_"+p1+"_induced")
        line_repressed = "%s\t%s\t%s\n"%(p1_bam, p0_bam, cl+"_"+p1+"_"+p0+"_repressed")
        step3file.write(line_induced)
        step3file.write(line_repressed)
        print line_induced,line_repressed
    step3file.close()
    if run:
        a = os.system("/home/elkon/bash_scripts/multiple_macs_two_samples_bamfiles.sh %s Hs 0 0"%(filename))
        os.system("rm ./macs2_out/*.bdg")
        return a
    return 0

def Step3_RNA_old(title, SRRIds, run = False, force = False):
    for srr in SRRIds:
        if not force and glob.glob("/home/gaga/data-scratch/idannuri/ChipSeq/Step3/FPKM/%s.FPKM"%srr)!=[]:
            continue
        line = "/home/elkon/tools/cufflinks-2.2.1/cufflinks --GTF /home/elkon/annots/GENCODE/Hs/hg19/v25/gencode.v25lift37.annotation.gtf accepted_hits.bam"
        print srr, line
        if run:
            os.chdir("/home/gaga/data-scratch/idannuri/ChipSeq/Step2/tophat2/%s"%srr)
            a = os.system(line)
            a = os.system('cp ./genes.fpkm_tracking /home/gaga/data-scratch/idannuri/ChipSeq/Step3/FPKM/%s.FPKM'%srr)
    return 0

def Step3_RNA(title, run = False, force = False):
    if not force and glob.glob("/home/gaga/data-scratch/idannuri/ChipSeq/Step3/FPKM/%s_rpkms*"%title)!=[]:
        print "Counts already processed for %s"%title
        return 0
    os.chdir("/home/gaga/data-scratch/idannuri/ChipSeq/Step2")
    line = "/home/elkon/bash_scripts/cnts_data_to_RPKM_qnorm_10qfloor_using_R.sh %s ./%s_HTseq_merged_levels.txt Hs"%(title,title)
    print title, line
    if run:
        a = os.system(line)
        f = "\t" + file("./%s_rpkms_data_qnorm.20q.floor.txt"%title).read()
        new = file("./%s_rpkms_data_qnorm.20q.floor.txt"%title,"wb")
        new.write(f)
        new.close()
        a = os.system("cp ./%s_rpkms_data_qnorm.20q.floor.txt ../Step3/FPKM/"%title)
    return 0

def GEOpipeline_Chipseq(title, filename, ID, run = False, force = False, skip1 = False, skip2 = False):
    """
    each line in the filename should contain 2 SRRIds and the cell line
    example : "SRR0001\tSRR0002\tMCF7\n"
    """
    pairs = file(filename).read()
    pairs = [x.split("\t") for x in pairs.splitlines()]
    SRRIds = list(set(sum_list([x[:2] for x in pairs])))
    print "running Step1"
    if not skip1:
        a = mult_Step1(SRRIds, ID, run, force)
    print "running Step2"
    if not skip2:
        a = Step2(title, SRRIds, ID, run, force)
    print "running Step3"
    a = Step3(pairs, ID, run)
    if run:
        return mult_generateFigure4a_Chipseq(filename)
    

def GEOpipeline_RNASeq(title, filename = "", ID = "", force = False, run = False, skip1 = False, skip2 = False):
    """
    each object in pairs should contain 2 SRRIds and the cell line
    example : ["SRR0001", "SRR0002", "MCF7"]
    """
    if filename == "":
        filename = "/home/gaga/data-scratch/idannuri/ChipSeq/pipeline_input/%s.txt"%title
    if ID == "":
        ID = title
    # if pairs are a string and not a list
    pairs = file(filename).read()
    pairs = [x.split("\t") for x in pairs.splitlines()]
    SRRIds = list(set(sum_list([x[:2] for x in pairs])))
    if not skip1:
        print "running Step1"
        a = mult_Step1(SRRIds, ID, run,force)
    if not skip2:
        print "running Step2"
        a = Step2(title, SRRIds, ID, run, force, data = "RNA")
    print "running Step3"
    a = Step3_RNA(title, run, force)
    if run:
        return mult_generateFigure5a_RNAseq(filename, FPKM = True)
