#!/usr/bin/python
import os,sys,re,string,math
import matplotlib as plt
import numpy as nm
import sqlite3
import random

cwd=os.path.realpath('./')
sys.path.append('%s/libsvm-3.16/python' % cwd)
import svmutil


#============= constant =======================
default_trainoptions='-s 0 -t 1 -b 1 -v 10'
""" -s:svm_type=C-SVC -t:kernel-type=RBF -b:probability-estimate  -v:cross validation n-fold"""
default_predicoptions='-b 1'
aa_table={
'GCA':'A','GCC':'A','GCG':'A','GCT':'A',                      # Alanine
'TGC':'C','TGT':'C',                                          # Cysteine
'GAC':'D','GAT':'D',                                          # Aspartic Acid
'GAA':'E','GAG':'E',                                          # Glutamic Acid
'TTC':'F','TTT':'F',                                          # Phenylalanine
'GGA':'G','GGC':'G','GGG':'G','GGT':'G',                      # Glycine
'CAC':'H','CAT':'H',                                          # Histidine
'ATA':'I','ATC':'I','ATT':'I',                                # Isoleucine
'AAA':'K','AAG':'K',                                          # Lysine
'CTA':'L','CTC':'L','CTG':'L','CTT':'L','TTA':'L','TTG':'L',  # Leucine
'ATG':'M',                                                    # Methionine
'AAC':'N','AAT':'N',                                          # Asparagine
'CCA':'P','CCC':'P','CCG':'P','CCT':'P',                      # Proline
'CAA':'Q','CAG':'Q',                                          # Glutamine
'CGA':'R','CGC':'R','CGG':'R','CGT':'R','AGA':'R','AGG':'R',  # Arginine
'TCA':'S','TCC':'S','TCG':'S','TCT':'S','AGC':'S','AGT':'S',  # Serine
'ACA':'T','ACC':'T','ACG':'T','ACT':'T',                      # Threonine
'GTA':'V','GTC':'V','GTG':'V','GTT':'V',                      # Valine
'TGG':'W',                                                    # Tryptophan
'TAC':'Y','TAT':'Y',                                          # Tyrosine
'TAA':'U','TAG':'U','TGA':'U'                                 # Stop
}

# amino acids are translated in bi-encoding,every amino acid is encoded as a tuple('...length=21..'),an instance contains $kmer tuples
bi_table={'A':'1000000000000000000000', 'C':'0100000000000000000000', 'D':'0010000000000000000000', 'E':'0001000000000000000000', 'F':'0000100000000000000000', 'G':'0000010000000000000000', 'H':'0000001000000000000000', 'I':'0000000100000000000000', 'K':'0000000010000000000000', 'L':'0000000001000000000000', 'M':'0000000000100000000000', 'N':'0000000000010000000000', 'P':'0000000000001000000000', 'Q':'0000000000000100000000', 'R':'0000000000000010000000', 'S':'0000000000000001000000', 'T':'0000000000000000100000', 'V':'0000000000000000010000', 'W':'0000000000000000001000', 'Y':'0000000000000000000100', 'U':'0000000000000000000010','X':'0000000000000000000001'}
# amino acids are translated in int-encoding,every amino acid is encoded as a single int.an instance contains $kmer ints
int_table={'A' :1, 'C' :2, 'D' :3, 'E' :4, 'F' :5, 'G' :6, 'H' :7, 'I' :8, 'K' :9, 'L' :10, 'M' :11, 'N' :12, 'P' :13, 'Q' :14, 'R' :15, 'S' :16, 'T' :17, 'V' :18, 'W' :19, 'Y' :20, 'U' :21, 'X':22}

mhcp_table={}
#----------------------

def read_fasta(file):
    f=open(f,'r')
    lines=[]
    while True:
        line=f.read()
        if re.match(line,"^>"):
            contiue
    lines.append(line)
    return lines

def read_mhcpep(filename):
    """ read 'mhcpep.txt' formated file,translate each records into a dict,return an iterator"""
    with open(filename,'r') as file:
        ID=re.compile('^>')
        MOLECULE=re.compile('^MHC MOLECULE')
        ACTIVITY=re.compile('^ACTIVITY')
        BINDING=re.compile('^BINDING')
        SEQUENCE=re.compile('^SEQUENCE')
        PREFERNCES=re.compile('^PREFERNCES')
        END=re.compile('...')
        while True:
            line=file.readline().strip()
            if not line:break
            if ID.match(line):
                id=line[1:len(line)-1]
            if MOLECULE.match(line):
                re.sub('\s+','',line)
                sublines=re.split("[:,\s]+",line)
                hla=re.sub('\(.*\)','',sublines[2])
                type=sublines[3]
                species=re.sub('[\(\)#]','',sublines[4])
            if ACTIVITY.match(line):
                if line == 'ACTIVITY: none#':
                    activity=-1
                if line == 'ACTIVITY: ?#':
                    activity=0
                if line == 'ACTIVITY: yes, ?#':
                    activity=1
                if line == 'ACTIVITY: yes, little#':
                    activity=2
                if line == 'ACTIVITY: yes, moderate#':
                    activity=3
                if line == 'ACTIVITY: yes, high#':
                    activity=4
            if BINDING.match(line):
                if line == 'BINDING: yes, ?#':
                    binding=0
                if line == 'BINDING: yes, little#':
                    binding=1
                if line == 'BINDING: yes, moderate#':
                    binding=2
                if line == 'BINDING: yes, high#':
                    binding=3
            if SEQUENCE.match(line):
                re.sub('[\s]','',line)
                sublines=re.split("[:,\(\)\s]+",line)
                sequence=sublines[1]
                sequence=re.sub('#','',sequence)
            if line.strip() == '...#':
                record={'id':id,'hla':hla,'type':type,'species':species,'activity':activity,'binding':binding,'sequence':sequence}
                yield record

def check_table(curs,tablename):
    query='select count(*) from sqlite_master where type=\'table\' and name=\'%s\'' % tablename
    curs.execute(query)
    value=curs.fetchone()[0]
    return value

def create_table():
    pass 
                
def storeDB_mhcpep(conn,record_iter):
    tablename='mhcpep'
    query='insert into mhcpep (id,hla,type,species,activity,binding,sequence) values (?,?,?,?,?,?,?)'
    for record in record_iter:
        values=[record['id'],record['hla'],record['type'],record['species'],record['activity'],record['binding'],record['sequence']]
        try:
            curs.execute(query,values)
        except Exception,error:
            print error
    
def split_kmer(sequence,k):
    "split random sequence to k-mer"
    seg=len(sequence) // k
    i=0
    while i < seg:
        kmer=sequence[i*k:k*(i+1)]
        i+=1
        yield kmer

def aa_translate(sequence,table):
    "translate amino acid to bi-code"
    bi_list=(table[i.upper()] for i in sequence)
    return bi_list

def divide_listlike(listlike,ratio):
    """ divide a listlike data into 2 parts: positive counts for $ratio of the entire index, return 2 lists of index """ 
    length=len(listlike)
#    positive_index=[i for i in range(length) if random.random() <=p ratio]
    positive_index=random.sample(range(length),int(length*ratio))
    negtive_index=[i for i in range(length) if i not in positive_index]
    return positive_index,negtive_index

def get_randomchunk(sequence,chunklength):
    """ cut a segment from a very long fasta file at random position"""
    totallength=len(sequence)
    if totallength>chunklength:
        start=int(random.uniform(0,totallength-chunklength-1))
        chunk=sequence[start:start+chunklength]
        return chunk
    else:
        raise Exception('chunk is too long, shrink it.')

def score_pssm():
    pass

def score_slide(sequence,k):
    """ slide the sequence,return all the k-mer"""
    l=len(sequence)
    if l>=k:
        return [sequence[i:i+k] for i in range(l-k)]
    else:
        return []

##-------------------- main loop --------------
if __name__=='__main__':
    here=os.path.realpath('./')
    #    conn=sqlite3.connect(here+'/mhcpep.db')
    #    curs=conn.cursor()
    #
    ###step1 : convert mhcpep.txt into sqltie3
    file=here+'/mhcpep/mhcpep.txt'
    records=read_mhcpep(file)
    #    if not check_table(curs,'mhcpep'):
    #        pass
    #    storeDB_mhcpep(curs,records)
    #    conn.commit()
    #    conn.close()
    ##step2.1 : get random negtive samples,random.aa is a fasta file with one-line title and one-line sequence.
    random.aa=here+'/random.aa'
    try:
        with open(random.aa,'r') as f:
            random_seq=f.readline() # pass the title
            random_seq=f.readline().strip() 
    except Exception,e:
        print e
    random_chunk=get_randomchunk(random_seq,200)
    random_kmers=split_kmer(random_chunk,9)
    negtive=[]
    for random_kmer in random_kmers:
        negtive.append( list(aa_translate(random_kmer,int_table)))
## step2.2 : get mhcp positive kmers
    hla='HLA-DR1'
    limit=10 # take at most 200 hlas
    hla_seq=[i['sequence'].strip('*') for i in records if i['hla'] == hla and i['type'] == 'CLASS-2']
    if len(hla_seq) > limit:
        hla_seq=hla_seq[0:limit-1]
    positive=[]
    for seq in hla_seq:
        hla_kmers=score_slide(seq,9)
        for hla_kmer in hla_kmers:
            positive.append(list(aa_translate(hla_kmer,int_table)))
## step 2.3 : differ the negtive samples

## step3 : train svm
#    positive_index,negtive_index = divide_listlike(trainingset,0.4)
#    labels=[1 if i in positive_index else -1 for i in range(len(trainingset))]
#    Prob=svmutil.svm_problem(labels,trainingset)
#    Pram=svmutil.svm_parameter('-c 4')
#    Train=svmutil.svm_train(Prob,Pram)
#    pwd=os.popen("pwd").read().strip()
#    svmutil.svm_save_model('%s/test.model' % pwd,Train)
#    testlines=read_fasta('')
#    p_label,p_acc,p_val=svmutil.svm_predict(y, x, m, '-b 1')
#    ACC, MSE, SCC = evaluations(y, p_label)

"""
model=svm_train(y,x,'-t 2 -c 5') / (problem) / (problem,param) all  is good.
     Y type must be int/double. param is generated by svm_parameter('training options'), problem is generated by svm_problem(y,x[isKernel=True]).

p_labels, p_acc, p_vals=svm_predict(y,x,model,'predict options')
     p_vals is a list of decision values or probability estimates (if '-b 1' is specified) . If k is the number of classes in training data,for decision values,each element  includes results of predictiong K(K-1)/2 binary-class SVMs. For classification,k=1 is a special case. Decision value [+1] is returned for each testing instance,instead of an empty list. For probabilityes ,each element contains K values indicating the probability that the testing instance is in each class.The order of classes is the same as the  'model.label' field in the mdoel structure.

svm_read_problem(y,x)
y,x = svm_load_model('../mhcp.txt')
svm_save_model()
ACC,MSE,SCC = evaluations(ty,pv): ty is a list of true values;pv is a list of predict values;SCC is squared correlation coefficien.

struct svm_parameter{
int svm_type;
int kernel_type;
int degree; /* for poly */
int gamma; /* for poly,rbf,sigmoid */
/* these are for training */
doble cache_size ; /* in MB */
double eps;
double C;
int nr_weight; /* fro C_SVC */
int *weight_label;
double *weight;
double nu; /* for NU_SVC,ONE_CLASS,NU_SVR */
double p; /* for EPSILON_SVR */
int shrinking; /* use the shirinking heuristics */
int probability; /* do probability estimates */
}
"""
