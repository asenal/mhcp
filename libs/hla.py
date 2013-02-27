#!/usr/bin/python
import os,sys,re,string,math
import matplotlib as plt
import numpy as nm
import sqlite3
import random
here=os.getcwd()
sys.path.append(here+'/libs/libsvm-3.12/python')
#sys.path.append(here+'/libsvm-3.16/python')
import svmutil
#------constant declaration
default_trainoptions='-s 0 -t 1 -b 1 -v 10'
""" -s:svm_type=C-SVC -t:kernel-type=RBF -b:probability-estimate  -v:cross validation n-fold"""

odefault_predicoptions='-b 1'

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
int_table={'A' :1, 'C' :2, 'D' :3, 'E' :4, 'F' :5, 'G' :6, 'H' :7, 'I' :8, 'K' :9, 'L' :10, 'M' :11, 'N' :12, 'P' :13, 'Q' :14, 'R' :15, 'S' :16, 'T' :17, 'V' :18, 'W' :19, 'Y' :20, 'U' :21,'X':22} # X represent for unknow or ambiguity
rint_table={'1' :'A','2' :'C', '3' :'D', '4' :'E', '5' :'F', '6' :'G', '7' :'H', '8' :'I', '9' :'K', '10' :'L', '11' :'M', '12' :'N', '13' :'P', '14' :'Q', '15' :'R', '16' :'S', '17' :'T', '18' :'V', '19' :'W', '20' :'Y', '21' :'U','22':'X'} # X represent for unknow or ambiguiyt

#------------subroutin defination
def read_fasta(filename):
    """read  fasta file,return a generator.An item in the generator is (title,seq)"""
    f=open(filename)
    print "read_fasta: open %s failed" % filename
    name, seq = None, []
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    f.close()


def read_mhcpep(filename):
    """ read 'mhcpep.txt' formated file,translate each records into a dict,return an iterator
<record_iter=fun(filename)>"""
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

def storeDB_mhcpep(conn,record_iter):
    """store records to database.
    <fun(conn,record_iter)>"""
    tablename='mhcpep'
    query='insert into mhcpep (id,hla,type,species,activity,binding,sequence) values (?,?,?,?,?,?,?)'
    for record in record_iter:
        values=[record['id'],record['hla'],record['type'],record['species'],record['activity'],record['binding'],record['sequence']]
        try:
            curs.execute(query,values)
        except Exception,error:
            print error
    
def split_kmer(sequence,k):
    """split random sequence to k-mer
    <fun(sequecne,k)>"""
    seg=len(sequence) // k
    i=0
    kmers=[]
    while i < seg:
        kmer=sequence[i*k:k*(i+1)]
        i+=1
        kmers.append(kmer)
    return kmers
#        yield kmer

def aa_translate(sequence,table):
    """Translate amino acid to bi-code
    <translated=fun(sequence,table)>"""
    bi_list=[table[i.upper()] for i in sequence]
    return bi_list

def divide_listlike(listlike,ratio):
    """ divide a listlike data into 2 parts: positive counts for $ratio of the entire index
    return 2 lists of index
    <postive_index,negtive_index=fun(L,ratio)> """ 
    length=len(listlike)
#    positive_index=[i for i in range(length) if random.random() <=p ratio]
    positive_index=random.sample(range(length),int(length*ratio))
    negtive_index=[i for i in range(length) if i not in positive_index]
    return positive_index,negtive_index

def get_randomchunk(sequence,chunklength):
    """ cut a segment from a very long fasta file at random position
    <chunnk=fun(sequence,chunklength)>"""
    totallength=len(sequence)
    if totallength>chunklength:
        start=int(random.uniform(0,totallength-chunklength-1))
        chunk=sequence[start:start+chunklength]
        return chunk
    else:
        raise Exception('chunk is too long, shrink it.')

def get_random_kmers(chunk_size,window_size):
    """get_randomchunk from a predefined random.aa then cut the chunk into non-overlap kmers
    <kmers=fun(chunk_size,window_size)>"""
    random.aa=here+'/data/random.aa'
    try:
        with open(random.aa,'r') as f:
            random_seq=f.readline() # pass the title
            random_seq=f.readline().strip() 
    except Exception,e:
        print e
    random_chunk=get_randomchunk(random_seq,chunk_size)
    random_kmers=split_kmer(random_chunk,window_size)
    return random_kmers

def trainset_convert(hla_kmers,random_kmers,ratio,trans_table):
    """input random_kmers & hla_kmers
    1. adjust the ratio of 2 sets
    2. allocate positive & negtive sample,label them
    3. bi-code the samples to SVM format.
    4. sample ratio% from each of the 2 sets to be trainset,the rest to be test set.If only trainset is wanted,set ratio to 1.
    <[train_data,train_label=fun(hla_kmers,random_kmers,ratio)]>
    """
    # step1:adjust ration
    intersection= set(hla_kmers) & set(random_kmers)
    hla_kmers=list(set(hla_kmers) - intersection)
    random_kmers=list(set(random_kmers) - intersection)
    if len(hla_kmers) > 3*len(random_kmers):
        hla_kmers=hla_kmers[0:3*len(random_kmers)-1]
    if len(random_kmers) > 3*len(hla_kmers):
        random_kmers=random_kmers[0:3*len(hla_kmers)-1]
        #   print 'hla:random kmers: %s : %s' % (len(hla_kmers),len(random_kmers))
    # step2:allocate sample,label the sample.
    # train:test=4:1, 
    train_positive_index,test_positive_index=divide_listlike(hla_kmers,ratio)                                                                                                                                                                                                                                                                                                                                                                                                                 
    #print 'trian-p-i:test-p-i: %s : %s' % (len(train_positive_index),len(test_positive_index))
    train_negtive_index,test_negtive_index=divide_listlike(random_kmers,ratio)
    #print 'trian-n-i:test-n-i: %s : %s' % (len(train_negtive_index),len(test_negtive_index))
                                                                         
    train_data=[aa_translate(hla_kmers[i],trans_table) for i in range(len(hla_kmers)) if i in train_positive_index]
    train_labels=[1]*len(train_positive_index)
    train_data.extend([aa_translate(random_kmers[i],trans_table) for i in range(len(random_kmers)) if i in train_negtive_index])
    train_labels.extend([-1]*len(train_negtive_index))
    
    test_data=[aa_translate(hla_kmers[i],trans_table) for i in range(len(hla_kmers)) if i in test_positive_index]
    test_labels=[1]*len(test_positive_index)
    test_data.extend([aa_translate(random_kmers[i],trans_table) for i in range(len(random_kmers)) if i in test_negtive_index])
    test_labels.extend([-1]*len(test_negtive_index))
    #    print 'trian-label:test-labels: %s : %s' % (len(train_labels),len(test_labels))
    return [train_data,train_labels,test_data,test_labels]

def score_pssm():
    pass

# 'score_slide' is identify with 'slide',now it do nothing but return all the kmers.
# should add a 'score-function' soon. asenal!!!!
def score_slide(sequence,k):
    """ slide the sequence,return all the k-mer"""
    l=len(sequence)
    if l>=k:
        return [sequence[i:i+k] for i in range(l-k)]
    else:
        return []

def slide(sequence,k):
    """ slide the sequence,return all the k-mer"""
    l=len(sequence)
    if l>=k:
        return [sequence[i:i+k] for i in range(l-k)]
    else:
        return []

#--------- analyse the kmer composition 
def kmer_composition(L):
    """L is a container of sequence,call split-kmer to split each sequence
    in L into 1-mer,2-mer,3-mer..., analyse the comosition of each k-mer group."""

#--------- svm subs---
# check this !!!!!!!!! asenal
#def report_pram(filename,parameter_object):
#    with open(filename,'w') as file:
#        print "PRAMETERS:"
#        for i in eval("%s._names" % parameter_object):
#            value=eval("%s.%s" % (parameter_object,i))
#            file.write("%s :%s" % (i,value))
##-------------------- main loop --------------
__all__=['aa_table','bi_table','int_table','read_mhcpep','check_table','create_table','storeDB_mhcpep','split_kmer','aa_translate','divide_listlike','get_randomchunk','score_pssm','score_slide','trainset_convert']
