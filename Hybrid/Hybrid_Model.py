import sys
import pickle
import os
import pandas as pd
import numpy as np
vcf_header = ['CHROM',  'POS', 'ID',  'REF' , 'ALT' ,'QUAL' ,'FILTER' ,'INFO' , 'FORMAT' , 'SAMPLE']


def hybrid(gatk_samp, dv_samp, file_name, output_file):
    #os.system('ml anaconda; conda activate bio')
    command = "ml anaconda; conda activate bio;rtg vcfeval -b {} -c {} --all-records --squash-ploidy --no-roc -t /scratch4/jluo2/drabiza1/ICI_data/Manuscript_Exome/hg19_SDF/ -o {}".format(gatk_samp, dv_samp, file_name)
    os.system(command)
    #import each of the fp and fn with the existing function (slighlt modified), run each model. Keep the ones that are predicted to be true. 
    gatk = Format_GATK(file_name+'/fn.vcf.gz')
    dv = Format_DV(file_name+'/fp.vcf.gz')
    both = pd.read_csv(file_name+'/tp-baseline.vcf.gz', sep = '\t', comment= '#', header = None)
    all_df = pd.concat([both,gatk,dv ])
    all_df.columns = vcf_header
    all_df = all_df.sort_values(['CHROM', 'POS'])
    
    head = pd.read_csv('FULL_GATK_Header.txt', sep ='\n', header = None)
    head = '\n'.join(list(head[0]))
    head = head.replace('46A', file_name)
    output_VCF = output_file+'.vcf'

    with open(output_VCF, 'w') as vcf:
        vcf.write(head+'\n')
    all_df.to_csv(output_VCF, sep="\t", mode='a', index=False, header = None)
    os.system('rm -r '+file_name)

def Format_DV(path):
    with open('DV_RF.pkl', 'rb') as f: 
        DV_clf = pickle.load(f)
    df = pd.read_csv(path, sep = '\t', comment= '#', header = None)
    Train = pd.DataFrame()
    Train['SNP_INDEL'] = df[3].apply(len) - df[4].apply(len)
    Train['Qual'] = df[5]
    Train['AQ'] = df[7].str.split(';', expand=True)[1].str.split('=', expand=True)[1]
    Train = pd.concat([Train,df[9].str.split(':', expand=True)[[0,1,3]].rename(columns={0:'GT',
                                                                                        1:'DP', 3:'GQ'}) ], axis =1)
    Train = pd.concat([Train,df[9].str.split(':', expand=True)[2].str.split(',', 
                                                                     expand=True).rename(columns={0:'RD', 1:'AD'})],
                       axis =1)
    Train['GT'] = Train['GT'].replace('0/1', 0).replace('0|1', 0).replace('1/1',
                    1).replace('1|1', 1).replace('1/0', 0).replace('1|0', 0).replace('./1', 1 )
    Train = Train.replace('.', 0)
    Train = Train.astype(float)
    Train['VAF'] = Train['AD']/Train['DP']
    Train['Caller'] = [1 for i in range(Train.shape[0])]
    Train = Train.fillna(0)
    results = DV_clf.predict(Train.to_numpy() )
    return df.loc[np.where(results==0)[0]]

def Format_GATK(path):
    with open('GATK_RF.pkl', 'rb') as f:
        GATK_clf = pickle.load(f)
    df = pd.read_csv(path, sep = '\t', comment= '#', header = None)
    Train = pd.DataFrame()
    Train['SNP_INDEL'] = df[3].apply(len) - df[4].apply(len)
    status = df[6]=='PASS'
    Train['Filter']= status.replace(True,0).replace(False,1)
    status = df[7].str.contains('NEGATIVE')
    Train['Negative'] = status.replace(True,0).replace(False,1)
    status = df[7].str.contains('POSITIVE')
    Train['Positive'] = status.replace(True,0).replace(False,1)
    Train['Qual'] = df[5]
    Train['SOR'] = df[7].apply(get_sor)
    Train['FS'] = df[7].apply(get_FS)
    Train['QD'] = df[7].apply(get_QD)
    Train['FS'] = df[7].apply(get_FS)
    Train['MQRANK'] = df[7].apply(get_MQRankSum)
    Train = pd.concat([Train,df[9].str.split(':', expand=True)[[0,2,3]].rename(columns={0:'GT',
                                                                                        2:'DP', 3:'GQ'}) ], axis =1)
    Train = pd.concat([Train,df[9].str.split(':', expand=True)[[1]][1].str.split(',',
                                                                     expand=True).rename(columns={0:'RD', 1:'AD'})],axis =1)
    Train['GT'] = Train['GT'].replace('0/1', 0).replace('0|1', 0).replace('1/1',
                    1).replace('1|1', 1).replace('1/0', 0).replace('1|0', '0')
    Train = Train.astype(float)
    Train['VAF'] = Train['AD']/Train['DP']
    Train['Caller'] = [1 for i in range(Train.shape[0])]
    Train = Train.fillna(0)
    #predict
    results = GATK_clf.predict(Train.to_numpy() )
    return df.loc[np.where(results==0)[0]]



def get_sor(x):
    return(float(x.split('SOR=')[1].split(';')[0]))
def get_FS(x):
    return(float(x.split('FS=')[1].split(';')[0]))
def get_MQ(x):
    return(float(x.split('MQ=')[1].split(';')[0]))
def get_QD(x):
    try:
        return(float(x.split('QD=')[1].split(';')[0]))
    except:
        return 0
def get_MQRankSum(x):
    try:
        return(float(x.split('MQRankSum=')[1].split(';')[0]))
    except:
        return 0
def get_GQ(x):
    try:
        return int(x.split(':')[-2])
    except:
        return 0
if __name__ == '__main__':
    # Map command line arguments to function arguments.
    hybrid(*sys.argv[1:])
