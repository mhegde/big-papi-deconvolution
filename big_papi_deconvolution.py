import pandas as pd
import csv, argparse, os
from itertools import izip
import gzip

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--read-file',
        type=str,
        help='Gzipped FASTQ file with construct reads')
    parser.add_argument('--barcode-file',
        type=str,
        help='Gzipped FASTQ file with sample barcodes')
    parser.add_argument('--ref-u6',
        type=str,
        help='Reference file for U6 guides;.csv file with no headers')
    parser.add_argument('--ref-h1',
        type=str,
        help='Reference file for H1 guides;.csv file with no headers') 
    parser.add_argument('--cond-file',
        type=str,
        help='Sample conditions file;.csv file with no headers')
    parser.add_argument('--outputfile',
        type=str,
        help='Path to outputfile')
    return parser

def get_ref_hash(ref_df):
    ref_hash = {}
    sp = list(ref_df['Construct Barcode'])
    for i,r in enumerate(sp):
        ref_hash[r] = ref_df['Construct IDs'][i]
    return ref_hash

def get_cond_hash(bcs_df):
    cond_hash = {}
    for i,b in enumerate(bcs_df['Barcode']):
        cond_hash[b]=bcs_df['Conditions'][i]+'_'+b
    return cond_hash

def get_count_hash(cond, ref):
    count_hash = {}
    for i,r in enumerate(ref):
        count_hash[r] = count_hash.get(r,{})
        for c in cond:
            count_hash[r][c] = 0
    return count_hash

def get_comb_ref(ref_u6,ref_h1):
    comb_ref = {}
    for i,r in ref_u6.iterrows():
        for j,r1 in ref_h1.iterrows():    
            id_new = r['Construct IDs']+';'+r1['Construct IDs']
            new_const = r['Construct Barcode']+';'+r1['Construct Barcode']
            comb_ref[new_const] = id_new
    return comb_ref

def get_output(output_hash, cond_hash, ref_hash):
    final_df=pd.DataFrame(output_hash).T.fillna(0)
    colnames = final_df.columns
    new_colnames = [cond_hash[x] for x in colnames]
    indexes = list(final_df.index)
    ref_col = [ref_hash[x] for x in indexes]
    final_df.columns = new_colnames
    final_df = final_df[sorted(final_df.columns)]
    final_df.insert(0, 'Reference', ref_col)
    return final_df      

if __name__ == '__main__':
    args = get_parser().parse_args()
    read_file = args.read_file
    barcode_file = args.barcode_file
    ref_u6 = pd.read_csv(args.ref_u6,header=None)
    ref_u6.columns = ['Construct Barcode','Construct IDs']
    ref_u6_hash = get_ref_hash(ref_u6)
    ref_u6_len = len(ref_u6_hash.keys()[0])
    ref_h1 = pd.read_csv(args.ref_h1,header=None)
    ref_h1.columns = ['Construct Barcode','Construct IDs']
    ref_h1_hash = get_ref_hash(ref_h1)
    ref_h1_len = len(ref_h1_hash.keys()[0])
    cond = pd.read_csv(args.cond_file,header=None)
    cond.columns = ['Barcode','Conditions']
    cond_hash = get_cond_hash(cond)
    comb_hash = get_comb_ref(ref_u6,ref_h1)
    count_hash = get_count_hash(cond_hash.keys(),comb_hash.keys())
    fixed_dist = 194
    search_str_1 = 'CACCG'
    search_str_2 = 'CTTAAAC'
    barcode_count=0
    search_str_1_count=0
    search_str_2_count=0
    u6_guide_count=0
    h1_guide_count=0
    count=-1
    with gzip.open(read_file,'r') as r, gzip.open(barcode_file,'r') as b:
        for read, bar in izip(r,b):
            count += 1
            if count%4 == 1:
                bar = bar[:-1]
                if bar in cond_hash.keys():
                    print 'Barcode found'
                    barcode_count+=1
                    print count/4
                    pos = read.find(search_str_1)
                    if pos != -1:
                        search_str_1_count+=1
                        u6_guide = read[pos+5:(pos+ref_u6_len+5)]
                        if u6_guide in ref_u6_hash.keys():
                            print 'U6 guide found'
                            u6_guide_count+=1
                            h1_pos = pos+25+fixed_dist
                            h1_guide = read[h1_pos:(h1_pos+ref_h1_len)]
                            pos2 = read.find(search_str_2)
                            h1_guide_test = read[pos2+7:(pos2+ref_h1_len+7)]
                            if h1_guide == h1_guide_test:
                                search_str_2_count+=1
                                if h1_guide in ref_h1_hash.keys():
                                    print 'H1 guide found '
                                    h1_guide_count+=1
                                    new_key = u6_guide+';'+h1_guide
                                    count_hash[new_key][bar]+=1
    print 'Total number of reads:'+str(count/4)
    print 'Number of barcode matches:'+str(barcode_count)
    print 'Number of CACCG matches:'+str(search_str_1_count)
    print 'Number of CTTAAAC matches:'+str(search_str_2_count)
    print 'Number of u6 guide matches:'+str(u6_guide_count)
    print 'Number of h1 guide matches:'+str(h1_guide_count)

    count_df = get_output(count_hash,cond_hash,comb_hash)
    count_df.to_csv(args.outputfile,sep='\t')

