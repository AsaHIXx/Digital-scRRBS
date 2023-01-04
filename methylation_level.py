import pandas as pd
analysis_file = [31, 59, 61, 64, 69]

def calc_meth_new(df):
#     for i in range(len(df)):
#         if df['Met'].iloc[i] >=3:
#             df['Met'].iloc[i] = 1
#         else:
#             df['Met'].iloc[i] =0
    df = df.loc[df['Type']=='CpG']
    
    df = df.loc[(df['Met']+df['UnMet'])>=3]
    df.loc[df[df.MetRate <0.9].index.tolist(),'Met_stat'] = 0
    df.loc[df[df.MetRate>=0.9].index.tolist(),'Met_stat'] = 1
    #df[df['Met']<3] =0
    met_num = sum(df['Met_stat'])
    total = len(df)
    if total == 0:
        return 0
    methyl_level = (round((met_num/total)*100, 2))
    return methyl_level

def methylation_level_calculation_new(df):
    df_intergenic = df[df['GeneType']=='intergenic']
    df_exonic = df[(df['GeneType'] == 'exonic')|(df['GeneType'] == 'UTR5')|(df['GeneType']=='UTR3')|(df['GeneType']=='UTR5;UTR3')]
    df_intronic = df[(df['GeneType'] == 'intronic')]
    df_splicing = df[(df['GeneType'] == 'splicing')]
    df_updown = df[(df['GeneType'] == 'upstream')|(df['GeneType'] == 'downstream')|(df['GeneType'] =='upstream;downstream')]
    intergenic_rate = calc_meth_new(df_intergenic)
    exonic_rate = calc_meth_new(df_exonic)
    intronic_rate =calc_meth_new(df_intronic)
    splicing_rate = calc_meth_new(df_splicing)
    updown_rate = calc_meth_new(df_updown)
    return [intergenic_rate,exonic_rate,intronic_rate,splicing_rate,updown_rate]
    
def vcf_input_generate_allChr(df_merge):
    ID = ['.' for i in range(len(df_merge['#Chr']))]
    score = [100 for i in range(len(df_merge['#Chr']))]
    filte = ['PASS' for i in range(len(df_merge['#Chr']))]
    INFO = ['DP=4' for i in range(len(df_merge['#Chr']))]
    #Meth_stat = [0 for i in range(len(df_merge['#Chr']))]
    vcf_input = pd.DataFrame({'#CHROM':df_merge['#Chr'], 'POS':df_merge['Pos'], 'ID':ID, 'REF':df_merge['Ref'], 'ALT':df_merge['Ref'], 'QUAL':score, 'FILTER':filte, 'INFO':INFO})
    return vcf_input    
for i in analysis_file:
    df = pd.read_csv('DMF-'+str(i)+'.single', sep='\t')
    