import pandas as pd
from Bio.Seq import Seq
from diff_match_patch import diff_match_patch
from collections import Counter
from os import listdir


def aminoacidchangefunc(position, initialbp, mutatedbp):
    '''Show amino acid change from reference to query.'''
    mutposintriplet = position%3
    aastartpos = position - mutposintriplet
    originaltriplet = pereference[aastartpos:aastartpos+3]
    mutatedtriplet = pereference[aastartpos:aastartpos+mutposintriplet] + mutatedbp + pereference[aastartpos + mutposintriplet+1:aastartpos+3]
    # print(originaltriplet)
    # print(mutatedtriplet)
    
    initialaa = str(Seq(originaltriplet).translate())
    mutatedaa = str(Seq(mutatedtriplet).translate())
    aaposition = int(aastartpos/3)
    return initialaa+str(aaposition)+mutatedaa
    
    
def differencestrings(ref,query):
    '''Compare reference to query and return list with mutations on nucleotide basis.'''
    dif = diff_match_patch()
    differencelist = dif.diff_main(ref,query)  # reference vs query
    
    mutations = []
    aachange = []
    positionref = 0
    positionquery = 0
    for index, stretch in enumerate(differencelist):
        if stretch[0] == 0:
            positionref+=len(stretch[1])
            positionquery+=len(stretch[1])
        elif stretch[0] == -1: 
            if differencelist[index+1][0] == 1:  # if the following tuple is insertion, convert deletion/insertion to replacement
                initialbp = str(stretch[1])
                mutatedbp = str(differencelist[index+1][1])
                mutation = str(positionref)+'_'+initialbp+'/'+mutatedbp
                mutations.append((positionref, 'Replacement', str(stretch[1]), str(differencelist[index+1][1]), mutation))
                aminoacidchange = aminoacidchangefunc(positionref, initialbp, mutatedbp)
                aachange.append(aminoacidchange)
                print(aminoacidchange)
                positionref+=len(stretch[1])
                positionquery+=len(stretch[1])
            else: #deletion in query
                mutation = str(positionref)+'_'+str(stretch[1])+'/'+'-'
                mutations.append((positionref, 'Deletion', str(stretch[1]), '-', mutation))
                positionref+=len(stretch[1])
                
        elif stretch[0] == 1: #insertion in query
            if not differencelist[index-1][0] == -1: # skip if previous was deletion since then it was converted to replacement
                mutation = str(positionref)+'_'+'-'+'/'+str(stretch[1])
                mutations.append((positionref, 'Insertion', '-', str(stretch[1]),mutation))
    mutationshortlist = [x[4] for x in mutations]
    replacementonly = [x[4] for x in mutations if x[1] == 'Replacement']
    return mutations, mutationshortlist, replacementonly, aachange


def fasta2df(infile):
    '''Read .fasta file and store it as dataframe. Remove spaces in fasta filename first'''
    seq_df = pd.DataFrame()
    new_fasta_file = ''
    with open(infile) as fasta_file:
        for line in fasta_file.readlines():
            if '>' in line:
                text = line.replace(' ','')
                name = line.replace(' ','')
            else:
                text = line
                seq_df.loc[name,'sequence'] = text.strip()  # strip is needed to remove \n at end of sequence
            new_fasta_file = new_fasta_file+text
    
    return seq_df


referencepath = 'Define_Path'
def list_files1(path):
    return [f for f in listdir(path) if '.fa' in f]

referencelist = list_files1(referencepath)
for reference in referencelist:  # loop through all evolutions/domains
    referencename = reference
    folder1name = '_'.join(referencename.split('_')[:2])
    editorname = referencename.split('_')[-1][:-3]
    fastapath = referencepath+folder1name+'\\'+editorname+'\\'
    fastaname = editorname+".fa" 
    
    pereference = pd.read_csv(referencepath+referencename).iloc[0,0]
    pereferencelength = len(pereference)
    
    fasta = fasta2df(fastapath+fastaname)
    fasta['length'] = fasta.sequence.apply(lambda x: len(x))
    
    # remove reads which deviate more than x bp from reference sequence:
    maxlengthdeviation = 10
    fasta = fasta[(fasta.length <pereferencelength+maxlengthdeviation) & (fasta.length >pereferencelength-maxlengthdeviation)]  
    
    mutationdf = pd.DataFrame()
    mutationsummary = []
    aasummary = []
    for index, query in fasta.iterrows():
        try:
            mutations, mutationshortlist, replacementonly, aachange = differencestrings(pereference,query['sequence'])
        except IndexError:
            print(query['sequence'])
            print()
        mutationsummary+=mutationshortlist
        aasummary+=aachange
        mutationdf.loc[index,'replacementmutations'] = str(replacementonly)
        mutationdf.loc[index,'mutationshortlist'] = str(mutationshortlist)
        mutationdf.loc[index,'mutations'] = str(mutations)
        mutationdf.loc[index,'aachangereplacement'] = str(aachange)
        print(replacementonly)
        print()
    
    aacount = Counter(aasummary)
    aacountdf = pd.DataFrame.from_dict(aacount, orient='index').reset_index()
    
    aacombinationcount = Counter(mutationdf['aachangereplacement'])
    aacombinationcountdf = pd.DataFrame.from_dict(aacombinationcount, orient='index').reset_index()
    
    mutationcount = Counter(mutationsummary)
    mutationcountdf = pd.DataFrame.from_dict(mutationcount, orient='index').reset_index()
    
    mutationcombinationcount = Counter(mutationdf['mutationshortlist'])
    mutationcombinationcountdf = pd.DataFrame.from_dict(mutationcombinationcount, orient='index').reset_index()
    
    
    aacountdf.to_csv(fastapath+'\\Define_path'+referencename[:-3]+'.csv')
    aacombinationcountdf.to_csv(fastapath+'\\Define_path'+referencename[:-3]+'.csv')
    mutationdf.to_csv(fastapath+'\\Define_path'+referencename[:-3]+'.csv')
    mutationcountdf.to_csv(fastapath+'\\Define_path'+referencename[:-3]+'.csv')
    mutationcombinationcountdf.to_csv(fastapath+'\\Define_path'+referencename[:-3]+'.csv')
