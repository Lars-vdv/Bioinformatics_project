import Bio
from Bio import SeqIO                                                  
prefixes = []                                                          
suffixes = []                                                          
handle = open(r'test_data\rosalind_grph.txt', 'r')                                 
for record in SeqIO.parse(handle, 'fasta'):                            
    count1 = 0                                                         
    count2 = 0                                                         
    prefix = [record.id]                                               
    suffix = [record.id]                                               
    pre = ''                                                           
    suf = ''                                                           
    for nt in record.seq:                                              
        if count1 < 3:                                                 
            pre += nt                                                  
            count1 += 1                                                
    prefix.append(pre)                                                 
    for tn in reversed(record.seq):                                    
        if count2 < 3:                                                 
            suf += tn                                                  
            count2 += 1                                                
    suffix.append(''.join(reversed(suf)))                              
    prefixes.append(prefix)                                            
    suffixes.append(suffix)                                            
handle.close()                                                         
                                                                       
for i, k in enumerate(suffixes): 
    from utilities import *                                      
    currentsf = suffixes[i][1]                                         
    currentid = suffixes[i][0]                                         
    for j, l in enumerate(prefixes):                                   
        if currentsf == prefixes[j][1] and currentid != prefixes[j][0]:
            write_file("test.txt", data=f"{currentid} {prefixes[j][0]}", mode='a')
            print(currentid, prefixes[j][0])                           



#write_file("test.py", data= mode='a')