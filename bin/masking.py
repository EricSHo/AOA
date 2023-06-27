# Import packages
import pandas as pd
import re
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

def convert(s):
    str1 = ''
    return(str1.join(s))
  
def rev_comp_st(seq): # Create reverse compliment
    seq = seq.replace('A', 't').replace('C', 'g').replace('T', 'a').replace('G', 'c')
    seq = seq.upper()
  
    # reverse strand 
    seq = seq[::-1]
  
    return seq

def create_kmers(s, k): # Create kmers
    kFreq = {}
  
    for i in range(0, len(s) - k+1):
      kmer = s[i:i+k]
  
      if kmer in kFreq:
        kFreq[kmer] += 1
      else:
        kFreq[kmer] = 1
    return kFreq

# Assume at the root directory AOA
path = "./data"

for bacterium in ['hInfluenzae','sAuresus','pGingivalis']:

  # STEP 1: Set up and read files
  # Set up
  sequence_RNAfold = [] # sequences from RNAfold
  deg_ids_RNAfold = []  # deg_ids from RNAfold
  sequence_mxfold = []  # sequences from mxfold
  deg_ids_mxfold = []   # deg_ids from mxfold
  rna_file = f'{path}/{bacterium}_eg_RNAfold.txt'
  mxfold_file = f'{path}/{bacterium}_eg_mxfold2.txt'
  
  # Read fasta file and store sequences in list 
  for seq_record in SeqIO.parse(open(rna_file),'fasta'): # RNAfold
    # Add record to list
    sequence_RNAfold.append(str(seq_record.seq))
    deg_ids_RNAfold.append(seq_record.id)
  
  for seq_record in SeqIO.parse(open(mxfold_file),'fasta'): # mxfold
    sequence_mxfold.append(str(seq_record.seq))  # sequences from mxfold
    deg_ids_mxfold.append(seq_record.id)
  
  # STEP 2: Separate lines by delimiter
  # Set up
  input_list_RNAfold = [] # Store input sequences RNAfold
  prediction_list_RNAfold = [] # Store secondary predicition RNAfold
  input_list_mxfold = [] # Store input sequences mxfold
  prediction_list_mxfold = [] # Store input sequences mxfold
  
  for i in range(len(sequence_RNAfold)): 
    count = 0
    s = sequence_RNAfold[i]
    if count == 0:
      input_seq = re.split('[.()]', s) # get input sequence
      input_list_RNAfold.append(input_seq[0])
  
      prediction = s[len(input_seq[0]):2*len(input_seq[0])] # get secondary prediction
      prediction_list_RNAfold.append(prediction)
  
      count += 1
  
  for i in range(len(sequence_mxfold)): 
    count = 0
    s = sequence_mxfold[i]
    if count == 0:
      input_seq = re.split('[.()]', s) # get input sequence
      input_list_mxfold.append(input_seq[0])
  
      prediction = s[len(input_seq[0]):2*len(input_seq[0])] # get secondary prediction
      prediction_list_mxfold.append(prediction)
  
      count += 1
  
  # STEP 3: Masking     FIX
  # Set up 
  masked_seq_RNAfold = [] # Store masked sequences
  masked_seq_mxfold = []
  
  for i in range(len(input_list_RNAfold)): 
    s_RNAfold = input_list_RNAfold[i].lower() # Everything starts off masked (aka lower case)
    p_RNAfold = prediction_list_RNAfold[i]
  
    # String to list
    list_sequence = list(s_RNAfold)
  
    # Inner for loop (to iterate through prediction list)
    for k in range(len(s_RNAfold)):
      if p_RNAfold[k] == '.':
        unmask = list_sequence[k].upper()
        list_sequence[k] = unmask
  
    # Append to masked sequence list, convert list back to string
    masked_seq_RNAfold.append(''.join(list_sequence))
  
  for i in range(len(input_list_mxfold)): 
    # Set up
    s_mxfold = input_list_mxfold[i].lower()
    p_mxfold = prediction_list_mxfold[i]
  
    # String to list
    list_sequence = list(s_mxfold)
  
    # Inner for loop (to iterate through prediction list)
    for k in range(len(s_mxfold)):
      if p_mxfold[k] == '.':
        unmask = list_sequence[k].upper()
        list_sequence[k] = unmask
  
    # Append to masked sequence list, convert list back to string
    masked_seq_mxfold.append(''.join(list_sequence))
  
  print(len(masked_seq_RNAfold))
  print(len(masked_seq_mxfold))
  
  # Write to file
  # File for RNAfold
  ofname = f"{path}/{bacterium}_eg_masked_RNAfold.txt"
  ofile = open(ofname, "w")
  
  for i in range(len(masked_seq_RNAfold)):
    ofile.write(">" + deg_ids_RNAfold[i] + "\n" +masked_seq_RNAfold[i] + "\n")
  
  ofile.close()
  
  # File for mxfold
  ofname = f"{path}/{bacterium}_eg_masked_mxfold.txt"
  ofile = open(ofname, "w")
  
  for i in range(len(masked_seq_mxfold)):
    ofile.write(">" + deg_ids_mxfold[i] + "\n" +masked_seq_mxfold[i] + "\n")
  
  ofile.close()
 
  ## Homology masking for humans
  # Step 1: Set up and read files
  human_file = f'{path}/{bacterium}_vs_humans.txt'
  
  blast_report = pd.read_csv(human_file, sep = '\t')
  blast_report.columns = ['qseqid','qlen','qstart','qend','sseqid','sstart','send','qseq','sseq','evalue','pident','nident','sstrand','length','stitle','staxid','ssciname','scomname','sskingdom']
  
  # A) Select rows that satisfy requirements (Subject (sstrand) is plus, and Subject sequence is RefSeqGene or accession prefix is NG_)
  ## Create empty dataframe to store report_new
  header = list(blast_report.columns)
  
  ## Scan 1: Scan for sstrand == plus
  temp = blast_report[blast_report['sstrand'] == 'plus'] # 25 fulfill requriement
  #print(temp)
  
  ## Scan 2: Scan for accession prefix is NG_ (under sseqid)
  filter = temp['sseqid'].str.contains('NG_') # 3 fulfill requirement (index = 7, 12, 28)
  report_dictionary = temp[filter]
  
  ### Generate dictionary 
  narrow_df = report_dictionary[['qseqid','qseq']]
  qseq_dict = {} #key = id (all unique), value = qseq (list)
  
  for id in range(len(narrow_df)): 
    temp_id = narrow_df['qseqid'].values[id]
    seq = narrow_df['qseq'].values[id]
    if temp_id not in qseq_dict:
      seq_list = []
      qseq_dict.update({temp_id : [seq]})
      continue
      qseq_dict[temp_id].append(seq)
  
  # Step 2: Read fasta file and store sequences in dictionary
  seq_dict_h = {}
  
  for seq_record in SeqIO.parse(open(f"{path}/{bacterium}_eg_original.fasta"),'fasta'):
    # Add record to dictionary
    seq_dict_h.update({seq_record.id : str(seq_record.seq)})
  
  # Step 3: Masking
  for hdr, seq in seq_dict_h.items(): # Iterate through original gene sequence
    gene_id = hdr[0:] # 1 through end
    if gene_id in qseq_dict:
      alignment_list = qseq_dict[gene_id]
      for alignment in alignment_list: # Find alignments 
        q_part = alignment.split('-') 
        pos_5 = seq.find(q_part[0]) # 5' end
        pos_3 = seq.find(q_part[-1]) # 3' end
        masked_range = pos_3 + len(q_part[-1])
        # Masking 
        homology = seq[pos_5 : masked_range]
        replace = re.compile(re.escape(homology), re.IGNORECASE)
        seq_dict_h[gene_id] = replace.sub(homology.lower(), seq_dict_h[gene_id])
          
  # Write to file 
  # File for humans
  headers_human = list(seq_dict_h.keys())
  sequences = []
  
  for i in seq_dict_h: # Fill sequence list
    sequences.append(seq_dict_h[i])
  
  ofname = f"{path}/{bacterium}_eg_masked_human.txt"
  ofile = open(ofname, "w")
  
  for i in range(len(seq_dict_h)):
    ofile.write(">" + headers_human[i] + "\n" +sequences[i] + "\n")
  
  ofile.close()   

  ## Homology masking for non-humans
  # Step 1: Set up and read files
  nh_file = f'{path}/{bacterium}_vs_nonhumans.txt'
  
  blast_report = pd.read_csv(nh_file, sep = '\t')
  blast_report.columns = ['qseqid','qlen','qstart','qend','sseqid','sstart','send','qseq','sseq',	'evalue','pident','nident','sstrand','length','stitle','staxid','ssciname','scomname','sskingdom']
  
  # A) Select rows that satisfy requirements (Subject (sstrand) is plus, and Subject sequence is RefSeqGene or accession prefix is NG_)
  ## Create empty dataframe to store report_new
  header = list(blast_report.columns) # Store headers
  
  ## Scan 1: Scan for sstrand == plus 
  temp = blast_report[blast_report['sstrand']=='plus']
  
  ## Scan 2: Scan to exclude different strains of P. gingivalis (under header stitle)
  filter = ~temp['stitle'].str.contains('influenzae|Haemophilus')
  report_dictionary = temp[filter]
  
  ### Generate dictionary 
  narrow_df = report_dictionary[['qseqid','qseq']]
  qseq_dict = {} #key = id (all unique), value = qseq (list)
  
  for id in range(len(narrow_df)): 
    temp_id = narrow_df['qseqid'].values[id]
    seq = narrow_df['qseq'].values[id]
    if temp_id not in qseq_dict:
      seq_list = []
      qseq_dict.update({temp_id : [seq]})
      continue
      qseq_dict[temp_id].append(seq)
  
  # Step 2: Read fasta file and store sequences in dictionary
  seq_dict_nh = {}
  
  fastaf = f"{path}/{bacterium}_eg_original.fasta"
  for seq_record in SeqIO.parse(open(fastaf),'fasta'):
    # Add record to dictionary
    seq_dict_nh.update({seq_record.id : str(seq_record.seq)})
  
  # Step 3: Masking
  for hdr, seq in seq_dict_nh.items(): # Iterate through original gene sequence
    gene_id = hdr[0:] # 1 through end
    if gene_id in qseq_dict:
      alignment_list = qseq_dict[gene_id]
      for alignment in alignment_list: # Find alignments 
        q_part = alignment.split('-') 
        pos_5 = seq.find(q_part[0]) # 5' end
        pos_3 = seq.find(q_part[-1]) # 3' end
        masked_range = pos_3 + len(q_part[-1])
        # Masking 
        homology = seq[pos_5 : masked_range]
        replace = re.compile(re.escape(homology), re.IGNORECASE)
        seq_dict_nh[gene_id] = replace.sub(homology.lower(), seq_dict_nh[gene_id])

  # File for nonhumans
  headers_nh = list(seq_dict_nh.keys())
  sequences = []
  
  for i in seq_dict_nh: # Fill sequence list
    sequences.append(seq_dict_nh[i])
  
  ofname = f"{path}/{bacterium}_eg_masked_nonhuman.txt"
  ofile = open(ofname, "w")
  
  for i in range(len(headers_nh)):
    ofile.write(">" + headers_nh[i] + "\n" +sequences[i] + "\n")
  
  ofile.close()    

  ## Combine masking
  # Set up
  sequence_RNAfold = [] # sequences from RNAfold
  deg_ids_RNAfold = []  # deg_ids from RNAfold
  
  # Read fasta file and store sequences in list 
  for seq_record in SeqIO.parse(open(f'{path}/{bacterium}_eg_masked_RNAfold.txt'),'fasta'): # RNAfold
    # Add record to list
    sequence_RNAfold.append(str(seq_record.seq))
    deg_ids_RNAfold.append(seq_record.id)
  
  # 2 sets of scans (1st, turn U -> T, 2nd turn u -> t)
  temp = [] # temporary list to store sequence
  converted_RNAfold = [] # store converted sequence (into DNA)
  for i in range(len(sequence_RNAfold)):
    x = sequence_RNAfold[i]
    new_sequence = x.replace('U', 'T')
    temp.append(new_sequence) 
  
  for i in range(len(temp)):
    x = temp[i]
    new_sequence = x.replace('u', 't')
    converted_RNAfold.append(new_sequence)
  
  print(deg_ids_RNAfold)
  
  # Create dictionaries and functions
  all_genes = {}
  rnafold = {}
  mxfold = {}
  blast_hs = {}
  blast_nonhs = {}
  
  all_genes = {rec.id : rec.seq for rec in SeqIO.parse(open(f'{path}/{bacterium}_eg_original.fasta'),'fasta')}
  mxfold = {rec.id : rec.seq for rec in SeqIO.parse(open(f'{path}/{bacterium}_eg_masked_mxfold.txt'),'fasta')}
  blast_hs = {rec.id : rec.seq for rec in SeqIO.parse(open(f'{path}/{bacterium}_eg_masked_human.txt'),'fasta')}
  blast_nonhs = {rec.id : rec.seq for rec in SeqIO.parse(open(f'{path}/{bacterium}_eg_masked_nonhuman.txt'),'fasta')}
  if (len(rnafold) == 0):
    for key in deg_ids_RNAfold:
        for value in converted_RNAfold:
            rnafold[key] = value
            converted_RNAfold.remove(value)
            break
 
  # Step 1: Combine RNAfold and mxfold results (AND)
    # Add missing genes into mxfold dictionary. Use "." in place of bases (len = len(gene))
  
    # Find difference (aka deg_ids that are not in mxfold)
    # Compare everything to all_genes, if not in rnafold or mxfold then directly copy seq
  gene = set(all_genes) - set(mxfold)
  
  for i in gene:
    gene_len = len(all_genes[i])
    sec_str = mxfold.get(i, "."*gene_len)
    mxfold.update({i : sec_str})
  
  
  # Compare sequences and make combined secondary structure dictionary (temp_dict)
  temp_dict = {}
  keys = list(all_genes.keys())
  
  for i in keys:
    temp_rnafold = rnafold[i]
    temp_mxfold = mxfold[i]
    deg_seq = all_genes[i]
    temp_seq = []
  
    if '.' in temp_mxfold: # Copy rnafold sequence if not in mxfold (represented as '.')
      #temp_dict.update({i: temp_rnafold})
      temp_dict[i] = temp_rnafold
      continue
  
    else:
      for pos in range(len(deg_seq)):
        if temp_rnafold[pos] in 'atgc' and temp_mxfold[pos] in 'atgc':
          temp_seq += temp_rnafold[pos]
          continue
        temp_seq += deg_seq[pos]
    temp_dict.update({i : convert(temp_seq)})

  # Step 2: Combine human and nonhuman homologies
  temp_dict_blast = {}
  keys = list(all_genes.keys())
  
  for i in keys:
    temp_hs = blast_hs[i]
    temp_nonhs = blast_nonhs[i]
    deg_seq = all_genes[i]
    temp_seq = []
  
    for pos in range(len(deg_seq)):
      if temp_hs[pos] in 'atgc':
        temp_seq += temp_hs[pos]
        continue
      if temp_nonhs[pos] in 'atgc':
        temp_seq += temp_nonhs[pos]
        continue
      temp_seq += deg_seq[pos] 
  
    temp_dict_blast.update({i : convert(temp_seq)})
  
  # Step 3: Merge the two to form a final file
  combined_mask_final = {}
  keys = list(all_genes.keys())
  
  for i in keys:
    temp_blast = temp_dict_blast[i]
    temp_structure = temp_dict[i]
    deg_seq = all_genes[i]
    temp_seq = []
  
    for pos in range(len(deg_seq)):
      if temp_blast[pos] in 'atgc':
        temp_seq += temp_blast[pos]
        continue
      if temp_structure[pos] in 'atgc':
        temp_seq += temp_structure[pos]
        continue
      temp_seq += deg_seq[pos]  
    combined_mask_final.update({i : convert(temp_seq)})

  # Write to fasta file
  headers = list(all_genes.keys())
  sequences = []
  
  for i in combined_mask_final: # Fill sequence list
    sequences.append(combined_mask_final[i])
  
  ofile = open(f'{path}/{bacterium}_eg_combined_masked_AND.txt', "w")
  
  for i in range(len(combined_mask_final)):
    ofile.write(">" + headers[i] + "\n" +sequences[i] + "\n")
  
  ofile.close()

  # Step 1: Save combined masking FASTA file as a dictionary
  seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(f'{path}/{bacterium}_eg_combined_masked_AND.txt'), 'fasta')} # file name = (baceria)_eg_combined_masked_AND.txt
  
  # Step 3: Generate kmers that only contain unmasked nucleotides 
  
  #k = input("k = ") # Values for k are user input, values that will be used are 19, 20, 21
  k = 20
  keys = list(seq_dict.keys())
  temp = []
  test_dict = {}
  genes = []
  
  for i in keys:
    seq_in = str(seq_dict[i])
    x = create_kmers(seq_in, int(k))
    temp.append(x)
  
  #print(temp)
  
  for i in range(len(temp)): # Outer loop to iterate through list of dictionaries
    temp_dict = temp[i]
    seq = list(temp_dict.keys())
  
    for k in range(len(seq)): # Inner loop to check if kmer is all unmasked (kmers are in the dictionaries)
      #print(seq[k])
      is_all_unmasked = seq[k].isupper() and seq[k].isalpha()
      if is_all_unmasked == True:
        header = keys[i] + '_'+str(list(temp_dict).index(seq[k]))
        #print(header)
        #print(seq[k])
        test_dict.update({header : seq[k]})
        genes.append(keys[i])
        
  print("number of kmers: " + str(len(test_dict)))
  
  # Get unique ids
  unique_genes = list(set(genes))
  
  print("number of unique genes: " + str(len(unique_genes)))
  
  # Step 4: Create dictionary with reverse compliment
  keys = list(test_dict.keys())
  reverse_comp = {}
  
  for i in keys:
    forward = test_dict[i]
  
    reverse = rev_comp_st(forward)
    reverse_comp.update({i : reverse})       

  # Write to fasta file
  headers = list(reverse_comp.keys())
  sequences = []
  
  for i in reverse_comp: # Fill sequence list
    sequences.append(reverse_comp[i])
  
  ofile = open(f'{path}/ASO_20mer_{bacterium}_AND.txt', "w")
  
  for i in range(len(reverse_comp)):
    ofile.write(">" + headers[i] + "\n" +sequences[i] + "\n")
  
  ofile.close()
  
