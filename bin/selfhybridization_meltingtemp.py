# Save combined masking FASTA file as a dictionary
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

# Define functions
def pair_type(b1: str, b2: str):
  if b1 == 'A' and b2 == 'T':
    return True
  if b1 == 'T' and b2 == 'A':
    return True
  if b1 == 'C' and b2 == 'G':
    return True
  if b1 == 'G' and b2 == 'C':
    return True
  return False

def mergeDictionary(dict_1, dict_2):
  dict_3 = {**dict_1, **dict_2}
  for key, value in dict_3.items():
    if key in dict_1 and key in dict_2:
      dict_3[key] = [value, dict_1[key]]
    return dict_3


# file name, e.g. data/ASO_20mer_hInfluenzae_AND.txt
aso_dict = {rec.id : rec.seq for rec in SeqIO.parse(path + input('file name = '), 'fasta')}

keys = list(aso_dict.keys())
SH_score = {}

# Loop 1: loop through ASO's
for i in keys:
  forward = aso_dict[i]
  reverse = aso_dict[i][::-1]
  max_shift = int(len(forward)/2)
  score_list = []
  pair_count = 0

  # Loop 2: loop through shifts
  for k in range(max_shift):
    count = 0
    f_temp = forward[k:] 
    r_temp = reverse[0:len(reverse)-k]

    # Loop 3: loop through positions in sequence
    for pos in range(len(f_temp)):
      pair_count += 1
      if pair_type(f_temp[pos], r_temp[pos]) == True:
        count += 1
    score_list.append(count)
  SH_score.update({i : max(score_list)})

keys = list(aso_dict.keys())
mt_score = {}
gc_score = {}

test = keys[:100]

for i in keys:
  gc_score.update({i : mt.Tm_GC(aso_dict[i], strict = False)})
  mt_score.update({i : mt.Tm_Wallace(aso_dict[i], strict = False)})


# Generate CSV file + Graphing
keys = list(aso_dict.keys())
deg_data = []
position_data = []
aso_data = []
sh_data = []
tm_data = []
gc_data = []

for i in keys:
  x = i.split('_')
  deg_data.append(x[0])
  position_data.append(x[1])
  aso_data.append(str(aso_dict[i]))
  sh_data.append(SH_score[i])
  tm_data.append(mt_score[i])
  gc_data.append(gc_score[i])

deg_df = pd.DataFrame(deg_data)
pos_df = pd.DataFrame(position_data)
aso_df = pd.DataFrame(aso_data)
sh_df = pd.DataFrame(sh_data)
tm_df = pd.DataFrame(tm_data)
gc_df = pd.DataFrame(gc_data)


df = pd.concat([deg_df, pos_df, aso_df, sh_df, tm_df, gc_df], axis = 1)
df.columns = ['DEG id', 'Position', 'ASO seq', 'SH score', 'TM score', 'GC score']

ofname = input("Name of output file e.g. hInfluenzae_long.csv :")
df.to_csv(ofname, index=False)

