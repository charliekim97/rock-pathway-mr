import pandas as pd
b=pd.read_csv("/sessions/determined-eager-davinci/mnt/outputs/rock2_mr/rock2_eqtlgen_converted.txt",sep="\t")
lead=b.iloc[0]; pk=b[b.SNP=="rs7581184"].iloc[0]
print("file rows blood:",len(b))
