from Bio import SeqIO

input_file="N.meningitidis_MC58.embl"
output_file="N.meningitidis_MC58.gbk"

count = SeqIO.convert(input_file, "embl", output_file, "genbank")
print(f"Convertidos {count} registros de EMBL a GenBank.")