pip install biopython

from Bio import Entrez
import os

Entrez.email = "pg55697@alunos.uminho.pt"

# Lista dos genes e os seus IDs correspondentes para o nucleótido e proteína
genes = ["dinB", "dnaA", "FE838_RS15320", "susG"]
nucleotide_id = ["NZ_CP040530.1", "NZ_CP040530.1", "NZ_CP040530.1", "NZ_CP040530.1"] 
protein_id = ["WP_008760075.1", "WP_008759760.1", "WP_011107229.1", "WP_011108937.1"]    

for i in range(len(genes)):
    if nucleotide_id[i] == "NA":
        print(f"No nucleotide information available for {genes[i]}.")
        continue 
    else:
        handle = Entrez.efetch(db="nucleotide", id=nucleotide_id[i], rettype="gb", retmode="text")
        with open(f'Informação_genes/{genes[i]}.gb', "w") as doc:
            doc.write(handle.read()) 
            handle.close()

for i in range(len(genes)):
        handle = Entrez.efetch(db="protein", id=protein_id[i], rettype="gb", retmode="text")
        
        with open(f'Informação_genes/{genes[i]}_protein.gb', "w") as doc:
            doc.write(handle.read())
    finally:
        handle.close()
        
#Anotação da proteína de cada gene
from Bio import SeqIO
from Bio import SeqFeature

for gene in genes:
    f=open(f'{gene}.txt','a')
    f.writelines('\n***INFORMACAO RELATIVA A PROTEINA***\n\n')
    record_prot=SeqIO.read(f"Informação_genes\{gene}_protein.gb",'genbank')
    id=record_prot.name
    seq=record_prot.seq
    tamanho=len(record_prot.seq)
    des=record_prot.description
    tipo_molecula=record_prot.annotations["molecule_type"]
    features=record_prot.features
    f.writelines(f'Tipo de molecula: {tipo_molecula}\nAccession number proteina: {id}\nTamanho proteina: {len(seq)} aa\n')
    f.writelines(f'Descricao da proteina: {des}\n')
    feat_site=[]
    for I,feature in enumerate(features):
        if feature.type=='Protein':
            peso_molecualr=feature.qualifiers["calculated_mol_wt"][0]
            f.writelines(f'Peso molecular: {peso_molecualr} Dalton\nNumero de features da proteina: {len(features)}\n')
        if feature.type=='Site':    
            feat_site.append(I)

    if len(feat_site)!=0:
        f.writelines('\nLocais de interesse da proteina:\n\n')
        for indice in feat_site:
            feature=features[indice]        
            localizacao=feature.location
            f.writelines(f'Site {indice+1}:\n \t-Localizacao: {localizacao}')
    
    f.writelines(f'\nSequencia proteina:\n')
    for I in range(0,len(seq),100):
        seq_=seq[I:I+100]
        f.writelines(f"{seq_}\n")

    f.close()

