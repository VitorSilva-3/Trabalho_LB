from Bio.Blast import NCBIXML 
from Bio.Blast import NCBIWWW 
from Bio import SeqIO
NCBIWWW.email = "PG55697@alunos.uminho.pt"
NCBIXML.email = "PG55697@alunos.uminho.pt"

genes=["dinB","dnaA","FE838_RS15320","susG"] # lista de genes em estudo

#O ciclo for realiza o blast para as proteínas codificadas pelos genes em estudo
#Foi utilizada a base de dados UniprotKB para a realização do blastp e para restringir a procura apenas a proteinas conhecidas e com as suas caracteristicas estudadas

for gene in genes:
    record=SeqIO.read(f'Informação_genes\{gene}_protein.gb','genbank')
    sequencia=record.seq
    result=NCBIWWW.qblast('blastp','swissprot',sequencia,hitlist_size=15,expect=6)
    save_result=open(f"{gene}_blast.txt","w")
    save_result.writelines(result.read())
    save_result.close()
    result.close()

#Escolha dos 6 melhores resultados
for gene in genes:
    record=SeqIO.read(f'Informação_genes/{gene}_protein.gb','genbank')
    
    id=record.id
    
    f=open(f'{gene}_blast_results.txt','w')

    result_handle=open(f'{gene}_blast.txt','r')
    record_blast=NCBIXML.read(result_handle)
    result_handle.close()

    matriz=record_blast.matrix
    gap_pen=record_blast.gap_penalties
    database=record_blast.database
    n_hits=len(record_blast.alignments)
    f.writelines('***RESULTADOS DO BLAST PARA A PROTEINA CODIFICADA PELO GENE {}***\n\n'.format(gene))
    f.write(f"Accession number (NCBI) da proteina codificada pelo gene {gene}: {id}\nNumero de hits: {n_hits}\n")
    f.writelines(f'Matriz usada no alinhamento: {matriz}\nPenalidades de espacamentos (abertura,propagacao): {gap_pen}\nBase de proteinas usada: {database}\n') 

    alinhamentos=record_blast.alignments

    if len(alinhamentos)>=5:
        f.write('\n---MELHORES ALINHAMENTOS---\n\n')
        for I in range(0,6):
                if I==0:
                    alinhamento=alinhamentos[I]
                    accession_match=alinhamento.accession 
                    hit_def=alinhamento.hit_def
                    f.write(f'Alinhameto {I+1}: Proteina codifica pelo gene {gene}\n')
                    f.writelines(f'\tAccession number (relativo a base de dados usado no blast): {accession_match}\n\tDefinicao: {hit_def}\n\n')

                else:
                    alinhamento=alinhamentos[I]
                    accession_match=alinhamento.accession 
                    hit_def=alinhamento.hit_def 
                    hit_len=alinhamento.length
                    f.write(f'Alinhamento {I+1}\n')
                    f.writelines(f'\tAccession number hit (relativo a base de dados usado no blast): {accession_match}\n\tDefinicao: {hit_def}\n\tTamanho do hit: {hit_len} aa\n')
                    hsps=alinhamento.hsps
                    f.write(f'\tNumero de hsps: {len(hsps)}\n')
                    for N,hsp in enumerate(hsps):
                        E_value=hsp.expect
                        Score=hsp.score
                        Length= hsp.align_length
                        identities= hsp.identities
                        positives= hsp.positives
                        gaps=hsp.gaps
                        f.writelines(f'\t\tNumero hsp: {N+1}\n\t\t\tE-value: {E_value}\n\t\t\tScore: {Score}\n\t\t\tTamanho hsp: {Length} aa\n')
                        f.writelines(f'\t\t\tNumero de aa identicos: {identities}\n\t\t\tNumero de aa identicos e positivos: {positives}\n\t\t\tNumero de gaps: {gaps}\n')           
                    f.write('\n')

    else:
        f.write('\n---MELHORES ALINHAMENTOS---\n\n')
        for I,alinhamento in enumerate(alinhamentos):
            if I==0:
                accession_match=alinhamento.accession 
                hit_def=alinhamento.hit_def
                f.write(f'Alinhameto {I+1}: Proteina codifica pelo gene {gene}\n')
                f.writelines(f'\tAccession number (relativo a base de dados usado no blast): {accession_match}\n\tDefinicao: {hit_def}\n\n')

            else:
                accession_match=alinhamento.accession 
                hit_def=alinhamento.hit_def 
                hit_len=alinhamento.length
                f.write(f'Alinhamento {I+1}\n')
                f.writelines(f'\tAccession number hit (relativo a base de dados usado no blast): {accession_match}\n\tDefinicao: {hit_def}\n\tTamanho do hit: {hit_len} aa\n')
                hsps=alinhamento.hsps
                f.write(f'\tNumero de hsps: {len(hsps)}\n')
                for N,hsp in enumerate(hsps):
                    E_value=hsp.expect
                    Score=hsp.score
                    Length= hsp.align_length
                    identities= hsp.identities
                    positives= hsp.positives
                    gaps=hsp.gaps
                    f.writelines(f'\t\tNumero hsp: {N+1}\n\t\t\tE-value: {E_value}\n\t\t\tScore: {Score}\n\t\t\tTamanho hsp: {Length} aa\n')
                    f.writelines(f'\t\t\tNumero de aa identicos: {identities}\n\t\t\tNumero de aa identicos e positivos: {positives}\n\t\t\tNumero de gaps: {gaps}\n')           
                f.write('\n')    
    
    f.close()

# Leitura dos resultados do Blast
for gene in genes:
    
    f=open(f'{gene}_blast_results.txt','a')

    result_handle=open(f'{gene}_blast.txt','r')
    record_blast=NCBIXML.read(result_handle)
    result_handle.close()
    
    alinhamentos=record_blast.alignments
    f.write('\n---ALINHAMENTOS---\n\n')

    for I,alinhamento in enumerate(alinhamentos):
        f.write(f'Numero alinhamento: {I+1}\n')
        accession_match=alinhamento.accession 
        hit_def=alinhamento.hit_def 
        hit_len=alinhamento.length
        n_hsps=len(alinhamento.hsps)
        hsps=alinhamento.hsps
        
        if I==0: 
            f.writelines(f'\tProteina codificada pelo gene accession number (relativo a base de dados usado no blast): {accession_match}\n')
            f.writelines(f'\tDefinicao: {hit_def}\n\tTamanho da proteina: {hit_len} aa\n\tNumero de hsps: {n_hsps}\n')
            
        else:
            f.writelines(f'\tHit accession number (relativo a base de dados usada no blast): {accession_match}\n')
            f.writelines(f'\tDefinicao: {hit_def}\n\tTamanho do hit: {hit_len} aa\n\tNumero de hsps: {n_hsps}\n')

        for N,hsp in enumerate(hsps):
            E_value=hsp.expect
            Score=hsp.score
            Length= hsp.align_length
            identities= hsp.identities
            positives= hsp.positives
            gaps=hsp.gaps
            query_seq=hsp.query
            hit_seq=hsp.sbjct
            match_seq=hsp.match
            f.writelines(f'\t\tNumero do hsp:{N+1}\n\t\t\tE-value: {E_value}\n\t\t\tScore: {Score}\n\t\t\tTamanho do hsp: {Length}\n')
            f.writelines(f'\t\t\tNumero de aa identicos: {identities}\n\t\t\tNumero de aa identicos e positivos: {positives}\n\t\t\tNumero de gaps: {gaps}\n')
            if I==0:
                f.write(f'\t\t\tSequencia da proteina:\n')
                for I in range(0,len(query_seq),100):
                    seq_=query_seq[I:I+100]
                    f.writelines(f"\t\t\t\t{seq_}\n")
                f.write('\n')
            else:
                f.write(f'\t\t\tSequencia do hit:\n')
                for I in range(0,len(hit_seq),100):
                    seq_=hit_seq[I:I+100]
                    f.writelines(f"\t\t\t\t{seq_}\n")
                f.write(f'\t\t\tSequencia do alinhamento:\n')
                for I in range(0,len(match_seq),100):
                    seq_=match_seq[I:I+100]
                    f.writelines(f"\t\t\t\t{seq_}\n")
                f.write('\n')
    f.close()