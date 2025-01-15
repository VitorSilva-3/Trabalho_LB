from Bio.Blast import NCBIXML
from Bio import SeqIO, Phylo, AlignIO
import subprocess
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator

# Defina seu e-mail para NCBI
email = "PG55697@alunos.uminho.pt"

genes = ['dinB', 'dnaA', 'FE838_RS15320', 'susG']
blast_files = [f"{gene}_blast.txt" for gene in genes]

for gene, blast_file in zip(genes, blast_files):
    # Load BLAST results from XML file
    blast_results = NCBIXML.parse(open(blast_file))
    
#Alinhamento múltiplo
    sequences = []
    id_counts = {}
    for result in blast_results:
        for alignment in result.alignments:
            for hsp in alignment.hsps:
                sequence_id = alignment.hit_id
                if sequence_id in id_counts:
                    id_counts[sequence_id] += 1
                    sequence_id = f"{sequence_id}_hsps{id_counts[sequence_id]}"
                else:
                    id_counts[sequence_id] = 1

                sequence = f">{sequence_id}\n{hsp.sbjct}\n"
                sequences.append(sequence)

    fasta_file_path = f"{gene}_aligned.fasta"
    with open(fasta_file_path, "w") as fasta_file:
        fasta_file.writelines(sequences)

#Construção da árvore filogenética
def run_clustal_omega(input_fasta, output_fasta, clustal_executable):
    """
    Executa o Clustal Omega para realizar alinhamento múltiplo.
    :param input_fasta: Caminho para o arquivo FASTA de entrada
    :param output_fasta: Caminho para o arquivo FASTA de saída
    :param clustal_executable: Caminho para o executável do Clustal Omega
    """
    try:
        subprocess.run(
            [clustal_executable, "-i", input_fasta, "-o", output_fasta, "--force", "--outfmt", "clu"],
            check=True
        )
        print(f"Alinhamento completo para {input_fasta}! Resultado salvo em: {output_fasta}")
    except subprocess.CalledProcessError as e:
        print(f"Erro ao executar o Clustal Omega para {input_fasta}: {e}")


def generate_phylogenetic_tree_upgma(alignment_file, tree_output_file):
    """
    Gera uma árvore filogenética usando o método UPGMA.
    :param alignment_file: Caminho para o arquivo de alinhamento no formato CLUSTAL
    :param tree_output_file: Caminho para salvar a árvore filogenética no formato Newick
    """
    try:
        # Carregar o alinhamento
        alignment = AlignIO.read(alignment_file, "clustal")
        print(f"Alinhamento carregado para: {alignment_file}")

        # Calcular a matriz de distâncias
        calculator = DistanceCalculator('identity')
        distance_matrix = calculator.get_distance(alignment)

        # Construir a árvore usando o método UPGMA
        constructor = DistanceTreeConstructor(calculator, method="upgma")
        tree = constructor.build_tree(alignment)

        # Salvar a árvore no formato Newick
        Phylo.write(tree, tree_output_file, "newick")
        print(f"Árvore filogenética (UPGMA) salva em: {tree_output_file}")

        # Exibir a árvore no console
        print("\nÁrvore filogenética (UPGMA):")
        Phylo.draw_ascii(tree)

    except Exception as e:
        print(f"Erro ao gerar a árvore filogenética (UPGMA) para {alignment_file}: {e}")


def align_and_generate_trees(genes, clustal_executable):
    """
    Realiza alinhamento múltiplo e gera árvores filogenéticas para uma lista de genes.
    :param genes: Lista de nomes dos genes
    :param clustal_executable: Caminho para o executável do Clustal Omega
    """
    for gene in genes:
        input_fasta = f"{gene}_homology.fa"  # Arquivo de entrada para cada gene
        output_fasta = f"{gene}_co_aligned.fasta"  # Arquivo de saída correspondente
        tree_output_file = f"{gene}__tree_upgma.newick"  # Arquivo para salvar a árvore filogenética

        # Executar Clustal Omega
        print(f"\nProcessando o gene: {gene}")
        run_clustal_omega(input_fasta, output_fasta, clustal_executable)

        # Gerar a árvore filogenética com UPGMA
        generate_phylogenetic_tree_upgma(output_fasta, tree_output_file)


if __name__ == "__main__":
    # Caminho do executável do Clustal Omega
    clustal_executable = r"C:\Users\filip\Bioinformática\LabBio\Trabalho\Repositório\clustal-omega\clustalo.exe"

    # Lista de genes
    genes = ['dinB', 'dnaA', 'FE838_RS15320', 'susG']

    # Realizar alinhamento múltiplo e gerar árvores filogenéticas
    align_and_generate_trees(genes, clustal_executable)
