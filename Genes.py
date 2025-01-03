{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "ec6e000e-0547-4571-81a3-a8fe619f004f",
   "metadata": {
    "scrolled": true
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Requirement already satisfied: biopython in c:\\users\\filip\\anaconda3\\lib\\site-packages (1.84)Note: you may need to restart the kernel to use updated packages.\n",
      "\n",
      "Requirement already satisfied: numpy in c:\\users\\filip\\anaconda3\\lib\\site-packages (from biopython) (1.26.4)\n"
     ]
    }
   ],
   "source": [
    "pip install biopython"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "id": "9ec732da-a94b-40ca-b85e-fa4342120972",
   "metadata": {},
   "outputs": [],
   "source": [
    "from Bio import Entrez\n",
    "import os\n",
    "\n",
    "Entrez.email = \"pg55697@alunos.uminho.pt\"\n",
    "\n",
    "# Lista dos genes e os seus IDs correspondentes para o nucleótido e proteína\n",
    "genes = [\"dinB\", \"dnaA\", \"FE838_RS15320\", \"susG\"]\n",
    "nucleotide_id = [\"NZ_CP040530.1\", \"NZ_CP040530.1\", \"NZ_CP040530.1\", \"NZ_CP040530.1\"] \n",
    "protein_id = [\"WP_008760075.1\", \"WP_008759760.1\", \"WP_011107229.1\", \"WP_011108937.1\"]    \n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "d81eec40-e43e-41d2-9616-a9f0de101ddc",
   "metadata": {},
   "outputs": [],
   "source": [
    "for i in range(len(genes)):\n",
    "    if nucleotide_id[i] == \"NA\":\n",
    "        print(f\"No nucleotide information available for {genes[i]}.\")\n",
    "        continue \n",
    "    else:\n",
    "        handle = Entrez.efetch(db=\"nucleotide\", id=nucleotide_id[i], rettype=\"gb\", retmode=\"text\")\n",
    "        with open(f'genes_information/{genes[i]}.gb', \"w\") as doc:\n",
    "            doc.write(handle.read()) \n",
    "            handle.close()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "id": "64a730a7-c860-4da8-84e6-a79887003ec5",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Successfully fetched and saved protein information for dinB.\n",
      "Successfully fetched and saved protein information for dnaA.\n",
      "Successfully fetched and saved protein information for FE838_RS15320.\n",
      "Successfully fetched and saved protein information for susG.\n"
     ]
    }
   ],
   "source": [
    "for i in range(len(genes)):\n",
    "    try:\n",
    "        handle = Entrez.efetch(db=\"protein\", id=protein_id[i], rettype=\"gb\", retmode=\"text\")\n",
    "        \n",
    "        with open(f'genes_information/{genes[i]}_protein.gb', \"w\") as doc:\n",
    "            doc.write(handle.read())\n",
    "        print(f\"Successfully fetched and saved protein information for {genes[i]}.\")\n",
    "    \n",
    "    except Exception as e:\n",
    "        print(f\"An error occurred while fetching data for {genes[i]}: {e}\")\n",
    "    \n",
    "    finally:\n",
    "        handle.close()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 16,
   "id": "f87a2166-eff6-4b66-95a1-bcd8f918c7a9",
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "<string>:8: SyntaxWarning: invalid escape sequence '\\{'\n",
      "<>:8: SyntaxWarning: invalid escape sequence '\\{'\n",
      "<string>:8: SyntaxWarning: invalid escape sequence '\\{'\n",
      "<>:8: SyntaxWarning: invalid escape sequence '\\{'\n",
      "C:\\Users\\filip\\AppData\\Local\\Temp\\ipykernel_21588\\3508622098.py:8: SyntaxWarning: invalid escape sequence '\\{'\n",
      "  record_prot=SeqIO.read(f\"genes_information\\{gene}_protein.gb\",'genbank')\n"
     ]
    }
   ],
   "source": [
    "#Anotação da proteína de cada gene\n",
    "from Bio import SeqIO\n",
    "from Bio import SeqFeature\n",
    "\n",
    "for gene in genes:\n",
    "    f=open(f'{gene}.txt','a')\n",
    "    f.writelines('\\n***INFORMACAO RELATIVA A PROTEINA***\\n\\n')\n",
    "    record_prot=SeqIO.read(f\"genes_information\\{gene}_protein.gb\",'genbank')\n",
    "    id=record_prot.name\n",
    "    seq=record_prot.seq\n",
    "    tamanho=len(record_prot.seq)\n",
    "    des=record_prot.description\n",
    "    tipo_molecula=record_prot.annotations[\"molecule_type\"]\n",
    "    features=record_prot.features\n",
    "    f.writelines(f'Tipo de molecula: {tipo_molecula}\\nAccession number proteina: {id}\\nTamanho proteina: {len(seq)} aa\\n')\n",
    "    f.writelines(f'Descricao da proteina: {des}\\n')\n",
    "    feat_site=[]\n",
    "    for I,feature in enumerate(features):\n",
    "        if feature.type=='Protein':\n",
    "            peso_molecualr=feature.qualifiers[\"calculated_mol_wt\"][0]\n",
    "            f.writelines(f'Peso molecular: {peso_molecualr} Dalton\\nNumero de features da proteina: {len(features)}\\n')\n",
    "        if feature.type=='Site':    \n",
    "            feat_site.append(I)\n",
    "\n",
    "    if len(feat_site)!=0:\n",
    "        f.writelines('\\nLocais de interesse da proteina:\\n\\n')\n",
    "        for indice in feat_site:\n",
    "            feature=features[indice]        \n",
    "            localizacao=feature.location\n",
    "            f.writelines(f'Site {indice+1}:\\n \\t-Localizacao: {localizacao}')\n",
    "    \n",
    "    f.writelines(f'\\nSequencia proteina:\\n')\n",
    "    for I in range(0,len(seq),100):\n",
    "        seq_=seq[I:I+100]\n",
    "        f.writelines(f\"{seq_}\\n\")\n",
    "\n",
    "    f.close()\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "86bad72f-22ec-487a-8e43-5c727e383aa2",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.12.4"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
