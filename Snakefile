import os
from glob import glob

# Directories
MAG_DIR = "data/45_long_reads"
UNCOMPRESSED_DIR = "data/45_long_reads_uncompressed"
REFORMAT_DIR = "data/45_long_reads_reformat"
CONTIGS_DB_DIR = "results/contigs_dbs"
NCBI_DIR = "results/ncbi"
PAN_DIR = "results/pangenome"
TREE_DIR = "results/tree"

# MAG names
MAG_NAMES = [os.path.basename(f).replace(".fasta.gz","") for f in glob(f"{MAG_DIR}/*.fasta.gz")]

# NCBI genomes will be downloaded into NCBI_DIR/genomes
NCBI_GENOMES_DIR = f"{NCBI_DIR}/genomes"
NCBI_REFORMAT_DIR = f"{NCBI_DIR}/reformat"
NCBI_CONTIGS_DIR = f"{NCBI_DIR}/contigs"

rule all:
    input:
        f"{PAN_DIR}/Mediterraneibacter-PAN/Mediterraneibacter-PAN-PAN.db.done"


###########################################
# Step 1: Process MAGs
###########################################

rule uncompress_fasta:
    input:
        f"{MAG_DIR}/{{genome}}.fasta.gz"
    output:
        f"{UNCOMPRESSED_DIR}/{{genome}}.fasta"
    conda:
        "/home/ajayakod/miniconda3/envs/anvio-8"
    shell:
        "gunzip -c {input} > {output}"

rule reformat_fasta:
    input:
        f"{UNCOMPRESSED_DIR}/{{genome}}.fasta"
    output:
        f"{REFORMAT_DIR}/{{genome}}.fa"
    conda:
        "/home/ajayakod/miniconda3/envs/anvio-8"
    shell:
        "anvi-script-reformat-fasta {input} -o {output} --simplify-names --prefix {wildcards.genome}"

rule contigs_db:
    input:
        f"{REFORMAT_DIR}/{{genome}}.fa"
    output:
        f"{CONTIGS_DB_DIR}/{{genome}}-contigs.db"
    conda:
        "/home/ajayakod/miniconda3/envs/anvio-8"
    shell:
        "anvi-gen-contigs-database -f {input} -o {output} -n {{genome}}"

rule run_hmms:
    input:
        f"{CONTIGS_DB_DIR}/{{genome}}-contigs.db"
    output:
        touch(f"{CONTIGS_DB_DIR}/{{genome}}.hmms.done")
    conda:
        "/home/ajayakod/miniconda3/envs/anvio-8"
    shell:
        "anvi-run-hmms -c {input}"


###########################################
# Step 2: Download and process NCBI genomes
###########################################

rule download_ncbi:
    output:
        directory(NCBI_GENOMES_DIR)
    conda:
        "/home/ajayakod/miniconda3/envs/anvio-8"
    shell:
        """
        mkdir -p {output}
        ncbi-genome-download --assembly-level complete \
                             bacteria \
                             --genus Mediterraneibacter \
                             --metadata {output}/metadata.txt \
                             -o {output}
        """

rule reformat_ncbi:
    input:
        NCBI_GENOMES_DIR
    output:
        directory(NCBI_REFORMAT_DIR)
    conda:
        "/home/ajayakod/miniconda3/envs/anvio-8"
    shell:
        """
        mkdir -p {output}
        for f in {input}/refseq/bacteria/*/*/*.fna.gz; do
            base=$(basename $f .fna.gz)
            gunzip -c $f | anvi-script-reformat-fasta -o {output}/$base.fa --simplify-names --prefix $base
        done
        """

rule ncbi_contigs_db:
    input:
        f"{NCBI_REFORMAT_DIR}/{{genome}}.fa"
    output:
        f"{NCBI_CONTIGS_DIR}/{{genome}}-contigs.db"
    conda:
        "/home/ajayakod/miniconda3/envs/anvio-8"
    shell:
        "anvi-gen-contigs-database -f {input} -o {output} -n {{genome}}"

rule ncbi_run_hmms:
    input:
        f"{NCBI_CONTIGS_DIR}/{{genome}}-contigs.db"
    output:
        touch(f"{NCBI_CONTIGS_DIR}/{{genome}}.hmms.done")
    conda:
        "/home/ajayakod/miniconda3/envs/anvio-8"
    shell:
        "anvi-run-hmms -c {input}"


###########################################
# Step 3: Create genomes storage and build pan-genome
###########################################

rule create_genomes_storage:
    input:
        contigs_dbs=expand(f"{CONTIGS_DB_DIR}/{{genome}}-contigs.db", genome=MAG_NAMES),
        hmms_done=expand(f"{CONTIGS_DB_DIR}/{{genome}}.hmms.done", genome=MAG_NAMES)
    output:
        f"{PAN_DIR}/Mediterraneibacter-GENOMES.db"
    params:
        genomes_txt=f"{PAN_DIR}/genomes.txt"
    conda:
        "/home/ajayakod/miniconda3/envs/anvio-8"
    shell:
        """
        mkdir -p {PAN_DIR}
        echo -e "name\\tcontigs_db_path" > {params.genomes_txt}
        for db in {input.contigs_dbs}; do
            genome=$(basename $db -contigs.db)
            full_path=$(realpath $db)
            echo -e "$genome\\t$full_path" >> {params.genomes_txt}
        done
        anvi-gen-genomes-storage -e {params.genomes_txt} -o {output}
        """

rule pangenome:
    input:
        genomes_storage=f"{PAN_DIR}/Mediterraneibacter-GENOMES.db"
    output:
        f"{PAN_DIR}/Mediterraneibacter-PAN/Mediterraneibacter-PAN-PAN.db"
    params:
        outdir=f"{PAN_DIR}/Mediterraneibacter-PAN"
    conda:
        "/home/ajayakod/miniconda3/envs/anvio-8"
    shell:
        "anvi-pan-genome -g {input.genomes_storage} -n Mediterraneibacter-PAN -o {params.outdir} --num-threads 8"


###########################################
# Step 4: Generate phylogenomic tree & import layer order
###########################################

rule build_sgc_tree:
    input:
        f"{PAN_DIR}/Mediterraneibacter-PAN/Mediterraneibacter-PAN-PAN.db"
    output:
        f"{TREE_DIR}/Mediterraneibacter_sgc_tree.newick"
    conda:
        "/home/ajayakod/miniconda3/envs/anvio-8"
    shell:
        """
        mkdir -p {TREE_DIR}
        anvi-get-sequences-for-hmm-hits -c {input} --hmm-source Bacteria_71 -o {TREE_DIR}/sgc_sequences.fa
        muscle -in {TREE_DIR}/sgc_sequences.fa -out {TREE_DIR}/sgc_sequences.aln
        fasttree -nt {TREE_DIR}/sgc_sequences.aln > {output}
        """

rule import_layer_order:
    input:
        db=f"{PAN_DIR}/Mediterraneibacter-PAN/Mediterraneibacter-PAN-PAN.db",
        tree=f"{TREE_DIR}/Mediterraneibacter_sgc_tree.newick"
    output:
        touch(f"{PAN_DIR}/Mediterraneibacter-PAN/layer_order.done")
    conda:
        "/home/ajayakod/miniconda3/envs/anvio-8"
    shell:
        """
        echo -e "item_name\tdata_type\tdata_value" > {PAN_DIR}/Mediterraneibacter-PAN/layer_order.txt
        echo -e "SCGs_Bayesian_Tree\tnewick\t$(cat {input.tree})" >> {PAN_DIR}/Mediterraneibacter-PAN/layer_order.txt
        anvi-import-misc-data {PAN_DIR}/Mediterraneibacter-PAN/layer_order.txt -p {input.db}
        touch {output}
        """

###########################################
# Step 5: Annotate Pfams
###########################################

rule run_pfams:
    input:
        db=f"{PAN_DIR}/Mediterraneibacter-PAN/Mediterraneibacter-PAN-PAN.db"
    output:
        touch(f"{PAN_DIR}/Mediterraneibacter-PAN/pfams.done")
    conda:
        "/home/ajayakod/miniconda3/envs/anvio-8"
    shell:
        "anvi-run-pfams -p {input.db} --num-threads 8 && touch {output}"

###########################################
# Step 6: Functional enrichment
###########################################

rule functional_enrichment:
    input:
        db=f"{PAN_DIR}/Mediterraneibacter-PAN/Mediterraneibacter-PAN-PAN.db"
    output:
        touch(f"{PAN_DIR}/Mediterraneibacter-PAN/enrichment.done")
    conda:
        "/home/ajayakod/miniconda3/envs/anvio-8"
    shell:
        "anvi-compute-functional-enrichment -p {input.db} -g {CONTIGS_DB_DIR}/*.db --category Source && touch {output}"

###########################################
# Step 8: Functional clustering of gene clusters
###########################################

rule cluster_gene_functions:
    input:
        db=f"{PAN_DIR}/Mediterraneibacter-PAN/Mediterraneibacter-PAN-PAN.db"
    output:
        f"{PAN_DIR}/Mediterraneibacter-PAN/functional_tree.newick"
    conda:
        "/home/ajayakod/miniconda3/envs/anvio-8"
    shell:
        "anvi-matrix-to-newick {input.db} --compute-functions -o {output}"

###########################################
# Step 9: Visualize pan-genome interactively
###########################################

rule display_pan:
    input:
        pan_db=f"{PAN_DIR}/Mediterraneibacter-PAN/Mediterraneibacter-PAN-PAN.db",
        genomes_storage=f"{PAN_DIR}/Mediterraneibacter-GENOMES.db"
    output:
        touch(f"{PAN_DIR}/Mediterraneibacter-PAN/display.done")
    conda:
        "/home/ajayakod/miniconda3/envs/anvio-8"
    shell:
        "anvi-display-pan -p {input.pan_db} -g {input.genomes_storage} --server-only --port 8080 --ip 0.0.0.0 && touch {output}"

###########################################
# Step 10: Export gene cluster tables
###########################################

rule export_gene_clusters:
    input:
        db=f"{PAN_DIR}/Mediterraneibacter-PAN/Mediterraneibacter-PAN-PAN.db"
    output:
        gene_clusters=f"{PAN_DIR}/Mediterraneibacter-PAN/gene_clusters.txt",
        frequencies=f"{PAN_DIR}/Mediterraneibacter-PAN/gene_cluster_frequencies.txt",
        presence_absence=f"{PAN_DIR}/Mediterraneibacter-PAN/gene_cluster_presence_absence.txt"
    conda:
        "/home/ajayakod/miniconda3/envs/anvio-8"
    shell:
        """
        anvi-export-table {input.db} --table gene_clusters -o {output.gene_clusters}
        anvi-export-table {input.db} --table gene_cluster_frequencies -o {output.frequencies}
        anvi-export-table {input.db} --table gene_cluster_presence_absence -o {output.presence_absence}
        """

###########################################
# Final completion marker
###########################################

rule final_completion:
    input:
        f"{PAN_DIR}/Mediterraneibacter-PAN/gene_clusters.txt",
        f"{PAN_DIR}/Mediterraneibacter-PAN/gene_cluster_frequencies.txt",
        f"{PAN_DIR}/Mediterraneibacter-PAN/gene_cluster_presence_absence.txt",
        f"{PAN_DIR}/Mediterraneibacter-PAN/display.done"
    output:
        touch(f"{PAN_DIR}/Mediterraneibacter-PAN/Mediterraneibacter-PAN-PAN.db.done")
    shell:
        "echo 'Pan-genome analysis completed successfully! Interactive visualization available at http://localhost:8080'"
