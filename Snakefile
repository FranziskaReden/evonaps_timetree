rule RetrieveNames:
    input:
        "data/TimeTre_v5_Final.nwk"
    output:
        "data/TimeTree5_leaves.txt"
    conda:
        "envs/timetree.yaml"
    shell:
        """
        python scripts/parse_timetree.py --tree {input}
        """

rule GetTaxIds:
    input:
        "data/TimeTree5_leaves.txt"
    output:
        tax_ids="data/TimeTree5_tax_ids.tsv",
        failed_tax_ids = "data/TimeTree5_tax_ids_failed.txt"
        lineages="data/TimeTree5_lineage.tsv"
    conda:
        "envs/parse_taxa.yaml"
    shell:
        """
        parse_taxon_name.py -nf {input} --prefix data/TimeTree5 --mode lenient -l full -q --score 80
        """

rule ReadLineages:
    input:
        "data/TimeTree5_lineage.tsv"
    output:
        linegaes_dict = "data/TimeTree5_lineages_unresolved.json",
        reversed_dict = "data/TimeTree5_lineages_reversed.json"
    conda:
        "envs/timetree.yaml"
    shell:
        """
        python scripts/parse_lineages.py --lineages {input}
        """

rule GetCladeAges:
    input:
        linegaes_dict = "data/TimeTree5_lineages_unresolved.json",
        reversed_dict = "data/TimeTree5_lineages_reversed.json"
    output:
        lineages_resolved = "data/TimeTree5_lineages_resolved.json",
        log = "data/retrieve_age.log",
        tree_resolved = "data/TimeTree5_renamed_resolved.nwk",
        tree_renamed = "data/TimeTree5_renamed.nwk"
    conda:
        "envs/timetree.yaml"
    shell:
        """"
        python scripts/get_mrca.py \
        --tree data/TimeTre_v5_Final.nwk \
        --lineages data/TimeTree5_lineages_unresolved.json \
        --tax_ids data/TimeTree5_tax_ids.tsv --prefix data/
        """"

rule GetEvoNAPS:
    input:
        config="config/EvoNAPS_credentials.cnf"
    output:
        aa_ali = "data/aa_alignments.tsv",
        aa_tax = "data/aa_alignments_taxonomy.tsv",
        dna_ali = "data/dna_alignments.tsv",
        dna_tax = "data/dna_alignments_taxonomy.tsv"
    conda:
        "envs/timetree.yaml"
    shell:
        """"
        python scripts/get_evonaps.py --config {input.config} --prefix ./data
        """"

rule GetEvoNAPSAges:
    input:
        lineages_resolved = "data/TimeTree5_lineages_resolved.json",
        aa_ali = "data/aa_alignments.tsv",
        aa_tax = "data/aa_alignments_taxonomy.tsv",
        dna_ali = "data/dna_alignments.tsv",
        dna_tax = "data/dna_alignments_taxonomy.tsv"
    output:
        ages="data/EvoNAPS_ages.tsv"
    shell:
        """"
        python scripts/get_evonaps_ages.py --prefix ./data
        """"

