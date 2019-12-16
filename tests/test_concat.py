from Concat import concat, concat_config

aln_its = 'data/Senecio_its/updt_aln.fasta'
aln_ets = 'data/Senecio_ets/updt_aln.fasta'
aln_schema = 'fasta'

tre_its = 'data/Senecio_its/updt_tre.tre'
tre_ets = 'data/Senecio_ets/updt_tre.tre'
tre_schema = 'newick'

workdir_comb = "tests/output/concat"

configfi = "data/concat.config"  # configuration file

config = concat_config.ConcatConfigObj(configfi, workdir_comb)
config.backbone = True
config.update_tree = False


spn_id_list = {"ITS": "data/Senecio_its/tipname_spn.csv",
            "ETS": "data/Senecio_ets/tipname_spn.csv"
                }


aln_list = {"ITS": aln_its,
            "ETS": aln_ets}

tree_list = {"ITS": tre_its,
            "ETS": tre_ets}

conc = concat.Concat(config, workdir_comb)
conc = conc.run_fresh(aln_dict=aln_list, spn_dict=spn_id_list, schema=aln_schema, tre_dict=tree_list, tre_schema=tre_schema)


