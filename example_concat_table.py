from Concat import concat, concat_config


aln_its = 'data/Senecio_its/updt_aln.fasta'
aln_ets = 'data/Senecio_ets/updt_aln.fasta'
aln_schema = 'fasta'

workdir_comb = "output/concat"

configfi = "data/concat.config"  # configuration file


genelist = {"ITS": aln_its,
            "ETS": aln_ets
            }



config = concat_config.ConcatConfigObj(configfi, workdir_comb)
conc = concat.Concat(config, workdir_comb)
conc = conc.write_comb_table(genelist, suggest=True)

