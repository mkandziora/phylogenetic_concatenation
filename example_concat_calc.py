from Concat import concat, concat_config


aln_its = 'data/Senecio_its/updt_aln.fasta'
aln_ets = 'data/Senecio_ets/updt_aln.fasta'

workdir_comb = "output/concat"

genelist = {"ITS": aln_its,
            "ETS": aln_ets
            }

configfi = "./data/concat.config"  # configuration file


config = concat_config.ConcatConfigObj(configfi, workdir_comb)

conc = concat.Concat(config, workdir_comb)

conc = conc.run_phylup(genelistdict=genelist, self_select=True)

