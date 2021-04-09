from Concat import concat, concat_config

workdir_its = 'data/Senecio_its'
workdir_ets = 'data/Senecio_ets'

workdir_comb = "tests/output/concat"

genelist = {"ITS": workdir_its,
            "ETS": workdir_ets}

configfi = "data/concat.config"  # configuration file


config = concat_config.ConcatConfigObj(configfi, workdir_comb)

conc = concat.Concat(config, workdir_comb)

conc = conc.run_phylup(genelistdict=genelist, self_select=False)

