# user settings class

import os
import sys
import configparser


def is_number(s):
    """test if string can be coerced to float"""
    try:
        float(s)
        return True
    except ValueError:
        return False


class ConcatConfigObj(object):
    def __init__(self, configfi, workdir):
        """
        Build a configuration class.

        During the initializing process the following self objects are generated:
            * **self.workdir**: working directory
            * **self.taxon_missingness** amount of missing data per otu across loci
            * **self.modeltest_criteria** BIC, AIC or AICc


        :param configfi: a configuration file in a specific format. The file needs to have a heading of the format:
                        [blast] and then somewhere below that heading as string, e.g. e_value_thresh = value
        :param workdir: the working directory
        """
        sys.stdout.write("Building concat config object\n")
        assert os.path.isfile(configfi), "file `%s` does not exists" % configfi

        self.workdir = workdir
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        assert self.workdir

        config = configparser.ConfigParser()
        config.read_file(open(configfi))

        # read in blast settings
        self.taxon_missingness = float(config["concat"]["taxon_missingness"])
        assert is_number(self.taxon_missingness), ("value `%s` does not exists" % self.taxon_missingness)

        self.subst_model_criteria = config['concat']['modeltest_criteria']
        assert self.subst_model_criteria in ['AIC', 'AICc', 'BIC'], ("self.modeltest_criteria `%s` "
                                                                   "is not AIC, AICc or BIC" % self.subst_model_criteria)

        self.backbone = config["concat"]["backbone"]
        if self.backbone == "True" or self.backbone == "true":
            self.backbone = True
        else:
            self.backbone = False

        self.ncbi_parser_nodes_fn = config["ncbi_parser"]["nodes_fn"]
        self.ncbi_parser_names_fn = config["ncbi_parser"]["names_fn"]
