"""
Concat: automatically concatenate single gene alignments.
Copyright (C) 2019  Martha Kandziora
martha.kandziora@yahoo.com

Package to automatically concatenate single gene alignments and to calculate concatenated phylogeny.
Process can be split into two:
1. it will be checked which taxa and sequences are available and printed to file: 'concattable_selfselect.txt'.
    This can be modified by the user, see example_concat_table.py
2. the concatenated alignment file is being written: example_concat_table.py

All classes and methods are distributed under the following license.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""

import os
import subprocess
import dendropy
import random
import shutil
import math

import pandas as pd
from Bio import AlignIO, SeqIO
from dendropy import DnaCharacterMatrix, Tree

from . import cd
from PhylUp.phylogenetic_helpers import estimate_number_threads_raxml, run_modeltest, \
    build_table_from_file, replace_uid_with_name, check_align


def trim(aln, workdir, taxon_missingness=0.8):  # slight modification from aln_updater trim
    """ It removes bases at the start and end of alignments, if they are represented by less than the value
    specified in the config file. E.g. 0.8 given in config.

    This means, that 80% of the sequences need to have a base present. This ensures, that not whole chromosomes
    get dragged in by cutting the ends of long sequences or that we have lots of missing data at the beginning
    or the end.
    """
    # print('in trim')
    i = 0
    seqlen = len(aln[i])
    for tax in aln:
        assert len(aln[tax]) == seqlen, "can't trim un-aligned inputs, moving on"
    start = 0
    stop = seqlen
    cutoff = len(aln) * taxon_missingness
    range_start = range(seqlen)
    start = get_base_position(aln, cutoff, range_start, start)
    range_end = range(seqlen-1, 0, -1)
    stop = get_base_position(aln, cutoff, range_end, stop)
    # here alignment gets shortened to start:stop
    for taxon in aln:
        aln[taxon] = aln[taxon][start:stop-1]
    # make sure that tre is presented in aln
    lfd = os.path.join(workdir, "logfile")
    with open(lfd, "a") as log:
        log.write("Original aln was trimmed from 1-{} to start {}, stop {}.\n".format(seqlen, start, stop))
    return aln


def get_base_position(aln, cutoff, range_val, baseposition):
    """
    Find baseposition where sufficient data is available.

    :param aln: alignment as dendropy object
    :param cutoff: min value that needs to be present in aln
    :param range_val: range to test for - different for start and stop
    :param baseposition: position for cutting
    :return: baseposition
    """
    for i in range_val:
        counts = 0
        for tax in aln:
            base = aln[tax][i].label
            if base in ['?', '-']:
                counts += 1
        if counts <= cutoff:
            baseposition = i
            break
    return baseposition


def add_alignedseq_to_table(aln, table):
    """
    Puts input sequences from alignment into the pandas table.

    :return:
    """
    queried_taxa = []
    for taxon, seq in aln.items():
        seq = seq.symbols_as_string().replace('?', '-')
        for index in table.index:
            tip_name = table.loc[index, 'accession']
            if taxon.label in tip_name:
                table.at[index, 'sseq'] = seq
                queried_taxa.append(tip_name)
    assert tip_name in queried_taxa, (tip_name, queried_taxa)
    return table


def clean_aln(aln_fn):
    """
    Replace '?' with '-' in alignment

    :param aln_fn:
    :return:
    """
    aln = DnaCharacterMatrix.get(path=aln_fn, schema='fasta')
    aln_string = aln.as_string('fasta')
    aln_string = aln_string.replace('?', '-')
    with open(aln_fn, "w") as aln_file:
        aln_file.write(aln_string)
    aln = DnaCharacterMatrix.get(path=aln_fn, schema='fasta')
    return aln

# todo-done: make sure backbone run functions
# todo-done: implement which substitution model for which locus
# todo-done: use start tree for raxml
# todo-done: user concat table option, instead of random per taxid


class Concat(object):
    """Combines several single runs into a concatenated alignment and calculates a phylogeny.

    Currently, the concatenation is done automatically, based on taxon ids.
    It is planned to provide an option where the user can provide a concatenation table.

    User need to make sure, that there are at least some overlapping lineages between the different datasets.
    Do not concatenate data from the same loci (if you want to expand an alignment, run the program first!).

    To build the class the following is needed:
        workdir_comb: the path to your directory where the data shall be stored

    During the initializing process the following self objects are generated:
        self.workdir: the path to your directory
        self.comb_table: pandas table, similar format to the single loci table
        self.locus_len: length of sequences/aln
        self.loci data: trees per locus
        self.concatenated_aln = concatenated alignment
        self.del_columns = list of columns that were deleted because of gap only characters
    """

    def __init__(self, config, workdir_comb):
        self.workdir = workdir_comb
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.config = config
        self.comb_table = None
        self.loci_data = {}
        self.concatenated_aln = None
        self.del_columns = []  # contains info which alignment bases were deleted, through gap only
        self.counter = 0
        self.max_locus = {}  # contains name of the locus with most sequences in concatenation
        with open(os.path.join(self.workdir, "logfile"), "w+") as log:
            log.write("{Concatenation:}\n")

    def get_len_aln(self):
        """
        Get the length of alignment per locus.

        :return:
        """
        locus_len = {}
        loci = set(self.comb_table['locus'])
        for i in loci:
            # aln_locus = DnaCharacterMatrix.get(path=os.path.join(self.workdir, '{}_nogap.fas'.format(i)),
            #                                    schema='fasta')
            # len_locus = len(aln_locus[1].symbols_as_string())
            len_locus = len(self.comb_table[self.comb_table['locus'] == i]['sseq'].values[0])
            locus_len[i] = len_locus
        return locus_len

    def run_fresh(self, aln_dict, spn_dict, schema='fasta', tre_dict=None, tre_schema=None, self_select=False):
        """
        This is the main method to concatenate different runs into a single alignment and tree.

        :param aln_dict:
        :param spn_dict:  dict with gene names as key and the corresponding spn_id_file
        :param schema: format of the alignment
        :param tre_dict: dict with gene names as key and the correspinding tree file name
        :param tre_schema: format of the tree
        :return:
        """
        # self.concatfile = user_concat_fn
        if not os.path.exists(os.path.join(self.workdir, "concat_full.raxml*")):
            for locus in aln_dict.keys():
                fn = aln_dict[locus]
                fn_new = os.path.join(self.workdir, '{}.fasta'.format(locus))
                aln = DnaCharacterMatrix.get(path=fn, schema=schema)
                aln.write(path=fn_new, schema=schema)
                aln = clean_aln(fn_new)
                aln = trim(aln, self.workdir, self.config.taxon_missingness)
                self.rm_gap_only(aln, "{}.fas".format(locus))
                if tre_dict is not None:
                    assert tre_schema is not None, 'Please provide format of input tree, for concatenation.\n'
                    self.loci_data[locus] = {'tre': tre_dict[locus], 'schema': tre_schema}
                table = build_table_from_file(spn_dict[locus], self.config)
                table = add_alignedseq_to_table(aln, table)
                table['locus'] = locus
                self.comb_table = pd.concat([self.comb_table, table], ignore_index=True)
        self.calc_concat(self_select)

    def run_phylup(self, genelistdict, self_select=False):
        """
        This is the main method to concatenate different runs into a single alignment and tree.

        :param self_select: T/F: True will read  in user supplied concatenation table.
        :param genelistdict:  dictionary with gene names as key and the corresponding workdir
        :return:
        """
        # self.concatfile = user_concat_fn
        if not os.path.exists(os.path.join(self.workdir, "concat_full.raxml*")):
            for item in genelistdict.keys():
                table = self.load_single_genes(genelistdict[item], item)
                table['locus'] = item
                aln_fn = os.path.join(self.workdir, "{}_nogap.fas".format(item))
                aln = DnaCharacterMatrix.get(path=aln_fn, schema='fasta')
                table = add_alignedseq_to_table(aln, table)
                self.comb_table = pd.concat([self.comb_table, table], ignore_index=True)
                print(self.comb_table)
                assert self.comb_table.index.is_unique, ('Index is not unique')
        self.calc_concat(self_select)

    def write_comb_table(self, genelistdict, suggest):
        """
        This is the main method to concatenate different runs into a single alignment and tree.

        :param genelistdict:  dictionary with gene names as key and the corresponding workdir
        :return:
        """
        # self.concatfile = user_concat_fn
        print(self.comb_table)
        print(suggest)
        for item in genelistdict.keys():
            table = self.load_single_genes(genelistdict[item], item)
            table['locus'] = item
            aln_fn = os.path.join(self.workdir, "{}_nogap.fas".format(item))
            aln = DnaCharacterMatrix.get(path=aln_fn, schema='fasta')
            table = add_alignedseq_to_table(aln, table)
            assert table.accession.is_unique,  table[table.duplicated(['accession'])]

            self.comb_table = pd.concat([self.comb_table, table], ignore_index=True)
            print(self.comb_table)

            assert self.comb_table.index.is_unique, ('Index is not unique')
        if suggest:
            self.make_comb_table(suggest=suggest)
        else:
            self.comb_table.to_csv(os.path.join(self.workdir, 'concattable_selfselect.txt'), index=False)

    def calc_concat(self, self_select):
        """
        Main function to concatenate different alignments and to calculate a concatenated tree.

        :return:
        """
        if not os.path.exists(os.path.join(self.workdir, "concat_full.raxml*")):
            self.get_len_aln()  #
            if self_select == False:
                self.make_comb_table(suggest=False)
            else:
                self.comb_table = pd.read_csv(os.path.join(self.workdir, 'concattable_selfselect.txt'))
                assert 'concat_id' in self.comb_table.columns, self.comb_table.columns
                loci_len = len(set(self.comb_table['locus']))
                for concat_id in self.comb_table['concat_id']:
                    # print(concat_id)
                    if not math.isnan(concat_id):
                        assert (self.comb_table['concat_id'] == concat_id).sum() <= loci_len, \
                            (len(self.comb_table.loc[self.comb_table['concat_id'] == concat_id, 'accession']),
                             loci_len, concat_id)
                self.comb_table['concat_id'] = self.comb_table['concat_id'].astype('Int64')
                self.counter = self.comb_table['concat_id'].max()
                print(self.comb_table.loc[self.comb_table['concat_id'].isnull() == True, 'accession'])
                print(self.comb_table.loc[self.comb_table['concat_id'].isna() == True, 'accession'])

                # assert self.comb_table['concat_id'].isnull().values.any(),
                #        self.comb_table['concat_id'].isnull().values.tolist()
                self.comb_table.loc[(self.comb_table['concat_id'].isnull()), 'status'] = -500

            self.write_concat_aln()

            self.concatenated_aln = clean_aln(os.path.join(self.workdir, 'concat.fasta'))
            self.rm_gap_only(self.concatenated_aln, "concat.fasta")
            replace_uid_with_name(os.path.join(self.config.workdir, 'concat.fasta'), self.comb_table, 'aln')
            self.write_partition_modeltest()  # used for estimate best subst model - needs DNA instead of subst model
            partition_file = os.path.abspath(os.path.join(self.workdir, 'partition'))
            run_modeltest('concat_nogap.fas', self.workdir, self.config.subst_model_criteria, partition_file)
            file_extension = self.config.subst_model_criteria.lower()
            shutil.copy(os.path.join(self.workdir, 'concat_nogap.fas.part.{}'.format(file_extension)),
                        os.path.join(self.workdir, 'partition'))
        num_threads = estimate_number_threads_raxml(self.workdir, 'concat_nogap.fas', 'DNA')

        print(self.config.update_tree)
        if self.config.backbone is False:
            self.calculate_bootstrap_ng(num_threads)
            if self.config.update_tree:
                replace_uid_with_name(os.path.join(self.config.workdir, 'concat_full.raxml.bestTree'),
                                      self.comb_table, 'tree')
                replace_uid_with_name(os.path.join(self.config.workdir, 'concat_full.raxml.bootstraps'),
                                      self.comb_table, 'tree')
                replace_uid_with_name(os.path.join(self.config.workdir, 'concat_full.raxml.support'),
                                      self.comb_table, 'tree')
        else:
            self.est_full_tree_ng(num_threads)
            if self.config.update_tree:
                replace_uid_with_name(os.path.join(self.config.workdir, 'RAxML_bestTree.backbone_concat'),
                                      self.comb_table, 'tree')

    def make_comb_table(self, suggest):
        """
        Adds the information to the table which data shall be combined.

        """
        print(suggest)
        loci = set(self.comb_table['locus'])
        for txid in set(self.comb_table['ncbi_txid']):
            data = self.comb_table[self.comb_table['ncbi_txid'] == txid]
            # add average per loci representatives
            num_per_loci = []
            for i in loci:
                data_locus = data[data['locus'] == i]
                num_per_loci.append(len(data_locus))
            avg_num = sum(num_per_loci) / len(num_per_loci)
            max_count = 0  # counter for self.counter
            for i in loci:
                data_locus = data[data['locus'] == i]
                count = 0
                if data_locus['status_note'].str.contains('unpublished').any():
                    unpublished = data_locus[data_locus['status_note'] == 'unpublished']
                    published = data_locus[data_locus['status_note'] != 'unpublished']

                    count = self.select_for_concat(avg_num, unpublished, count)
                    count = self.select_for_concat(avg_num, published, count)
                else:
                    count = self.select_for_concat(avg_num, data_locus, count)

                if max_count < count:
                    max_count = count
            self.counter += max_count

        max_seqs = 0
        max_locus = None
        for locus in set(self.comb_table['locus']):
            per_locus = self.comb_table[self.comb_table['locus'] == locus]
            # print(per_locus)
            amnt_seq_concat_per_locus = len(per_locus[per_locus['concat_id'] != -500])
            if amnt_seq_concat_per_locus > max_seqs:
                max_seqs = amnt_seq_concat_per_locus
                max_locus = locus
        self.max_locus = max_locus
        if suggest:
            self.comb_table.to_csv(os.path.join(self.workdir, 'concattable_selfselect.txt'), index=False)
        else:
            self.comb_table.to_csv(os.path.join(self.workdir, 'concattable.txt'), index=False)

    def select_for_concat(self, avg_num, data_locus, count):
        """
        Select which sequences shall be used for concatenation.

        :param avg_num:
        :param data_locus:
        :param count:
        :return:
        """

        full_val = int(avg_num)
        if full_val == 0:
            full_val = 1
        if len(data_locus) >= full_val and len(data_locus) > 0:
            seq_to_concat = data_locus.sample(full_val)
            dropped = data_locus.drop(seq_to_concat.index.tolist())  # remove segs that are used
            for idx in dropped.index:
                self.comb_table.loc[idx, 'status'] = -500
                self.comb_table.loc[idx, 'status_note'] = 'too many seqs for locus'
        elif len(data_locus) < int(avg_num):
            seq_to_concat = data_locus
        elif int(avg_num) <= 0.5:  # if sample is very small - not many seq for sp avail
            if len(data_locus) >= 1:
                seq_to_concat = data_locus.sample(1)
            else:
                seq_to_concat = pd.DataFrame()  # is needed to overwrite old seq_to_concat
        for idx in seq_to_concat.index:
            count += 1
            self.comb_table.loc[idx, 'concat_id'] = self.counter + count
            self.comb_table.loc[idx, 'status'] = 500
            self.comb_table.loc[idx, 'status_note'] = 'used in concat'

        return count

    def write_concat_aln(self):
        """
        Write alignment of concatenated data - as file and for self.concatenated_aln.

        :return: self.concatenated_aln
        """
        #print('write concat aln')
        loci = set(self.comb_table['locus'])

        print('check len in file')
        for locus in loci:
            print(locus)
            seq_locus = self.comb_table[self.comb_table['locus'] == locus]
            first = True

            for i in range(1, self.counter + 1):
                val = seq_locus[seq_locus['concat_id'] == i]
                if i in seq_locus['concat_id'].values:
                    seq = val['sseq'].values[0]
                    seq_len = len(seq)
                    if first == True:
                        len_seqs = seq_len
                        first = False
                    else:
                        name = val['accession'].values[0]
                        assert len_seqs == seq_len, (name, len_seqs, seq_len)

        with open(os.path.join(self.workdir, 'concat.fasta'), 'w') as the_file:
            for i in range(1, self.counter+1):
                cid = '>concat_{}\n'.format(i)
                row = ''
                comb_seqs = self.comb_table[self.comb_table['concat_id'] == i]
                if not comb_seqs.empty:
                    print(cid)
                    for locus in loci:
                        seq_locus = self.comb_table[self.comb_table['locus'] == locus]
                        locus_len = self.get_len_aln()
                        # aln_locus = DnaCharacterMatrix.get(path=os.path.join(self.workdir, '{}_nogap.fas'.format(locus)),
                        #                                    schema='fasta')
                        if locus in comb_seqs['locus'].values:
                            # print(seq_locus[seq_locus['concat_id'] == i])
                            # print(seq_locus[seq_locus['concat_id'] == i]['accession'].values)
                            #if comb_seqs['locus'].str.contains(locus).any():
                            # acc = comb_seqs[comb_seqs['locus'] == locus]['accession'].values[0]
                            # found = False
                            # for taxon, seq in aln_locus.items():
                            #     if acc == taxon:
                            #         seq = "{}".format(seq.symbols_as_string())
                            #         found = True
                            #
                            # #assert found == True
                            # print(locus)
                            # print(seq)
                            seq = comb_seqs[comb_seqs['locus'] == locus]['sseq'].values[0]
                        else:
                            print('no seq')
                            seq = '-'*locus_len[locus]
                            # seq_len = len(seq)
                        row = "{}{}".format(row, seq)
                        assert len(seq) == locus_len[locus], (len(seq), locus_len[locus], locus)
                    row = cid + row + '\n'
                    the_file.write(row)
        self.concatenated_aln = DnaCharacterMatrix.get(path=os.path.join(self.workdir, 'concat.fasta'), schema='fasta')
        seq_len = check_align(self.concatenated_aln)

        assert type(seq_len) == int, (type(seq_len), seq_len)

    def rm_gap_only(self, input_aln, fn="concat.fasta", mformat="fasta"):
        """
        Remove gap only char from input aln and writes out and reloads modified aln.

        Partial copy from
        https://stackoverflow.com/questions/28220301/python-remove-special-column-from-multiple-sequence-alignment
        """
        fn_begin = fn.split(".")[0]
        self.del_columns = []
        input_aln.write(path=os.path.join(self.workdir, fn), schema="fasta")
        aln = AlignIO.read(os.path.join(self.workdir, fn), mformat)
        n_row = float(len(aln))
        n_col = float(len(aln[0]))
        i = 0
        while i < n_col:
            ct = 0
            while i + ct < n_col and aln[:, i + ct].count('-') / n_row == 1:
                ct += 1
                self.del_columns.append(i + ct)
            if ct > 0:  # delete columns [i:i+ct]
                if i == 0:
                    aln = aln[:, ct:]
                elif i + ct == n_col:
                    aln = aln[:, :i]
                else:
                    aln = aln[:, :i] + aln[:, i + ct:]
                n_col -= ct  # seq. ct positions shorter
            else:  # nothing to delete, proceed
                i += 1
        with open(os.path.join(self.workdir, "rm_gap.txt"), "a+") as del_file:
            for item in self.del_columns:
                del_file.write("{}\n".format(item))
        SeqIO.write(aln, os.path.join(self.workdir, "{}_nogap.fas".format(fn_begin)), mformat)
        input_aln = dendropy.DnaCharacterMatrix.get(
            file=open(os.path.join(self.workdir, "{}_nogap.fas".format(fn_begin))), schema=mformat)
        return input_aln

    def write_partition_modeltest(self):
        """
        Write the partitioning file for RAxML.

        Takes the info from rm_gap_only to reduce the partition by the columns that have been removed.
        """
        print("write_partition")
        count = 0
        shortend_len1 = 0
        loci = set(self.comb_table['locus'])
        for locus in loci:
            locus_len = self.get_len_aln()
            org_len_gene = locus_len[locus]
            if count == 0:  # used to delimit end of first loci - rest done under else
                # subtract removed columns (rm_gap_only) from len_gene
                # count number of cols which are smaller than len_gene
                rm_col_a = []
                for num in self.del_columns:
                    if num <= org_len_gene:
                        rm_col_a.append(num)
                len_gene0 = org_len_gene
                shortend_len0 = org_len_gene - len(rm_col_a)
                # part_len0 = shortend_len0
                with open(os.path.join(self.workdir, "partition"), "w") as partition:
                    partition.write("DNA, {} = 1-{}\n".format(locus, shortend_len0))
                count = 1
                part_len1 = len_gene0
                end = shortend_len0
            else:
                start = end + 1
                # subtract removed columns from len_gene
                # count number of cols which are smaller than len_gene, must be done with original col length (rm_col_a)
                rm_col = []
                for num in self.del_columns:
                    # same as if num > part_len1 and num <= (part_len1 + org_len_gene):
                    if part_len1 < num <= (part_len1 + org_len_gene):
                        rm_col.append(num)
                shortend_len1 = shortend_len1 + org_len_gene - len(rm_col)
                end = shortend_len0 + shortend_len1
                part_len1 = part_len1 + org_len_gene
                with open(os.path.join(self.workdir, "partition"), "a") as partition:
                    partition.write("DNA, {} = {}-{}\n".format(locus, start, end))

    def load_single_genes(self, workdir_singlerun, genename):
        """
        Load PhylUp class objects and make a single dict per run.

        Removes abandoned nodes and gap only columns first.

        :param workdir_singlerun: directory of single gene run
        :param genename: string, name for locus provided by user
        :return: table
        """
        print("load_single_genes: {}".format(genename))
        fn = os.path.join(workdir_singlerun, 'updt_aln.fasta')
        aln = DnaCharacterMatrix.get(path=fn, schema="fasta")
        seq_len = check_align(aln)
        assert type(seq_len) == int, (type(seq_len), seq_len)
        aln.write(path=fn, schema='fasta')
        aln = clean_aln(fn)
        table = pd.read_csv(os.path.join(workdir_singlerun, 'table.updated'))
        table = table[table['status'] >= 0]
        # remove gap only char from input aln and return modified aln
        aln = trim(aln, self.workdir, self.config.taxon_missingness)
        self.rm_gap_only(aln, "{}.fas".format(genename))
        return table

    # def write_labelled(self, tree_file):
    #     """
    #     write labelled  tree
    #
    #     :param tree_file: file name of tree file
    #     :return:
    #     """
    #     tr_fn = os.path.join(self.workdir, tree_file)
    #     with open(tr_fn, "r") as label_new:
    #         new_tree = label_new.read()
    #         concat_ids = self.comb_table['concat_id'].dropna().unique()
    #         for item in concat_ids:
    #             with open("{}_relabel".format(tr_fn), "wt") as fout:
    #                 current_id = "concat_{}:".format(str(item).split('.')[0])
    #                 assert len(self.comb_table.loc[self.comb_table['concat_id'] == item, 'ncbi_txn'].unique()) == 1, \
    #                     self.comb_table.loc[self.comb_table['concat_id'] == item, 'ncbi_txn']
    #                 spn = self.comb_table.loc[self.comb_table['concat_id'] == item,
    #                                           'ncbi_txn'].to_list()[0].replace(" ", "_")
    #                 new_id = "{}_{}".format(spn, current_id)
    #                 new_tree = new_tree.replace(current_id, new_id)
    #                 fout.write(new_tree)

    def get_starting_tree(self):
        """
        Get a starting tree from the single runs - its the one with most sequences represented in the concatenated
        dataset.  It will prune all tips not present in concatenation.
        
        :return:
        """
        print('get starting tree')
        starting_tree_fn = self.loci_data[self.max_locus]['tre']
        tre_schema = self.loci_data[self.max_locus]['schema']
        if starting_tree_fn is not "None" or starting_tree_fn is not None:
            tree = Tree.get(path=starting_tree_fn, schema=tre_schema, preserve_underscores=True)
            seq_present = self.comb_table[self.comb_table['locus'] == self.max_locus]
            seq_present = seq_present[seq_present['status'] != -500]
            tree_label = list(taxon.label for taxon in tree.taxon_namespace)
            for tipname in tree_label:
                if tipname not in seq_present['accession'].to_list():
                    tree.prune_taxa([tipname])
                    tree.prune_taxa_with_labels([tipname])
                    tree.prune_taxa_with_labels([tipname])
            tree.write(path=os.path.join(self.workdir, 'tmp.tre'),
                       schema='newick', unquoted_underscores=True, suppress_rooting=True)
            with open(os.path.join(self.workdir, 'tmp.tre'), "r") as label_new:
                new_tree = label_new.read()
                with open(os.path.join(self.workdir, "starting.tre"), "wt") as fout:
                    for tipname in tree.taxon_namespace:
                        if tipname.label in seq_present['accession'].to_list():
                            row = seq_present[seq_present['accession'] == tipname.label]
                            conc_id = str(seq_present.loc[row.index, 'concat_id'].values[0]).split('.')[0]
                            concat_id = "concat_{}".format(conc_id)
                            if tipname.label in new_tree:
                                new_tree = new_tree.replace(tipname.label, concat_id)
                    fout.write(new_tree)

    def est_full_tree_ng(self, num_threads=2):
        """Full raxml run from the placement tree as starting tree.
        """
        print("run full tree")
        seed = random.randint(1, 21)
        num_threads = str(num_threads)
        if self.config.backbone is True:
            self.get_starting_tree()
            namefix = 'backbone_concat'
        else:
            self.get_starting_tree()
            namefix = 'concat'
        if os.path.exists('starting.tre'):
            starting_tree_fn = 'starting.tre'
        else:
            starting_tree_fn = None

        with cd(self.workdir):
            if os.path.exists("concat_nogap.fas.reduced") and os.path.exists("partition.reduced"):
                aln = "concat_nogap.fas.reduced"
                partition = "partition.reduced"
            else:
                aln = "concat_nogap.fas"
                partition = "partition"
            if self.config.update_tree is True:
                if starting_tree_fn is None:
                    print("no starting tree")
                    subprocess.call(["raxml-ng-mpi",
                                     "--threads", "{}".format(num_threads),
                                     "--msa", aln,
                                     '--model', partition,
                                     "--seed", "{}".format(seed),
                                     "--prefix", namefix])
                elif self.config.backbone is not True:
                    print("no backbone")
                    subprocess.call(["raxml-ng-mpi",
                                     "--threads", "{}".format(num_threads),
                                     "--msa", aln,
                                     '--model', partition,
                                     "--tree", starting_tree_fn,
                                     "--seed", "{}".format(seed),
                                     "--prefix", namefix])
                else:
                    print("backbone")
                    subprocess.call(["raxml-ng-mpi",
                                     "--threads", "{}".format(num_threads),
                                     "--msa", aln,
                                     '--model', partition,
                                     "--tree-constraint", starting_tree_fn,
                                     "--seed", "{}".format(seed),
                                     "--prefix", namefix])
            else:
                todo = 'To update the data run the following command in your working directory.'
                if starting_tree_fn is None:
                    cmd1 = "raxml-ng-mpi --msa {} --model {} --seed {} --threads {} " \
                           "--prefix {}".format(aln, partition, seed, num_threads, namefix)
                elif self.config.backbone is not True:
                    cmd1 = "raxml-ng-mpi --msa {} --model {} --tree {} --seed {} --threads {} " \
                           "--prefix {}".format(aln, partition, starting_tree_fn, seed, num_threads, namefix)
                else:
                    cmd1 = "raxml-ng-mpi --msa {} --model {} --tree-constraint {} --seed {} --threads {} " \
                           "--prefix {}".format(aln, partition, starting_tree_fn, seed, num_threads, namefix)

                print(todo)
                print(cmd1)

                lfd = os.path.join(self.workdir, "logfile")
                with open(lfd, "a") as log:
                    log.write("{}\n".format(todo))
                    log.write("{}\n".format(cmd1))

    def calculate_bootstrap_ng(self, num_threads=2):
        """Calculate bootstrap and consensus trees using raxml-ng.
        """
        print("calc bootstrap")
        with cd(self.workdir):
            if os.path.exists("concat_nogap.fas.reduced"):
                aln_fn = "concat_nogap.fas.reduced"
            else:
                aln_fn = "concat_nogap.fas"
            if os.path.exists("starting.tre"):
                starting_fn = "starting.tre"
            else:
                starting_fn = None
            partition = "partition"
            ntasks = os.environ.get('SLURM_NTASKS_PER_NODE')
            nnodes = os.environ.get("SLURM_JOB_NUM_NODES")
            # env_var = int(nnodes) * int(ntasks)
            seed = str(random.randint(1, 21))
            num_threads = str(num_threads)
            mpi = False

            if nnodes is not None and ntasks is not None:
                env_var = int(nnodes) * int(ntasks)
                mpi = True
            if self.config.update_tree is True:
                if starting_fn is None:
                    if mpi:
                        print("run with mpi")
                        subprocess.call(["mpiexec", "-n", "{}".format(env_var), "raxml-ng-mpi", '--parse',
                                         '--all', "--msa", "{}".format(aln_fn), '--model', partition, '--bs-trees',
                                         'autoMRE', '--seed', seed, "--threads", "{}".format(num_threads),
                                         "--prefix", "concat_full"])
                    else:
                        print('run without mpi')
                        subprocess.call(["raxml-ng-mpi", '--all', "--msa", aln_fn, '--model', partition, '--bs-trees',
                                         'autoMRE', '--seed', seed, "--threads", num_threads, "--prefix",
                                         'concat_full'])
                else:
                    if mpi:
                        print("run with mpi")
                        subprocess.call(["mpiexec", "-n", "{}".format(env_var), "raxml-ng-mpi", '--parse',
                                         '--all', "--msa", "{}".format(aln_fn), '--model', partition, '--bs-trees',
                                         'autoMRE', '--seed', seed, "--threads", "{}".format(num_threads),
                                         '--tree', starting_fn, "--prefix", "concat_full"])
                    else:
                        print('run without mpi')
                        subprocess.call(["raxml-ng-mpi", '--all', "--msa", aln_fn, '--tree', starting_fn, '--model',
                                         partition,  '--bs-trees', 'autoMRE', '--seed', seed, "--threads", num_threads,
                                         "--prefix", 'concat_full'])
                #subprocess.call(["raxml-ng-mpi", '--support', '--tree', 'concat_full.raxml.bestTree', '--bs-trees',
                #                 'concat_full.raxml.bootstraps', "--prefix", 'support'])
                subprocess.call(["raxml-ng-mpi", '--consense', 'MRE', '--tree', 'concat_full.raxml.bootstraps',
                                 "--prefix", 'consMRE'])
                subprocess.call(["raxml-ng-mpi", '--consense', 'STRICT', '--tree', 'concat_full.raxml.bootstraps',
                                 "--prefix", 'consSTRICT'])
                subprocess.call(["raxml-ng-mpi", '--consense', 'MR', '--tree', 'concat_full.raxml.bootstraps',
                                 "--prefix", 'consMR'])
            else:
                todo = 'To update the data run the following command in your working directory.'
                if starting_fn is None:
                    cmd1 = "raxml-ng-mpi --all --msa {} --model {} --bs-trees autoMRE --seed {} --threads {} " \
                           "--prefix concat_full".format(aln_fn, partition, seed, num_threads)
                else:
                    cmd1 = "raxml-ng-mpi --all --msa {} --tree {} --model {} --bs-trees autoMRE --seed {} " \
                           "--threads {} --prefix concat_full".format(aln_fn, starting_fn, partition, seed, num_threads)
                cmd2 = "raxml-ng-mpi --consense MRE --tree concat_full.raxml.bootstraps --prefix consMRE"
                cmd3 = "raxml-ng-mpi --consense STRICT --tree concat_full.raxml.bootstraps --prefix consSTRICT"
                cmd4 = "raxml-ng-mpi --consense MR --tree concat_full.raxml.bootstraps --prefix consMR"

                print(todo)
                print(cmd1)
                print(cmd2)
                print(cmd3)
                print(cmd4)

                lfd = os.path.join(self.workdir, "logfile")
                with open(lfd, "a") as log:
                    log.write("{}\n".format(todo))
                    log.write("{}\n".format(cmd1))
                    log.write("{}\n".format(cmd2))
                    log.write("{}\n".format(cmd3))
                    log.write("{}\n".format(cmd4))

