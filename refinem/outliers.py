###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import os
import sys
import ast
import itertools
import logging
from collections import defaultdict

from scipy.stats import pearsonr
from numpy import (mean as np_mean)

import biolib.seq_io as seq_io
from biolib.common import find_nearest, alphanumeric_sort, remove_extension
from biolib.genomic_signature import GenomicSignature


class Outliers():
    """Identify scaffolds with divergent or compatible genomic characteristics."""

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger()

        self.min_required_coverage = 0.01

    def remove_outliers(self, genome_file, outlier_file, out_genome):
        """Remove sequences specified as outliers.

        Any scaffolds lists in the first column of
        the outlier file are removed from the specified
        genome.

        Parameters
        ----------
        genome_file : str
            Fasta file of binned scaffolds.
        outlier_file : str
            File specifying outlying scaffolds.
        out_genome : str
            Name of output genome.
        """

        genome_seqs = seq_io.read(genome_file)

        # remove scaffolds
        with open(outlier_file) as f:
            f.readline()

            for line in f:
                line_split = line.split('\t')
                scaffold_id = line_split[0]
                genome_seqs.pop(scaffold_id, None)

        # save modified bin
        seq_io.write_fasta(genome_seqs, out_genome)

    def add_compatible_unique(self, scaffold_file, genome_file, compatible_file, out_genome):
        """Add sequences specified as compatible.

        Only sequences specified exactly once in the
        compatibility file are added.

        Parameters
        ----------
        scaffold_file : str
            Fasta file containing scaffolds to add.
        genome_file : str
            Fasta file of binned scaffolds.
        compatible_file : str
            File specifying compatible scaffolds.
        out_genome : str
            Name of output genome.
        """

        cur_bin_id = remove_extension(genome_file)

        # determine scaffolds compatible with genome
        scaffold_ids = []
        bin_ids = {}
        with open(compatible_file) as f:
            f.readline()

            for line in f:
                line_split = line.split('\t')
                scaffold_id = line_split[0]
                bin_id = line_split[1].strip()

                scaffold_ids.append(scaffold_id)
                bin_ids[scaffold_id] = bin_id

        compatible_scaffolds = set()
        for scaffold_id, bin_id in bin_ids.iteritems():
            if scaffold_ids.count(scaffold_id) == 1 and bin_id == cur_bin_id:
                compatible_scaffolds.add(scaffold_id)

        # add compatible sequences to genome
        genome_seqs = seq_io.read(genome_file)
        for seq_id, seq in seq_io.read_seq(scaffold_file):
            if seq_id in compatible_scaffolds:
                genome_seqs[seq_id] = seq

        # save modified bin
        seq_io.write_fasta(genome_seqs, out_genome)

    def add_compatible_closest(self, scaffold_file, genome_file, compatible_file, out_genome):
        """Add sequences specified as compatible.

        A sequences is added to a bin if and only if it is
        closest to that bin in GC, tetranuclotide, and
        coverage space.

        Parameters
        ----------
        scaffold_file : str
            Fasta file containing scaffolds to add.
        genome_file : str
            Fasta file of binned scaffolds.
        compatible_file : str
            File specifying compatible scaffolds.
        out_genome : str
            Name of output genome.
        """

        cur_bin_id = remove_extension(genome_file)

        # determine statistics for each potentially compatible scaffold
        scaffold_ids = defaultdict(dict)
        with open(compatible_file) as f:
            headers = [x.strip() for x in f.readline().split('\t')]
            scaffold_gc_index = headers.index('Scaffold GC')
            genome_gc_index = headers.index('Mean genome GC')
            td_dist_index = headers.index('Scaffold TD')
            scaffold_cov_index = headers.index('Mean scaffold coverage')
            genome_cov_index = headers.index('Mean genome coverage')

            for line in f:
                line_split = line.split('\t')
                scaffold_id = line_split[0]
                bin_id = line_split[1].strip()

                scaffold_gc = float(line_split[scaffold_gc_index])
                genome_gc = float(line_split[genome_gc_index])
                gc_dist = abs(scaffold_gc - genome_gc)

                td_dist = float(line_split[td_dist_index])

                scaffold_cov = float(line_split[scaffold_cov_index])
                genome_cov = float(line_split[genome_cov_index])
                cov_dist = abs(scaffold_cov - genome_cov)

                scaffold_ids[scaffold_id][bin_id] = [gc_dist, td_dist, cov_dist]

        # determine scaffolds that are closest to a single bin
        # in terms of GC, tetranucleotide distance, and coverage
        compatible_scaffolds = set()
        for scaffold_id, bin_stats in scaffold_ids.iteritems():
            best_gc = [1e9, None]
            best_td = [1e9, None]
            best_cov = [1e9, None]
            for bin_id, stats in bin_stats.iteritems():
                gc, td, cov = stats
                if gc < best_gc[0]:
                    best_gc = [gc, bin_id]
                if td < best_td[0]:
                    best_td = [td, bin_id]
                if cov < best_cov[0]:
                    best_cov = [cov, bin_id]

            # check if scaffold is closest to a single bin
            if (best_gc[1] == best_td[1] == best_cov[1]) and best_gc[1] == cur_bin_id:
                compatible_scaffolds.add(scaffold_id)

        # add compatible sequences to genome
        genome_seqs = seq_io.read(genome_file)
        for seq_id, seq in seq_io.read_seq(scaffold_file):
            if seq_id in compatible_scaffolds:
                genome_seqs[seq_id] = seq

        # save modified bin
        seq_io.write_fasta(genome_seqs, out_genome)

    def identify(self, scaffold_stats, genome_stats,
                        gc_per, td_per,
                        cov_corr, cov_perc,
                        report_type, output_file):
        """Identify scaffolds with divergent genomic characteristics.

        Outliers are identified independently based on GC content,
        tetranucleotide signatures, coverage profile correlation, and
        mean absolute percent error of coverage profile. The coverage correlation
        check is ignored if the coverage profile consists of a single value.

        Parameters
        ----------
        scaffold_stats : ScaffoldStats
            Statistics for individual scaffolds.
        genome_stats : GenomeStats
            Statistics for individual genomes.
        gc_per : int.
            Percentile for identifying GC outliers
        td_per : int
            Percentile for identifying TD outliers.
        cov_corr : int
            Correlation for identifying divergent coverage profiles.
        cov_perc : int
            Mean absolute percent error for identifying divergent coverage profiles.
        report_type : str
            Report scaffolds that are outliers in 'all' or 'any' distribution.
        output_file : str
            Name of output file.
        """

        # read reference distributions from file
        self.logger.info('  Reading reference distributions.')
        self.gc_dist = self._read_distribution('gc_dist')
        self.td_dist = self._read_distribution('td_dist')

        # identify outliers in each genome
        fout = open(output_file, 'w')
        fout.write('Scaffold id\tGenome id\tScaffold length (bp)\tOutlying distributions')
        fout.write('\tScaffold GC\tMean genome GC\tLower GC bound (%s%%)\tUpper GC bound (%s%%)' % (gc_per, gc_per))
        fout.write('\tScaffold TD\tMean genome TD\tUpper TD bound (%s%%)' % td_per)
        fout.write('\tMean scaffold coverage\tMean genome coverage\tCoverage correlation\tMean coverage error\n')

        genomic_signature = GenomicSignature(0)

        processed_genomes = 0
        for genome_id, scaffold_ids in scaffold_stats.scaffolds_in_genome.iteritems():
            processed_genomes += 1

            sys.stdout.write('    Finding outliers in %d of %d (%.1f%%) genomes.\r' % (processed_genomes,
                                                                                     scaffold_stats.num_genomes(),
                                                                                     processed_genomes * 100.0 / scaffold_stats.num_genomes()))
            sys.stdout.flush()

            # find keys into GC and TD distributions
            # gc -> [mean GC][scaffold length][percentile]
            # td -> [scaffold length][percentile]
            gs = genome_stats[genome_id]
            closest_gc = find_nearest(self.gc_dist.keys(), gs.mean_gc / 100.0)
            sample_seq_len = self.gc_dist[closest_gc].keys()[0]
            d = self.gc_dist[closest_gc][sample_seq_len]
            gc_lower_bound_key = find_nearest(d.keys(), (100 - gc_per) / 2.0)
            gc_upper_bound_key = find_nearest(d.keys(), (100 + gc_per) / 2.0)

            td_bound_key = find_nearest(self.td_dist[self.td_dist.keys()[0]].keys(), td_per)

            for scaffold_id in scaffold_ids:
                stats = scaffold_stats.stats[scaffold_id]

                # find GC and TD bounds
                closest_seq_len = find_nearest(self.gc_dist[closest_gc].keys(), stats.length)
                gc_lower_bound = self.gc_dist[closest_gc][closest_seq_len][gc_lower_bound_key]
                gc_upper_bound = self.gc_dist[closest_gc][closest_seq_len][gc_upper_bound_key]

                closest_seq_len = find_nearest(self.td_dist.keys(), stats.length)
                td_bound = self.td_dist[closest_seq_len][td_bound_key]

                # find changes from mean
                delta_gc = (stats.gc - gs.mean_gc) / 100.0
                delta_td = genomic_signature.manhattan(stats.signature, gs.mean_signature)

                # determine if scaffold is an outlier
                outlying_dists = []
                if delta_gc < gc_lower_bound or delta_gc > gc_upper_bound:
                    outlying_dists.append('GC')

                if delta_td > td_bound:
                    outlying_dists.append('TD')

                corr_r = 1.0
                if len(gs.mean_coverage) > 1:
                    corr_r, _corr_p = pearsonr(gs.mean_coverage, stats.coverage)
                    if  corr_r < cov_corr:
                        outlying_dists.append('COV_CORR')

                mean_cp = []
                for cov_genome, cov_scaffold in itertools.izip(gs.mean_coverage, stats.coverage):
                    if cov_genome >= self.min_required_coverage:
                        mean_cp.append(abs(cov_scaffold - cov_genome) * 100.0 / cov_genome)

                if len(mean_cp) == 0:
                    # genome has zero coverage which is general
                    # will indicate something is wrong
                    mean_cp = -1
                    outlying_dists.append('COV_PERC')
                else:
                    mean_cp = np_mean(mean_cp)
                    if mean_cp > cov_perc:
                        outlying_dists.append('COV_PERC')

                # report outliers
                if (report_type == 'any' and len(outlying_dists) >= 1) or (report_type == 'all' and len(outlying_dists) >= 3):
                    fout.write('%s\t%s\t%s\t%s' % (scaffold_id, genome_id, stats.length, ','.join(outlying_dists)))
                    fout.write('\t%.2f\t%.2f\t%.2f\t%.2f' % (stats.gc, gs.mean_gc, gs.mean_gc + gc_lower_bound * 100, gs.mean_gc + gc_upper_bound * 100))
                    fout.write('\t%.3f\t%.3f\t%.3f' % (delta_td, gs.mean_td, td_bound))
                    fout.write('\t%.2f\t%.2f\t%.2f\t%.2f' % (np_mean(stats.coverage), np_mean(gs.mean_coverage), corr_r, mean_cp))
                    fout.write('\n')

        sys.stdout.write('\n')
        fout.close()

    def compatible(self, scaffolds_of_interest,
                        scaffold_stats,
                        genome_stats,
                        gc_per, td_per,
                        cov_corr, cov_perc,
                        report_type, output_file):
        """Identify scaffolds with compatible genomic characteristics.

        Compatible scaffolds are identified based on GC content,
        tetranucleotide signatures, coverage profile correlation, and
        mean absolute percent error of coverage profile. The coverage correlation
        check is ignored if the coverage profile consists of a single value.

        Parameters
        ----------
        scaffolds_of_interest : d[scaffold_id] -> [no. genes, perc. genes with homology]
            Scaffolds to consider for compatibility.
        scaffold_stats : ScaffoldStats
            Statistics for individual scaffolds to check.
        genome_stats : GenomeStats
            Statistics for individual genomes.
        gc_per : int
            Percentile for identifying GC outliers.
        td_per : int
            Percentile for identifying TD outliers.
        cov_corr : int
            Correlation for identifying divergent coverage profiles.
        cov_perc : int
            Mean absolute percent error for identifying divergent coverage profiles.
        report_type : str
            Report scaffolds that are outliers in 'all' or 'any' distribution.
        output_file : str
            Name of output file.
        """

        # read reference distributions from file
        self.logger.info('')
        self.logger.info('  Reading reference distributions.')
        self.gc_dist = self._read_distribution('gc_dist')
        self.td_dist = self._read_distribution('td_dist')

        # identify compatible scaffolds in each genome
        fout = open(output_file, 'w')
        fout.write('Scaffold id\tGenome id\tScaffold length (bp)\tCompatible distributions')
        fout.write('\tScaffold GC\tMean genome GC\tLower GC bound (%s%%)\tUpper GC bound (%s%%)' % (gc_per, gc_per))
        fout.write('\tScaffold TD\tMean genome TD\tUpper TD bound (%s%%)' % td_per)
        fout.write('\tMean scaffold coverage\tMean genome coverage\tCoverage correlation\tMean coverage error')
        fout.write('\t# genes\t% genes with homology\n')

        genomic_signature = GenomicSignature(0)

        self.logger.info('  Identifying scaffolds compatible with bins.')
        processed_scaffolds = 0
        for scaffold_id, ss in scaffold_stats.stats.iteritems():
            processed_scaffolds += 1
            sys.stdout.write('    Processed %d of %d (%.1f%%) scaffolds.\r' % (processed_scaffolds,
                                                                         len(scaffold_stats.stats),
                                                                         processed_scaffolds * 100.0 / len(scaffold_stats.stats)))
            sys.stdout.flush()

            if scaffold_id not in scaffolds_of_interest:
                continue

            for genome_id, gs in genome_stats.iteritems():
                # find keys into GC and TD distributions
                # gc -> [mean GC][scaffold length][percentile]
                # td -> [scaffold length][percentile]
                closest_gc = find_nearest(self.gc_dist.keys(), gs.mean_gc / 100.0)
                sample_seq_len = self.gc_dist[closest_gc].keys()[0]
                d = self.gc_dist[closest_gc][sample_seq_len]
                gc_lower_bound_key = find_nearest(d.keys(), (100 - gc_per) / 2.0)
                gc_upper_bound_key = find_nearest(d.keys(), (100 + gc_per) / 2.0)

                td_bound_key = find_nearest(self.td_dist[self.td_dist.keys()[0]].keys(), td_per)

                # find GC and TD bounds
                closest_seq_len = find_nearest(self.gc_dist[closest_gc].keys(), ss.length)
                gc_lower_bound = self.gc_dist[closest_gc][closest_seq_len][gc_lower_bound_key]
                gc_upper_bound = self.gc_dist[closest_gc][closest_seq_len][gc_upper_bound_key]

                closest_seq_len = find_nearest(self.td_dist.keys(), ss.length)
                td_bound = self.td_dist[closest_seq_len][td_bound_key]

                # find changes from mean
                delta_gc = (ss.gc - gs.mean_gc) / 100.0
                delta_td = genomic_signature.manhattan(ss.signature, gs.mean_signature)

                # determine if scaffold compatible
                compatible_dists = []
                if delta_gc >= gc_lower_bound and delta_gc <= gc_upper_bound:
                    compatible_dists.append('GC')

                if delta_td <= td_bound:
                    compatible_dists.append('TD')

                corr_r = 1.0
                if len(gs.mean_coverage) > 1:
                    corr_r, _corr_p = pearsonr(gs.mean_coverage, ss.coverage)
                    if  corr_r >= cov_corr:
                        compatible_dists.append('COV_CORR')

                mean_cp = []
                for cov_genome, cov_scaffold in itertools.izip(gs.mean_coverage, ss.coverage):
                    if cov_genome >= self.min_required_coverage:
                        mean_cp.append(abs(cov_genome - cov_scaffold) * 100.0 / cov_genome)

                mean_cp = np_mean(mean_cp)
                if mean_cp <= cov_perc:
                    compatible_dists.append('COV_PERC')

                # report compatible scaffolds
                if (report_type == 'any' and len(compatible_dists) >= 1) or (report_type == 'all' and len(compatible_dists) >= 3):
                    fout.write('%s\t%s\t%s\t%s' % (scaffold_id, genome_id, ss.length, ','.join(compatible_dists)))
                    fout.write('\t%.2f\t%.2f\t%.2f\t%.2f' % (ss.gc, gs.mean_gc, gs.mean_gc + gc_lower_bound * 100, gs.mean_gc + gc_upper_bound * 100))
                    fout.write('\t%.3f\t%.3f\t%.3f' % (delta_td, gs.mean_td, td_bound))
                    fout.write('\t%.2f\t%.2f\t%.2f\t%.2f' % (np_mean(ss.coverage), np_mean(gs.mean_coverage), corr_r, mean_cp))
                    fout.write('\t%d\t%.1f' % (scaffolds_of_interest[scaffold_id][0], scaffolds_of_interest[scaffold_id][1]))
                    fout.write('\n')

        sys.stdout.write('\n')
        fout.close()

    def create_html_index(self, plot_dir, genome_plots):
        """Create HTML index for navigating outlier plots.

        Parameters
        ----------
        plot_dir : str
          Directory containing plots.
        genome_plots : d[genome_id] -> [(plot_type, plot_filename), ...]
          Hash indicating the plot types and filenames for each genome of interest.
        """

        sorted_genome_ids = alphanumeric_sort(genome_plots.keys())

        starting_plot_filename = genome_plots[sorted_genome_ids[0]][0][1]
        starting_plot_str = sorted_genome_ids[0] + '<br>' + genome_plots[sorted_genome_ids[0]][0][0]

        fout = open(os.path.join(plot_dir, 'index.html'), 'w')
        fout.write('<html>\n')
        fout.write('<head>')
        fout.write('<title>RefineM outlier plots</title>\n')
        fout.write('</head>\n')
        fout.write('<frameset cols="15%,85%">\n')
        fout.write('<frame src="plot_menu.html" name="menu">\n')
        fout.write('<frame src="%s" name="plot">\n' % starting_plot_filename)
        fout.write('</frameset>\n')
        fout.write('</html>\n')
        fout.close()

        fout = open(os.path.join(plot_dir, 'plot_menu.html'), 'w')
        fout.write('<html>\n')
        fout.write('<script>\n')
        fout.write('    function change_title(name) {\n')
        fout.write('        document.getElementById("active_plot").innerHTML = name;\n')
        fout.write('    }\n')
        fout.write('</script>\n\n')
        fout.write('<style>\n')
        fout.write('ul {\n')
        fout.write('margin-top: 0px;\n')
        fout.write('margin-bottom: 12px;\n')
        fout.write('}\n')
        fout.write('</style>\n\n')
        fout.write('<body>\n')
        fout.write('<div><b>Active plot:</b>\n')
        fout.write('<div id="active_plot">%s</div>\n' % starting_plot_str)
        fout.write('</div>\n')
        fout.write('<br>\n')
        fout.write('<div><b>Plots:</b></div>\n')
        for genome_id in sorted_genome_ids:
            fout.write('<i>  %s:</i>\n' % genome_id)
            fout.write('    <ul>\n')
            for (plot_type, plot_filename) in genome_plots[genome_id]:
                fout.write('    <li><a href="%s" target="plot" onclick="change_title(\'%s\');">%s</a><br></li>\n' % (plot_filename,
                                                                                                                     genome_id + '<br>' + plot_type,
                                                                                                                     plot_type))
            fout.write('    </ul>\n')
        fout.write('</body>\n')
        fout.write('</html>\n')
        fout.close()

    def _read_distribution(self, prefix):
        """Read distribution file.

        Parameters
        ----------
        prefix : str
            Prefix of distibution to read (gc_dist or td_dist)

        Returns
        -------
        dict : d[percentile] -> critical value
            Critical values at integer percentiles.
        """

        dist_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'distributions', prefix + '.txt')

        with open(dist_file, 'r') as f:
            s = f.read()
            d = ast.literal_eval(s)

        return d
