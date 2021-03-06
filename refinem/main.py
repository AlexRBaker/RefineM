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
import logging
import ast #for safe eval function
from collections import defaultdict

from refinem.scaffold_stats import ScaffoldStats
from refinem.genome_stats import GenomeStats
from refinem.gene_profile import GeneProfile
from refinem.bin_comparer import BinComparer
from refinem.reference import Reference
from refinem.unbinned import Unbinned
from refinem.coverage import Coverage
from refinem.tetranucleotide import Tetranucleotide
from refinem.outliers import Outliers
from refinem.cluster import Cluster
from refinem.plots.gc_plots import GcPlots
from refinem.plots.td_plots import TdPlots
from refinem.plots.cov_perc_plots import CovPercPlots
from refinem.plots.cov_corr_plots import CovCorrPlots
from refinem.plots.distribution_plots import DistributionPlots
from refinem.plots.gc_cov_plot import GcCovPlot
from refinem.plots.tetra_pca_plot import TetraPcaPlot
from refinem.plots.combined_plots import CombinedPlots
from refinem.singlegenome import WindowGen

import biolib.seq_io as seq_io
import biolib.genome_tk as genome_tk
from biolib.common import (make_sure_path_exists,
                           check_dir_exists,
                           check_file_exists,
                           query_yes_no)
from biolib.external.prodigal import Prodigal
from biolib.misc.time_keeper import TimeKeeper
from biolib.external.execute import check_dependencies


class OptionsParser():
    def __init__(self):
        """Initialization"""
        self.logger = logging.getLogger()
        self.time_keeper = TimeKeeper()
        
    #~ def item_eval(item):
        #~ try:
            #~ return ast.literal_eval(item)
        #~ except ValueError:
            #~ return item

    def _genome_files(self, genome_dir, genome_ext):
        """Identify genomes files.

        Parameters
        ----------
        genome_dir : str
            Directory containing genomes of interest.
        genome_ext : str
            Extension of genome files.

        Returns
        -------
        list
            Path to genome files.
        """

        check_dir_exists(genome_dir)

        genome_files = []
        for f in os.listdir(genome_dir):
            if f.endswith(genome_ext):
                genome_files.append(os.path.join(genome_dir, f))

        if not genome_files:
            self.logger.warning('  [Warning] No genomes found. Check the --genome_ext or --protein_ext flag used to identify genomes.')
            sys.exit()

        return genome_files

    def _check_nuclotide_seqs(self, seq_files):
        """Check if files contain sequences in nucleotide space.

        Parameters
        ----------
        seq_files : iterable
            Sequence files to check.

        Returns
        -------
        boolean
            True if files can be treated as containing nucleotide sequences.
        """

        for seq_file in seq_files:
            if not seq_io.is_nucleotide(seq_file):
                print('Expected all files to contain sequences in nucleotide space.')
                print('File %s appears like it may contain amino acids sequences.' % seq_file)

                yes_response = query_yes_no('Do all files contain only nucleotide sequences?', default='no')
                if not yes_response:
                    return False

        return True

    def _check_protein_seqs(self, seq_files):
        """Check if files contain sequences in amino acid space.

        Parameters
        ----------
        seq_files : iterable
            Sequence files to check.

        Returns
        -------
        boolean
            True if files can be treated as containing amino acid sequences.
        """

        for seq_file in seq_files:
            if not seq_io.is_protein(seq_file):
                print('Expected all files to contain sequences in amino acid space.')
                print('File %s appears like it may contain nucleotide sequences.' % seq_file)

                yes_response = query_yes_no('Do all files contain only amino acid sequences?', default='no')
                if not yes_response:
                    return False

        return True

    def scaffold_stats(self, options):
        """Scaffold statistics command"""
        print options
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [RefineM - scaffold_stats] Calculating statistics for scaffolds.')
        self.logger.info('*******************************************************************************')

        check_file_exists(options.scaffold_file)

        if not self._check_nuclotide_seqs([options.scaffold_file]):
            self.logger.warning('[Warning] Scaffold file must contain nucleotide sequences.')
            sys.exit()

        genome_files = self._genome_files(options.genome_nt_dir, options.genome_ext)
        if not self._check_nuclotide_seqs(genome_files):
            self.logger.warning('[Warning] All files must contain nucleotide sequences.')
            sys.exit()

        make_sure_path_exists(options.output_dir)

        # get coverage information
        if not options.coverage_file:
            if not options.bam_files:
                self.logger.warning('\n  [Warning] One or more BAM files must be specified in order to calculate coverage profiles.')
                coverage_file = None
            else:
                coverage = Coverage(options.cpus)
                coverage_file = os.path.join(options.output_dir, 'coverage.tsv')
                coverage.run(options.bam_files, coverage_file, options.cov_all_reads, options.cov_min_align, options.cov_max_edit_dist)
                self.logger.info('')
                self.logger.info('  Coverage profiles written to: %s' % coverage_file)
        else:
            coverage_file = options.coverage_file

        # get tetranucleotide signatures - ALEX - IMPORTANT FOR MY STUFF 
        if not options.tetra_file:
            self.logger.info('')
            tetra = Tetranucleotide(options.cpus)
            tetra_file = os.path.join(options.output_dir, 'tetra.tsv')
            signatures = tetra.run(options.scaffold_file)
            tetra.write(signatures, tetra_file)
            self.logger.info('  Tetranucleotide signatures written to: %s' % tetra_file)
        else:
            tetra_file = options.tetra_file

        # write out scaffold statistics
        stats_output = os.path.join(options.output_dir, 'scaffold_stats.tsv')
        stats = ScaffoldStats(options.cpus)
        stats.run(options.scaffold_file, genome_files, tetra_file, coverage_file, stats_output)

        self.logger.info('  Scaffold statistic written to: %s' % stats_output)

        self.time_keeper.print_time_stamp()

    def genome_stats(self, options):
        """Genomes statistics command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [RefineM - genome_stats] Calculating statistics for genomes.')
        self.logger.info('*******************************************************************************')

        check_file_exists(options.scaffold_stats_file)

        self.logger.info('')
        self.logger.info('  Reading scaffold statistics.')
        scaffold_stats = ScaffoldStats(options.cpus)
        scaffold_stats.read(options.scaffold_stats_file)

        genome_stats = GenomeStats()
        genome_stats.run(scaffold_stats)
        genome_stats.write(options.output_file)

        self.logger.info('  Genome statistic written to: %s' % options.output_file)

        self.time_keeper.print_time_stamp()

    def gene_profile(self, options):
        """Call genes command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [RefineM - gene_profile] Generating taxonomic profiles from genes.')
        self.logger.info('*******************************************************************************')

        make_sure_path_exists(options.output_dir)
        check_file_exists(options.scaffold_stats_file)
        check_file_exists(options.taxonomy_file)
        check_file_exists(options.db_file)

        gene_files = self._genome_files(options.genome_prot_dir, options.protein_ext)
        if not self._check_protein_seqs(gene_files):
            self.logger.warning('[Warning] All files must contain amino acid sequences.')
            sys.exit()

        # build gene profile
        gene_profile = GeneProfile(options.cpus, options.output_dir)
        gene_profile.run(gene_files,
                            options.scaffold_stats_file,
                            options.db_file,
                            options.taxonomy_file,
                            options.per_to_classify,
                            options.evalue,
                            options.per_identity)

        self.logger.info('')
        self.logger.info('  Results written to: %s' % options.output_dir)

        self.time_keeper.print_time_stamp()

    def outliers(self, options):
        """Outlier command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [RefineM - outliers] Identifying scaffolds with divergent characteristics.')
        self.logger.info('*******************************************************************************')

        check_file_exists(options.scaffold_stats_file)
        make_sure_path_exists(options.output_dir)

        self.logger.info('')
        self.logger.info('  Reading scaffold statistics.')
        scaffold_stats = ScaffoldStats()
        scaffold_stats.read(options.scaffold_stats_file)

        genome_stats = GenomeStats()
        genome_stats = genome_stats.run(scaffold_stats)

        # identify outliers
        outliers = Outliers()
        outlier_file = os.path.join(options.output_dir, 'outliers.tsv')
        outliers.identify(scaffold_stats, genome_stats,
                                      options.gc_perc, options.td_perc,
                                      options.cov_corr, options.cov_perc,
                                      options.report_type, outlier_file)
        self.logger.info('  Outlier information written to: ' + outlier_file)

        # create outlier plots
        self.logger.info('')

        highlight_scaffolds_ids = {}
        if options.highlight_file:
            for line in open(options.highlight_file):
                line_split = line.strip().split('\t')
                if len(line_split) > 1:
                    highlight_scaffolds_ids[line_split[0]] = [float(x.strip()) / 255.0 for x in line_split[1].split(',')]
                else:
                    highlight_scaffolds_ids[line_split[0]] = [1.0, 0, 0]

        link_scaffold_ids = []
        if options.links_file:
            with open(options.links_file) as links_file:
                for line in links_file:
                    #print line.strip().split('\t')
                    link_scaffold_ids.append([ast.literal_eval(item) if i not in (0,2) else item for i,item in enumerate((line.strip().split('\t')))])
            #link_scaffold_ids.append(line.strip().split('\t') for line in open(options.links_file))
            
        #print list(link_scaffold_ids[0])
        
        # create plots
        genomes_processed = 0
        plot_dir = os.path.join(options.output_dir, 'plots')
        make_sure_path_exists(plot_dir)
        genome_plots = defaultdict(list)
        for genome_id, gs in genome_stats.iteritems():
            genomes_processed += 1

            sys.stdout.write('  Plotting scaffold distribution for %d of %d (%.1f%%) genomes.\r' %
                                                                                            (genomes_processed,
                                                                                             len(genome_stats),
                                                                                             genomes_processed * 100.0 / len(genome_stats)))
            sys.stdout.flush()

            genome_scaffold_stats = {}
            for scaffold_id in scaffold_stats.scaffolds_in_genome[genome_id]:
                genome_scaffold_stats[scaffold_id] = scaffold_stats.stats[scaffold_id]

            if options.individual_plots:
                #~ # GC plot
                #~ gc_plots = GcPlots(options)
                #~ gc_plots.plot(genome_scaffold_stats, highlight_scaffolds_ids, link_scaffold_ids, gs.mean_gc, outliers.gc_dist, [options.gc_perc])
#~ 
                #~ output_plot = os.path.join(plot_dir, genome_id + '.gc_plots.' + options.image_type)
                #~ gc_plots.save_plot(output_plot, dpi=options.dpi)
                #~ gc_plots.save_html(os.path.join(plot_dir, genome_id + '.gc_plots.html'))

                # TD plot
                td_plots = TdPlots(options)
                td_plots.plot(genome_scaffold_stats, highlight_scaffolds_ids, link_scaffold_ids, gs.mean_signature, outliers.td_dist, [options.td_perc])

                output_plot = os.path.join(plot_dir, genome_id + '.td_plots.' + options.image_type)
                td_plots.save_plot(output_plot, dpi=options.dpi)
                td_plots.save_html(os.path.join(plot_dir, genome_id + '.td_plots.html'))

                #~ # mean absolute deviation of coverage profiles
                #~ cov_perc_plots = CovPercPlots(options)
                #~ cov_perc_plots.plot(genome_scaffold_stats, highlight_scaffolds_ids, link_scaffold_ids, gs.mean_coverage, [options.cov_perc])
#~ 
                #~ output_plot = os.path.join(plot_dir, genome_id + '.cov_perc.' + options.image_type)
                #~ cov_perc_plots.save_plot(output_plot, dpi=options.dpi)
                #~ cov_perc_plots.save_html(os.path.join(plot_dir, genome_id + '.cov_perc.html'))
#~ 
                #~ # coverage correlation plots
                #~ if len(gs.mean_coverage) > 1:
                    #~ cov_corr_plots = CovCorrPlots(options)
                    #~ cov_corr_plots.plot(genome_scaffold_stats, highlight_scaffolds_ids, gs.mean_coverage, [options.cov_corr])
#~ 
                    #~ output_plot = os.path.join(plot_dir, genome_id + '.cov_corr.' + options.image_type)
                    #~ cov_corr_plots.save_plot(output_plot, dpi=options.dpi)
                    #~ cov_corr_plots.save_html(os.path.join(plot_dir, genome_id + '.cov_corr.html'))

            #~ # combined distribution, GC vs. coverage, and tetranucleotide signature plots
            #~ combined_plots = CombinedPlots(options)
            #~ combined_plots.plot(genome_scaffold_stats,
                            #~ highlight_scaffolds_ids, link_scaffold_ids, gs,
                            #~ outliers.gc_dist, outliers.td_dist,
                            #~ options.gc_perc, options.td_perc, options.cov_perc)
#~ 
            #~ output_plot = os.path.join(plot_dir, genome_id + '.combined.' + options.image_type)
            #~ combined_plots.save_plot(output_plot, dpi=options.dpi)
            #~ combined_plots.save_html(os.path.join(plot_dir, genome_id + '.combined.html'))
#~ 
            #~ genome_plots[genome_id].append(('Combined', genome_id + '.combined.html'))
#~ 
            #~ # combined plot of distributions
            #~ dist_plots = DistributionPlots(options)
            #~ dist_plots.plot(genome_scaffold_stats,
                            #~ highlight_scaffolds_ids,
                            #~ link_scaffold_ids,
                            #~ gs,
                            #~ outliers.gc_dist, outliers.td_dist,
                            #~ options.gc_perc, options.td_perc, options.cov_perc)
#~ 
            #~ output_plot = os.path.join(plot_dir, genome_id + '.dist_plot.' + options.image_type)
            #~ dist_plots.save_plot(output_plot, dpi=options.dpi)
            #~ dist_plots.save_html(os.path.join(plot_dir, genome_id + '.dist_plot.html'))
#~ 
            #~ genome_plots[genome_id].append(('Distributions', genome_id + '.dist_plot.html'))
#~ 
            #~ # GC vs. coverage plot
            #~ gc_cov_plot = GcCovPlot(options)
            #~ gc_cov_plot.plot(genome_scaffold_stats,
                             #~ highlight_scaffolds_ids, link_scaffold_ids,
                             #~ gs.mean_gc, gs.mean_coverage)
#~ 
            #~ output_plot = os.path.join(plot_dir, genome_id + '.gc_coverge.' + options.image_type)
            #~ gc_cov_plot.save_plot(output_plot, dpi=options.dpi)
            #~ gc_cov_plot.save_html(os.path.join(plot_dir, genome_id + '.gc_coverge.html'))
#~ 
            #~ genome_plots[genome_id].append(('GC vs. coverage', genome_id + '.gc_coverge.html'))

            # tetranucleotide signature PCA plot
            tetra = TetraPcaPlot(options)
            tetra.plot(genome_scaffold_stats, highlight_scaffolds_ids, link_scaffold_ids)

            output_plot = os.path.join(plot_dir, genome_id + '.tetra_pca.' + options.image_type)
            tetra.save_plot(output_plot, dpi=options.dpi)
            tetra.save_html(os.path.join(plot_dir, genome_id + '.tetra_pca.html'))

            genome_plots[genome_id].append(('Tetra PCA', genome_id + '.tetra_pca.html'))

        sys.stdout.write('\n')

        outliers.create_html_index(plot_dir, genome_plots)

        self.logger.info('  Outlier plots written to: ' + plot_dir)

        self.time_keeper.print_time_stamp()

    def cluster(self, options):
        """Cluster command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [RefineM - cluster] Partitioning bin into clusters.')
        self.logger.info('*******************************************************************************')

        check_file_exists(options.scaffold_stats_file)
        check_file_exists(options.genome_file)
        make_sure_path_exists(options.output_dir)

        self.logger.info('')
        self.logger.info('  Reading scaffold statistics.')
        scaffold_stats = ScaffoldStats()
        scaffold_stats.read(options.scaffold_stats_file)

        cluster = Cluster(options.cpus)
        cluster.run(scaffold_stats,
                    options.num_clusters,
                    options.num_components,
                    options.K,
                    options.no_coverage,
                    options.no_pca,
                    options.iterations,
                    options.genome_file,
                    options.output_dir)

        self.logger.info('')
        self.logger.info('  Partitioned sequences written to: ' + options.output_dir)

        self.time_keeper.print_time_stamp()

    def reference(self, options):
        """Reference command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info('[RefineM - reference] Identifying scaffolds similar to specific genome(s).')
        self.logger.info('*******************************************************************************')

        check_file_exists(options.scaffold_prot_file)
        check_file_exists(options.scaffold_stats_file)
        make_sure_path_exists(options.output_dir)

        ref_gene_files = self._genome_files(options.ref_genome_prot_dir, options.protein_ext)
        if not self._check_protein_seqs(ref_gene_files):
            self.logger.warning('[Warning] All files must contain amino acid sequences.')
            sys.exit()

        reference = Reference(options.cpus, options.output_dir)
        reference_out = reference.run(options.scaffold_prot_file,
                                        options.scaffold_stats_file,
                                        ref_gene_files,
                                        options.db_file,
                                        options.evalue,
                                        options.per_identity)

        self.logger.info('')
        self.logger.info('  Results written to: ' + reference_out)

        self.time_keeper.print_time_stamp()

    def compatible(self, options):
        """Compatible command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info('[RefineM - compatible] Identify scaffolds with compatible genomic statistics.')
        self.logger.info('*******************************************************************************')

        check_file_exists(options.reference_file)
        check_file_exists(options.scaffold_stats_file)
        make_sure_path_exists(options.output_dir)

        # read scaffold statistics and calculate genome stats
        self.logger.info('')
        self.logger.info('  Reading scaffold statistics.')
        scaffold_stats = ScaffoldStats()
        scaffold_stats.read(options.scaffold_stats_file)

        genome_stats = GenomeStats()
        genome_stats = genome_stats.run(scaffold_stats)

        # identify putative homologs to reference genomes
        reference = Reference(1, None)
        putative_homologs = reference.homology_check(options.reference_file,
                                                         options.min_genes,
                                                         float(options.perc_genes))

        # identify scaffolds compatible with bins
        outliers = Outliers()
        output_file = os.path.join(options.output_dir, 'compatible.tsv')
        outliers.compatible(putative_homologs, scaffold_stats, genome_stats,
                                      options.gc_perc, options.td_perc,
                                      options.cov_corr, options.cov_perc,
                                      options.report_type, output_file)

        self.logger.info('')
        self.logger.info('  Results written to: ' + output_file)

        self.time_keeper.print_time_stamp()

    def modify(self, options):
        """Modify command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [RefineM - modify] Modifying scaffolds in genome.')
        self.logger.info('*******************************************************************************')

        make_sure_path_exists(os.path.dirname(options.output_genome))

        if not (options.add or options.remove or options.outlier_file or options.compatible_file):
            self.logger.warning('  [Warning] No modification to bin requested.\n')
            sys.exit()

        if (options.add or options.remove) and (options.outlier_file or options.compatible_file):
            self.logger.warning("  [Warning] The 'outlier_file' and 'compatible_file' options cannot be specified with 'add' or 'remove'.\n")
            sys.exit()

        if options.outlier_file and options.compatible_file:
            self.logger.warning("  [Warning] The 'outlier_file' and 'compatible_file' options cannot be specified at the same time.\n")
            sys.exit()

        failed_to_add = []
        failed_to_remove = []
        if options.add or options.remove:
            failed_to_add, failed_to_remove = genome_tk.modify(options.genome_file,
                                                               options.scaffold_file,
                                                               options.add,
                                                               options.remove,
                                                               options.output_genome)
        elif options.outlier_file:
            outliers = Outliers()
            outliers.remove_outliers(options.genome_file, options.outlier_file, options.output_genome)
        elif options.compatible_file:
            outliers = Outliers()
            if options.unique_only:
                outliers.add_compatible_unique(options.scaffold_file, options.genome_file, options.compatible_file, options.output_genome)
            else:
                outliers.add_compatible_closest(options.scaffold_file, options.genome_file, options.compatible_file, options.output_genome)

        if failed_to_add:
            self.logger.warning('  [Warning] Failed to add the following sequence(s):')
            for seq_id in failed_to_add:
                self.logger.warning('    %s' % seq_id)

        if failed_to_remove:
            self.logger.warning('  [Warning] Failed to remove the following sequence(s):')
            for seq_id in failed_to_remove:
                self.logger.warning('    %s' % seq_id)

        self.logger.info('')
        self.logger.info('  Modified genome written to: ' + options.output_genome)

        self.time_keeper.print_time_stamp()

    def call_genes(self, options):
        """Call genes command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [RefineM - call_genes] Identifying genes within genomes.')
        self.logger.info('*******************************************************************************')

        check_dir_exists(options.genome_nt_dir)
        make_sure_path_exists(options.output_dir)

        genome_files = self._genome_files(options.genome_nt_dir, options.genome_ext)
        if not self._check_nuclotide_seqs(genome_files):
            self.logger.warning('[Warning] All files must contain nucleotide sequences.')
            sys.exit()

        # call genes in genomes
        prodigal = Prodigal(options.cpus)
        prodigal.run(genome_files, options.output_dir)
        self.logger.info('  Genes in genomes written to: %s' % options.output_dir)

        # call genes in unbinned scaffolds
        if options.unbinned_file:
            unbinned_output_dir = os.path.join(options.output_dir, 'unbinned')
            prodigal.run([options.unbinned_file], unbinned_output_dir, meta=True)
            self.logger.info('  Genes in unbinned scaffolds written to: %s' % unbinned_output_dir)

        self.time_keeper.print_time_stamp()

    def unique(self, options):
        """Unique command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info('[RefineM - unique] Ensuring sequences are assigned to a single genome.')
        self.logger.info('*******************************************************************************')

        genome_files = self._genome_files(options.genome_nt_dir, options.genome_ext)
        if not self._check_nuclotide_seqs(genome_files):
            self.logger.warning('[Warning] All files must contain nucleotide sequences.')
            sys.exit()

        duplicates = genome_tk.unique(genome_files)

        self.logger.info('')
        if len(duplicates) == 0:
            self.logger.info('  Pass: All sequences were identified exactly once.')
        else:
            self.logger.info('  Fail: One or more sequences were observed multiple times.')

            genome_ids = sorted(duplicates.keys())
            for i in xrange(0, len(genome_ids)):
                genome_idA = genome_ids[i]

                for j in xrange(i, len(genome_ids)):
                    genome_idB = genome_ids[j]

                    dup_seq_ids = duplicates[genome_idA][genome_idB]
                    if len(dup_seq_ids) == 0:
                        continue

                    self.logger.info('')
                    if genome_idA == genome_idB:
                        self.logger.info('  There are %d sequences present more than once in %s:' % (len(dup_seq_ids), genome_idA))
                    else:
                        self.logger.info('  There are %d sequences shared between %s and %s:' % (len(dup_seq_ids), genome_idA, genome_idB))

                    for seq_id in dup_seq_ids:
                        self.logger.info('    %s' % seq_id)

        self.time_keeper.print_time_stamp()

    def bin_compare(self, options):
        """Bin compare command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info('[RefineM - bin_compare] Comparing two sets of genomes.')
        self.logger.info('*******************************************************************************')

        check_dir_exists(options.genome_nt_dir1)
        check_dir_exists(options.genome_nt_dir2)

        genomes_files1 = self._genome_files(options.genome_nt_dir1, options.genome_ext1)
        if not self._check_nuclotide_seqs(genomes_files1):
            self.logger.warning('[Warning] All files must contain nucleotide sequences.')
            sys.exit()

        genomes_files2 = self._genome_files(options.genome_nt_dir2, options.genome_ext2)
        if not self._check_nuclotide_seqs(genomes_files2):
            self.logger.warning('[Warning] All files must contain nucleotide sequences.')
            sys.exit()

        bin_comparer = BinComparer()
        bin_comparer.run(genomes_files1, genomes_files2, options.scaffold_file, options.output_file)

        self.logger.info('')
        self.logger.info('  Detailed bin comparison written to: ' + options.output_file)

        self.time_keeper.print_time_stamp()

    def unbinned(self, options):
        """Unbinned Command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [RefineM - unbinned] Identify unbinned scaffolds.')
        self.logger.info('*******************************************************************************')

        check_dir_exists(options.genome_nt_dir)

        genomes_files = self._genome_files(options.genome_nt_dir, options.genome_ext)
        if not self._check_nuclotide_seqs(genomes_files):
            self.logger.warning('[Warning] All files must contain nucleotide sequences.')
            sys.exit()

        unbinned = Unbinned()
        unbinned_seqs = unbinned.run(genomes_files, options.scaffold_file, options.min_seq_len)

        seq_io.write_fasta(unbinned_seqs, options.output_file)

        self.logger.info('')
        self.logger.info('  Unbinned scaffolds written to: ' + options.output_file)

        self.time_keeper.print_time_stamp()
        
        
    def tetra_compare(self, options):
        """Tetranucleotide comparison command"""
        self.logger.info('')
        self.logger.info('*******************************************************************************')
        self.logger.info(' [RefineM - tetra_compare] compare tetranucleotide frequencies')
        self.logger.info('*******************************************************************************')
        
        check_file_exists(options.scaffold_file)

        if not self._check_nuclotide_seqs([options.scaffold_file]):
            self.logger.warning('[Warning] Scaffold file must contain nucleotide sequences.')
            sys.exit()

        genome_files = self._genome_files(options.genome_nt_dir, options.genome_ext)
        
        if not self._check_nuclotide_seqs(genome_files):
            self.logger.warning('[Warning] All files must contain nucleotide sequences.')
            sys.exit()

        make_sure_path_exists(options.output_dir)
        
        windows=WindowGen(options.cpus)
        windows_file, links_file=windows.write_windows(options.scaffold_file,options.output_dir,options.window_size,options.gap_size)
        
        options.scaffold_file=windows_file
        print options.scaffold_file
        options.genome_nt_dir=os.path.split(windows_file)[0] #Expects one genome - the scaffolds file
        print options.genome_nt_dir
        options.links_file=links_file
        print options.links_file
        
        self.scaffold_stats(options)
        
        options.scaffold_stats_file=os.path.join(options.output_dir, 'scaffold_stats.tsv')
        print options.scaffold_stats_file
        
        self.outliers(options)
        

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""

        logging.basicConfig(format='', level=logging.INFO)

        check_dependencies(('diamond', 'ktImportText'))

        if(options.subparser_name == 'scaffold_stats'):
            print options
            self.scaffold_stats(options)
        elif(options.subparser_name == 'genome_stats'):
            self.genome_stats(options)
        elif(options.subparser_name == 'gene_profile'):
            self.gene_profile(options)
        elif(options.subparser_name == 'outliers'):
            self.outliers(options)
        elif(options.subparser_name == 'cluster'):
            self.cluster(options)
        elif(options.subparser_name == 'reference'):
            self.reference(options)
        elif(options.subparser_name == 'compatible'):
            self.compatible(options)
        elif(options.subparser_name == 'unique'):
            self.unique(options)
        elif(options.subparser_name == 'bin_compare'):
            self.bin_compare(options)
        elif(options.subparser_name == 'modify'):
            self.modify(options)
        elif(options.subparser_name == 'call_genes'):
            self.call_genes(options)
        elif(options.subparser_name == 'unbinned'):
            self.unbinned(options)
        elif (options.subparser_name == 'tetra_compare'):
            self.tetra_compare(options)
        else:
            self.logger.error('  [Error] Unknown RefineM command: ' + options.subparser_name + '\n')
            sys.exit()

        return 0
