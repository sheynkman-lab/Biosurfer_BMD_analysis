# %% Library imports
import csv
import multiprocessing as mp
import os
import re
import sys
from collections import Counter
from functools import reduce
from itertools import chain, groupby, product
from operator import attrgetter, itemgetter
from typing import TYPE_CHECKING, NamedTuple
from warnings import warn

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from biosurfer.analysis.sqtl import (
    get_cblocks_attributed_to_transcript_event,
    get_transcript_events_associated_with_junction,
    split_transcripts_on_junction_usage, get_pblocks_related_to_junction)
from biosurfer.core.alignments import (CodonAlignment, ProteinAlignment,
                                       TranscriptAlignment)
from biosurfer.core.constants import OTHER_EXCLUSIVE, SQANTI, Strand
from biosurfer.core.database import Database
from biosurfer.core.helpers import ExceptionLogger
from biosurfer.core.models.biomolecules import (GencodeTranscript, Gene,
                                                Junction, PacBioTranscript,
                                                Transcript)
from biosurfer.plots.plotting import IsoformPlot, mpatches
from IPython.display import display
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec
from more_itertools import first, only
from scipy.sparse import coo_matrix
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql.expression import and_, not_, or_, select
from tqdm import tqdm

if TYPE_CHECKING:
    from biosurfer.core.alignments import (CodonAlignmentBlock,
                                           ProteinAlignmentBlock)
    from biosurfer.core.models.features import ProteinFeature
    from biosurfer.core.splice_events import BasicTranscriptEvent

plt.switch_backend('agg')

data_dir = '../data/bone_nmd'
output_dir = '../output/bone_nmd_latest/'

# sqtlite database instance for 'containing' transcripts
db = Database('bone_cds_latest')
# sqtlite database instance for 'lacking' transcripts
db_not = Database('bone_cds_not_latest')

# %%
# Reading expression table for both lacking and containing isoforms
print('Loading normalized isoform expression table...')

expression = pd.read_csv(f'{data_dir}/NMD_mix_abund.tsv', sep=' ', index_col='rn')
timepoints = ['t0a_norm', 't0c_norm', 't2a_norm', 't2b_norm', 't2c_norm', 't4a_norm', 't4b_norm', 't4c_norm', 't10a_norm', 't10b_norm', 't10c_norm']
expression['average'] = expression[timepoints].mean(axis=1)

# %%
# Reading colocalized sQTL table

print('Loading colocalized sQTL table...')
sqtls: pd.DataFrame = pd.read_csv(f'{data_dir}/sqtl_coloc.tsv', sep=' ', nrows=None)
sqtls[['chr', 'start', 'end', 'cluster', 'gene_id']] = sqtls['phenotype_id'].str.split(':', expand=True)
sqtls['gene_id_stem'] = sqtls['gene_id'].str.split('.').str.get(0)
sqtls[['start', 'end']] = sqtls[['start', 'end']].astype(int)

# %%
mapped_gene_count = 0
gene_id_mapper = dict()
with db.get_session() as session:
    for gene_id_stem in sqtls['gene_id_stem'].unique():
        gene_id = session.execute(
            select(Gene.accession).
            where(Gene.accession.startswith(gene_id_stem))
        ).scalar()
        if gene_id:
            gene_id_mapper[gene_id_stem] = gene_id
with db_not.get_session() as session_not:
    for gene_id_stem in sqtls['gene_id_stem'].unique():
        gene_id = session_not.execute(
            select(Gene.accession).
            where(Gene.accession.startswith(gene_id_stem))
        ).scalar()
        if gene_id:
            gene_id_mapper[gene_id_stem] = gene_id

# %%
# Helper functions
junc_colors = {
    (True, False): '#FEEFF4',
    (False, True): '#F0F0FD',
    (True, True): '#F8EEFD'
}

def sortkey(transcript: 'Transcript'):
    return -abundance(transcript), getattr(transcript, 'sqanti', SQANTI.OTHER)

def abundance(transcript: 'Transcript') -> float:
    try:
        # print("Transcript name-abundance", transcript, expression.average[transcript.accession])
        return expression.average[transcript.accession]
    except KeyError:
        return 0.0

class CBlockEventTuple(NamedTuple):
    cblock: 'CodonAlignmentBlock'
    event: 'BasicTranscriptEvent'
    anchor: 'Transcript'
    other: 'Transcript'

class BlockTuple(NamedTuple):
    pblock: 'ProteinAlignmentBlock'
    cblock: 'CodonAlignmentBlock'
    events: tuple['BasicTranscriptEvent', ...]
    anchor: 'Transcript'
    other: 'Transcript'

class FeatureTuple(NamedTuple):
    feature: 'ProteinFeature'
    cblock: 'CodonAlignmentBlock'
    affects_whole: bool
    anchor: 'Transcript'
    other: 'Transcript'

def get_augmented_sqtl_record(row):
    chr, gene_id, start, end, cluster, pval, slope, maf = row
    try:
        gene_id = gene_id_mapper[gene_id]
    except KeyError:
        print("Gene from sqtl coloc missing in expression table", gene_id)
        return None
    
    with db.get_session() as session, db_not.get_session() as session_not:
        gene = Gene.from_accession(session, gene_id)
        gene_not = Gene.from_accession(session_not, gene_id)
        if not gene:
            print("Gene from sqtl coloc missing in db cds", gene)
            return None
            # Modify single base
        if not gene_not:
            print("Gene from sqtl coloc missing in db cds not", gene_not)
            return None

        start += 1
        end -= 1
        strand = gene.strand
        if strand is Strand.MINUS:
            start, end = end, start
        junc = Junction.from_coordinates(chr, strand, start, end)
        # junc_not = Junction.from_coordinates(chr, strand_not, start, end)
        junc_info = {
            'chr': chr,
            'strand': str(strand),
            'donor': start,
            'acceptor': end,
            'gene': gene.name, 
            'cluster': cluster,
            'pval': pval,
            'slope': slope,
            'maf': maf
        }
        junc_info_not = {
            'chr': chr,
            'strand': str(strand),
            'donor': start,
            'acceptor': end,
            'gene': gene_not.name,
            'pval': pval,
            'slope': slope,
            'maf': maf
        }

        junc_info.update(junc_info_not)
        #pb_transcripts = {tx for tx in gene.transcripts if tx.accession in expression.index and any(orf.has_stop_codon for orf in tx.orfs)}
        # Since transcripts are protein coding already, no need for finding stop codin in orf + all transcripts are containing

        # Splitting containing and lacking using two database CDS.gtf and CDS_not_exact.gtf

        pb_transcripts = {tx for tx in gene.transcripts if tx.accession in expression.index}
        if not pb_transcripts:
            print("Gene from sqtl coloc missing pb transcripts in db cds", gene)
            return None

        pb_transcripts_not = {tx for tx in gene_not.transcripts if tx.accession in expression.index}
        if not pb_transcripts_not:
            print("Gene from sqtl coloc missing pb transcripts in db cds not", gene_not)
            # return None

        # all_pb_transcripts = pb_transcripts.union(pb_transcripts_not)
        # lacking, containing = split_transcripts_on_junction_usage(junc, filter(attrgetter('orfs'), all_pb_transcripts))
        containing = pb_transcripts
        print("PacBio tx containing :", containing)
        containing = sorted(containing, key=sortkey)
        lacking = pb_transcripts_not
        lacking = sorted(lacking, key=sortkey)
        print("PacBio tx in lacking :", lacking)
        junc_info['containing'] = len(containing)
        junc_info['lacking'] = len(lacking)

        if not containing:
            return None
        
        gc_transcripts = sorted((tx for tx in gene.transcripts if isinstance(tx, GencodeTranscript)), key=attrgetter('appris'), reverse=True)
        anchor_tx = gc_transcripts[0]

        pairs = list(product(lacking, containing))

        abundance_denom = sum(abundance(tx) for tx in containing) * sum(abundance(tx) for tx in lacking)
        weights = {
            (tx1, tx2): abundance(tx1) * abundance(tx2) / abundance_denom
            for tx1, tx2 in pairs
        }
        if weights and abs(error := sum(weights.values()) - 1) > 2**-8:
            warn(f'Weights add up to 1{error:+.3e}')

        top_pair, junc_info['highest_pair_weight'] = max(weights.items(), key=itemgetter(1), default=(None, None))

        # calculate fraction of pairs where one isoform is NMD and other is not
        junc_info['NMD'] = None

        for tx1, tx2 in pairs:
            if tx1.primary_orf is None:
                print("tx1-tx2 orf", tx1)
            elif tx2.primary_orf is None:
                print("tx1-tx2 orf", tx2)

        with ExceptionLogger(f'Error for {gene.name} {junc}'):
            nmd_sum = 0
            for tx1, tx2 in pairs:
                try:
                    # if (tx1.primary_orf is not None and tx2.primary_orf is not None):
                        # ORFs should have a stop codon, region should be div by 3
                        # ORFs with stop codons at least 50 bp upstream of the last splice site in the mature transcript
                        # (i.e. the beginning of the last exon) are considered candidates for nonsense-mediated decay (NMD)
                    
                    nmd_value = float(tx1.primary_orf.nmd ^ tx2.primary_orf.nmd) * weights[tx1, tx2]
                    # print("NMD vals + weight:",tx1.primary_orf.nmd, tx2.primary_orf.nmd, weights[tx1, tx2], nmd_value)
                    nmd_sum += nmd_value
                except:
                    print("No primary orf found for: ", tx2)
            junc_info['NMD'] = nmd_sum

        tx_aln_to_events: dict['TranscriptAlignment', set['BasicTranscriptEvent']] = dict()
        for anchor, other in pairs:
            with ExceptionLogger(f'Error for {anchor} and {other}'):
                anchor.nucleotides, other.nucleotides
                print("Other transcript donor and acceptor :", other, other.start, other.stop)
                tx_aln = TranscriptAlignment.from_transcripts(anchor, other)

                events = set(get_transcript_events_associated_with_junction(junc, tx_aln))
                tx_aln_to_events[tx_aln] = events
        cblock_event_tuples = [
            CBlockEventTuple(cblock, event, tx_aln.anchor, tx_aln.other)
            for tx_aln, events in tx_aln_to_events.items()
            for event in events
            for cblock in get_cblocks_attributed_to_transcript_event(event, CodonAlignment.from_proteins(tx_aln.anchor.protein, tx_aln.other.protein))
        ]
        cblock_event_tuples.sort(key=lambda cet: (cet.anchor.name, cet.other.name, cet.cblock, cet.event.start, cet.event.stop))
        block_tuples = [
            BlockTuple(
                ProteinAlignment.from_proteins(anchor.protein, other.protein).cblock_to_pblock[cblock],
                cblock,
                tuple(cet.event for cet in group),
                anchor, other
            )
            for (anchor, other, cblock), group in groupby(cblock_event_tuples, key=lambda cet: (cet.anchor, cet.other, cet.cblock))
        ]

        feature_tuples = list(chain(
            (
                FeatureTuple(
                    feature,
                    bt.cblock,
                    bt.cblock.anchor_range.start <= feature.protein_start - 1 and feature.protein_stop <= bt.cblock.anchor_range.stop,
                    bt.anchor,
                    bt.other
                )
                for bt in block_tuples if (bt.anchor, bt.other) == top_pair
                for feature in bt.anchor.protein.features if feature.protein_start - 1 < bt.cblock.anchor_range.stop and bt.cblock.anchor_range.start < feature.protein_stop
            ),
            (
                FeatureTuple(
                    feature,
                    bt.cblock,
                    bt.cblock.other_range.start <= feature.protein_start - 1 and feature.protein_stop <= bt.cblock.other_range.stop,
                    bt.anchor,
                    bt.other
                )
                for bt in block_tuples if (bt.anchor, bt.other) == top_pair
                for feature in bt.other.protein.features if feature.protein_start - 1 < bt.cblock.other_range.stop and bt.cblock.other_range.start < feature.protein_stop
            )
        ))
        # Calculate the average delta_aa for a set of block_tuples by grouping them based on their anchor and other
        # properties, and then weighting their delta_aa values based on their relative importance
        delta_aa_weights = tuple(zip(*(
            (sum(bt.cblock.delta_length for bt in group), weights[pair])
            for pair, group in groupby(block_tuples, key=lambda bt: (bt.anchor, bt.other))
        )))
        if delta_aa_weights:
            junc_info['avg_delta_aa'] = np.average(delta_aa_weights[0], weights=delta_aa_weights[1])
        else:
            junc_info['avg_delta_aa'] = None

        # calculate event frequencies
        event_freqs = reduce(
            Counter.__add__,
            (
                Counter({k: v*weight for k, v in counts.items()})
                for counts, weight in (
                    (Counter(events), weights[tx_aln.anchor, tx_aln.other])
                    for tx_aln, events in tx_aln_to_events.items()
                )
            ),
            Counter()
        )
    
    return junc_info

rows = [(row.chr, row.gene_id_stem, row.start, row.end, row.cluster, row.pval_nominal, row.slope, row.maf) for row in sqtls.itertuples()]
records = []

tx_no_porf = []
with tqdm(desc='Analyzing sQTL junctions', total=sqtls.shape[0], file=sys.stdout, unit='junctions') as t:
    with tqdm(desc='Annotated junctions', file=sys.stdout, unit='junctions', miniters=1) as t2:
        for row in rows:
            print(row)
            result = get_augmented_sqtl_record(row)
            t.update()
            if result:
                records.append(result)
                t2.update()

sqtls_augmented = pd.DataFrame.from_records(records)
display(sqtls_augmented)
print(f'Annotated {sqtls_augmented.shape[0]} sQTL junctions')
sqtls_augmented.to_csv(f'{output_dir}/coloc_sqtls_annotated.tsv', sep='\t', index=False)

# %%
