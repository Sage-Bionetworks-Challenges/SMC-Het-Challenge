#!/usr/bin/env python2
from __future__ import print_function

# Requires PyVCF. To install: pip3 install pyvcf
import vcf
import argparse
from collections import namedtuple, defaultdict, OrderedDict
import random

class VariantParser(object):
  def __init__(self):
    # Child classes must give the following variables sensible values in
    # constructor so that list_variants() works subsequently.
    self._cnvs = None
    self._vcf_filename = None

  def list_variants(self):
    variants = self._filter(self._vcf_filename)
    variants_and_reads = []
    for variant in variants:
      ref_reads, total_reads = self._calc_read_counts(variant)
      variants_and_reads.append((variant, ref_reads, total_reads))
    return variants_and_reads

  def _calc_read_counts(self, variant):
    raise Exception('Not implemented -- use child class')

  def _parse_vcf(self, vcf_filename):
    vcfr = vcf.Reader(filename=vcf_filename)
    records = list(vcfr)
    return records

  def _is_good_chrom(self, chrom):
    chrom = chrom.lower()
    if chrom.startswith('chrun') or chrom.endswith('_random'):
      # Variant unmapped ('chrUn') or mapped to fragmented chromosome
      # ('_random'), so ignore.
      return False

    # Mitochondrial are weird, so ignore them.
    if chrom == 'chrm':
      return False
    # Sex chromosomes difficult to deal with, as expected frequency depends on
    # whether patient is male or female, so ignore them for now.
    #
    # TODO: properly deal with sex chromsomes.
    if chrom in ('chrx', 'chry'):
      return False

    return True

  def _does_variant_pass_filters(self, variant):
    if variant.FILTER is None:
      return True
    if len(variant.FILTER) > 0:
      # Variant failed one or more filters.
      return False
    return True

  def _filter(self, vcf_filename):
    variants = []

    all_variants = self._parse_vcf(vcf_filename)

    for variant in all_variants:
      if not self._is_good_chrom(variant.CHROM):
        continue
      if not self._does_variant_pass_filters(variant):
        continue

      variants.append(variant)

    variants.sort(key = lambda v: (v.CHROM, v.POS))
    return variants

class MutectParser(VariantParser):
  '''
  Works with Battenberg-formatted CNV calls.
  '''
  def __init__(self, vcf_filename):
    self._vcf_filename = vcf_filename

  def _calc_read_counts(self, variant):
    # Currently hardcodes tumour sample as the second column.
    # Might not always be true
    ref_reads = variant.samples[1]['AD'][0]
    variant_reads = variant.samples[1]['AD'][1]
    total_reads = ref_reads + variant_reads

    return (ref_reads, total_reads)

class SangerParser(VariantParser):
  '''
  Works with Battenberg-formatted CNV calls.
  '''
  def __init__(self, vcf_filename):
    self._vcf_filename = vcf_filename

  def _find_ref_and_variant_nt(self, variant):
    # Try most probable genotype called by CaVEMan. If can't find variant nt in
    # it, try the second most probable genotype.
    genotypes = [variant.INFO['TG'][0], variant.INFO['SG'][0]]
    variant_set = set()

    while len(variant_set) == 0:
      if len(genotypes) == 0:
        raise Exception('No more genotypes to find variant_nt in for %s' % variant)
      genotype = genotypes.pop(0)
      normal_genotype, tumor_genotype = genotype.split('/')
      # TODO: We throw out hetero germline. Deal with this later.
      if normal_genotype[0] != normal_genotype[1]:
        print('Ignoring heterozygous normal genotype %s' % normal_genotype, file=sys.stderr)
      reference_nt = normal_genotype[0]
      variant_set = set(tumor_genotype) - set(reference_nt)

    variant_nt = variant_set.pop()
    return (reference_nt, variant_nt)

  def _calc_read_counts(self, variant):
    normal = variant.genotype('NORMAL')
    tumor = variant.genotype('TUMOUR')

    reference_nt, variant_nt = self._find_ref_and_variant_nt(variant)
    tumor_reads = {
      'forward': {
        'A': int(tumor['FAZ']),
        'C': int(tumor['FCZ']),
        'G': int(tumor['FGZ']),
        'T': int(tumor['FTZ']),
      },
      'reverse': {
        'A': int(tumor['RAZ']),
        'C': int(tumor['RCZ']),
        'G': int(tumor['RGZ']),
        'T': int(tumor['RTZ']),
      },
    }

    ref_reads = tumor_reads['forward'][reference_nt] + tumor_reads['reverse'][reference_nt]
    # For now, variant reads are defined as only the non-reference nucleotide in
    # the inferred tumor SNP. We ignore reads of a third or fourth base.
    variant_reads = tumor_reads['forward'][variant_nt] + tumor_reads['reverse'][variant_nt]
    total_reads = ref_reads + variant_reads

    return (ref_reads, total_reads)

class OncoscanParser(VariantParser):
  '''
  Works with prostate data from Oncoscan.
  '''
  def __init__(self, vcf_filename, oncoscan_filename):
    self._vcf_filename = vcf_filename
    self._cnvs = self._parse_oncoscan(oncoscan_filename)

  def _calc_read_counts(self, variant):
    tumor = variant.genotype('TUMOR')
    ref_forward, ref_reverse, alt_forward, alt_reverse = [int(e) for e in tumor['DP4']]
    ref_reads = ref_forward + ref_reverse
    total_reads = ref_reads + alt_forward + alt_reverse
    return (ref_reads, total_reads)

  def _parse_oncoscan(self, oncoscan_filename):
    cnvs = defaultdict(list)
    with open(oncoscan_filename) as cnv_file:
      for line in cnv_file:
        chr, start, end, delta = line.strip().split()
        cnvs[chr].append((int(start), int(end)))
    return cnvs

class BattenbergParser(object):
  def __init__(self, bb_filename):
    self._bb_filename = bb_filename

  def _is_normal_cn(self, cnv):
    return cnv['frac'] == cnv['nmaj'] == cnv['nmin'] == 1

  def _parse_float(self, val):
    try:
      return float(val)
    except ValueError:
      return None

  def _parse_int(self, val):
    try:
      return int(val)
    except ValueError:
      return None

  def _compute_cn(self, cnv1, cnv2):
    cn1 = (cnv1['nmaj'] + cnv1['nmin']) * cnv1['frac']
    if cnv2:
      cn2 = (cnv2['nmaj'] + cnv2['nmin']) * cnv2['frac']
    else:
      cn2 = 0
    total_cn = cn1 + cn2
    return total_cn

  def parse(self):
    cn_normal = defaultdict(list)
    cn_abnormal = defaultdict(list)
    pval_threshold = 0.05

    with open(self._bb_filename) as bbf:
      header = bbf.next()
      for line in bbf:
        fields = line.strip().split()
        chrom = fields[1].lower()
        start = int(fields[2])
        end = int(fields[3])
        pval = float(fields[5])

        cnv1 = OrderedDict()
        cnv1['start'] = start
        cnv1['end'] = end
        cnv1['nmaj'] = int(fields[8])
        cnv1['nmin'] = int(fields[9])
        cnv1['frac'] = float(fields[10])

        cnv2 = None
        # Stefan's comment on p values: The p-values correspond "to whether a
        # segment should be clonal or subclonal copynumber. We first fit a
        # clonal copynumber profile for the whole sample and then perform a
        # simple two-sided t-test twhere the null hypothesis is: A particular
        # segment is clonal. And the alternative: It is subclonal."
        #
        # Thus: if t-test falls below significance threshold, we push cnv1 to
        # clonal frequency.
        if pval <= pval_threshold:
          cnv2 = OrderedDict()
          cnv2['start'] = start
          cnv2['end'] = end
          cnv2['nmaj'] = int(fields[11])
          cnv2['nmin'] = int(fields[12])
          cnv2['frac'] = float(fields[13])
        else:
          cnv1['frac'] = 1.0

        if cnv2 is None:
          if self._is_normal_cn(cnv1):
            cn_normal[chrom].append(cnv1)
          else:
            cn_abnormal[chrom].append(cnv1)
        else:
          # cnv1 must be abnormal if cnv2 exists -- if cnv1 was normal, cnv2
          # would be "NA". cnv2 may be normal or abnormal -- it could, for
          # example, be normal, indicaating that the CNV (represented by the
          # abnormal cnv1) occurred in only some cells.
          assert not self._is_normal_cn(cnv1)
          for cnv in (cnv1, cnv2):
            if not (cnv['nmin'] == cnv['nmaj'] == 1):
              cn_abnormal[chrom].append(cnv)

    return (cn_normal, cn_abnormal)

class CnvFormatter(object):
  def _find_overlapping_variants(self, chrom, cnv, variants):
    overlapping = []

    start = cnv['start']
    end = cnv['end']
    for variant in variants:
      if chrom.lower() == variant['chrom'].lower():
        if start <= variant['pos'] <= end:
          overlapping.append(variant['ssm_id'])
    return overlapping

  def _calc_ref_reads(self, pop_frac, cellularity, total_reads):
    tumor_cells_frac = cellularity * pop_frac
    vaf = tumor_cells_frac / 2
    ref_reads = int((1 - vaf) * total_reads)
    return ref_reads

  def _calc_total_reads(self, pop_frac, cellularity, locus_start, locus_end, new_cn, read_depth, read_length):
    # Proportion of all cells carrying CNV.
    p = cellularity * pop_frac
    if new_cn == 2:
      # If no net change in copy number -- e.g., because (major, minor) went
      # from (1, 1) to (2, 0) -- force the delta_cn to be 1.
      delta_cn = 1.
      no_net_change = True
    else:
      delta_cn = float(new_cn - 2)
      no_net_change = False

    region_length = locus_end - locus_start + 1
    fn = (read_depth * region_length) / read_length
    d = (delta_cn**2 / 4) * (fn * p * (2 - p)) / (1 + (delta_cn  * p) / 2)

    if no_net_change:
      # If no net change in CN occurred, the estimate was just based on BAFs,
      # meaning we have lower confidence in it. Indicate this lack of
      # confidence via d by multiplying it by (read length / distance between
      # common SNPs), with the "distance between common SNPs" taken to be 1000 bp.
      d *= (read_length / 1000.)

    # Cap at 1e6 * read_depth.
    return int(round(min(d, 1e6 * read_depth)))

  def _format_overlapping_variants(self, variants, maj_cn, min_cn):
      variants = [(ssm_id, str(min_cn), str(maj_cn)) for ssm_id in variants]
      return variants

  def _format_cnvs(self, cnvs, variants, cellularity, read_depth):
    # This is just a guess. Setting this properly (e.g., from the BAM files?) would be preferable.
    read_length = 100
    print('Estimated read depth: %s' % read_depth)

    for chrom, chrom_cnvs in cnvs.items():
      for cnv in chrom_cnvs:
        overlapping_variants = self._find_overlapping_variants(chrom, cnv, variants)
        total_reads = self._calc_total_reads(
          cnv['frac'],
          cellularity,
          cnv['start'],
          cnv['end'],
          cnv['nmaj'] + cnv['nmin'],
          read_depth,
          read_length
        )
        yield {
          'frac': cnv['frac'],
          'ref_reads': self._calc_ref_reads(cnv['frac'], cellularity, total_reads),
          'total_reads': total_reads,
          'overlapping_variants': self._format_overlapping_variants(overlapping_variants, cnv['nmaj'], cnv['nmin']),
        }

  def _merge_variants(self, cnv1, cnv2):
    cnv1_variant_names = set([v[0] for v in cnv1['overlapping_variants']])
    for variant in cnv2['overlapping_variants']:
      variant_name = variant[0]
      if variant_name not in cnv1_variant_names:
        cnv1['overlapping_variants'].append(variant)
      else:
        # If variant already in cnv1's list, ignore it. This should only occur
        # if two subclonal CNVs have close to 0.5 frequency each. In this case,
        # we lose information about major/minor status of the cnv2 relative to
        # its SSMs.
        print('%s already in %s' % (variant, cnv1['cnv_id']))

  def format_and_merge_cnvs(self, cnvs, variants, cellularity, read_depth):
    formatted = list(self._format_cnvs(cnvs, variants, cellularity, read_depth))
    formatted.sort(key = lambda f: f['ref_reads'])
    delta = 0.001

    merged, formatted = formatted[:1], formatted[1:]
    merged[0]['cnv_id'] = 'c0'
    counter = 1

    for current in formatted:
      last = merged[-1]
      last_frac = float(last['ref_reads']) / last['total_reads']
      current_frac = float(current['ref_reads']) / current['total_reads']

      # Only merge CNVs if they're clonal. If they're subclonal, leave them
      # free to move around the tree.
      if current['frac'] == last['frac'] == 1.0 and abs(last_frac - current_frac) <= delta:
        # Merge the CNVs.
        last['ref_reads'] += current['ref_reads']
        last['total_reads'] += current['total_reads']
        self._merge_variants(last, current)
      else:
        # Do not merge the CNVs.
        current['cnv_id'] = 'c%s' % counter
        merged.append(current)
        counter += 1

    return merged

class VariantFormatter(object):
  def __init__(self):
    self._counter = 0

  def _split_types(self, genotype):
    types = [int(e) for e in genotype.split('/')]
    if len(types) != 2:
      raise Exception('Not diploid: %s' % types)
    return types

  def _calc_ref_freq(self, ref_genotype, error_rate):
    types = self._split_types(ref_genotype)
    num_ref = len([t for t in types if t == 0])
    freq = (num_ref / 2) - error_rate
    if freq < 0:
      freq = 0.0
    if freq > 1:
      raise Exception('Nonsensical frequency: %s' % freq)
    return freq

  def format_variants(self, variant_list):
    for variant, ref_reads, total_reads in variant_list:
      ssm_id = 's%s' % self._counter
      variant_name = '%s_%s' % (variant.CHROM, variant.POS)

      # TODO: switch back to using calc_ref_freq() when we no longer want mu_r
      # and mu_v fixed.
      # This is mu_r in PhyloWGS.
      expected_ref_freq = 0.999
      if variant.CHROM.lower() in ('chrx', 'chry', 'chrm'):
        # Haploid, so should only see non-variants when sequencing error
        # occurred. Note that chrY and chrM are always haploid; chrX is haploid
        # only in men, so script must know sex of patient to choose correct
        # value. Currently, I just assume that all data comes from men.
        #
        # This is mu_v in PhyloWGS.
        expected_var_freq = 0.001
      else:
        # Diploid, so should see variants in (0.5 - error_rate) proportion of
        # reads.
        #
        # This is mu_v in PhyloWGS.
        expected_var_freq = 0.499

      yield {
        'ssm_id': ssm_id,
        'chrom': variant.CHROM,
        'pos': variant.POS,
        'variant_name': variant_name,
        'ref_reads': ref_reads,
        'total_reads': total_reads,
        'expected_ref_freq': expected_ref_freq,
        'expected_var_freq': expected_var_freq,
      }
      self._counter += 1

def restricted_float(x):
  x = float(x)
  if x < 0.0 or x > 1.0:
    raise argparse.ArgumentTypeError('%r not in range [0.0, 1.0]' % x)
  return x

class VariantAndCnvGroup(object):
  def __init__(self):
    self._variants = None
    self._normal_cn = None
    self._abnormal_cn = None

  def add_variants(self, variants_and_reads):
    self._variants_and_reads  = variants_and_reads

  def add_cnvs(self, normal_cn, abnormal_cn):
    self._normal_cn, self._abnormal_cn = normal_cn, abnormal_cn

  def has_cnvs(self):
    return self._abnormal_cn is not None

  def retain_only_variants_in_normal_cn_regions(self):
    filtered = []

    if self._normal_cn is None:
      raise Exception('Normal CN regions not yet provided')
    for cn_chrom, norm_cn_regions in self._normal_cn.items():
      for norm_cn in norm_cn_regions:
        for variant, ref_reads, total_reads in self._variants_and_reads:
          if variant.CHROM.lower() == cn_chrom.lower() and norm_cn['start'] <= variant.POS <= norm_cn['end']:
            filtered.append((variant, ref_reads, total_reads))

    # Use a set, as the same variant may overlap multiple CNVs, and so get
    # added multiple times.
    abnormal_variants = set()
    for cn_chrom, abnorm_cn_regions in self._abnormal_cn.items():
      for abnorm_cn in abnorm_cn_regions:
        for variant, ref_reads, total_reads in self._variants_and_reads:
          if variant.CHROM.lower() == cn_chrom.lower() and abnorm_cn['start'] <= variant.POS <= abnorm_cn['end']:
            abnormal_variants.add(variant)
    print('normal_cn_variants=%s abnormal_cn_variants=%s omitted_variants=%s total_variants=%s' % (
      len(filtered),
      len(abnormal_variants),
      len(self._variants_and_reads) - len(filtered) - len(abnormal_variants),
      len(self._variants_and_reads),
    ))

    self._variants_and_reads = filtered

  def subsample_variants(self, sample_size):
    random.shuffle(self._variants_and_reads)
    self._variants_and_reads = self._variants_and_reads[:sample_size]
    self._variants_and_reads.sort(key = lambda v: (v[0].CHROM, v[0].POS))

  def write_variants(self, outfn):
    formatter = VariantFormatter()
    self._variants = list(formatter.format_variants(self._variants_and_reads))

    with open(outfn, 'w') as outf:
      print('\t'.join(('id', 'gene', 'a', 'd', 'mu_r', 'mu_v')), file=outf)
      for variant in self._variants:
        vals = (
          'ssm_id',
          'variant_name',
          'ref_reads',
          'total_reads',
          'expected_ref_freq',
          'expected_var_freq',
        )
        vals = [variant[k] for k in vals]
        print('\t'.join([str(v) for v in vals]), file=outf)

  def _estimate_read_depth(self):
    read_sum = 0
    for variant, ref_reads, total_reads in self._variants_and_reads:
      read_sum += total_reads
    return float(read_sum) / len(self._variants_and_reads)

  def write_cnvs(self, cellularity, outfn):
    with open(outfn, 'w') as outf:
      print('\t'.join(('cnv', 'a', 'd', 'ssms')), file=outf)
      formatter = CnvFormatter()
      for cnv in formatter.format_and_merge_cnvs(self._abnormal_cn, self._variants, cellularity, self._estimate_read_depth()):
        overlapping = [','.join(o) for o in cnv['overlapping_variants']]
        vals = (
          cnv['cnv_id'],
          str(cnv['ref_reads']),
          str(cnv['total_reads']),
          ';'.join(overlapping),
        )
        print('\t'.join(vals), file=outf)

def main():
  parser = argparse.ArgumentParser(
    description='Create SSM input file for PhyloWGS from VCF and CNV data',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  # TODO: re-enable --error-rate argument once we decide that values of mu_r
  # and mu_v shouldn't always be fixed.
  #parser.add_argument('--error-rate', dest='error_rate', type=float, default=0.001,
    #help='Expected error rate of sequencing platform')
  parser.add_argument('-s', '--sample-size', dest='sample_size', type=int,
    help='Subsample reported number of SSMs')
  parser.add_argument('-b', '--battenberg', dest='battenberg',
    help='Path to Battenberg-formatted list of CNVs')
  parser.add_argument('--only-normal-cn', dest='only_normal_cn', action='store_true', default=False,
      help='Only output variants lying in normal CN regions. Do not output CNV data directly.')
  parser.add_argument('--output-cnvs', dest='output_cnvs', default='cnv_data.txt',
    help='Output destination for CNVs')
  parser.add_argument('--output-variants', dest='output_variants', default='ssm_data.txt',
    help='Output destination for variants')
  parser.add_argument('-c', '--cellularity', dest='cellularity', type=restricted_float, default=1.0,
    help='Fraction of sample that is cancerous rather than somatic')
  parser.add_argument('-v', '--variant-type', dest='input_type', required=True, choices=('sanger', 'oncoscan', 'mutect'),
      help='Type of VCF file')
  parser.add_argument('vcf_file')
  args = parser.parse_args()

  # Fix random seed to ensure same set of SSMs chosen on each invocation.
  random.seed(1)

  if args.input_type == 'sanger':
    grouper = VariantAndCnvGroup()
    variant_parser = SangerParser(args.vcf_file)
    variants_and_reads = variant_parser.list_variants()
    grouper.add_variants(variants_and_reads)

    if args.battenberg:
      cnv_parser = BattenbergParser(args.battenberg)
      normal_cn, abnormal_cn = cnv_parser.parse()
      grouper.add_cnvs(normal_cn, abnormal_cn)

    if args.only_normal_cn:
      grouper.retain_only_variants_in_normal_cn_regions()

    if args.sample_size:
      grouper.subsample_variants(args.sample_size)

    grouper.write_variants(args.output_variants)
    if not args.only_normal_cn and grouper.has_cnvs():
      grouper.write_cnvs(args.cellularity, args.output_cnvs)
  
  elif args.input_type == 'mutect':
    grouper = VariantAndCnvGroup()
    variant_parser = MutectParser(args.vcf_file)
    variants_and_reads = variant_parser.list_variants()
    grouper.add_variants(variants_and_reads)

    if args.battenberg:
      cnv_parser = BattenbergParser(args.battenberg)
      normal_cn, abnormal_cn = cnv_parser.parse()
      grouper.add_cnvs(normal_cn, abnormal_cn)

    if args.only_normal_cn:
      grouper.retain_only_variants_in_normal_cn_regions()

    if args.sample_size:
      grouper.subsample_variants(args.sample_size)

    grouper.write_variants(args.output_variants)
    if not args.only_normal_cn and grouper.has_cnvs():
      grouper.write_cnvs(args.cellularity, args.output_cnvs)
      
  elif args.input_type == 'oncoscan':
    raise Exception('Various changes necessary until this works again')
    #parser = OncoscanParser(args.vcf_file, args.cnv_file)


if __name__ == '__main__':
  main()