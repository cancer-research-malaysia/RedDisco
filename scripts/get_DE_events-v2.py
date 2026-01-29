#!/usr/bin/env python
"""
Differential RNA Editing Analysis

Identifies differentially edited sites between control and disease samples
using Mann-Whitney U test or linear model with coverage stability analysis.

Metadata File Format (CSV):
    SampleName,Condition,/absolute/path/to/outTable.txt
    
    Example:
    Sample1,CTRL,/data/Sample1/editing/outTable_12345.txt
    Sample2,CTRL,/data/Sample2/results/outTable_67890.txt
    Patient_A,DIS,/data/Patient_A/output/outTable_99999.txt

Usage:
    python get_DE_events.py -input_file metadata.csv [OPTIONS]

Options:
    -input_file FILE         Metadata CSV file (required)
    -c INT                   Minimum coverage (default: 10)
    -f FLOAT                 Minimum editing frequency (default: 0.1)
    -mts FLOAT               Min percentage of samples per category (default: 50.0)
    -sig yes/no              Return only significant events (default: no)
    -cpval INT               P-value correction: 0=none, 1=Bonferroni, 2=BH (default: 0)
    -linear                  Enable linear model with stability analysis
    -h, --help               Show this help message
"""

import os
import sys
import argparse
from scipy import stats
from scipy.stats import wilcoxon, mannwhitneyu, fisher_exact
import numpy as np
import pandas as pd
import math

# ============================================================================
# ARGUMENT PARSING
# ============================================================================

parser = argparse.ArgumentParser(
    description='Differential RNA editing analysis between conditions',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__
)

parser.add_argument("-c", action='store', dest='min_coverage',
                    type=int, default=10, help='Minimum coverage-q30 (default: 10)')
parser.add_argument("-cpval", action='store', dest='pvalue_correction',
                    type=int, default=0, help='P-value correction: 0=none, 1=Bonferroni, 2=BH (default: 0)')
parser.add_argument("-input_file", action='store', dest='samples_informations_file',
                    type=str, required=True, help='CSV file with format: SampleName,Condition,/path/to/outTable.txt')
parser.add_argument("-f", action='store', dest='min_edit_frequency',
                    type=float, default=0.1, help='Minimum editing frequency (default: 0.1)')
parser.add_argument("-mts", action='store', dest='min_sample_testing',
                    type=float, default=50.0, help="Min percentage of each sample category (default: 50.0)")
parser.add_argument("-sig", action='store', dest='only_significant',
                    type=str, default='no', help='Return only significant editing events (default: no)')
parser.add_argument("-linear", action='store_true', dest='enable_linear_model',
                    help='Enable linear model with coverage stability analysis')

args = parser.parse_args()

min_coverage = args.min_coverage
min_edit_frequency = args.min_edit_frequency
min_sample_testing = args.min_sample_testing
only_significants = args.only_significant
pvalue_correction = args.pvalue_correction
samples_informations_file = args.samples_informations_file
enable_linear_model = args.enable_linear_model

# ============================================================================
# VALIDATE INPUT FILE
# ============================================================================

if not os.path.exists(samples_informations_file):
    sys.exit('ERROR: Metadata file not found: %s' % samples_informations_file)

print 'Reading metadata file: %s' % samples_informations_file

# ============================================================================
# LINEAR MODEL FUNCTION (UNCHANGED ALGORITHM)
# ============================================================================

def call_differential_editing_sites(config_file):
    """Linear model with coverage stability analysis - ALGORITHM UNCHANGED"""
    stability_value = 0.03
    min_disease_people = 5
    min_control_people = 5
    min_disease_people_5_cov = 10
    min_control_people_5_cov = 10
    editing_file = './temp.csv'
    output_file = './editing_sites.with_stats_converted_disease.csv'
    
    # Read in files
    editing_table = pd.read_csv(editing_file, sep='\t')
    config_table = pd.read_csv(config_file, sep=',', skiprows=1, header=None)
    all_people = config_table[0]
    disease_people = config_table[0][config_table[1] == "DIS"].reset_index(drop=True)
    control_people = config_table[0][config_table[1] == "CTRL"].reset_index(drop=True)

    # Get editing and coverage tables
    edit_level_table = editing_table[all_people]

    def get_editing_levels_for_cov_table(i):
        info = i.astype(str).str.split(pat="\\^")
        editing_levels = info.apply(lambda x: float('nan') if x[0] == "nan" else x[2])
        return editing_levels
    
    cov_table = edit_level_table.apply(get_editing_levels_for_cov_table)
    cov_table = cov_table.apply(lambda x: pd.to_numeric(x))

    def get_editing_levels(i):
        info = i.astype(str).str.split(pat="\\^")
        editing_levels = info.apply(lambda x: float('nan') if x[0] == "nan" else x[0])
        return editing_levels
    
    edit_level_table = edit_level_table.apply(get_editing_levels)
    edit_level_table = edit_level_table.apply(lambda x: pd.to_numeric(x))

    # Initialize arrays
    coverage_threshold_used = np.repeat(0., edit_level_table.shape[0])
    stability_based_on = np.repeat(0., edit_level_table.shape[0])
    stable_mean_disease_editing_level = np.repeat(0., edit_level_table.shape[0])
    stable_std_dev_disease_editing_level = np.repeat(0., edit_level_table.shape[0])
    stable_mean_control_editing_level = np.repeat(0., edit_level_table.shape[0])
    stable_std_dev_control_editing_level = np.repeat(0., edit_level_table.shape[0])
    stable_number_disease_with_at_least_min_coverage = np.repeat(0., edit_level_table.shape[0])
    stable_number_disease_nonzero_editing_and_min_coverage = np.repeat(0., edit_level_table.shape[0])
    stable_disease_prevalence = np.repeat(0., edit_level_table.shape[0])
    stable_number_control_with_at_least_min_coverage = np.repeat(0., edit_level_table.shape[0])
    stable_number_control_nonzero_editing_and_min_coverage = np.repeat(0., edit_level_table.shape[0])
    stable_control_prevalence = np.repeat(0., edit_level_table.shape[0])
    stable_total_number_individuals_nonzero_editing_and_min_coverage = np.repeat(0., edit_level_table.shape[0])
    stable_mann_whitney_p_value = np.repeat(0., edit_level_table.shape[0])
    stable_editing_level_effect_size = np.repeat(0., edit_level_table.shape[0])
    stable_frequency_fishers_p_value = np.repeat(0., edit_level_table.shape[0])
    stable_frequency_OR = np.repeat(0., edit_level_table.shape[0])
    stable_prevalence_effect_size = np.repeat(0., edit_level_table.shape[0])

    # Process each row
    for i in range(0, edit_level_table.shape[0]):
        print i
        disease_edit_row = edit_level_table.loc[i, disease_people]
        control_edit_row = edit_level_table.loc[i, control_people]
        disease_cov_row = cov_table.loc[i, disease_people]
        control_cov_row = cov_table.loc[i, control_people]
        
        # Find coverage threshold
        number_disease_20_cov = disease_cov_row[disease_cov_row >= 20].count()
        number_control_20_cov = control_cov_row[control_cov_row >= 20].count()
        number_disease_15_cov = disease_cov_row[disease_cov_row >= 15].count()
        number_control_15_cov = control_cov_row[control_cov_row >= 15].count()
        number_disease_10_cov = disease_cov_row[disease_cov_row >= 10].count()
        number_control_10_cov = control_cov_row[control_cov_row >= 10].count()
        number_disease_5_cov = disease_cov_row[disease_cov_row >= 5].count()
        number_control_5_cov = control_cov_row[control_cov_row >= 5].count()
        
        if number_disease_20_cov >= min_disease_people and number_control_20_cov >= min_control_people:
            stability_based_on[i] = 20
        elif number_disease_15_cov >= min_disease_people and number_control_15_cov >= min_control_people:
            stability_based_on[i] = 15
        elif number_disease_10_cov >= min_disease_people and number_control_10_cov >= min_control_people:
            stability_based_on[i] = 10
        elif number_disease_5_cov >= min_disease_people_5_cov and number_control_5_cov >= min_control_people_5_cov:
            stability_based_on[i] = 5
        else:
            stability_based_on[i] = float('nan')

        if np.isnan(stability_based_on[i]):
            coverage_threshold_used[i] = 5
        else:
            current_stability_cov = stability_based_on[i]
            stability_disease_mean = disease_edit_row[disease_cov_row >= current_stability_cov].mean()
            stability_control_mean = control_edit_row[control_cov_row >= current_stability_cov].mean()
            
            for j in np.arange(5, stability_based_on[i] + 1e-4, 5):
                disease_mean = disease_edit_row[disease_cov_row >= j].mean()
                control_mean = control_edit_row[control_cov_row >= j].mean()
                if np.absolute(disease_mean - stability_disease_mean) <= stability_value and np.absolute(control_mean - stability_control_mean) <= stability_value:
                    coverage_threshold_used[i] = j
                    break
        
        # Calculate statistics
        stable_min_cov = coverage_threshold_used[i]
        disease_adju_edit_row = disease_edit_row[np.logical_and(np.logical_and((~np.isnan(disease_edit_row)), (~np.isnan(disease_cov_row))), (disease_cov_row >= stable_min_cov))]
        disease_adju_cov_row = disease_cov_row[np.logical_and((~np.isnan(disease_cov_row)), (disease_cov_row >= stable_min_cov))]
        control_adju_edit_row = control_edit_row[np.logical_and(np.logical_and((~np.isnan(control_edit_row)), (~np.isnan(control_cov_row))), (control_cov_row >= stable_min_cov))]
        control_adju_cov_row = control_cov_row[np.logical_and((~np.isnan(control_cov_row)), (control_cov_row >= stable_min_cov))]
        
        stable_mean_disease_editing_level[i] = disease_adju_edit_row.mean()
        stable_std_dev_disease_editing_level[i] = disease_adju_edit_row.std()
        stable_mean_control_editing_level[i] = control_adju_edit_row.mean()
        stable_std_dev_control_editing_level[i] = control_adju_edit_row.std()
        stable_number_disease_with_at_least_min_coverage[i] = disease_adju_cov_row[disease_adju_cov_row >= stable_min_cov].count()
        stable_number_disease_nonzero_editing_and_min_coverage[i] = disease_adju_cov_row[(~np.isnan(disease_adju_cov_row)) & (disease_adju_cov_row >= stable_min_cov) & (disease_adju_edit_row > 0)].count()
        stable_disease_prevalence[i] = stable_number_disease_nonzero_editing_and_min_coverage[i] / stable_number_disease_with_at_least_min_coverage[i]
        stable_number_control_with_at_least_min_coverage[i] = control_adju_cov_row[control_adju_cov_row >= stable_min_cov].count()
        stable_number_control_nonzero_editing_and_min_coverage[i] = control_adju_cov_row[(~np.isnan(control_adju_cov_row)) & (control_adju_cov_row >= stable_min_cov) & (control_adju_edit_row > 0)].count()
        stable_control_prevalence[i] = stable_number_control_nonzero_editing_and_min_coverage[i] / stable_number_control_with_at_least_min_coverage[i]
        stable_total_number_individuals_nonzero_editing_and_min_coverage[i] = (stable_number_disease_nonzero_editing_and_min_coverage[i] + stable_number_control_nonzero_editing_and_min_coverage[i]).sum()
        
        if (len(disease_adju_edit_row) >= 1) & (len(control_adju_edit_row) >= 1):
            if (np.all(disease_adju_edit_row.values == control_adju_edit_row.values)):
                stable_mann_whitney_p_value[i] = float('nan')
            else:
                temp, stable_mann_whitney_p_value[i] = mannwhitneyu(disease_adju_edit_row, control_adju_edit_row, alternative='two-sided')
        else:
            stable_mann_whitney_p_value[i] = float('nan')
        
        stable_editing_level_effect_size[i] = np.absolute(stable_mean_disease_editing_level[i] - stable_mean_control_editing_level[i])
        fisher_matrix = np.matrix([[stable_number_disease_nonzero_editing_and_min_coverage[i], stable_number_disease_with_at_least_min_coverage[i] - stable_number_disease_nonzero_editing_and_min_coverage[i]], [stable_number_control_nonzero_editing_and_min_coverage[i], stable_number_control_with_at_least_min_coverage[i] - stable_number_control_nonzero_editing_and_min_coverage[i]]])
        stable_frequency_OR[i], stable_frequency_fishers_p_value[i] = fisher_exact(fisher_matrix)
        stable_prevalence_effect_size[i] = np.absolute(stable_disease_prevalence[i] - stable_control_prevalence[i])

    # Build output table
    header_info = editing_table[['chromosome', 'position', 'type_editing']]
    stats_table = pd.DataFrame(coverage_threshold_used)
    stats_table = stats_table.rename(columns={stats_table.columns[0]: 'coverage_threshold_used'})
    stats_table['stability_based_on'] = pd.DataFrame(stability_based_on)
    stats_table['stable_mean_disease_editing_level'] = pd.DataFrame(stable_mean_disease_editing_level)
    stats_table['stable_std_dev_disease_editing_level'] = pd.DataFrame(stable_std_dev_disease_editing_level)
    stats_table['stable_mean_control_editing_level'] = pd.DataFrame(stable_mean_control_editing_level)
    stats_table['stable_std_dev_control_editing_level'] = pd.DataFrame(stable_std_dev_control_editing_level)
    stats_table['stable_number_disease_with_at_least_min_coverage'] = pd.DataFrame(stable_number_disease_with_at_least_min_coverage)
    stats_table['stable_number_disease_nonzero_editing_and_min_coverage'] = pd.DataFrame(stable_number_disease_nonzero_editing_and_min_coverage)
    stats_table['stable_disease_prevalence'] = pd.DataFrame(stable_disease_prevalence)
    stats_table['stable_number_control_with_at_least_min_coverage'] = pd.DataFrame(stable_number_control_with_at_least_min_coverage)
    stats_table['stable_number_control_nonzero_editing_and_min_coverage'] = pd.DataFrame(stable_number_control_nonzero_editing_and_min_coverage)
    stats_table['stable_control_prevalence'] = pd.DataFrame(stable_control_prevalence)
    stats_table['stable_total_number_individuals_nonzero_editing_and_min_coverage'] = pd.DataFrame(stable_total_number_individuals_nonzero_editing_and_min_coverage)
    stats_table['stable_mann_whitney_p_value'] = pd.DataFrame(stable_mann_whitney_p_value)
    stats_table['stable_editing_level_effect_size'] = pd.DataFrame(stable_editing_level_effect_size)
    stats_table['stable_frequency_fishers_p_value'] = pd.DataFrame(stable_frequency_fishers_p_value)
    stats_table['stable_frequency_OR'] = pd.DataFrame(stable_frequency_OR)
    stats_table['stable_prevalence_effect_size'] = pd.DataFrame(stable_prevalence_effect_size)

    full_table = pd.concat([header_info, stats_table, editing_table[all_people]], axis=1)
    full_table.to_csv(output_file, sep='\t', index=False)

    print "Linear model job completed\n"

# ============================================================================
# UTILITY FUNCTIONS (UNCHANGED ALGORITHMS)
# ============================================================================

def Set_Chr_Nr(Chr):
    """Sort by chromosome"""
    if Chr:
        New = Chr.lstrip('chr').split('_')[0]
        if New == 'X':
            New = 23
        elif New == 'Y':
            New = 24
        elif New == 'M':
            New = 25
        else:
            New = int(New)
    else:
        New = 0
    return New

def Sample_percentage(row):
    """Percentage of samples from each type"""
    percentage = (len(filter(lambda x: x != '-', row)) / float(len(row))) * 100
    return round(percentage)

def Sample_count(row):
    """Number of samples from each type"""
    count = len(filter(lambda x: x != '-', row))
    return count

def get_bh(pvalue, siglevel):
    """Benjamini-Hochberg correction"""
    pvalue.sort()
    x = 1
    y = 0
    p = 0
    for i in pvalue:
        nf = i[0] * len(pvalue)
        fdr = nf / x
        if fdr <= siglevel:
            i[1].append('True')
            p = i[0]
            y += 1
        else:
            i[1].append('False')
        x += 1
    return pvalue, y, p

def get_b(pvalue, siglevel):
    """Bonferroni correction"""
    pvalue.sort()
    y = 0
    pp = 1.0
    for i in pvalue:
        p = i[0] * len(pvalue)
        if p <= siglevel:
            i[1].append('True')
            y += 1
            if p < pp:
                pp = p
        else:
            i[1].append('False')
    return pvalue, y, pp

def only_sig(row_a, row):
    """Returns only significant events"""
    if (row_a[-1] != '-' and row_a[-1] != 0.0 and row_a[-1] <= 0.05):
        row = row[0].split('_') + row[2:]
        row.insert(2, 'A.to.G')
        print '\t'.join(map(str, row))

def tuple_replace(i):
    if type(i) == tuple:
        return i[0]
    else:
        return i

def tuple_replace_bis(k):
    if type(k) == tuple:
        return k[1]
    else:
        return k

def remove_underscore(lis):
    lis = lis[:lis.index('_')]
    return lis

# ============================================================================
# LOAD METADATA AND VALIDATE PATHS
# ============================================================================

sample_informations = {}
sample_paths = {}

print 'Validating sample files...'
with open(samples_informations_file, 'r') as f:
    line_num = 0
    for line in f:
        line_num += 1
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        
        parts = map(str.strip, line.split(','))
        
        if len(parts) != 3:
            sys.exit('ERROR: Line %d in metadata file must have 3 columns: SampleName,Condition,Path\nGot: %s' % (line_num, line))
        
        sample_name, condition, filepath = parts
        
        # Validate condition
        if condition.upper() not in ['CTRL', 'DIS']:
            sys.exit('ERROR: Line %d - Condition must be CTRL or DIS, got: %s' % (line_num, condition))
        
        # Validate file path
        if not os.path.exists(filepath):
            sys.exit('ERROR: Line %d - File not found: %s' % (line_num, filepath))
        
        sample_informations[sample_name] = condition.upper()
        sample_paths[sample_name] = filepath
        print '  [OK] %s (%s): %s' % (sample_name, condition.upper(), filepath)

if len(sample_informations) == 0:
    sys.exit('ERROR: No valid samples found in metadata file')

num_ctrl = len([x for x in sample_informations.values() if x == 'CTRL'])
num_dis = len([x for x in sample_informations.values() if x == 'DIS'])

print '\nSummary:'
print '  Total samples: %d' % len(sample_informations)
print '  Control samples: %d' % num_ctrl
print '  Disease samples: %d' % num_dis
print ''

if num_ctrl == 0 or num_dis == 0:
    sys.exit('ERROR: Need at least one CTRL and one DIS sample')

# ============================================================================
# LOAD EDITING DATA FROM OUTTABLES
# ============================================================================

print 'Loading editing sites from outTable files...'

all_available_sites = []
sample_edited_sites = {}

for sample_name in sorted(sample_informations.keys()):
    table_path = sample_paths[sample_name]
    sites_found = 0
    
    with open(table_path, 'r') as a:
        for line in a:
            if line.startswith('chr'):
                s = map(str.strip, line.split("\t"))
                if s[7] == 'AG':
                    site = s[0] + "_" + s[1]
                    freq = s[8]
                    coverage = s[4]
                    freq_gnum_cov = '%s^%s^%s' % (s[8], eval(s[6])[2], s[4])
                    
                    if site not in all_available_sites:
                        all_available_sites.append(site)
                    
                    if (int(coverage) >= min_coverage) and (float(freq) >= min_edit_frequency):
                        sample_edited_sites.setdefault((sample_name, site), []).append((freq, freq_gnum_cov))
                        sites_found += 1
    
    print '  %s: %d AG editing sites (coverage>=%d, freq>=%.2f)' % (sample_name, sites_found, min_coverage, min_edit_frequency)

print '\nTotal unique editing sites across all samples: %d\n' % len(all_available_sites)

# ============================================================================
# PREPARE OUTPUT COLUMNS
# ============================================================================

table_columns = map(lambda x: x + '_' + sample_informations[x], sorted(sample_informations.keys()))
disease = [i for i in table_columns if i.upper().find('DIS') != -1]
controls = [i for i in table_columns if i.upper().find('CTRL') != -1]

# ============================================================================
# LINEAR MODEL MODE
# ============================================================================

if enable_linear_model:
    print 'Running linear model analysis...\n'
    outtable = ''
    header = ['chromosome', 'position', 'type_editing'] + map(remove_underscore, controls) + map(remove_underscore, disease)
    outtable += '\t'.join(header)
    outtable += '\n'
    
    for chrom in sorted(all_available_sites, key=lambda x: Set_Chr_Nr(x)):
        row = [chrom]
        for col in header[2:]:
            row.append(sample_edited_sites.get((col.split('_')[0], chrom), ['-'])[0])
        ctrls = zip(*(zip(controls, row[1:])))[1]
        dss = zip(*(zip(disease, row[len(ctrls) + 1:])))[1]
        ctrls_freq = map(tuple_replace, ctrls)
        dss_freq = map(tuple_replace, dss)
        row.append(str([Sample_count(ctrls), Sample_count(dss)]))

        row_b = map(tuple_replace_bis, row)
        row_b = row_b[0].split('_') + row_b[2:]
        row_b.insert(2, 'A.to.G')
        final_list = row_b[:-1]
        outtable += '\t'.join(map(str, final_list)).replace('-', 'NA')
        outtable += '\n'

    with open('temp.csv', 'w') as t:
        t.write(outtable)

    call_differential_editing_sites(samples_informations_file)

# ============================================================================
# MANN-WHITNEY MODE (DEFAULT)
# ============================================================================

else:
    print 'Running Mann-Whitney U test analysis...\n'
    
    header = ['chromosome', 'position', 'type_editing'] + controls + disease + ['[num_controls/num_disease]'] + ['delta_diff'] + ['pvalue (Mannwhitney)']

    if pvalue_correction == 1:
        header += ['pvalue Bonferroni corrected']
    if pvalue_correction == 2:
        header += ['pvalue BH corrected']

    print '\t'.join(header)

    for chrom in sorted(all_available_sites, key=lambda x: Set_Chr_Nr(x)):
        row = [chrom]
        for col in header[3:header.index('[num_controls/num_disease]')]:
            row.append(sample_edited_sites.get((col.split('_')[0], chrom), ['-'])[0])
        ctrls = zip(*(zip(controls, row[1:])))[1]
        dss = zip(*(zip(disease, row[len(ctrls) + 1:])))[1]
        ctrls_freq = map(tuple_replace, ctrls)
        dss_freq = map(tuple_replace, dss)
        row.append(str([Sample_count(ctrls), Sample_count(dss)]))
        
        if (Sample_percentage(ctrls) >= min_sample_testing) and (Sample_percentage(dss) >= min_sample_testing):
            ctrls_mean = sum(map(float, filter(lambda x: x != '-', ctrls_freq))) / len(filter(lambda x: x != '-', ctrls_freq))
            dss_mean = sum(map(float, filter(lambda x: x != '-', dss_freq))) / len(filter(lambda x: x != '-', dss_freq))
            delta_diff = abs(ctrls_mean - dss_mean)
            
            # Check if all values are identical (Mann-Whitney will fail)
            all_values = filter(lambda x: x != '-', ctrls_freq + dss_freq)
            if len(set(all_values)) == 1:
                # All values are identical, no difference
                row.append(round(delta_diff, 3))
                row.append('-')
                if pvalue_correction == 1:
                    row.append('-')
                elif pvalue_correction == 2:
                    row.append('-')
            else:
                try:
                    pvalue = stats.mannwhitneyu(ctrls_freq, dss_freq, alternative='two-sided')
                    row.append(round(delta_diff, 3))
                    row.append(str(round(pvalue[1], 3)))
                    correction_argmnt = [(pvalue[1], ctrls_freq + dss_freq)]

                    if pvalue_correction == 1:
                        row.append(round(get_b(correction_argmnt, 0.05)[-1], 6))
                    elif pvalue_correction == 2:
                        row.append(round(get_bh(correction_argmnt, 0.05)[-1], 6))
                except ValueError:
                    # Fallback for any other Mann-Whitney errors
                    row.append(round(delta_diff, 3))
                    row.append('-')
                    if pvalue_correction == 1:
                        row.append('-')
                    elif pvalue_correction == 2:
                        row.append('-')
        else:
            if pvalue_correction == 0:
                row += ['-', '-']
            else:
                row += ['-', '-', '-']
        
        row_a = map(tuple_replace, row)
        row_b = map(tuple_replace_bis, row)
        
        if pvalue_correction != 0 and only_significants == 'yes':
            only_sig(row_a, row_b)
        else:
            row_b = row_b[0].split('_') + row_b[2:]
            row_b.insert(2, 'A.to.G')
            print '\t'.join(map(str, row_b))

