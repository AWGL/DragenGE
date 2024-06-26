
import csv
import sys
import glob

ped = sys.argv[1]
seq_id = sys.argv[2]

fam_dict = {}
no_fam_int = 0

# We need the metrics files (Dragen v3.7) and target bed coverage metrics files (Dragen v3.10.8) for excluding low coverage samples
metrics_files = glob.glob('*/*.mapping_metrics.csv')
target_files = glob.glob('*/*target_bed_coverage_metrics.csv')

min_depth = 5
sample_dict = {}

# Check if we need to be using metrics or target bed metrics file for depth, using just the first file
if 'Average sequenced coverage over target region' in open(metrics_files[0]).read():
	for coverage_file in metrics_files:
		with open(coverage_file) as csvfile:
			spamreader = csv.reader(csvfile, delimiter=',')
			for row in spamreader:
				key = row[2]
				value = row[3]
				if key == 'Average sequenced coverage over target region':
					if float(value) > min_depth:

						if value == 'NA':

							value = 0.0

						sample_id = coverage_file.split('/')[0]
						sample_dict[sample_id] = sample_id
						break
else:
	for coverage_file in target_files:
		with open(coverage_file) as csvfile:
			spamreader = csv.reader(csvfile, delimiter=',')
			for row in spamreader:
				key = row[2]
				value = row[3]
				if key == 'Average alignment coverage over target region':
					
					if value == 'NA':

						value = 0.0
					if float(value) > min_depth:
						sample_id = coverage_file.split('/')[0]
						sample_dict[sample_id] = sample_id
						break

with open(ped) as csvfile:
	spamreader = csv.reader(csvfile, delimiter='\t')
	for row in spamreader:
		fam_id = row[0]
		sample_name = row[1]
		if fam_id == 0 or fam_id == '0':
			fam_dict[f'singleton_{no_fam_int}'] = [sample_name]
			no_fam_int += 1
		else:
			if fam_id not in fam_dict:
				fam_dict[fam_id] = [sample_name]
			else:
				fam_dict[fam_id].append(sample_name)

for key in fam_dict:
	out_file = f'{key}_for_sv.family'
	family_rows = []
	for sample in fam_dict[key]:
		if sample in sample_dict:
			family_rows.append(f'--bam-input {sample}/{seq_id}_{sample}.bam \\')
	if len(family_rows) > 0:
		with open(out_file, 'w') as csvfile:
			spamwriter = csv.writer(csvfile, delimiter='\t', lineterminator='\n')
			for row in family_rows:
				spamwriter.writerow([row])

