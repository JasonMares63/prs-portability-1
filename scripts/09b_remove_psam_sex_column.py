import csv


original = open('data/ukb_populations/EUR_all.psam', 'r')
filtered = open('data/ukb_populations/EUR_all_nosex.psam', 'w')

reader = csv.reader(original, delimiter = '\t')
writer = csv.writer(filtered, delimiter = '\t')
    
for line in reader:
    writer.writerow(line[0:2])

original.close()
filtered.close()

# Filter the covariates file to remove age * sex, age^2, age^2 * sex
original = open('data/ukb_filtered/covar_all_samples.covar', 'r')
filtered = open('data/ukb_filtered/covar_all_samples_filtered.covar', 'w')

reader = csv.DictReader(original, delimiter = '\t')
fieldnames = reader.fieldnames.copy()
[fieldnames.remove(field) for field in ['age_sq', 'age_sex', 'age_sq_sex']]
writer = csv.DictWriter(filtered, delimiter = '\t', fieldnames = fieldnames)
writer.writeheader()

#FID IID sex age age_sq age_sex age_sq_sex PC1_AVG PC2_AVG PC3_AVG PC4_AVG PC5_AVG PC6_AVG PC7_AVG PC8_AVG PC9_AVG PC10_AVG PC11_AVG PC12_AVG PC13_AVG PC14_AVG PC15_AVG PC16_AVG PC17_AVG PC18_AVG PC19_AVG PC20_AVG    
for line in reader:
    [line.pop(field) for field in ['age_sq', 'age_sex', 'age_sq_sex']]
    writer.writerow(line)

original.close()
filtered.close()


