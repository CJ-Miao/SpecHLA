#!/usr/bin/env python3
"""
Convert allelefrequencies.net per-population HLA data into SpecHLA format.
Reads A.tsv, B.tsv, C.tsv, DPA1.tsv, DPB1.tsv, DQA1.tsv, DQB1.tsv, DRB1.tsv
and produces HLA_FREQ_HLA_I_II.txt with Caucasian/Black/Asian aggregated frequencies.
"""
import os
import glob
import shutil
from datetime import datetime

def classify_population(pop):
    """Classify a population name into Caucasian, Black, or Asian."""
    p = pop.lower()

    african_regions = ['african american', 'african-american']
    african_countries = [
        'algeria', 'angola', 'benin', 'botswana', 'burkina faso', 'burundi',
        'cameroun', 'cameroon', 'cape verde', 'central african', 'chad',
        'congo', "cote d'ivoire", 'djibouti', 'egypt', 'equatorial guinea',
        'eritrea', 'ethiopia', 'gabon', 'gambia', 'ghana', 'guinea',
        'guinea-bissau', 'ivory coast', 'kenya', 'lesotho', 'liberia',
        'libya', 'madagascar', 'malawi', 'mali', 'mauritania', 'mauritius',
        'morocco', 'mozambique', 'namibia', 'niger', 'nigeria',
        'rwanda', 'senegal', 'seychelles', 'sierra leone', 'somalia',
        'south africa', 'sudan', 'swaziland', 'tanzania', 'togo',
        'tunisia', 'uganda', 'zaire', 'zambia', 'zimbabwe',
    ]
    african_keywords = [
        'pygmy', 'bantu', 'zulu', 'xhosa', 'sotho', 'tswana',
        'masai', 'maasai', 'fulani', 'mossi', 'rimaibe',
        'bamileke', 'beti', 'saa', 'sawa', 'yaounde',
        'baka', 'bakola', 'mende', 'temne', 'mandinka', 'wolof',
        'black', 'quilombol', 'african',
    ]

    asian_countries = [
        'afghanistan', 'bangladesh', 'bhutan', 'brunei', 'cambodia',
        'china', 'india', 'indonesia', 'japan', 'korea', 'laos',
        'malaysia', 'maldives', 'myanmar', 'nepal', 'pakistan',
        'philippines', 'singapore', 'sri lanka', 'taiwan', 'thailand',
        'timor', 'vietnam',
    ]
    asian_regions = [
        'south east asia', 'south-east asia', 'southeast asia',
        'east asia', 'central asia', 'south asia', 'middle east',
        'persian', 'caucasus',
    ]
    # Middle East → Asian
    middle_east = [
        'iran', 'israel', 'jordan', 'lebanon', 'oman', 'saudi',
        'gaza', 'palestin', 'druze', 'yemen', 'iraq',
    ]
    # Asian additional keywords
    asian_extra = [
        'chinese', 'japanese', 'korean', 'vietnamese', 'thai',
        'filipino', 'indonesian', 'malaysian', 'indian', 'pakistani',
        'bangladeshi', 'burmese', 'khmer', 'sinhalese', 'tamils',
        'dravidian', 'hazara', 'pashtun', 'punjabi', 'gujarati',
        'marathi', 'bengali', 'sinhala',
        # Middle East specific
        'parsi', 'zoroastrian', 'azeri', 'ashkenazi', 'jewish',
        'jews', 'kurd', 'kurdish',
        # Mongolia
        'mongolia', 'mongol', 'buriat', 'hoton', 'khalkha', 'oold',
        'tarialan',
        # Borneo / Southeast Asia
        'borneo', 'bandjarmasin', 'okinawa',
    ]

    # European
    european_countries = [
        'albania', 'andorra', 'armenia', 'austria', 'azerbaijan',
        'belarus', 'belgium', 'bosnia', 'bulgaria', 'croatia', 'cyprus',
        'czech', 'denmark', 'estonia', 'finland', 'france', 'georgia',
        'germany', 'greece', 'hungary', 'iceland', 'ireland', 'italy',
        'kosovo', 'latvia', 'liechtenstein', 'lithuania', 'luxembourg',
        'macedonia', 'malta', 'moldova', 'monaco', 'montenegro',
        'netherlands', 'norway', 'poland', 'portugal', 'romania',
        'russia', 'serbia', 'slovakia', 'slovenia', 'spain', 'sweden',
        'switzerland', 'turkey', 'ukraine', 'united kingdom', 'uk',
        'wales', 'scotland', 'england', 'britain', 'british',
        'azores', 'madeira', 'faroe', 'vatican',
    ]
    european_keywords = [
        'caucasian', 'european', 'roman', 'basque', 'sami',
        'romani', 'gypsy', 'roma',
    ]

    # Latin American (default Caucasian unless mixed)
    latin_american = [
        'argentina', 'bolivia', 'brazil', 'chile', 'colombia',
        'costa rica', 'cuba', 'dominican', 'ecuador', 'el salvador',
        'guatemala', 'honduras', 'mexico', 'nicaragua', 'panama',
        'paraguay', 'peru', 'puerto rico', 'uruguay', 'venezuela',
        'jamaica', 'martinique',
    ]

    # Skip: Indigenous American
    indigenous_american = [
        'mapuche', 'toba', 'wichi', 'mataco', 'chiriguano',
        'tehuelche', 'kaingang', 'guarani', 'xavante',
        'cree', 'chipewyan', 'athabaskan', 'penutian',
        'aborigine', 'aboriginal', 'indigenous american',
        'native american', 'inuit', 'maya', 'aztec', 'inc',
        'quechua', 'aymara', 'navajo', 'apache', 'sioux',
        'pima', 'yupik', 'aleut', 'ojibwa',
    ]

    # Skip: Oceania/Pacific
    oceania = [
        'papua', 'cook island', 'tonga', 'nauru', 'kiribati', 'niue',
        'tokelau', 'new caledonia', 'polynesia', 'melanesia',
        'micronesia', 'new zealand', 'maori', 'fiji', 'samoa',
        'american samoa',
    ]

    # Black first
    for kw in african_regions + african_countries + african_keywords:
        if kw in p:
            return 'Black'

    # Asian (countries + regions + keywords + Middle East)
    for kw in asian_countries + asian_regions + asian_extra:
        if kw in p:
            return 'Asian'

    # Caucasian (Europe)
    for kw in european_countries + european_keywords:
        if kw in p:
            return 'Caucasian'

    # Latin American (skip mixed)
    for kw in latin_american:
        if kw in p:
            if 'pard' in p or 'mixed' in p or 'mestizo' in p or 'hispanic' in p:
                return None
            return 'Caucasian'

    # Skip: Indigenous American / Oceania
    for kw in indigenous_american:
        if kw in p:
            return None
    for kw in oceania:
        if kw in p:
            return None

    # USA-specific
    if 'usa' in p:
        if 'asian' in p:
            return 'Asian'
        if 'hispanic' in p or 'mestizo' in p or 'chicano' in p:
            return None
        if 'hawaiian' in p or 'pacific islander' in p:
            return None
        if 'native' in p or 'alaska' in p:
            return None
        # USA general (Olmsted, San Diego etc.) → Caucasian
        return 'Caucasian'

    return None


def main():
    input_dir = '/home/miaocj/docker_dir/amplicon/allelefrequencies/allelefrequencies.net/hla'
    output_file = '/home/miaocj/docker_dir/biosoft/SpecHLA/db/HLA_FREQ_HLA_I_II.txt'
    old_file = '/home/miaocj/docker_dir/biosoft/SpecHLA/db/old_HLA_FREQ_HLA_I_II.txt'

    # Rename old file
    if os.path.exists(output_file):
        shutil.move(output_file, old_file)
        print(f'Moved old file to: {old_file}')

    files = sorted(glob.glob(os.path.join(input_dir, '*.tsv')))

    # {allele: {group: [(freq, weight), ...]}}
    data = {}

    for fpath in files:
        with open(fpath, 'r') as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue
                allele = parts[0]
                pop = parts[1]
                freq_str = parts[3]
                n_str = parts[4]

                if not freq_str or freq_str.strip() == '' or freq_str == '0':
                    continue
                try:
                    freq = float(freq_str)
                except ValueError:
                    continue
                try:
                    n = int(n_str) if n_str and n_str.strip() else None
                except ValueError:
                    n = None

                group = classify_population(pop)
                if group is None:
                    continue

                if allele not in data:
                    data[allele] = {'Caucasian': [], 'Black': [], 'Asian': []}

                weight = n if n and n > 0 else 1
                data[allele][group].append((freq, weight))

    # Aggregate weighted average per allele per group
    results = []
    for allele in sorted(data.keys()):
        groups = data[allele]
        freqs = {}
        for group_name in ['Caucasian', 'Black', 'Asian']:
            entries = groups[group_name]
            if not entries:
                freqs[group_name] = 0
            else:
                total_w = sum(w for _, w in entries)
                weighted_avg = sum(f * w for f, w in entries) / total_w if total_w > 0 else 0
                freqs[group_name] = weighted_avg
        results.append((allele, freqs['Caucasian'], freqs['Black'], freqs['Asian']))

    # Sort by gene order
    gene_order = {'A*': 1, 'B*': 2, 'C*': 3, 'DPA1*': 4, 'DPB1*': 5,
                  'DQA1*': 6, 'DQB1*': 7, 'DRB1*': 8}

    def sort_key(item):
        allele = item[0]
        prefix = '*'.join(allele.split('*')[:1]) + '*'
        order = gene_order.get(prefix, 99)
        return (order, allele)

    results.sort(key=sort_key)

    # Write output
    with open(output_file, 'w') as f:
        f.write('\t'.join(['Allele', 'Caucasian', 'Black', 'Asian']) + '\n')
        for allele, cauc, black, asian in results:
            f.write(f'{allele}\t{cauc}\t{black}\t{asian}\n')

    print(f'Wrote {len(results)} alleles to {output_file}')
    genes = {}
    for allele, _, _, _ in results:
        gene = '*'.join(allele.split('*')[:1])
        genes[gene] = genes.get(gene, 0) + 1
    for gene in sorted(genes.keys()):
        print(f'  {gene}: {genes[gene]} alleles')


if __name__ == '__main__':
    main()
