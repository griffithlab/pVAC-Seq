import argparse
import sys
import re
import csv

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('pvacseq binding_filter')
    parser.add_argument('input_file', type=argparse.FileType('r'),
                        help="Combined parsed epitope file")
    parser.add_argument('output_file', type=argparse.FileType('w'),
                        help="Output .tsv file containing list of filtered " +
                        "epitopes based on binding affinity")
    parser.add_argument('-c', '--minimum-fold-change', type=int,
                        help="Minimum fold change between mutant binding " +
                        "score and wild-type score. The default is 0, which " +
                        "filters no results, but 1 is often a sensible " +
                        "default (requiring that binding is better to the MT " +
                        "than WT)",
                        default=0)
    parser.add_argument('-b', '--binding-threshold', type=int,
                        help="Report only epitopes where the mutant allele " +
                        "has ic50 binding scores below this value; default 500",
                        default=500)
    parser.add_argument('-m', '--top-score-metric',
                        choices=['lowest', 'median'],
                        default='median',
                        help="The ic50 scoring metric to use when filtering epitopes by binding-threshold. " +
                        "lowest: Best MT Score - lowest MT ic50 binding score of all chosen prediction methods. " +
                        "median: Median MT Score All Methods - median MT ic50 binding score of all chosen prediction methods. " +
                        "Default: median")

    args = parser.parse_args(args_input)

    prediction = {}
    fieldnames = []

    reader = csv.DictReader(args.input_file, delimiter='\t')
    fieldnames = reader.fieldnames

    for entry in reader:
        name = entry['Gene Name']
        if args.top_score_metric == 'median':
            score = float(entry['Median MT Score All Methods'])
        elif args.top_score_metric == 'lowest':
            score = float(entry['Best MT Score'])
        fold_change = sys.maxsize if entry['Fold Change'] == 'NA' else float(entry['Fold Change'])

        if score > args.binding_threshold or fold_change < args.minimum_fold_change:
            continue

        if (name not in prediction or
                score < prediction[name]['SCORE']):
            prediction[name] = {
                'GENES' : [entry],
                'SCORE' : score
            }
        elif score == prediction[name]['SCORE']:
            prediction[name]['GENES'].append(entry)
    args.input_file.close()

    writer = csv.DictWriter(
        args.output_file,
        fieldnames,
        delimiter = '\t',
        lineterminator = '\n'
    )

    writer.writeheader()

    writer.writerows(
        entry
        for gene in sorted(prediction)
        for entry in prediction[gene]['GENES']
    )

    args.output_file.close()



if __name__ == "__main__":
    main()
