import argparse
import requests
import sys
import re

def main(args_input = sys.argv[1:]):
    parser = argparse.ArgumentParser('pvacseq call_iedb')
    parser.add_argument('input_file', type=argparse.FileType('r'),
                        help="Input FASTA file")
    parser.add_argument('output_file', type=argparse.FileType('w'),
                        help="Output file from iedb")
    parser.add_argument('method',
                        choices=['netmhcpan', 'ann', 'smmpmbec', 'smm', 'netmhccons', 'pickpocket'],
                        help="The iedb analysis method to use")
    parser.add_argument('allele',
                        help="Allele for which to make prediction")
    parser.add_argument('epitope_length', type=int, choices=[8,9,10,11,12,13,14,15],
                        help="Length of subpeptides (epitopes) to predict")
    args = parser.parse_args(args_input)

    #Insert * in the allele name
    allele = re.sub(r'(\w*-[\w|\d])(.*)', r'\1*\2', args.allele)

    data = {
        'sequence_text': args.input_file.read(),
        'method':        args.method,
        'allele':        allele,
        'length':        args.epitope_length,
    }

    response = requests.post('http://tools-api.iedb.org/tools_api/mhci/', data=data)
    if response.status_code != 200:
        sys.exit("Error posting request to IEDB.\n%s" % response.text)
    args.output_file.write(response.text)

    args.input_file.close()
    args.output_file.close()

if __name__ == "__main__":
    main()